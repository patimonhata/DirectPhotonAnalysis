#include <TFile.h>
#include <TKey.h>
#include <TObject.h>
#include <TTree.h>

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <memory>
#include <regex>
#include <string>
#include <vector>

namespace
{
bool has_clean_layout(TFile& file)
{
  bool has_event_tree = false;
  bool has_metadata = false;
  bool ok = file.GetListOfKeys()->GetSize() == 2;
  TIter next(file.GetListOfKeys());
  while (auto* key = dynamic_cast<TKey*>(next()))
  {
    const std::string name = key->GetName();
    const bool is_tree = std::string(key->GetClassName()) == "TTree";
    has_event_tree |= name == "event_tree" && is_tree;
    has_metadata |= name == "metadata" && is_tree;
    ok &= (name == "event_tree" || name == "metadata") && is_tree;
  }
  return ok && has_event_tree && has_metadata;
}

bool clone_tree(TTree& source, TFile& output, const char* name)
{
  output.cd();
  std::unique_ptr<TTree> clone(source.CloneTree(-1));
  if (!clone)
  {
    return false;
  }
  clone->SetName(name);
  clone->SetDirectory(&output);
  const bool written = clone->Write(name, TObject::kOverwrite) > 0;
  clone->SetDirectory(nullptr);
  return written;
}
}

int copy_analysis_trees(const char* input_path, const char* output_path)
{
  if (!input_path || !output_path)
  {
    return 1;
  }
  std::error_code input_error;
  std::error_code output_error;
  const auto input_canonical = std::filesystem::weakly_canonical(input_path, input_error);
  const auto output_canonical = std::filesystem::weakly_canonical(output_path, output_error);
  if (!input_error && !output_error && input_canonical == output_canonical)
  {
    std::cerr << "Input and output must differ" << std::endl;
    return 2;
  }
  if (std::filesystem::exists(output_path))
  {
    std::cerr << "Refusing to overwrite output: " << output_path << std::endl;
    return 3;
  }
  const std::filesystem::path output_fs(output_path);
  if (!output_fs.parent_path().empty())
  {
    std::filesystem::create_directories(output_fs.parent_path());
  }

  std::unique_ptr<TFile> input(TFile::Open(input_path, "READ"));
  if (!input || input->IsZombie())
  {
    std::cerr << "Cannot open input: " << input_path << std::endl;
    return 4;
  }
  TTree* event_tree = input->Get<TTree>("event_tree");
  TTree* metadata = input->Get<TTree>("metadata");
  if (!event_tree || !metadata)
  {
    std::cerr << "Missing event_tree or metadata in " << input_path << std::endl;
    return 5;
  }
  const Long64_t expected_events = event_tree->GetEntries();
  const Long64_t expected_metadata = metadata->GetEntries();

  std::unique_ptr<TFile> output(TFile::Open(output_path, "RECREATE"));
  if (!output || output->IsZombie())
  {
    return 6;
  }
  if (!clone_tree(*event_tree, *output, "event_tree") ||
      !clone_tree(*metadata, *output, "metadata"))
  {
    output->Close();
    return 7;
  }
  output->Close();
  input->Close();

  std::unique_ptr<TFile> check(TFile::Open(output_path, "READ"));
  TTree* check_event_tree = check ? check->Get<TTree>("event_tree") : nullptr;
  TTree* check_metadata = check ? check->Get<TTree>("metadata") : nullptr;
  const bool ok = check && !check->IsZombie() && has_clean_layout(*check) &&
      check_event_tree && check_metadata &&
      check_event_tree->GetEntries() == expected_events &&
      check_metadata->GetEntries() == expected_metadata &&
      (expected_metadata == 0 || check_metadata->GetEntry(0) > 0);
  if (!ok)
  {
    std::cerr << "Clean-copy validation failed: " << output_path << std::endl;
    return 8;
  }
  return 0;
}

int copy_analysis_directory(const char* input_directory, const char* output_directory)
{
  if (!input_directory || !output_directory)
  {
    return 1;
  }
  const std::filesystem::path input_dir(input_directory);
  const std::filesystem::path output_dir(output_directory);
  if (!std::filesystem::is_directory(input_dir))
  {
    std::cerr << "Input directory does not exist: " << input_dir << std::endl;
    return 2;
  }
  std::filesystem::create_directories(output_dir);

  const std::regex selected(R"(^photon_analysis_tree_[0-9]{6}(_with_bdt|_scored)\.root$)");
  std::vector<std::filesystem::path> inputs;
  for (const auto& entry : std::filesystem::directory_iterator(input_dir))
  {
    if (entry.is_regular_file() &&
        std::regex_match(entry.path().filename().string(), selected))
    {
      inputs.push_back(entry.path());
    }
  }
  std::sort(inputs.begin(), inputs.end());
  if (inputs.empty())
  {
    std::cerr << "No BDT/scored ROOT files found in " << input_dir << std::endl;
    return 3;
  }

  std::size_t copied = 0;
  for (const auto& input : inputs)
  {
    const std::filesystem::path output = output_dir / input.filename();
    const int status = copy_analysis_trees(input.c_str(), output.c_str());
    if (status != 0)
    {
      std::cerr << "Failed after " << copied << " files: " << input
                << " (status " << status << ")" << std::endl;
      return 10 + status;
    }
    ++copied;
    if (copied % 100 == 0 || copied == inputs.size())
    {
      std::cout << "copy_analysis_directory - copied " << copied << "/"
                << inputs.size() << std::endl;
    }
  }
  return 0;
}
