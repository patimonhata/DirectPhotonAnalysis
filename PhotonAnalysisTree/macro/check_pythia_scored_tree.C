#include <TFile.h>
#include <TKey.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace
{
template <class T>
bool bind_pythia_scored(TTree* tree, const char* name, T* address)
{
  if (!tree->GetBranch(name))
  {
    std::cerr << "check_pythia_scored_tree - missing branch: " << name << std::endl;
    return false;
  }
  return tree->SetBranchAddress(name, address) >= 0;
}

bool has_only_expected_trees(TFile* file)
{
  bool has_event_tree = false;
  bool has_metadata = false;
  bool valid = file && file->GetListOfKeys()->GetSize() == 2;
  TIter next(file ? file->GetListOfKeys() : nullptr);
  while (auto* key = dynamic_cast<TKey*>(next()))
  {
    const std::string name = key->GetName();
    const bool is_tree = std::string(key->GetClassName()) == "TTree";
    has_event_tree |= name == "event_tree" && is_tree;
    has_metadata |= name == "metadata" && is_tree;
    valid &= (name == "event_tree" || name == "metadata") && is_tree;
  }
  return valid && has_event_tree && has_metadata;
}

bool valid_score_vectors(const std::vector<float>* scores,
                         const std::vector<unsigned char>* valid,
                         std::size_t expected,
                         const std::vector<unsigned char>* required_shower_valid,
                         Long64_t& valid_count)
{
  if (!scores || !valid || scores->size() != expected || valid->size() != expected)
  {
    return false;
  }
  for (std::size_t cluster = 0; cluster < expected; ++cluster)
  {
    if ((*valid)[cluster] > 1U ||
        ((*valid)[cluster] && !std::isfinite((*scores)[cluster])) ||
        (required_shower_valid &&
         (required_shower_valid->size() != expected ||
          ((*valid)[cluster] && !(*required_shower_valid)[cluster]))))
    {
      return false;
    }
    valid_count += (*valid)[cluster] ? 1 : 0;
  }
  return true;
}
}

int check_pythia_scored_tree(const char* input_path)
{
  if (!input_path)
  {
    return 1;
  }
  std::unique_ptr<TFile> input(TFile::Open(input_path, "READ"));
  if (!input || input->IsZombie() || !has_only_expected_trees(input.get()))
  {
    std::cerr << "check_pythia_scored_tree - expected only event_tree and metadata TTrees" << std::endl;
    return 2;
  }
  TTree* tree = input->Get<TTree>("event_tree");
  TTree* metadata = input->Get<TTree>("metadata");
  if (!tree || !metadata || metadata->GetEntries() != 1)
  {
    return 2;
  }

  int schema_version = 0;
  std::string* sample_type = nullptr;
  UInt_t metadata_source_file_id = 0;
  ULong64_t n_events_processed = 0;
  ULong64_t n_events_written = 0;
  bool metadata_ok = true;
  metadata_ok &= bind_pythia_scored(metadata, "schema_version", &schema_version);
  metadata_ok &= bind_pythia_scored(metadata, "sample_type", &sample_type);
  metadata_ok &= bind_pythia_scored(metadata, "source_file_id", &metadata_source_file_id);
  metadata_ok &= bind_pythia_scored(metadata, "n_events_processed", &n_events_processed);
  metadata_ok &= bind_pythia_scored(metadata, "n_events_written", &n_events_written);
  metadata_ok &= metadata->GetEntry(0) > 0;
  metadata_ok &= schema_version == 4 && sample_type && *sample_type == "pythia";
  metadata_ok &= n_events_written == static_cast<ULong64_t>(tree->GetEntries());
  metadata_ok &= n_events_processed >= n_events_written;
  if (!metadata_ok)
  {
    std::cerr << "check_pythia_scored_tree - invalid Pythia metadata" << std::endl;
    return 3;
  }
  metadata->ResetBranchAddresses();

  const char* forbidden_nosplit_scores[] = {
      "nosplit_cluster_bdt_base_v3E_score",
      "nosplit_cluster_bdt_base_v3E_valid",
      "nosplit_cluster_p_gamma",
      "nosplit_cluster_p_gamma_valid"};
  for (const char* name : forbidden_nosplit_scores)
  {
    if (tree->GetBranch(name))
    {
      std::cerr << "check_pythia_scored_tree - unexpected no-split score branch: "
                << name << std::endl;
      return 4;
    }
  }

  UInt_t source_file_id = 0;
  UInt_t split_ncluster = 0;
  std::vector<unsigned char>* shower_valid = nullptr;
  std::vector<float>* base_bdt_score = nullptr;
  std::vector<unsigned char>* base_bdt_valid = nullptr;
  std::vector<float>* ppg15_bdt_score = nullptr;
  std::vector<unsigned char>* ppg15_bdt_valid = nullptr;
  std::vector<float>* gamma_score = nullptr;
  std::vector<unsigned char>* gamma_valid = nullptr;
  bool ok = true;
  ok &= bind_pythia_scored(tree, "source_file_id", &source_file_id);
  ok &= bind_pythia_scored(tree, "split_ncluster", &split_ncluster);
  ok &= bind_pythia_scored(tree, "split_cluster_shower_valid", &shower_valid);
  ok &= bind_pythia_scored(tree, "split_cluster_bdt_base_v3E_score", &base_bdt_score);
  ok &= bind_pythia_scored(tree, "split_cluster_bdt_base_v3E_valid", &base_bdt_valid);
  ok &= bind_pythia_scored(tree, "split_cluster_bdt_ppg15v1_score", &ppg15_bdt_score);
  ok &= bind_pythia_scored(tree, "split_cluster_bdt_ppg15v1_valid", &ppg15_bdt_valid);
  ok &= bind_pythia_scored(tree, "split_cluster_p_gamma", &gamma_score);
  ok &= bind_pythia_scored(tree, "split_cluster_p_gamma_valid", &gamma_valid);
  if (!ok)
  {
    return 5;
  }
  const char* score_branches[] = {
      "split_cluster_bdt_base_v3E_score",
      "split_cluster_bdt_base_v3E_valid",
      "split_cluster_bdt_ppg15v1_score",
      "split_cluster_bdt_ppg15v1_valid",
      "split_cluster_p_gamma",
      "split_cluster_p_gamma_valid"};
  for (const char* name : score_branches)
  {
    TBranch* branch = tree->GetBranch(name);
    if (!branch || branch->GetEntries() != tree->GetEntries())
    {
      std::cerr << "check_pythia_scored_tree - branch entry count mismatch: "
                << name << std::endl;
      return 6;
    }
  }

  Long64_t valid_base_bdt = 0;
  Long64_t valid_ppg15_bdt = 0;
  Long64_t valid_gamma = 0;
  for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry)
  {
    if (tree->GetEntry(entry) <= 0 ||
        source_file_id != metadata_source_file_id ||
        !shower_valid || shower_valid->size() != split_ncluster ||
        !valid_score_vectors(base_bdt_score, base_bdt_valid, split_ncluster,
                             shower_valid, valid_base_bdt) ||
        !valid_score_vectors(ppg15_bdt_score, ppg15_bdt_valid, split_ncluster,
                             shower_valid, valid_ppg15_bdt) ||
        !valid_score_vectors(gamma_score, gamma_valid, split_ncluster,
                             nullptr, valid_gamma))
    {
      std::cerr << "check_pythia_scored_tree - malformed event " << entry << std::endl;
      return 7;
    }
  }

  std::cout << "check_pythia_scored_tree - events/valid_base_bdt/"
               "valid_ppg15_bdt/valid_gamma = "
            << tree->GetEntries() << "/" << valid_base_bdt << "/"
            << valid_ppg15_bdt << "/" << valid_gamma << std::endl;
  std::cout << "check_pythia_scored_tree - metadata source/processed/written = "
            << metadata_source_file_id << "/" << n_events_processed << "/"
            << n_events_written << std::endl;
  return 0;
}
