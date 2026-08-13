#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TSystem.h>
#include <TTree.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include <vector>

namespace
{
template <class T>
bool bind_branch(TTree* tree, const char* name, T* address)
{
  if (!tree->GetBranch(name))
  {
    std::cerr << "AccumulatePythiaTruthPtSpectra - missing branch: " << name << std::endl;
    return false;
  }
  tree->SetBranchStatus(name, true);
  return tree->SetBranchAddress(name, address) >= 0;
}

template <class T>
bool aligned(const std::vector<T>* values, const std::size_t expected)
{
  return values && values->size() == expected;
}

bool in_acceptance(const float eta, const double max_abs_eta)
{
  return max_abs_eta < 0.0 || std::abs(static_cast<double>(eta)) < max_abs_eta;
}

std::string join_path(const std::string& directory, const std::string& basename)
{
  if (directory.empty() || directory == ".")
  {
    return basename;
  }
  return directory.back() == '/' ? directory + basename : directory + "/" + basename;
}

bool valid_suffix(const std::string& suffix)
{
  return suffix.size() > 5 && suffix.compare(suffix.size() - 5, 5, ".root") == 0 &&
      suffix.find('/') == std::string::npos && suffix.find('\\') == std::string::npos;
}
}

int AccumulatePythiaTruthPtSpectra(
    const std::string manifest_path,
    const std::string input_directory,
    const Long64_t manifest_begin,
    const Long64_t manifest_end,
    const std::string output_file,
    const int n_bins = 100,
    const double pt_max = 20.0,
    const double max_abs_eta = 0.7,
    const bool use_event_weight = false)
{
  if (manifest_path.empty() || input_directory.empty() || output_file.empty() ||
      manifest_begin < 0 || manifest_end <= manifest_begin || n_bins <= 0 ||
      !std::isfinite(pt_max) || pt_max <= 0.0 || !std::isfinite(max_abs_eta))
  {
    std::cerr << "AccumulatePythiaTruthPtSpectra - invalid argument" << std::endl;
    return 1;
  }

  std::ifstream manifest(manifest_path);
  if (!manifest)
  {
    std::cerr << "AccumulatePythiaTruthPtSpectra - cannot open manifest: "
              << manifest_path << std::endl;
    return 1;
  }

  std::vector<std::string> suffixes;
  std::set<std::string> unique_suffixes;
  std::string line;
  Long64_t row = 0;
  while (std::getline(manifest, line))
  {
    if (!line.empty() && line.back() == '\r')
    {
      line.pop_back();
    }
    if (row >= manifest_begin && row < manifest_end)
    {
      if (!valid_suffix(line) || !unique_suffixes.insert(line).second)
      {
        std::cerr << "AccumulatePythiaTruthPtSpectra - invalid or duplicate suffix at row "
                  << row << ": " << line << std::endl;
        return 2;
      }
      suffixes.push_back(line);
    }
    ++row;
    if (row >= manifest_end)
    {
      break;
    }
  }
  const Long64_t requested_files = manifest_end - manifest_begin;
  if (static_cast<Long64_t>(suffixes.size()) != requested_files)
  {
    std::cerr << "AccumulatePythiaTruthPtSpectra - manifest has fewer rows than requested: "
              << suffixes.size() << "/" << requested_files << std::endl;
    return 2;
  }

  TChain tree("event_tree");
  for (const std::string& suffix : suffixes)
  {
    const std::string tag = suffix.substr(0, suffix.size() - 5);
    const std::string path = join_path(
        input_directory, "pythia_truth_spectrum_tree_" + tag + ".root");
    if (gSystem->AccessPathName(path.c_str()) || tree.Add(path.c_str()) != 1)
    {
      std::cerr << "AccumulatePythiaTruthPtSpectra - cannot add input: " << path << std::endl;
      return 2;
    }
  }
  const Long64_t n_entries = tree.GetEntries();
  if (n_entries <= 0)
  {
    std::cerr << "AccumulatePythiaTruthPtSpectra - no event_tree entries" << std::endl;
    return 2;
  }
  tree.SetBranchStatus("*", false);

  UChar_t event_weight_valid = 0U;
  double event_weight = 1.0;
  UInt_t truth_photon_n = 0U;
  std::vector<unsigned char>* photon_kinematics_valid = nullptr;
  std::vector<float>* photon_pt = nullptr;
  std::vector<float>* photon_eta = nullptr;
  std::vector<unsigned char>* photon_classification_valid = nullptr;
  std::vector<int>* photon_category = nullptr;
  std::vector<unsigned char>* photon_copy_chain_valid = nullptr;
  std::vector<int>* photon_origin_parent_count = nullptr;
  std::vector<int>* photon_origin_parent_pdg = nullptr;
  UInt_t truth_pi0_n = 0U;
  std::vector<unsigned char>* pi0_kinematics_valid = nullptr;
  std::vector<float>* pi0_pt = nullptr;
  std::vector<float>* pi0_eta = nullptr;
  UInt_t truth_pi0_decay_photon_n = 0U;
  std::vector<unsigned char>* pi0_decay_photon_kinematics_valid = nullptr;
  std::vector<float>* pi0_decay_photon_pt = nullptr;
  std::vector<float>* pi0_decay_photon_eta = nullptr;

  bool ok = true;
  if (use_event_weight)
  {
    ok &= bind_branch(&tree, "event_weight_valid", &event_weight_valid);
    ok &= bind_branch(&tree, "event_weight", &event_weight);
  }
  ok &= bind_branch(&tree, "truth_photon_n", &truth_photon_n);
  ok &= bind_branch(&tree, "truth_photon_kinematics_valid", &photon_kinematics_valid);
  ok &= bind_branch(&tree, "truth_photon_pt", &photon_pt);
  ok &= bind_branch(&tree, "truth_photon_eta", &photon_eta);
  ok &= bind_branch(&tree, "truth_photon_classification_valid", &photon_classification_valid);
  ok &= bind_branch(&tree, "truth_photon_category", &photon_category);
  ok &= bind_branch(&tree, "truth_photon_copy_chain_valid", &photon_copy_chain_valid);
  ok &= bind_branch(&tree, "truth_photon_origin_parent_count", &photon_origin_parent_count);
  ok &= bind_branch(&tree, "truth_photon_origin_parent_pdg", &photon_origin_parent_pdg);
  ok &= bind_branch(&tree, "truth_pi0_n", &truth_pi0_n);
  ok &= bind_branch(&tree, "truth_pi0_kinematics_valid", &pi0_kinematics_valid);
  ok &= bind_branch(&tree, "truth_pi0_pt", &pi0_pt);
  ok &= bind_branch(&tree, "truth_pi0_eta", &pi0_eta);
  ok &= bind_branch(&tree, "truth_pi0_decay_photon_n", &truth_pi0_decay_photon_n);
  ok &= bind_branch(&tree, "truth_pi0_decay_photon_kinematics_valid", &pi0_decay_photon_kinematics_valid);
  ok &= bind_branch(&tree, "truth_pi0_decay_photon_pt", &pi0_decay_photon_pt);
  ok &= bind_branch(&tree, "truth_pi0_decay_photon_eta", &pi0_decay_photon_eta);
  if (!ok)
  {
    return 3;
  }

  TH1::AddDirectory(false);
  TH1D prompt_photon("h_prompt_photon_truth_pt_raw", "", n_bins, 0.0, pt_max);
  TH1D pi0("h_pi0_truth_pt_raw", "", n_bins, 0.0, pt_max);
  TH1D pi0_decay_photon("h_pi0_decay_photon_truth_pt_raw", "", n_bins, 0.0, pt_max);
  TH1D hepmc_pi0_decay_photon(
      "h_hepmc_pi0_decay_photon_truth_pt_raw", "", n_bins, 0.0, pt_max);
  TH1D g4_pi0_decay_photon(
      "h_g4_pi0_decay_photon_truth_pt_raw", "", n_bins, 0.0, pt_max);
  for (TH1D* histogram : {&prompt_photon, &pi0, &pi0_decay_photon,
                           &hepmc_pi0_decay_photon, &g4_pi0_decay_photon})
  {
    histogram->Sumw2();
  }

  ULong64_t n_prompt_photon = 0ULL;
  ULong64_t n_pi0 = 0ULL;
  ULong64_t n_pi0_decay_photon = 0ULL;
  ULong64_t n_hepmc_pi0_decay_photon = 0ULL;
  ULong64_t n_g4_pi0_decay_photon = 0ULL;
  ULong64_t malformed_events = 0ULL;
  ULong64_t invalid_weight_events = 0ULL;
  for (Long64_t entry = 0; entry < n_entries; ++entry)
  {
    if (tree.GetEntry(entry) <= 0)
    {
      ++malformed_events;
      continue;
    }
    const std::size_t n_photon = truth_photon_n;
    const std::size_t n_pi0_event = truth_pi0_n;
    const std::size_t n_pi0_decay_photon_event = truth_pi0_decay_photon_n;
    const bool photon_aligned =
        aligned(photon_kinematics_valid, n_photon) && aligned(photon_pt, n_photon) &&
        aligned(photon_eta, n_photon) && aligned(photon_classification_valid, n_photon) &&
        aligned(photon_category, n_photon) && aligned(photon_copy_chain_valid, n_photon) &&
        aligned(photon_origin_parent_count, n_photon) &&
        aligned(photon_origin_parent_pdg, n_photon);
    const bool pi0_aligned = aligned(pi0_kinematics_valid, n_pi0_event) &&
        aligned(pi0_pt, n_pi0_event) && aligned(pi0_eta, n_pi0_event);
    const bool pi0_decay_photon_aligned =
        aligned(pi0_decay_photon_kinematics_valid, n_pi0_decay_photon_event) &&
        aligned(pi0_decay_photon_pt, n_pi0_decay_photon_event) &&
        aligned(pi0_decay_photon_eta, n_pi0_decay_photon_event);
    if (!photon_aligned || !pi0_aligned || !pi0_decay_photon_aligned)
    {
      ++malformed_events;
      continue;
    }
    if (use_event_weight && (!event_weight_valid || !std::isfinite(event_weight)))
    {
      ++invalid_weight_events;
      continue;
    }
    const double weight = use_event_weight ? event_weight : 1.0;

    for (std::size_t particle = 0; particle < n_photon; ++particle)
    {
      if (!(*photon_kinematics_valid)[particle] ||
          !in_acceptance((*photon_eta)[particle], max_abs_eta))
      {
        continue;
      }
      if ((*photon_classification_valid)[particle] &&
          ((*photon_category)[particle] == 1 || (*photon_category)[particle] == 2))
      {
        prompt_photon.Fill((*photon_pt)[particle], weight);
        ++n_prompt_photon;
      }
      if ((*photon_copy_chain_valid)[particle] &&
          (*photon_origin_parent_count)[particle] == 1 &&
          (*photon_origin_parent_pdg)[particle] == 111)
      {
        hepmc_pi0_decay_photon.Fill((*photon_pt)[particle], weight);
        pi0_decay_photon.Fill((*photon_pt)[particle], weight);
        ++n_hepmc_pi0_decay_photon;
        ++n_pi0_decay_photon;
      }
    }
    for (std::size_t particle = 0; particle < n_pi0_event; ++particle)
    {
      if ((*pi0_kinematics_valid)[particle] &&
          in_acceptance((*pi0_eta)[particle], max_abs_eta))
      {
        pi0.Fill((*pi0_pt)[particle], weight);
        ++n_pi0;
      }
    }
    for (std::size_t particle = 0; particle < n_pi0_decay_photon_event; ++particle)
    {
      if ((*pi0_decay_photon_kinematics_valid)[particle] &&
          in_acceptance((*pi0_decay_photon_eta)[particle], max_abs_eta))
      {
        g4_pi0_decay_photon.Fill((*pi0_decay_photon_pt)[particle], weight);
        pi0_decay_photon.Fill((*pi0_decay_photon_pt)[particle], weight);
        ++n_g4_pi0_decay_photon;
        ++n_pi0_decay_photon;
      }
    }
    if ((entry + 1) % 100000 == 0 || entry + 1 == n_entries)
    {
      std::cout << "AccumulatePythiaTruthPtSpectra - processed " << entry + 1
                << "/" << n_entries << " events" << std::endl;
    }
  }
  if (malformed_events != 0ULL || invalid_weight_events != 0ULL)
  {
    std::cerr << "AccumulatePythiaTruthPtSpectra - malformed/invalid-weight events = "
              << malformed_events << "/" << invalid_weight_events << std::endl;
    return 4;
  }

  TFile output(output_file.c_str(), "RECREATE");
  if (output.IsZombie())
  {
    std::cerr << "AccumulatePythiaTruthPtSpectra - cannot create output: "
              << output_file << std::endl;
    return 5;
  }
  prompt_photon.Write();
  pi0.Write();
  pi0_decay_photon.Write();
  hepmc_pi0_decay_photon.Write();
  g4_pi0_decay_photon.Write();

  Int_t schema_version = 2;
  std::string photon_selection = "prompt_category_1_or_2";
  std::string pi0_decay_photon_selection =
      "hepmc_final_photon_with_valid_single_pi0_origin_plus_"
      "g4_immediate_photon_daughter_of_signal_primary_pi0";
  Long64_t files_added = static_cast<Long64_t>(suffixes.size());
  std::string stored_manifest_path = manifest_path;
  std::string stored_input_directory = input_directory;
  Long64_t stored_manifest_begin = manifest_begin;
  Long64_t stored_manifest_end = manifest_end;
  Long64_t events_processed = n_entries;
  Int_t stored_n_bins = n_bins;
  double stored_pt_max = pt_max;
  double stored_max_abs_eta = max_abs_eta;
  UChar_t weighted = use_event_weight ? 1U : 0U;
  UChar_t bin_width_normalized = 0U;
  std::string first_suffix = suffixes.front();
  std::string last_suffix = suffixes.back();
  TTree metadata("metadata", "Pythia truth pT partial metadata");
  metadata.Branch("schema_version", &schema_version);
  metadata.Branch("photon_selection", &photon_selection);
  metadata.Branch("pi0_decay_photon_selection", &pi0_decay_photon_selection);
  metadata.Branch("manifest_path", &stored_manifest_path);
  metadata.Branch("input_directory", &stored_input_directory);
  metadata.Branch("manifest_begin", &stored_manifest_begin);
  metadata.Branch("manifest_end", &stored_manifest_end);
  metadata.Branch("files_added", &files_added);
  metadata.Branch("first_suffix", &first_suffix);
  metadata.Branch("last_suffix", &last_suffix);
  metadata.Branch("events_processed", &events_processed);
  metadata.Branch("n_bins", &stored_n_bins);
  metadata.Branch("pt_max", &stored_pt_max);
  metadata.Branch("max_abs_eta", &stored_max_abs_eta);
  metadata.Branch("use_event_weight", &weighted);
  metadata.Branch("bin_width_normalized", &bin_width_normalized);
  metadata.Branch("prompt_photon_count", &n_prompt_photon);
  metadata.Branch("pi0_count", &n_pi0);
  metadata.Branch("pi0_decay_photon_count", &n_pi0_decay_photon);
  metadata.Branch("hepmc_pi0_decay_photon_count", &n_hepmc_pi0_decay_photon);
  metadata.Branch("g4_pi0_decay_photon_count", &n_g4_pi0_decay_photon);
  metadata.Branch("malformed_event_count", &malformed_events);
  metadata.Branch("invalid_weight_event_count", &invalid_weight_events);
  metadata.Fill();
  metadata.Write();
  output.Close();
  if (output.TestBit(TFile::kWriteError))
  {
    std::cerr << "AccumulatePythiaTruthPtSpectra - write error: " << output_file << std::endl;
    return 5;
  }

  std::cout << "AccumulatePythiaTruthPtSpectra - files/events/prompt/pi0/"
               "pi0 decay (total/HepMC/G4) = "
            << files_added << "/" << n_entries << "/" << n_prompt_photon << "/"
            << n_pi0 << "/" << n_pi0_decay_photon << "/"
            << n_hepmc_pi0_decay_photon << "/" << n_g4_pi0_decay_photon << std::endl;
  return 0;
}
