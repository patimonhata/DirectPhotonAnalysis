#include <TFile.h>
#include <TH1D.h>
#include <TKey.h>
#include <TTree.h>

#include <cmath>
#include <iostream>
#include <set>
#include <algorithm>
#include <string>

namespace
{
template <class T>
bool bind_branch(TTree* tree, const char* name, T* address)
{
  if (!tree->GetBranch(name))
  {
    std::cerr << "check_pythia_truth_pt_partial - missing metadata branch: "
              << name << std::endl;
    return false;
  }
  return tree->SetBranchAddress(name, address) >= 0;
}

bool valid_histogram(
    const TH1D* histogram,
    const int n_bins,
    const double pt_max,
    const ULong64_t expected_entries)
{
  if (!histogram || histogram->GetNbinsX() != n_bins ||
      std::abs(histogram->GetXaxis()->GetXmin()) > 1e-12 ||
      std::abs(histogram->GetXaxis()->GetXmax() - pt_max) > 1e-12 ||
      histogram->GetSumw2N() == 0 ||
      std::abs(histogram->GetEntries() - static_cast<double>(expected_entries)) > 0.5)
  {
    return false;
  }
  for (int bin = 0; bin <= n_bins + 1; ++bin)
  {
    if (!std::isfinite(histogram->GetBinContent(bin)) ||
        !std::isfinite(histogram->GetBinError(bin)))
    {
      return false;
    }
  }
  return true;
}

bool histogram_sum_matches(
    const TH1D* total,
    const TH1D* hepmc,
    const TH1D* g4)
{
  if (!total || !hepmc || !g4 ||
      std::abs(total->GetEntries() - hepmc->GetEntries() - g4->GetEntries()) > 0.5)
  {
    return false;
  }
  for (int bin = 0; bin <= total->GetNbinsX() + 1; ++bin)
  {
    const double expected_content =
        hepmc->GetBinContent(bin) + g4->GetBinContent(bin);
    const double expected_error2 =
        std::pow(hepmc->GetBinError(bin), 2) + std::pow(g4->GetBinError(bin), 2);
    const double content_scale = std::max(1.0, std::abs(expected_content));
    const double error2_scale = std::max(1.0, std::abs(expected_error2));
    if (std::abs(total->GetBinContent(bin) - expected_content) >
            1e-10 * content_scale ||
        std::abs(std::pow(total->GetBinError(bin), 2) - expected_error2) >
            1e-10 * error2_scale)
    {
      return false;
    }
  }
  return true;
}
}

int check_pythia_truth_pt_partial(const std::string input_file)
{
  TFile input(input_file.c_str(), "READ");
  if (input.IsZombie())
  {
    std::cerr << "check_pythia_truth_pt_partial - cannot open " << input_file << std::endl;
    return 1;
  }

  const std::set<std::string> expected_keys = {
      "h_prompt_photon_truth_pt_raw",
      "h_pi0_truth_pt_raw",
      "h_pi0_decay_photon_truth_pt_raw",
      "h_hepmc_pi0_decay_photon_truth_pt_raw",
      "h_g4_pi0_decay_photon_truth_pt_raw",
      "metadata"};
  std::set<std::string> observed_keys;
  TIter next_key(input.GetListOfKeys());
  while (TKey* key = dynamic_cast<TKey*>(next_key()))
  {
    observed_keys.insert(key->GetName());
  }
  if (observed_keys != expected_keys)
  {
    std::cerr << "check_pythia_truth_pt_partial - unexpected top-level keys" << std::endl;
    return 2;
  }

  TTree* metadata = nullptr;
  TH1D* prompt_photon = nullptr;
  TH1D* pi0 = nullptr;
  TH1D* pi0_decay_photon = nullptr;
  TH1D* hepmc_pi0_decay_photon = nullptr;
  TH1D* g4_pi0_decay_photon = nullptr;
  input.GetObject("metadata", metadata);
  input.GetObject("h_prompt_photon_truth_pt_raw", prompt_photon);
  input.GetObject("h_pi0_truth_pt_raw", pi0);
  input.GetObject("h_pi0_decay_photon_truth_pt_raw", pi0_decay_photon);
  input.GetObject("h_hepmc_pi0_decay_photon_truth_pt_raw", hepmc_pi0_decay_photon);
  input.GetObject("h_g4_pi0_decay_photon_truth_pt_raw", g4_pi0_decay_photon);
  if (!metadata || metadata->GetEntries() != 1)
  {
    std::cerr << "check_pythia_truth_pt_partial - invalid metadata tree" << std::endl;
    return 2;
  }

  Int_t schema_version = 0;
  std::string* photon_selection = nullptr;
  std::string* pi0_decay_photon_selection = nullptr;
  std::string* manifest_path = nullptr;
  std::string* input_directory = nullptr;
  Long64_t manifest_begin = -1;
  Long64_t manifest_end = -1;
  Long64_t files_added = 0;
  std::string* first_suffix = nullptr;
  std::string* last_suffix = nullptr;
  Long64_t events_processed = 0;
  Int_t n_bins = 0;
  double pt_max = 0.0;
  double max_abs_eta = 0.0;
  UChar_t use_event_weight = 0U;
  UChar_t bin_width_normalized = 1U;
  ULong64_t prompt_photon_count = 0ULL;
  ULong64_t pi0_count = 0ULL;
  ULong64_t pi0_decay_photon_count = 0ULL;
  ULong64_t hepmc_pi0_decay_photon_count = 0ULL;
  ULong64_t g4_pi0_decay_photon_count = 0ULL;
  ULong64_t malformed_event_count = 0ULL;
  ULong64_t invalid_weight_event_count = 0ULL;

  bool ok = true;
  ok &= bind_branch(metadata, "schema_version", &schema_version);
  ok &= bind_branch(metadata, "photon_selection", &photon_selection);
  ok &= bind_branch(metadata, "pi0_decay_photon_selection", &pi0_decay_photon_selection);
  ok &= bind_branch(metadata, "manifest_path", &manifest_path);
  ok &= bind_branch(metadata, "input_directory", &input_directory);
  ok &= bind_branch(metadata, "manifest_begin", &manifest_begin);
  ok &= bind_branch(metadata, "manifest_end", &manifest_end);
  ok &= bind_branch(metadata, "files_added", &files_added);
  ok &= bind_branch(metadata, "first_suffix", &first_suffix);
  ok &= bind_branch(metadata, "last_suffix", &last_suffix);
  ok &= bind_branch(metadata, "events_processed", &events_processed);
  ok &= bind_branch(metadata, "n_bins", &n_bins);
  ok &= bind_branch(metadata, "pt_max", &pt_max);
  ok &= bind_branch(metadata, "max_abs_eta", &max_abs_eta);
  ok &= bind_branch(metadata, "use_event_weight", &use_event_weight);
  ok &= bind_branch(metadata, "bin_width_normalized", &bin_width_normalized);
  ok &= bind_branch(metadata, "prompt_photon_count", &prompt_photon_count);
  ok &= bind_branch(metadata, "pi0_count", &pi0_count);
  ok &= bind_branch(metadata, "pi0_decay_photon_count", &pi0_decay_photon_count);
  ok &= bind_branch(metadata, "hepmc_pi0_decay_photon_count", &hepmc_pi0_decay_photon_count);
  ok &= bind_branch(metadata, "g4_pi0_decay_photon_count", &g4_pi0_decay_photon_count);
  ok &= bind_branch(metadata, "malformed_event_count", &malformed_event_count);
  ok &= bind_branch(metadata, "invalid_weight_event_count", &invalid_weight_event_count);
  if (!ok || metadata->GetEntry(0) <= 0)
  {
    return 3;
  }

  const bool valid_metadata = schema_version == 2 && photon_selection &&
      *photon_selection == "prompt_category_1_or_2" &&
      pi0_decay_photon_selection &&
      *pi0_decay_photon_selection ==
          "hepmc_final_photon_with_valid_single_pi0_origin_plus_"
          "g4_immediate_photon_daughter_of_signal_primary_pi0" &&
      manifest_path &&
      !manifest_path->empty() && input_directory && !input_directory->empty() &&
      manifest_begin >= 0 && manifest_end > manifest_begin &&
      files_added == manifest_end - manifest_begin && first_suffix &&
      !first_suffix->empty() && last_suffix && !last_suffix->empty() &&
      events_processed > 0 && n_bins > 0 && std::isfinite(pt_max) && pt_max > 0.0 &&
      std::isfinite(max_abs_eta) && use_event_weight <= 1U &&
      bin_width_normalized == 0U &&
      pi0_decay_photon_count ==
          hepmc_pi0_decay_photon_count + g4_pi0_decay_photon_count &&
      malformed_event_count == 0ULL && invalid_weight_event_count == 0ULL;
  if (!valid_metadata)
  {
    std::cerr << "check_pythia_truth_pt_partial - invalid metadata values" << std::endl;
    return 4;
  }
  if (!valid_histogram(prompt_photon, n_bins, pt_max, prompt_photon_count) ||
      !valid_histogram(pi0, n_bins, pt_max, pi0_count) ||
      !valid_histogram(pi0_decay_photon, n_bins, pt_max, pi0_decay_photon_count) ||
      !valid_histogram(hepmc_pi0_decay_photon, n_bins, pt_max,
          hepmc_pi0_decay_photon_count) ||
      !valid_histogram(g4_pi0_decay_photon, n_bins, pt_max,
          g4_pi0_decay_photon_count) ||
      !histogram_sum_matches(
          pi0_decay_photon, hepmc_pi0_decay_photon, g4_pi0_decay_photon))
  {
    std::cerr << "check_pythia_truth_pt_partial - invalid histogram" << std::endl;
    return 5;
  }

  std::cout << "check_pythia_truth_pt_partial - range/files/events/prompt/pi0/"
               "pi0 decay (total/HepMC/G4) = ["
            << manifest_begin << ":" << manifest_end << "]/" << files_added << "/"
            << events_processed << "/" << prompt_photon_count << "/" << pi0_count
            << "/" << pi0_decay_photon_count << "/"
            << hepmc_pi0_decay_photon_count << "/" << g4_pi0_decay_photon_count
            << std::endl;
  return 0;
}
