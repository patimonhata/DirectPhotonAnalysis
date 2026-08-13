#include <TFile.h>
#include <TH1D.h>
#include <TKey.h>
#include <TTree.h>

#include <cmath>
#include <iostream>
#include <set>
#include <string>

namespace
{
template <class T>
bool bind(TTree* tree, const char* name, T* address)
{
  if (!tree->GetBranch(name))
  {
    std::cerr << "check_pythia_cluster_et_partial - missing branch: "
              << name << std::endl;
    return false;
  }
  return tree->SetBranchAddress(name, address) >= 0;
}

bool valid_histogram(const TH1D* histogram, int n_bins, double et_max,
                     unsigned long long expected_entries)
{
  if (!histogram || histogram->GetNbinsX() != n_bins ||
      std::abs(histogram->GetXaxis()->GetXmin()) > 1e-12 ||
      std::abs(histogram->GetXaxis()->GetXmax() - et_max) > 1e-12 ||
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
}

int check_pythia_cluster_et_partial(const std::string input_file)
{
  TFile input(input_file.c_str(), "READ");
  if (input.IsZombie())
  {
    return 1;
  }
  const std::set<std::string> expected_keys = {
      "h_prompt_cluster_et_raw", "h_pi0_cluster_et_raw",
      "h_pi0_separated_cluster_et_raw", "h_pi0_merged_cluster_et_raw",
      "h_pi0_missing_cluster_et_raw", "metadata"};
  std::set<std::string> observed_keys;
  TIter next(input.GetListOfKeys());
  while (TKey* key = dynamic_cast<TKey*>(next()))
  {
    observed_keys.insert(key->GetName());
  }
  if (observed_keys != expected_keys)
  {
    std::cerr << "check_pythia_cluster_et_partial - unexpected top-level keys"
              << std::endl;
    return 2;
  }

  TTree* metadata = nullptr;
  TH1D* prompt = nullptr;
  TH1D* pi0 = nullptr;
  TH1D* separated = nullptr;
  TH1D* merged = nullptr;
  TH1D* missing = nullptr;
  input.GetObject("metadata", metadata);
  input.GetObject("h_prompt_cluster_et_raw", prompt);
  input.GetObject("h_pi0_cluster_et_raw", pi0);
  input.GetObject("h_pi0_separated_cluster_et_raw", separated);
  input.GetObject("h_pi0_merged_cluster_et_raw", merged);
  input.GetObject("h_pi0_missing_cluster_et_raw", missing);
  if (!metadata || metadata->GetEntries() != 1)
  {
    return 2;
  }

  int schema_version = 0;
  std::string* manifest_path = nullptr;
  std::string* first_suffix = nullptr;
  std::string* last_suffix = nullptr;
  std::string* cluster_collection = nullptr;
  std::string* prompt_selection = nullptr;
  std::string* pi0_selection = nullptr;
  std::string* topology_priority = nullptr;
  std::string* projection_scheme = nullptr;
  long long manifest_begin = -1;
  long long manifest_end = -1;
  int signal_embedding_id = 0;
  int n_bins = 0;
  double et_max = 0.0;
  double truth_eta_max = 0.0;
  double cluster_eta_max = 0.0;
  double min_cluster_energy = 0.0;
  double dominant_fraction_min = 0.0;
  double pi0_contributor_fraction_min = 0.0;
  double separated_delta_r_cut = 0.0;
  double merged_delta_r_cut = 0.0;
  double response_min = 0.0;
  double response_max = 0.0;
  unsigned char bin_width_normalized = 1U;
  unsigned long long events_processed = 0;
  unsigned long long events_written = 0;
  unsigned long long events_invalid = 0;
  unsigned long long prompt_count = 0;
  unsigned long long pi0_count = 0;
  unsigned long long pi0_g4_count = 0;
  unsigned long long pi0_generator_count = 0;
  unsigned long long candidate_g4_count = 0;
  unsigned long long candidate_generator_count = 0;
  unsigned long long projection_failure_count = 0;
  unsigned long long separated_count = 0;
  unsigned long long merged_count = 0;
  unsigned long long missing_count = 0;
  unsigned long long none_count = 0;
  unsigned long long ambiguous_count = 0;
  unsigned long long separated_fill_count = 0;
  unsigned long long merged_fill_count = 0;
  unsigned long long missing_fill_count = 0;

  bool ok = true;
  ok &= bind(metadata, "schema_version", &schema_version);
  ok &= bind(metadata, "manifest_path", &manifest_path);
  ok &= bind(metadata, "manifest_begin", &manifest_begin);
  ok &= bind(metadata, "manifest_end", &manifest_end);
  ok &= bind(metadata, "first_suffix", &first_suffix);
  ok &= bind(metadata, "last_suffix", &last_suffix);
  ok &= bind(metadata, "cluster_collection", &cluster_collection);
  ok &= bind(metadata, "prompt_selection", &prompt_selection);
  ok &= bind(metadata, "pi0_selection", &pi0_selection);
  ok &= bind(metadata, "topology_priority", &topology_priority);
  ok &= bind(metadata, "projection_scheme", &projection_scheme);
  ok &= bind(metadata, "signal_embedding_id", &signal_embedding_id);
  ok &= bind(metadata, "n_bins", &n_bins);
  ok &= bind(metadata, "et_max", &et_max);
  ok &= bind(metadata, "truth_eta_max", &truth_eta_max);
  ok &= bind(metadata, "cluster_eta_max", &cluster_eta_max);
  ok &= bind(metadata, "min_cluster_energy", &min_cluster_energy);
  ok &= bind(metadata, "dominant_fraction_min", &dominant_fraction_min);
  ok &= bind(metadata, "pi0_contributor_fraction_min", &pi0_contributor_fraction_min);
  ok &= bind(metadata, "separated_delta_r_cut", &separated_delta_r_cut);
  ok &= bind(metadata, "merged_delta_r_cut", &merged_delta_r_cut);
  ok &= bind(metadata, "response_min", &response_min);
  ok &= bind(metadata, "response_max", &response_max);
  ok &= bind(metadata, "bin_width_normalized", &bin_width_normalized);
  ok &= bind(metadata, "events_processed", &events_processed);
  ok &= bind(metadata, "events_written", &events_written);
  ok &= bind(metadata, "events_invalid", &events_invalid);
  ok &= bind(metadata, "prompt_cluster_count", &prompt_count);
  ok &= bind(metadata, "pi0_cluster_count", &pi0_count);
  ok &= bind(metadata, "pi0_cluster_g4_decay_count", &pi0_g4_count);
  ok &= bind(metadata, "pi0_cluster_generator_decay_count", &pi0_generator_count);
  ok &= bind(metadata, "pi0_candidate_g4_decay_count", &candidate_g4_count);
  ok &= bind(metadata, "pi0_candidate_generator_decay_count", &candidate_generator_count);
  ok &= bind(metadata, "pi0_projection_failure_count", &projection_failure_count);
  ok &= bind(metadata, "pi0_separated_count", &separated_count);
  ok &= bind(metadata, "pi0_merged_count", &merged_count);
  ok &= bind(metadata, "pi0_missing_count", &missing_count);
  ok &= bind(metadata, "pi0_none_count", &none_count);
  ok &= bind(metadata, "pi0_ambiguous_count", &ambiguous_count);
  ok &= bind(metadata, "pi0_separated_cluster_fill_count", &separated_fill_count);
  ok &= bind(metadata, "pi0_merged_cluster_fill_count", &merged_fill_count);
  ok &= bind(metadata, "pi0_missing_cluster_fill_count", &missing_fill_count);
  if (!ok || metadata->GetEntry(0) <= 0)
  {
    return 3;
  }

  const bool valid_metadata = schema_version == 1 && manifest_path &&
      !manifest_path->empty() && manifest_begin >= 0 && manifest_end > manifest_begin &&
      first_suffix && !first_suffix->empty() && last_suffix && !last_suffix->empty() &&
      cluster_collection && *cluster_collection == "split" && prompt_selection &&
      *prompt_selection == "dominant_prompt_category_1_or_2" && pi0_selection &&
      topology_priority && *topology_priority ==
          "separated_then_merged_then_missing_then_none" &&
      projection_scheme && *projection_scheme ==
          "g4_photon_vertex_and_momentum_to_cemc_cylinder" &&
      signal_embedding_id == 1 && n_bins > 0 && et_max > 0.0 &&
      truth_eta_max > 0.0 && cluster_eta_max > 0.0 && min_cluster_energy >= 0.0 &&
      dominant_fraction_min >= 0.0 && dominant_fraction_min <= 1.0 &&
      pi0_contributor_fraction_min >= 0.0 && pi0_contributor_fraction_min <= 1.0 &&
      separated_delta_r_cut > 0.0 && merged_delta_r_cut >= separated_delta_r_cut &&
      response_min >= 0.0 && response_max > response_min &&
      bin_width_normalized == 0U && events_processed > 0 && events_invalid == 0 &&
      events_written == events_processed && pi0_count == pi0_g4_count + pi0_generator_count &&
      separated_fill_count == 2ULL * separated_count &&
      merged_fill_count == merged_count && missing_fill_count == missing_count &&
      candidate_g4_count + candidate_generator_count == separated_count + merged_count +
          missing_count + none_count + ambiguous_count + projection_failure_count;
  if (!valid_metadata)
  {
    std::cerr << "check_pythia_cluster_et_partial - invalid metadata values"
              << std::endl;
    return 4;
  }
  if (!valid_histogram(prompt, n_bins, et_max, prompt_count) ||
      !valid_histogram(pi0, n_bins, et_max, pi0_count) ||
      !valid_histogram(separated, n_bins, et_max, separated_fill_count) ||
      !valid_histogram(merged, n_bins, et_max, merged_fill_count) ||
      !valid_histogram(missing, n_bins, et_max, missing_fill_count))
  {
    std::cerr << "check_pythia_cluster_et_partial - invalid histogram" << std::endl;
    return 5;
  }
  std::cout << "check_pythia_cluster_et_partial - range/events/prompt/pi0/sep/merged/missing = ["
            << manifest_begin << ":" << manifest_end << "]/" << events_processed
            << "/" << prompt_count << "/" << pi0_count << "/" << separated_fill_count
            << "/" << merged_fill_count << "/" << missing_fill_count << std::endl;
  return 0;
}
