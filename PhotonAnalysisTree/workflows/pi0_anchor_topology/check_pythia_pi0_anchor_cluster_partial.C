#include <TFile.h>
#include <TH1D.h>
#include <TKey.h>
#include <TTree.h>

#include <array>
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
    std::cerr
        << "check_pythia_pi0_anchor_cluster_partial - missing branch: "
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
      std::abs(histogram->GetEntries() -
               static_cast<double>(expected_entries)) > 0.5)
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

int check_pythia_pi0_anchor_cluster_partial(const std::string input_file)
{
  TFile input(input_file.c_str(), "READ");
  if (input.IsZombie())
  {
    return 1;
  }

  const std::set<std::string> expected_keys = {
      "h_prompt_cluster_et_raw", "h_pi0_anchor_cluster_et_raw",
      "h_pi0_anchor_separated_cluster_et_raw",
      "h_pi0_anchor_merged_cluster_et_raw",
      "h_pi0_anchor_missing_cluster_et_raw",
      "h_pi0_anchor_other_cluster_et_raw", "metadata"};
  std::set<std::string> observed_keys;
  TIter next(input.GetListOfKeys());
  while (TKey* key = dynamic_cast<TKey*>(next()))
  {
    observed_keys.insert(key->GetName());
  }
  if (observed_keys != expected_keys)
  {
    std::cerr
        << "check_pythia_pi0_anchor_cluster_partial - unexpected top-level keys"
        << std::endl;
    return 2;
  }

  TTree* metadata = nullptr;
  input.GetObject("metadata", metadata);
  if (!metadata || metadata->GetEntries() != 1)
  {
    return 2;
  }

  int schema_version = 0;
  int signal_embedding_id = 0;
  int n_bins = 0;
  int matcher_version = 0;
  long long manifest_begin = -1;
  long long manifest_end = -1;
  double et_max = 0.0;
  double truth_eta_max = 0.0;
  double anchor_cluster_eta_max = 0.0;
  double partner_cluster_eta_max = 0.0;
  double min_cluster_energy = 0.0;
  double dominant_fraction_min = 0.0;
  double anchor_pi0_fraction_min = 0.0;
  double min_energy_contribution_fraction = 0.0;
  double min_photon_energy_recovery = 0.0;
  unsigned char bin_width_normalized = 1U;
  std::string* manifest_path = nullptr;
  std::string* first_suffix = nullptr;
  std::string* last_suffix = nullptr;
  std::string* cluster_collection = nullptr;
  std::string* classification_unit = nullptr;
  std::string* pi0_selection = nullptr;
  std::string* partner_selection = nullptr;
  std::string* topology_definition = nullptr;
  std::string* topology_priority = nullptr;
  std::string* response_policy = nullptr;
  std::string* photon_recovery_policy = nullptr;
  unsigned long long events_processed = 0;
  unsigned long long events_written = 0;
  unsigned long long events_invalid = 0;
  unsigned long long cluster_considered = 0;
  unsigned long long cluster_invalid_truth = 0;
  unsigned long long prompt_count = 0;
  unsigned long long candidate_g4 = 0;
  unsigned long long candidate_generator = 0;
  unsigned long long malformed_daughters = 0;
  unsigned long long anchor_count = 0;
  unsigned long long anchor_g4 = 0;
  unsigned long long anchor_generator = 0;
  unsigned long long ambiguous_main = 0;
  unsigned long long energy_match_invalid = 0;
  unsigned long long separated_count = 0;
  unsigned long long merged_count = 0;
  unsigned long long missing_count = 0;
  unsigned long long other_count = 0;

  bool ok = true;
  ok &= bind(metadata, "schema_version", &schema_version);
  ok &= bind(metadata, "manifest_path", &manifest_path);
  ok &= bind(metadata, "manifest_begin", &manifest_begin);
  ok &= bind(metadata, "manifest_end", &manifest_end);
  ok &= bind(metadata, "first_suffix", &first_suffix);
  ok &= bind(metadata, "last_suffix", &last_suffix);
  ok &= bind(metadata, "cluster_collection", &cluster_collection);
  ok &= bind(metadata, "classification_unit", &classification_unit);
  ok &= bind(metadata, "pi0_selection", &pi0_selection);
  ok &= bind(metadata, "partner_selection", &partner_selection);
  ok &= bind(metadata, "topology_definition", &topology_definition);
  ok &= bind(metadata, "topology_priority", &topology_priority);
  ok &= bind(metadata, "response_policy", &response_policy);
  ok &= bind(metadata, "photon_recovery_policy",
             &photon_recovery_policy);
  ok &= bind(metadata, "signal_embedding_id", &signal_embedding_id);
  ok &= bind(metadata, "n_bins", &n_bins);
  ok &= bind(metadata, "et_max", &et_max);
  ok &= bind(metadata, "truth_eta_max", &truth_eta_max);
  ok &= bind(metadata, "anchor_cluster_eta_max", &anchor_cluster_eta_max);
  ok &= bind(metadata, "partner_cluster_eta_max", &partner_cluster_eta_max);
  ok &= bind(metadata, "min_cluster_energy", &min_cluster_energy);
  ok &= bind(metadata, "dominant_fraction_min", &dominant_fraction_min);
  ok &= bind(metadata, "anchor_pi0_fraction_min",
             &anchor_pi0_fraction_min);
  ok &= bind(metadata, "min_energy_contribution_fraction",
             &min_energy_contribution_fraction);
  ok &= bind(metadata, "min_photon_energy_recovery",
             &min_photon_energy_recovery);
  ok &= bind(metadata, "pi0_truth_matching_algorithm_version",
             &matcher_version);
  ok &= bind(metadata, "bin_width_normalized", &bin_width_normalized);
  ok &= bind(metadata, "events_processed", &events_processed);
  ok &= bind(metadata, "events_written", &events_written);
  ok &= bind(metadata, "events_invalid", &events_invalid);
  ok &= bind(metadata, "cluster_considered_count", &cluster_considered);
  ok &= bind(metadata, "cluster_invalid_truth_count",
             &cluster_invalid_truth);
  ok &= bind(metadata, "prompt_cluster_count", &prompt_count);
  ok &= bind(metadata, "pi0_candidate_g4_decay_count", &candidate_g4);
  ok &= bind(metadata, "pi0_candidate_generator_decay_count",
             &candidate_generator);
  ok &= bind(metadata, "pi0_malformed_daughters_count",
             &malformed_daughters);
  ok &= bind(metadata, "anchor_cluster_count", &anchor_count);
  ok &= bind(metadata, "anchor_g4_decay_count", &anchor_g4);
  ok &= bind(metadata, "anchor_generator_decay_count", &anchor_generator);
  ok &= bind(metadata, "anchor_ambiguous_main_count", &ambiguous_main);
  ok &= bind(metadata, "energy_match_invalid_count", &energy_match_invalid);
  ok &= bind(metadata, "separated_count", &separated_count);
  ok &= bind(metadata, "merged_count", &merged_count);
  ok &= bind(metadata, "missing_count", &missing_count);
  ok &= bind(metadata, "other_count", &other_count);
  if (!ok || metadata->GetEntry(0) <= 0)
  {
    return 3;
  }

  const bool valid_metadata =
      schema_version == 3 && manifest_path && !manifest_path->empty() &&
      manifest_begin >= 0 && manifest_end > manifest_begin &&
      first_suffix && !first_suffix->empty() &&
      last_suffix && !last_suffix->empty() &&
      cluster_collection && *cluster_collection == "split" &&
      classification_unit &&
      *classification_unit ==
          "every_cluster_with_selected_pi0_as_grouped_main_contributor" &&
      pi0_selection &&
      *pi0_selection ==
          "signal_g4_primary_pi0_or_generator_pi0_with_exactly_two_g4_photons" &&
      partner_selection &&
      *partner_selection ==
          "same_energy_cut_as_anchor_partner_eta_cut_configurable" &&
      topology_definition &&
      *topology_definition ==
          "anchor_membership_in_recovered_direct_daughter_maximum_deposit_clusters" &&
      topology_priority &&
      *topology_priority ==
          "ambiguous_main_to_other_then_merged_then_separated_then_missing_then_other" &&
      response_policy &&
      *response_policy == "not_used_for_classification" &&
      photon_recovery_policy &&
      *photon_recovery_policy ==
          "cluster_energy_times_gamma_deposit_fraction_over_truth_energy_threshold" &&
      signal_embedding_id > 0 && n_bins > 0 &&
      std::isfinite(et_max) && et_max > 0.0 &&
      std::isfinite(truth_eta_max) && truth_eta_max > 0.0 &&
      std::isfinite(anchor_cluster_eta_max) &&
      anchor_cluster_eta_max > 0.0 &&
      std::isfinite(partner_cluster_eta_max) &&
      std::isfinite(min_cluster_energy) && min_cluster_energy >= 0.0 &&
      dominant_fraction_min >= 0.0 && dominant_fraction_min <= 1.0 &&
      anchor_pi0_fraction_min >= 0.0 &&
      anchor_pi0_fraction_min <= 1.0 &&
      min_energy_contribution_fraction >= 0.0 &&
      min_energy_contribution_fraction < 1.0 &&
      std::isfinite(min_photon_energy_recovery) &&
      min_photon_energy_recovery >= 0.0 &&
      min_photon_energy_recovery <= 1.0 &&
      matcher_version > 0 && bin_width_normalized == 0U &&
      events_processed > 0 &&
      events_written + events_invalid == events_processed &&
      cluster_invalid_truth <= cluster_considered &&
      anchor_count == anchor_g4 + anchor_generator &&
      anchor_count ==
          separated_count + merged_count + missing_count + other_count &&
      ambiguous_main <= other_count;
  if (!valid_metadata)
  {
    std::cerr
        << "check_pythia_pi0_anchor_cluster_partial - invalid metadata values"
        << std::endl;
    return 4;
  }

  const std::array<const char*, 6> names = {
      "h_prompt_cluster_et_raw", "h_pi0_anchor_cluster_et_raw",
      "h_pi0_anchor_separated_cluster_et_raw",
      "h_pi0_anchor_merged_cluster_et_raw",
      "h_pi0_anchor_missing_cluster_et_raw",
      "h_pi0_anchor_other_cluster_et_raw"};
  const std::array<unsigned long long, 6> counts = {
      prompt_count, anchor_count, separated_count, merged_count,
      missing_count, other_count};
  std::array<TH1D*, 6> histograms{};
  for (std::size_t index = 0; index < names.size(); ++index)
  {
    input.GetObject(names[index], histograms[index]);
    if (!valid_histogram(
            histograms[index], n_bins, et_max, counts[index]))
    {
      std::cerr
          << "check_pythia_pi0_anchor_cluster_partial - invalid histogram: "
          << names[index] << std::endl;
      return 5;
    }
  }

  for (int bin = 0; bin <= n_bins + 1; ++bin)
  {
    const double categories =
        histograms[2]->GetBinContent(bin) +
        histograms[3]->GetBinContent(bin) +
        histograms[4]->GetBinContent(bin) +
        histograms[5]->GetBinContent(bin);
    if (std::abs(histograms[1]->GetBinContent(bin) - categories) > 1e-9)
    {
      std::cerr
          << "check_pythia_pi0_anchor_cluster_partial - bin partition failure"
          << std::endl;
      return 5;
    }
  }

  std::cout
      << "check_pythia_pi0_anchor_cluster_partial - range/events/anchor"
      << "/separated/merged/missing/other = ["
      << manifest_begin << ":" << manifest_end << "]/"
      << events_processed << "/" << anchor_count << "/"
      << separated_count << "/" << merged_count << "/"
      << missing_count << "/" << other_count << std::endl;
  return 0;
}
