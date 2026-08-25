#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
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
    std::cerr << "check_pythia_pi0_cluster_multiplicity_partial - missing branch: "
              << name << std::endl;
    return false;
  }
  return tree->SetBranchAddress(name, address) >= 0;
}

bool finite_histogram(const TH1* histogram)
{
  if (!histogram || histogram->GetSumw2N() == 0)
  {
    return false;
  }
  for (int x = 0; x <= histogram->GetNbinsX() + 1; ++x)
  {
    if (histogram->GetDimension() == 1)
    {
      if (!std::isfinite(histogram->GetBinContent(x)) ||
          !std::isfinite(histogram->GetBinError(x)))
      {
        return false;
      }
      continue;
    }
    for (int y = 0; y <= histogram->GetNbinsY() + 1; ++y)
    {
      const int bin = histogram->GetBin(x, y);
      if (!std::isfinite(histogram->GetBinContent(bin)) ||
          !std::isfinite(histogram->GetBinError(bin)))
      {
        return false;
      }
    }
  }
  return true;
}

bool entries_equal(const TH1* histogram, unsigned long long expected)
{
  return histogram &&
      std::abs(histogram->GetEntries() - static_cast<double>(expected)) < 0.5;
}

std::set<std::string> expected_keys()
{
  std::set<std::string> result = {
      "metadata", "h_pi0_maximum_compatible_fraction_raw",
      "h_pi0_second_compatible_fraction_raw",
      "h_pi0_compatible_fraction_vs_cluster_energy_raw"};
  const std::array<std::string, 4> thresholds = {
      "0p0", "0p1", "0p3", "0p5"};
  const std::array<std::string, 2> pathways = {
      "g4_primary", "generator"};
  for (const std::string& threshold : thresholds)
  {
    result.insert("h_pi0_cluster_multiplicity_fmin_" + threshold + "_raw");
    result.insert("h_pi0_cluster_multiplicity_vs_truth_pt_fmin_" +
                  threshold + "_raw");
    for (const std::string& pathway : pathways)
    {
      result.insert("h_pi0_cluster_multiplicity_fmin_" + threshold + "_" +
                    pathway + "_raw");
    }
  }
  return result;
}
}

int check_pythia_pi0_cluster_multiplicity_partial(const std::string input_file)
{
  TFile input(input_file.c_str(), "READ");
  if (input.IsZombie())
  {
    return 1;
  }

  std::set<std::string> observed_keys;
  TIter next(input.GetListOfKeys());
  while (TKey* key = dynamic_cast<TKey*>(next()))
  {
    observed_keys.insert(key->GetName());
  }
  if (observed_keys != expected_keys())
  {
    std::cerr << "check_pythia_pi0_cluster_multiplicity_partial - unexpected top-level keys"
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
  int truth_matcher_algorithm_version = 0;
  int pt_bins = 0;
  int multiplicity_max = 0;
  int cluster_energy_bins = 0;
  double pt_max = 0.0;
  double cluster_energy_max = 0.0;
  double truth_eta_max = 0.0;
  unsigned char cluster_energy_cut_applied = 1U;
  std::array<double, 4> thresholds{};
  long long manifest_begin = -1;
  long long manifest_end = -1;
  std::string* manifest_path = nullptr;
  std::string* first_suffix = nullptr;
  std::string* last_suffix = nullptr;
  std::string* cluster_collection = nullptr;
  std::string* pi0_selection = nullptr;
  std::string* cluster_selection = nullptr;
  std::string* fraction_definition = nullptr;
  std::string* zero_threshold_definition = nullptr;
  unsigned long long events_processed = 0;
  unsigned long long events_written = 0;
  unsigned long long events_invalid = 0;
  unsigned long long cluster_considered = 0;
  unsigned long long cluster_invalid_truth = 0;
  unsigned long long candidate_count = 0;
  unsigned long long candidate_g4_primary = 0;
  unsigned long long candidate_generator = 0;
  unsigned long long malformed_daughters = 0;
  unsigned long long pair_evaluated = 0;
  unsigned long long pair_positive = 0;

  bool ok = true;
  ok &= bind(metadata, "schema_version", &schema_version);
  ok &= bind(metadata, "manifest_path", &manifest_path);
  ok &= bind(metadata, "manifest_begin", &manifest_begin);
  ok &= bind(metadata, "manifest_end", &manifest_end);
  ok &= bind(metadata, "first_suffix", &first_suffix);
  ok &= bind(metadata, "last_suffix", &last_suffix);
  ok &= bind(metadata, "cluster_collection", &cluster_collection);
  ok &= bind(metadata, "pi0_selection", &pi0_selection);
  ok &= bind(metadata, "cluster_selection", &cluster_selection);
  ok &= bind(metadata, "fraction_definition", &fraction_definition);
  ok &= bind(metadata, "zero_threshold_definition", &zero_threshold_definition);
  ok &= bind(metadata, "signal_embedding_id", &signal_embedding_id);
  ok &= bind(metadata, "truth_matcher_algorithm_version",
             &truth_matcher_algorithm_version);
  ok &= bind(metadata, "pt_bins", &pt_bins);
  ok &= bind(metadata, "pt_max", &pt_max);
  ok &= bind(metadata, "multiplicity_max", &multiplicity_max);
  ok &= bind(metadata, "cluster_energy_bins", &cluster_energy_bins);
  ok &= bind(metadata, "cluster_energy_max", &cluster_energy_max);
  ok &= bind(metadata, "truth_eta_max", &truth_eta_max);
  ok &= bind(metadata, "cluster_energy_cut_applied",
             &cluster_energy_cut_applied);
  ok &= bind(metadata, "fraction_thresholds", thresholds.data());
  ok &= bind(metadata, "events_processed", &events_processed);
  ok &= bind(metadata, "events_written", &events_written);
  ok &= bind(metadata, "events_invalid", &events_invalid);
  ok &= bind(metadata, "cluster_considered_count", &cluster_considered);
  ok &= bind(metadata, "cluster_invalid_truth_count", &cluster_invalid_truth);
  ok &= bind(metadata, "pi0_candidate_count", &candidate_count);
  ok &= bind(metadata, "pi0_candidate_g4_primary_count",
             &candidate_g4_primary);
  ok &= bind(metadata, "pi0_candidate_generator_count", &candidate_generator);
  ok &= bind(metadata, "pi0_malformed_daughters_count", &malformed_daughters);
  ok &= bind(metadata, "pi0_cluster_pair_evaluated_count", &pair_evaluated);
  ok &= bind(metadata, "pi0_cluster_pair_positive_count", &pair_positive);
  if (!ok || metadata->GetEntry(0) <= 0)
  {
    return 2;
  }

  const std::array<double, 4> expected_thresholds = {0.0, 0.1, 0.3, 0.5};
  const bool valid_metadata =
      schema_version == 3 && manifest_path && !manifest_path->empty() &&
      manifest_begin >= 0 && manifest_end > manifest_begin &&
      first_suffix && !first_suffix->empty() && last_suffix &&
      !last_suffix->empty() && cluster_collection &&
      *cluster_collection == "split" && pi0_selection &&
      !pi0_selection->empty() && cluster_selection &&
      *cluster_selection ==
          "finite_cluster_kinematics_without_eta_or_energy_threshold" &&
      fraction_definition &&
      !fraction_definition->empty() && zero_threshold_definition &&
      *zero_threshold_definition == "strictly_positive_fraction" &&
      signal_embedding_id > 0 && truth_matcher_algorithm_version > 0 &&
      pt_bins > 0 && pt_max > 0.0 && multiplicity_max > 0 &&
      cluster_energy_bins > 0 && cluster_energy_max > 0.0 &&
      truth_eta_max > 0.0 &&
      cluster_energy_cut_applied == 0U && thresholds == expected_thresholds &&
      events_written + events_invalid == events_processed &&
      cluster_invalid_truth <= cluster_considered &&
      candidate_count == candidate_g4_primary + candidate_generator &&
      pair_positive <= pair_evaluated;
  if (!valid_metadata)
  {
    std::cerr << "check_pythia_pi0_cluster_multiplicity_partial - invalid metadata"
              << std::endl;
    return 3;
  }

  const std::array<std::string, 4> threshold_tags = {
      "0p0", "0p1", "0p3", "0p5"};
  const std::array<std::string, 2> pathway_tags = {
      "g4_primary", "generator"};
  const std::array<unsigned long long, 2> pathway_entries = {
      candidate_g4_primary, candidate_generator};
  for (const std::string& threshold : threshold_tags)
  {
    TH1D* multiplicity = nullptr;
    TH2D* versus_pt = nullptr;
    input.GetObject(
        ("h_pi0_cluster_multiplicity_fmin_" + threshold + "_raw").c_str(),
        multiplicity);
    input.GetObject(
        ("h_pi0_cluster_multiplicity_vs_truth_pt_fmin_" + threshold +
         "_raw").c_str(),
        versus_pt);
    if (!finite_histogram(multiplicity) || !finite_histogram(versus_pt) ||
        !entries_equal(multiplicity, candidate_count) ||
        !entries_equal(versus_pt, candidate_count) ||
        multiplicity->GetNbinsX() != multiplicity_max + 1 ||
        versus_pt->GetNbinsX() != pt_bins ||
        versus_pt->GetNbinsY() != multiplicity_max + 1)
    {
      return 4;
    }
    for (std::size_t pathway = 0; pathway < pathway_tags.size(); ++pathway)
    {
      TH1D* histogram = nullptr;
      input.GetObject(
          ("h_pi0_cluster_multiplicity_fmin_" + threshold + "_" +
           pathway_tags[pathway] + "_raw").c_str(),
          histogram);
      if (!finite_histogram(histogram) ||
          !entries_equal(histogram, pathway_entries[pathway]) ||
          histogram->GetNbinsX() != multiplicity_max + 1)
      {
        return 4;
      }
    }
  }

  TH1D* maximum = nullptr;
  TH1D* second = nullptr;
  TH2D* fraction_vs_energy = nullptr;
  input.GetObject("h_pi0_maximum_compatible_fraction_raw", maximum);
  input.GetObject("h_pi0_second_compatible_fraction_raw", second);
  input.GetObject("h_pi0_compatible_fraction_vs_cluster_energy_raw", fraction_vs_energy);
  if (!finite_histogram(maximum) || !finite_histogram(second) ||
      !finite_histogram(fraction_vs_energy) ||
      !entries_equal(maximum, candidate_count) ||
      !entries_equal(second, candidate_count) ||
      !entries_equal(fraction_vs_energy, pair_positive) ||
      fraction_vs_energy->GetNbinsX() != cluster_energy_bins)
  {
    return 4;
  }

  std::cout << "check_pythia_pi0_cluster_multiplicity_partial - valid: "
            << input_file << ", events/pi0/positive pairs = "
            << events_written << "/" << candidate_count << "/"
            << pair_positive << std::endl;
  return 0;
}
