#include "../../macro/Utilities/sPhenixStyle.C"

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <THStack.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TObjArray.h>
#include <TPad.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace
{
constexpr std::size_t kHistogramCount = 14;
constexpr std::size_t kCategoryCount = 11;
constexpr std::array<std::size_t, kCategoryCount> kCategoryHistogramIndices = {2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13};
constexpr std::size_t kSummaryCategoryCount = 5;
constexpr std::array<std::size_t, kSummaryCategoryCount> kSummaryCategoryHistogramIndices = {2, 3, 4, 5, 13};
constexpr std::size_t kSummarySpectrumCount = 7;
constexpr std::array<std::size_t, kSummarySpectrumCount> kSummarySpectrumHistogramIndices = {0, 1, 2, 3, 4, 5, 13};
constexpr int kCanvasWidth = 1100;
constexpr int kCanvasHeight = 900;
constexpr double kPlotAreaTop = 0.62;
constexpr double kPlotLeftMargin = 0.13;
constexpr double kPlotRightMargin = 0.04;
constexpr double kPlotBottomMargin = 0.16;
constexpr double kPlotTopMargin = 0.04;
constexpr double kAnnotationX = 0.06;
constexpr double kLegendX = 0.42;
constexpr double kTextSize = 0.026;
constexpr double kAxisTextSize = 0.05;
constexpr double kXTitleOffset = 1.4;
constexpr double kYTitleOffset = 1.15;

void style_plot_axes(TAxis* x_axis, TAxis* y_axis) {
  x_axis->SetLabelSize(kAxisTextSize);
  x_axis->SetTitleSize(kAxisTextSize);
  x_axis->SetTitleOffset(kXTitleOffset);
  x_axis->CenterTitle();
  y_axis->SetLabelSize(kAxisTextSize);
  y_axis->SetTitleSize(kAxisTextSize);
  y_axis->SetTitleOffset(kYTitleOffset);
  y_axis->CenterTitle();
}

std::unique_ptr<TPad> make_plot_pad(const std::string& name, bool log_y = false) {
  auto pad = std::make_unique<TPad>(name.c_str(), "", 0.0, 0.0, 1.0, kPlotAreaTop);
  pad->SetLeftMargin(kPlotLeftMargin);
  pad->SetRightMargin(kPlotRightMargin);
  pad->SetBottomMargin(kPlotBottomMargin);
  pad->SetTopMargin(kPlotTopMargin);
  pad->SetTicks(1, 1);
  pad->SetLogy(log_y);
  pad->Draw();
  pad->cd();
  return pad;
}

struct PartialMetadata {
  std::string path;
  int schema_version = 0;
  std::string manifest_path;
  std::string cluster_collection;
  std::string tower_geom_node;
  std::string classification_unit;
  std::string pi0_selection;
  std::string partner_selection;
  std::string topology_definition;
  std::string topology_priority;
  std::string missing_category_priority;
  std::string response_policy;
  std::string photon_recovery_policy;
  std::string vertex_selection;
  long long manifest_begin = -1;
  long long manifest_end = -1;
  int signal_embedding_id = 0;
  int n_bins = 0;
  int matcher_version = 0;
  int topology_version = 0;
  double et_max = 0.0;
  double anchor_cluster_eta_max = 0.0;
  double partner_cluster_eta_max = 0.0;
  double cemc_acceptance_eta_max = 0.0;
  double pre_cemc_interaction_radius = 0.0;
  double min_cluster_energy = 0.0;
  double dominant_fraction_min = 0.0;
  double anchor_pi0_fraction_min = 0.0;
  double min_energy_contribution_fraction = 0.0;
  double min_photon_energy_recovery = 0.0;
  double min_direct_match_cluster_energy_coverage = 0.0;
  double missing_diagnostic_max_delta_r = 0.0;
  bool enable_missing_diagnostics = false;
  double max_abs_vertex_z = 0.0;
  unsigned char bin_width_normalized = 1U;
  unsigned long long events_processed = 0;
  unsigned long long events_written = 0;
  unsigned long long events_invalid = 0;
  unsigned long long events_vertex_rejected = 0;
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
  unsigned long long single_contaminated_count = 0;
  unsigned long long missing_count = 0;
  unsigned long long missing_energy_threshold_count = 0;
  unsigned long long missing_displaced_partner_cluster_count = 0;
  unsigned long long missing_acceptance_count = 0;
  unsigned long long missing_no_cemc_deposit_count = 0;
  unsigned long long missing_unclustered_deposit_count = 0;
  unsigned long long missing_match_incomplete_count = 0;
  unsigned long long missing_other_count = 0;
  unsigned long long other_count = 0;
};

template <class T>
bool bind(TTree* tree, const char* name, T* address) {
  return tree->GetBranch(name) && tree->SetBranchAddress(name, address) >= 0;
}

bool same_double(double left, double right) {
  return std::abs(left - right) <= 1e-12 * std::max({1.0, std::abs(left), std::abs(right)});
}

bool read_metadata(const std::string& path, PartialMetadata& value) {
  TFile input(path.c_str(), "READ");
  TTree* tree = nullptr;
  input.GetObject("metadata", tree);
  if (input.IsZombie() || !tree || tree->GetEntries() != 1) {
    return false;
  }

  std::string* manifest_path = nullptr;
  std::string* cluster_collection = nullptr;
  std::string* tower_geom_node = nullptr;
  std::string* classification_unit = nullptr;
  std::string* pi0_selection = nullptr;
  std::string* partner_selection = nullptr;
  std::string* topology_definition = nullptr;
  std::string* topology_priority = nullptr;
  std::string* missing_category_priority = nullptr;
  std::string* response_policy = nullptr;
  std::string* photon_recovery_policy = nullptr;
  std::string* vertex_selection = nullptr;
  bool ok = true;
  ok &= bind(tree, "schema_version", &value.schema_version);
  ok &= bind(tree, "manifest_path", &manifest_path);
  ok &= bind(tree, "manifest_begin", &value.manifest_begin);
  ok &= bind(tree, "manifest_end", &value.manifest_end);
  ok &= bind(tree, "cluster_collection", &cluster_collection);
  ok &= bind(tree, "tower_geom_node", &tower_geom_node);
  ok &= bind(tree, "classification_unit", &classification_unit);
  ok &= bind(tree, "pi0_selection", &pi0_selection);
  ok &= bind(tree, "partner_selection", &partner_selection);
  ok &= bind(tree, "topology_definition", &topology_definition);
  ok &= bind(tree, "topology_priority", &topology_priority);
  ok &= bind(tree, "missing_category_priority", &missing_category_priority);
  ok &= bind(tree, "response_policy", &response_policy);
  ok &= bind(tree, "photon_recovery_policy", &photon_recovery_policy);
  ok &= bind(tree, "vertex_selection", &vertex_selection);
  ok &= bind(tree, "signal_embedding_id", &value.signal_embedding_id);
  ok &= bind(tree, "n_bins", &value.n_bins);
  ok &= bind(tree, "et_max", &value.et_max);
  ok &= bind(tree, "anchor_cluster_eta_max", &value.anchor_cluster_eta_max);
  ok &= bind(tree, "partner_cluster_eta_max", &value.partner_cluster_eta_max);
  ok &= bind(tree, "cemc_acceptance_eta_max", &value.cemc_acceptance_eta_max);
  ok &= bind(tree, "pre_cemc_interaction_radius", &value.pre_cemc_interaction_radius);
  ok &= bind(tree, "min_cluster_energy", &value.min_cluster_energy);
  ok &= bind(tree, "dominant_fraction_min", &value.dominant_fraction_min);
  ok &= bind(tree, "anchor_pi0_fraction_min", &value.anchor_pi0_fraction_min);
  ok &= bind(tree, "min_energy_contribution_fraction", &value.min_energy_contribution_fraction);
  ok &= bind(tree, "min_photon_energy_recovery", &value.min_photon_energy_recovery);
  ok &= bind(tree, "min_direct_match_cluster_energy_coverage", &value.min_direct_match_cluster_energy_coverage);
  ok &= bind(tree, "missing_diagnostic_max_delta_r", &value.missing_diagnostic_max_delta_r);
  ok &= bind(tree, "enable_missing_diagnostics", &value.enable_missing_diagnostics);
  ok &= bind(tree, "max_abs_vertex_z", &value.max_abs_vertex_z);
  ok &= bind(tree, "pi0_truth_matching_algorithm_version", &value.matcher_version);
  ok &= bind(tree, "pi0_topology_algorithm_version", &value.topology_version);
  ok &= bind(tree, "bin_width_normalized", &value.bin_width_normalized);
  ok &= bind(tree, "events_processed", &value.events_processed);
  ok &= bind(tree, "events_written", &value.events_written);
  ok &= bind(tree, "events_invalid", &value.events_invalid);
  ok &= bind(tree, "events_vertex_rejected", &value.events_vertex_rejected);
  ok &= bind(tree, "cluster_considered_count", &value.cluster_considered);
  ok &= bind(tree, "cluster_invalid_truth_count", &value.cluster_invalid_truth);
  ok &= bind(tree, "prompt_cluster_count", &value.prompt_count);
  ok &= bind(tree, "pi0_candidate_g4_decay_count", &value.candidate_g4);
  ok &= bind(tree, "pi0_candidate_generator_decay_count", &value.candidate_generator);
  ok &= bind(tree, "pi0_malformed_daughters_count", &value.malformed_daughters);
  ok &= bind(tree, "anchor_cluster_count", &value.anchor_count);
  ok &= bind(tree, "anchor_g4_decay_count", &value.anchor_g4);
  ok &= bind(tree, "anchor_generator_decay_count", &value.anchor_generator);
  ok &= bind(tree, "anchor_ambiguous_main_count", &value.ambiguous_main);
  ok &= bind(tree, "energy_match_invalid_count", &value.energy_match_invalid);
  ok &= bind(tree, "separated_count", &value.separated_count);
  ok &= bind(tree, "merged_count", &value.merged_count);
  ok &= bind(tree, "single_contaminated_count", &value.single_contaminated_count);
  ok &= bind(tree, "missing_count", &value.missing_count);
  ok &= bind(tree, "missing_energy_threshold_count", &value.missing_energy_threshold_count);
  ok &= bind(tree, "missing_displaced_partner_cluster_count", &value.missing_displaced_partner_cluster_count);
  ok &= bind(tree, "missing_acceptance_count", &value.missing_acceptance_count);
  ok &= bind(tree, "missing_no_cemc_deposit_count", &value.missing_no_cemc_deposit_count);
  ok &= bind(tree, "missing_unclustered_deposit_count", &value.missing_unclustered_deposit_count);
  ok &= bind(tree, "missing_match_incomplete_count", &value.missing_match_incomplete_count);
  ok &= bind(tree, "missing_other_count", &value.missing_other_count);
  ok &= bind(tree, "other_count", &value.other_count);
  if (!ok || tree->GetEntry(0) <= 0 || !manifest_path ||
      !cluster_collection || !tower_geom_node || !classification_unit || !pi0_selection ||
      !partner_selection || !topology_definition || !topology_priority ||
      !missing_category_priority ||
      !response_policy || !photon_recovery_policy || !vertex_selection)
  {
    return false;
  }

  value.path = path;
  value.manifest_path = *manifest_path;
  value.cluster_collection = *cluster_collection;
  value.tower_geom_node = *tower_geom_node;
  value.classification_unit = *classification_unit;
  value.pi0_selection = *pi0_selection;
  value.partner_selection = *partner_selection;
  value.topology_definition = *topology_definition;
  value.topology_priority = *topology_priority;
  value.missing_category_priority = *missing_category_priority;
  value.response_policy = *response_policy;
  value.photon_recovery_policy = *photon_recovery_policy;
  value.vertex_selection = *vertex_selection;
  return true;
}

bool valid_metadata(const PartialMetadata& value) {
  return value.schema_version == 8 &&
      !value.manifest_path.empty() &&
      value.manifest_begin >= 0 &&
      value.manifest_end > value.manifest_begin &&
      value.cluster_collection == "split" &&
      !value.tower_geom_node.empty() &&
      value.classification_unit == "every_cluster_with_selected_pi0_as_grouped_main_contributor" &&
      value.pi0_selection == "signal_g4_primary_pi0_or_generator_pi0_with_exactly_two_g4_photons" &&
      value.partner_selection == "same_energy_cut_as_anchor_partner_eta_cut_configurable" &&
      value.topology_definition == "anchor_membership_in_recovered_direct_daughter_maximum_deposit_clusters_with_single_contaminated_pre_cemc_split" &&
      value.topology_priority == "ambiguous_main_to_other_then_single_contaminated_then_merged_then_separated_then_missing_then_other" &&
      value.missing_category_priority == "projection_then_acceptance_then_cemc_deposit_then_threshold_near_or_displaced_then_recovery_then_match_incomplete_then_unclustered_then_other" &&
      value.response_policy == "not_used_for_classification" &&
      value.photon_recovery_policy == "cluster_energy_times_gamma_deposit_fraction_over_truth_energy_threshold" &&
      value.vertex_selection == "signal_hepmc_collision_vertex_abs_z_lt_max" &&
      value.signal_embedding_id > 0 && value.n_bins > 0 &&
      value.et_max > 0.0 &&
      value.anchor_cluster_eta_max > 0.0 &&
      std::isfinite(value.partner_cluster_eta_max) &&
      std::isfinite(value.cemc_acceptance_eta_max) && value.cemc_acceptance_eta_max > 0.0 &&
      std::isfinite(value.pre_cemc_interaction_radius) && value.pre_cemc_interaction_radius > 0.0 &&
      value.min_cluster_energy >= 0.0 &&
      value.dominant_fraction_min >= 0.0 &&
      value.dominant_fraction_min <= 1.0 &&
      value.anchor_pi0_fraction_min >= 0.0 &&
      value.anchor_pi0_fraction_min <= 1.0 &&
      value.min_energy_contribution_fraction >= 0.0 &&
      value.min_energy_contribution_fraction < 1.0 &&
      std::isfinite(value.min_photon_energy_recovery) &&
      value.min_photon_energy_recovery >= 0.0 &&
      value.min_photon_energy_recovery <= 1.0 &&
      value.min_direct_match_cluster_energy_coverage >= 0.0 &&
      value.min_direct_match_cluster_energy_coverage <= 1.0 &&
      std::isfinite(value.missing_diagnostic_max_delta_r) && value.missing_diagnostic_max_delta_r > 0.0 &&
      std::isfinite(value.max_abs_vertex_z) && value.max_abs_vertex_z > 0.0 &&
      value.matcher_version > 0 && value.topology_version > 0 &&
      value.bin_width_normalized == 0U &&
      value.events_processed > 0 &&
      value.events_written + value.events_invalid +
          value.events_vertex_rejected == value.events_processed &&
      value.cluster_invalid_truth <= value.cluster_considered &&
      value.anchor_count == value.anchor_g4 + value.anchor_generator &&
      value.anchor_count == value.separated_count + value.merged_count + value.single_contaminated_count + value.missing_count + value.other_count &&
      value.missing_count == value.missing_energy_threshold_count + value.missing_displaced_partner_cluster_count + value.missing_acceptance_count +
          value.missing_no_cemc_deposit_count + value.missing_unclustered_deposit_count + value.missing_match_incomplete_count + value.missing_other_count &&
      (value.enable_missing_diagnostics || (value.missing_energy_threshold_count == 0 && value.missing_displaced_partner_cluster_count == 0 &&
          value.missing_no_cemc_deposit_count == 0 && value.missing_unclustered_deposit_count == 0 && value.missing_match_incomplete_count == 0)) &&
      value.ambiguous_main <= value.other_count;
}

bool compatible(const PartialMetadata& value, const PartialMetadata& reference) {
  return valid_metadata(value) &&
      value.schema_version == reference.schema_version &&
      value.manifest_path == reference.manifest_path &&
      value.cluster_collection == reference.cluster_collection &&
      value.tower_geom_node == reference.tower_geom_node &&
      value.classification_unit == reference.classification_unit &&
      value.pi0_selection == reference.pi0_selection &&
      value.partner_selection == reference.partner_selection &&
      value.topology_definition == reference.topology_definition &&
      value.topology_priority == reference.topology_priority &&
      value.missing_category_priority == reference.missing_category_priority &&
      value.response_policy == reference.response_policy &&
      value.photon_recovery_policy == reference.photon_recovery_policy &&
      value.vertex_selection == reference.vertex_selection &&
      value.signal_embedding_id == reference.signal_embedding_id &&
      value.n_bins == reference.n_bins &&
      value.matcher_version == reference.matcher_version &&
      value.topology_version == reference.topology_version &&
      same_double(value.et_max, reference.et_max) &&
      same_double(value.anchor_cluster_eta_max, reference.anchor_cluster_eta_max) &&
      same_double(value.partner_cluster_eta_max, reference.partner_cluster_eta_max) &&
      same_double(value.cemc_acceptance_eta_max, reference.cemc_acceptance_eta_max) &&
      same_double(value.pre_cemc_interaction_radius, reference.pre_cemc_interaction_radius) &&
      same_double(value.min_cluster_energy, reference.min_cluster_energy) &&
      same_double(value.dominant_fraction_min, reference.dominant_fraction_min) &&
      same_double(value.anchor_pi0_fraction_min, reference.anchor_pi0_fraction_min) &&
      same_double(value.min_energy_contribution_fraction, reference.min_energy_contribution_fraction) &&
      same_double(value.min_photon_energy_recovery, reference.min_photon_energy_recovery) &&
      same_double(value.min_direct_match_cluster_energy_coverage, reference.min_direct_match_cluster_energy_coverage) &&
      same_double(value.missing_diagnostic_max_delta_r, reference.missing_diagnostic_max_delta_r) &&
      value.enable_missing_diagnostics == reference.enable_missing_diagnostics &&
      same_double(value.max_abs_vertex_z, reference.max_abs_vertex_z);
}

bool valid_histogram(const TH1D* histogram, const PartialMetadata& metadata, unsigned long long expected_entries) {
  if (!histogram || histogram->GetNbinsX() != metadata.n_bins ||
      std::abs(histogram->GetXaxis()->GetXmin()) > 1e-12 ||
      !same_double(histogram->GetXaxis()->GetXmax(), metadata.et_max) ||
      histogram->GetSumw2N() == 0 ||
      std::abs(histogram->GetEntries() -
               static_cast<double>(expected_entries)) > 0.5)
  {
    return false;
  }
  for (int bin = 0; bin <= histogram->GetNbinsX() + 1; ++bin) {
    if (!std::isfinite(histogram->GetBinContent(bin)) || !std::isfinite(histogram->GetBinError(bin))) {
      return false;
    }
  }
  return true;
}

bool make_output_directory(const std::string& output_base) {
  const std::size_t slash = output_base.find_last_of('/');
  if (slash == std::string::npos) {
    return true;
  }
  const std::string directory = output_base.substr(0, slash);
  return directory.empty() ||
      !gSystem->AccessPathName(directory.c_str()) ||
      gSystem->mkdir(directory.c_str(), true) == 0;
}

double smallest_positive(const std::array<std::unique_ptr<TH1D>, kHistogramCount>& histograms) {
  double result = std::numeric_limits<double>::infinity();
  for (const auto& histogram : histograms) {
    for (int bin = 1; bin <= histogram->GetNbinsX(); ++bin) {
      const double value = histogram->GetBinContent(bin);
      if (value > 0.0 && value < result) {
        result = value;
      }
    }
  }
  return std::isfinite(result) ? result : 0.0;
}
}

int FinalizePythiaPi0AnchorClusterSpectra(
    const std::string partial_pattern = "output/pi0_anchor_topology_partial/eta07_zvtx60_full_partner_fgamma0p0_recovery0p2_clusterenergy/partial_*.root",
    const std::string output_base = "output/plots/pi0_anchor_topology/minimum_bias/with_vertex_cut60",
    const long long expected_manifest_begin = 0,
    const long long expected_manifest_end = -1,
    const std::string sample_label = "Pythia8 p+p MB")
{
  if (partial_pattern.empty() || output_base.empty() || sample_label.empty() ||
      expected_manifest_begin < 0 || expected_manifest_end < -1 ||
      (expected_manifest_end >= 0 &&
       expected_manifest_end <= expected_manifest_begin))
  {
    return 1;
  }

  TChain chain("metadata");
  const int matched = chain.Add(partial_pattern.c_str());
  const TObjArray* files = chain.GetListOfFiles();
  if (matched <= 0 || !files || files->GetEntries() <= 0) {
    std::cerr << "FinalizePythiaPi0AnchorClusterSpectra - no partials matched" << std::endl;
    return 2;
  }

  std::vector<PartialMetadata> partials;
  std::set<std::string> unique_paths;
  for (int index = 0; index < files->GetEntries(); ++index) {
    const TObject* element = files->At(index);
    const std::string path = element ? element->GetTitle() : "";
    PartialMetadata metadata;
    if (path.empty() || !unique_paths.insert(path).second || !read_metadata(path, metadata)) {
      std::cerr << "FinalizePythiaPi0AnchorClusterSpectra - invalid partial: " << path << std::endl;
      return 3;
    }
    partials.push_back(metadata);
  }

  std::sort(partials.begin(), partials.end(),
      [](const auto& left, const auto& right) {
        return left.manifest_begin < right.manifest_begin;
      });
  const PartialMetadata& reference = partials.front();
  long long next_begin = expected_manifest_begin;
  for (const PartialMetadata& partial : partials) {
    if (!compatible(partial, reference) || partial.manifest_begin != next_begin) {
      std::cerr << "FinalizePythiaPi0AnchorClusterSpectra - incompatible or noncontiguous partial: "
                << partial.path << ", expected begin " << next_begin
                << std::endl;
      return 4;
    }
    next_begin = partial.manifest_end;
  }
  if (expected_manifest_end >= 0 && next_begin != expected_manifest_end) {
    return 4;
  }

  const std::array<std::string, kHistogramCount> raw_names = {
      "h_prompt_cluster_et_raw", "h_pi0_anchor_cluster_et_raw",
      "h_pi0_anchor_separated_cluster_et_raw",
      "h_pi0_anchor_merged_cluster_et_raw",
      "h_pi0_anchor_single_contaminated_cluster_et_raw",
      "h_pi0_anchor_missing_cluster_et_raw",
      "h_pi0_anchor_missing_energy_threshold_cluster_et_raw",
      "h_pi0_anchor_missing_displaced_partner_cluster_et_raw",
      "h_pi0_anchor_missing_acceptance_cluster_et_raw",
      "h_pi0_anchor_missing_no_cemc_deposit_cluster_et_raw",
      "h_pi0_anchor_missing_unclustered_deposit_cluster_et_raw",
      "h_pi0_anchor_missing_match_incomplete_cluster_et_raw",
      "h_pi0_anchor_missing_other_cluster_et_raw",
      "h_pi0_anchor_other_cluster_et_raw"};
  const std::array<std::string, kHistogramCount> density_names = {
      "h_prompt_cluster_et_density",
      "h_pi0_anchor_cluster_et_density",
      "h_pi0_anchor_separated_cluster_et_density",
      "h_pi0_anchor_merged_cluster_et_density",
      "h_pi0_anchor_single_contaminated_cluster_et_density",
      "h_pi0_anchor_missing_cluster_et_density",
      "h_pi0_anchor_missing_energy_threshold_cluster_et_density",
      "h_pi0_anchor_missing_displaced_partner_cluster_et_density",
      "h_pi0_anchor_missing_acceptance_cluster_et_density",
      "h_pi0_anchor_missing_no_cemc_deposit_cluster_et_density",
      "h_pi0_anchor_missing_unclustered_deposit_cluster_et_density",
      "h_pi0_anchor_missing_match_incomplete_cluster_et_density",
      "h_pi0_anchor_missing_other_cluster_et_density",
      "h_pi0_anchor_other_cluster_et_density"};

  std::array<std::unique_ptr<TH1D>, kHistogramCount> raw;
  for (std::size_t index = 0; index < raw.size(); ++index) {
    raw[index] = std::make_unique<TH1D>(raw_names[index].c_str(), "", reference.n_bins, 0.0, reference.et_max);
    raw[index]->Sumw2();
  }

  PartialMetadata total = reference;
  total.events_processed = total.events_written = total.events_invalid = 0;
  total.events_vertex_rejected = 0;
  total.cluster_considered = total.cluster_invalid_truth = 0;
  total.prompt_count = total.candidate_g4 = total.candidate_generator = 0;
  total.malformed_daughters = total.anchor_count = 0;
  total.anchor_g4 = total.anchor_generator = total.ambiguous_main = 0;
  total.energy_match_invalid = 0;
  total.separated_count = total.merged_count = total.single_contaminated_count = 0;
  total.missing_count = total.missing_energy_threshold_count = total.missing_displaced_partner_cluster_count = 0;
  total.missing_acceptance_count = total.missing_no_cemc_deposit_count = total.missing_unclustered_deposit_count = 0;
  total.missing_match_incomplete_count = total.missing_other_count = total.other_count = 0;

  for (const PartialMetadata& partial : partials) {
    TFile input(partial.path.c_str(), "READ");
    const std::array<unsigned long long, kHistogramCount> counts = {
        partial.prompt_count, partial.anchor_count,
        partial.separated_count, partial.merged_count, partial.single_contaminated_count,
        partial.missing_count, partial.missing_energy_threshold_count, partial.missing_displaced_partner_cluster_count,
        partial.missing_acceptance_count, partial.missing_no_cemc_deposit_count, partial.missing_unclustered_deposit_count,
        partial.missing_match_incomplete_count, partial.missing_other_count, partial.other_count};
    std::array<TH1D*, kHistogramCount> partial_histograms{};
    for (std::size_t index = 0; index < raw.size(); ++index) {
      input.GetObject(raw_names[index].c_str(), partial_histograms[index]);
      if (input.IsZombie() ||
          !valid_histogram(partial_histograms[index], partial, counts[index]) ||
          !raw[index]->Add(partial_histograms[index]))
      {
        std::cerr
            << "FinalizePythiaPi0AnchorClusterSpectra - invalid histogram in "
            << partial.path << std::endl;
        return 5;
      }
    }
    for (int bin = 0; bin <= partial.n_bins + 1; ++bin) {
      double missing_categories = 0.0;
      for (std::size_t index = 6; index <= 12; ++index) missing_categories += partial_histograms[index]->GetBinContent(bin);
      const double categories = partial_histograms[2]->GetBinContent(bin) + partial_histograms[3]->GetBinContent(bin) +
          partial_histograms[4]->GetBinContent(bin) + missing_categories + partial_histograms[13]->GetBinContent(bin);
      if (std::abs(partial_histograms[5]->GetBinContent(bin) - missing_categories) > 1e-9 ||
          std::abs(partial_histograms[1]->GetBinContent(bin) - categories) > 1e-9) {
        return 5;
      }
    }

    total.events_processed += partial.events_processed;
    total.events_written += partial.events_written;
    total.events_invalid += partial.events_invalid;
    total.events_vertex_rejected += partial.events_vertex_rejected;
    total.cluster_considered += partial.cluster_considered;
    total.cluster_invalid_truth += partial.cluster_invalid_truth;
    total.prompt_count += partial.prompt_count;
    total.candidate_g4 += partial.candidate_g4;
    total.candidate_generator += partial.candidate_generator;
    total.malformed_daughters += partial.malformed_daughters;
    total.anchor_count += partial.anchor_count;
    total.anchor_g4 += partial.anchor_g4;
    total.anchor_generator += partial.anchor_generator;
    total.ambiguous_main += partial.ambiguous_main;
    total.energy_match_invalid += partial.energy_match_invalid;
    total.separated_count += partial.separated_count;
    total.merged_count += partial.merged_count;
    total.single_contaminated_count += partial.single_contaminated_count;
    total.missing_count += partial.missing_count;
    total.missing_energy_threshold_count += partial.missing_energy_threshold_count;
    total.missing_displaced_partner_cluster_count += partial.missing_displaced_partner_cluster_count;
    total.missing_acceptance_count += partial.missing_acceptance_count;
    total.missing_no_cemc_deposit_count += partial.missing_no_cemc_deposit_count;
    total.missing_unclustered_deposit_count += partial.missing_unclustered_deposit_count;
    total.missing_match_incomplete_count += partial.missing_match_incomplete_count;
    total.missing_other_count += partial.missing_other_count;
    total.other_count += partial.other_count;
  }

  for (int bin = 0; bin <= reference.n_bins + 1; ++bin) {
    double missing_categories = 0.0;
    for (std::size_t index = 6; index <= 12; ++index) missing_categories += raw[index]->GetBinContent(bin);
    const double categories = raw[2]->GetBinContent(bin) + raw[3]->GetBinContent(bin) +
        raw[4]->GetBinContent(bin) + missing_categories + raw[13]->GetBinContent(bin);
    if (std::abs(raw[5]->GetBinContent(bin) - missing_categories) > 1e-9 ||
        std::abs(raw[1]->GetBinContent(bin) - categories) > 1e-9) {
      return 5;
    }
  }

  std::array<std::unique_ptr<TH1D>, kHistogramCount> density;
  const std::array<int, kHistogramCount> colors = {
      kRed + 1, kBlue + 1, kAzure + 7, kMagenta + 1, kCyan + 2, kGreen + 2,
      kOrange + 7, kPink + 7, kViolet + 1, kYellow + 2, kSpring + 5, kBlue + 3, kGreen + 3, kGray + 2};
  for (std::size_t index = 0; index < raw.size(); ++index) {
    density[index].reset(static_cast<TH1D*>(raw[index]->Clone(density_names[index].c_str())));
    density[index]->SetDirectory(nullptr);
    density[index]->Scale(1.0, "width");
    density[index]->SetStats(false);
    density[index]->SetLineColor(colors[index]);
    density[index]->SetLineWidth(index < 2 ? 3 : 2);
    density[index]->GetXaxis()->SetTitle("Anchor Cluster E_{T} [GeV]");
    density[index]->GetYaxis()->SetTitle("Clusters / GeV");
    style_plot_axes(density[index]->GetXaxis(), density[index]->GetYaxis());

    raw[index]->SetStats(false);
    raw[index]->SetLineColor(colors[index]);
    raw[index]->SetLineWidth(index < 2 ? 3 : 2);
    raw[index]->GetXaxis()->SetTitle("Anchor Cluster E_{T} [GeV]");
    raw[index]->GetYaxis()->SetTitle("Counts");
    style_plot_axes(raw[index]->GetXaxis(), raw[index]->GetYaxis());
  }

  const std::array<std::string, kCategoryCount> detailed_fraction_names = {
      "h_pi0_anchor_separated_fraction",
      "h_pi0_anchor_merged_fraction",
      "h_pi0_anchor_single_contaminated_fraction",
      "h_pi0_anchor_missing_energy_threshold_fraction",
      "h_pi0_anchor_missing_displaced_partner_cluster_fraction",
      "h_pi0_anchor_missing_acceptance_fraction",
      "h_pi0_anchor_missing_no_cemc_deposit_fraction",
      "h_pi0_anchor_missing_unclustered_deposit_fraction",
      "h_pi0_anchor_missing_match_incomplete_fraction",
      "h_pi0_anchor_missing_other_fraction",
      "h_pi0_anchor_other_fraction"};
  std::array<std::unique_ptr<TH1D>, kCategoryCount> detailed_fractions;
  for (std::size_t index = 0; index < detailed_fractions.size(); ++index) {
    const std::size_t histogram_index = kCategoryHistogramIndices[index];
    detailed_fractions[index].reset(static_cast<TH1D*>(raw[histogram_index]->Clone(detailed_fraction_names[index].c_str())));
    detailed_fractions[index]->SetDirectory(nullptr);
    detailed_fractions[index]->Divide(raw[histogram_index].get(), raw[1].get(), 1.0, 1.0, "B");
    detailed_fractions[index]->SetStats(false);
    detailed_fractions[index]->SetLineColor(colors[histogram_index]);
    detailed_fractions[index]->SetMarkerColor(colors[histogram_index]);
    detailed_fractions[index]->SetMarkerStyle(20 + static_cast<int>(index));
    detailed_fractions[index]->SetMarkerSize(0.9);
    detailed_fractions[index]->SetLineWidth(2);
    detailed_fractions[index]->GetXaxis()->SetTitle("Anchor Cluster E_{T} [GeV]");
    detailed_fractions[index]->GetYaxis()->SetTitle("Fraction");
    style_plot_axes(detailed_fractions[index]->GetXaxis(), detailed_fractions[index]->GetYaxis());
  }

  std::array<std::unique_ptr<TH1D>, kCategoryCount> detailed_stacked_fractions;
  for (std::size_t index = 0; index < detailed_stacked_fractions.size(); ++index) {
    detailed_stacked_fractions[index].reset(static_cast<TH1D*>(detailed_fractions[index]->Clone(
        (detailed_fraction_names[index] + "_stack_component").c_str())));
    detailed_stacked_fractions[index]->SetDirectory(nullptr);
    detailed_stacked_fractions[index]->SetFillColor(colors[kCategoryHistogramIndices[index]]);
    detailed_stacked_fractions[index]->SetLineColor(kBlack);
    detailed_stacked_fractions[index]->SetLineWidth(1);
    detailed_stacked_fractions[index]->SetMarkerStyle(0);
  }

  const std::array<std::string, kSummaryCategoryCount> summary_fraction_names = {
      "h_pi0_anchor_summary_separated_fraction",
      "h_pi0_anchor_summary_merged_fraction",
      "h_pi0_anchor_summary_single_contaminated_fraction",
      "h_pi0_anchor_summary_missing_fraction",
      "h_pi0_anchor_summary_other_fraction"};
  std::array<std::unique_ptr<TH1D>, kSummaryCategoryCount> summary_fractions;
  for (std::size_t index = 0; index < summary_fractions.size(); ++index) {
    const std::size_t histogram_index = kSummaryCategoryHistogramIndices[index];
    summary_fractions[index].reset(static_cast<TH1D*>(raw[histogram_index]->Clone(summary_fraction_names[index].c_str())));
    summary_fractions[index]->SetDirectory(nullptr);
    summary_fractions[index]->Divide(raw[histogram_index].get(), raw[1].get(), 1.0, 1.0, "B");
    summary_fractions[index]->SetStats(false);
    summary_fractions[index]->SetLineColor(colors[histogram_index]);
    summary_fractions[index]->SetMarkerColor(colors[histogram_index]);
    summary_fractions[index]->SetMarkerStyle(20 + static_cast<int>(index));
    summary_fractions[index]->SetMarkerSize(0.9);
    summary_fractions[index]->SetLineWidth(2);
    summary_fractions[index]->GetXaxis()->SetTitle("Anchor Cluster E_{T} [GeV]");
    summary_fractions[index]->GetYaxis()->SetTitle("Fraction");
    style_plot_axes(summary_fractions[index]->GetXaxis(), summary_fractions[index]->GetYaxis());
  }

  std::array<std::unique_ptr<TH1D>, kSummaryCategoryCount> summary_stacked_fractions;
  for (std::size_t index = 0; index < summary_stacked_fractions.size(); ++index) {
    summary_stacked_fractions[index].reset(static_cast<TH1D*>(summary_fractions[index]->Clone(
        (summary_fraction_names[index] + "_stack_component").c_str())));
    summary_stacked_fractions[index]->SetDirectory(nullptr);
    summary_stacked_fractions[index]->SetFillColor(colors[kSummaryCategoryHistogramIndices[index]]);
    summary_stacked_fractions[index]->SetLineColor(kBlack);
    summary_stacked_fractions[index]->SetLineWidth(1);
    summary_stacked_fractions[index]->SetMarkerStyle(0);
  }

  if (!make_output_directory(output_base)) {
    return 6;
  }
  SetsPhenixStyle();

  TCanvas spectrum_canvas("c_pythia_pi0_anchor_cluster_et", "Pythia pi0 anchor cluster spectra", kCanvasWidth, kCanvasHeight);
  spectrum_canvas.SetCanvasSize(kCanvasWidth, kCanvasHeight);
  auto spectrum_plot_pad = make_plot_pad("pythia_pi0_anchor_cluster_et_plot", true);
  double maximum = 0.0;
  for (const auto& histogram : raw) {
    maximum = std::max(maximum, histogram->GetMaximum());
  }
  const double minimum = smallest_positive(raw);
  raw[0]->SetMinimum(minimum > 0.0 ? 0.5 * minimum : 0.5);
  raw[0]->SetMaximum(maximum > 0.0 ? 2.0 * maximum : 1.0);
  raw[0]->Draw("HIST");
  for (std::size_t index = 1; index < raw.size(); ++index) {
    raw[index]->Draw("HIST SAME");
  }
  spectrum_canvas.cd();
  TLegend spectrum_legend(0.50, 0.62, 0.95, 0.97);
  spectrum_legend.SetBorderSize(0);
  spectrum_legend.SetFillStyle(0);
  spectrum_legend.SetTextSize(0.016);
  spectrum_legend.AddEntry(raw[0].get(), "Prompt-#gamma cluster", "l");
  spectrum_legend.AddEntry(raw[1].get(), "#pi^{0}-main anchor", "l");
  spectrum_legend.AddEntry(raw[2].get(), "Separated", "l");
  spectrum_legend.AddEntry(raw[3].get(), "Merged", "l");
  spectrum_legend.AddEntry(raw[4].get(), "Single contaminated", "l");
  spectrum_legend.AddEntry(raw[5].get(), "Missing (total)", "l");
  spectrum_legend.AddEntry(raw[6].get(), "Missing: energy threshold", "l");
  spectrum_legend.AddEntry(raw[7].get(), "Missing: displaced partner cluster", "l");
  spectrum_legend.AddEntry(raw[8].get(), "Missing: acceptance", "l");
  spectrum_legend.AddEntry(raw[9].get(), "Missing: no CEMC deposit", "l");
  spectrum_legend.AddEntry(raw[10].get(), "Missing: unclustered deposit", "l");
  spectrum_legend.AddEntry(raw[11].get(), "Missing: match incomplete", "l");
  spectrum_legend.AddEntry(raw[12].get(), "Missing: other", "l");
  spectrum_legend.AddEntry(raw[13].get(), "Other", "l");
  spectrum_legend.Draw();

  TLatex spectrum_label;
  spectrum_label.SetNDC();
  spectrum_label.SetTextAlign(13);
  spectrum_label.SetTextSize(kTextSize);
  spectrum_label.DrawLatex(kAnnotationX, 0.96, "#it{#bf{sPHENIX}} Internal");
  spectrum_label.DrawLatex(kAnnotationX, 0.91, sample_label.c_str());
  spectrum_label.DrawLatex(kAnnotationX, 0.86, "Category unit: #pi^{0}-main anchor cluster");
  std::ostringstream anchor_label;
  anchor_label << "Anchor |#eta| < " << reference.anchor_cluster_eta_max;
  spectrum_label.DrawLatex(kAnnotationX, 0.81, anchor_label.str().c_str());
  std::ostringstream min_cluster_energy_label;
  min_cluster_energy_label << "E_{min}^{cluster} = " << reference.min_cluster_energy << " GeV";
  spectrum_label.DrawLatex(kAnnotationX, 0.76, min_cluster_energy_label.str().c_str());
  std::ostringstream recovery_label;
  recovery_label << "E_{clus} f_{dep}^{#gamma}/E_{truth}^{#gamma} #geq " << reference.min_photon_energy_recovery;
  spectrum_label.DrawLatex(kAnnotationX, 0.71, recovery_label.str().c_str());
  std::ostringstream vertex_label;
  vertex_label << "|z_{vtx}^{truth}| < " << reference.max_abs_vertex_z << " cm";
  spectrum_label.DrawLatex(kAnnotationX, 0.66, vertex_label.str().c_str());
  spectrum_plot_pad->cd();
  spectrum_plot_pad->RedrawAxis();
  spectrum_canvas.cd();
  spectrum_canvas.SaveAs((output_base + "_detailed.pdf").c_str());

  double summary_spectrum_maximum = 0.0;
  double summary_spectrum_minimum = std::numeric_limits<double>::infinity();
  for (const std::size_t histogram_index : kSummarySpectrumHistogramIndices) {
    summary_spectrum_maximum = std::max(summary_spectrum_maximum, raw[histogram_index]->GetMaximum());
    for (int bin = 1; bin <= raw[histogram_index]->GetNbinsX(); ++bin) {
      const double value = raw[histogram_index]->GetBinContent(bin);
      if (value > 0.0) summary_spectrum_minimum = std::min(summary_spectrum_minimum, value);
    }
  }
  if (!std::isfinite(summary_spectrum_minimum)) summary_spectrum_minimum = 0.0;

  TCanvas summary_spectrum_canvas("c_pythia_pi0_anchor_cluster_et_summary", "Summary Pythia pi0 anchor cluster spectra", kCanvasWidth, kCanvasHeight);
  summary_spectrum_canvas.SetCanvasSize(kCanvasWidth, kCanvasHeight);
  auto summary_spectrum_plot_pad = make_plot_pad("pythia_pi0_anchor_cluster_et_summary_plot", true);
  raw[0]->SetMinimum(summary_spectrum_minimum > 0.0 ? 0.5 * summary_spectrum_minimum : 0.5);
  raw[0]->SetMaximum(summary_spectrum_maximum > 0.0 ? 2.0 * summary_spectrum_maximum : 1.0);
  raw[0]->Draw("HIST");
  for (std::size_t index = 1; index < kSummarySpectrumHistogramIndices.size(); ++index) raw[kSummarySpectrumHistogramIndices[index]]->Draw("HIST SAME");
  summary_spectrum_canvas.cd();
  TLegend summary_spectrum_legend(0.55, 0.65, 0.94, 0.95);
  summary_spectrum_legend.SetBorderSize(0);
  summary_spectrum_legend.SetFillStyle(0);
  summary_spectrum_legend.SetTextSize(kTextSize);
  summary_spectrum_legend.AddEntry(raw[0].get(), "Prompt-#gamma cluster", "l");
  summary_spectrum_legend.AddEntry(raw[1].get(), "#pi^{0}-main anchor", "l");
  summary_spectrum_legend.AddEntry(raw[2].get(), "Separated", "l");
  summary_spectrum_legend.AddEntry(raw[3].get(), "Merged", "l");
  summary_spectrum_legend.AddEntry(raw[4].get(), "Single contaminated", "l");
  summary_spectrum_legend.AddEntry(raw[5].get(), "Missing", "l");
  summary_spectrum_legend.AddEntry(raw[13].get(), "Other", "l");
  summary_spectrum_legend.Draw();
  spectrum_label.DrawLatex(kAnnotationX, 0.96, "#it{#bf{sPHENIX}} Internal");
  spectrum_label.DrawLatex(kAnnotationX, 0.91, sample_label.c_str());
  spectrum_label.DrawLatex(kAnnotationX, 0.86, "Category unit: #pi^{0}-main anchor cluster");
  spectrum_label.DrawLatex(kAnnotationX, 0.81, anchor_label.str().c_str());
  spectrum_label.DrawLatex(kAnnotationX, 0.76, min_cluster_energy_label.str().c_str());
  spectrum_label.DrawLatex(kAnnotationX, 0.71, recovery_label.str().c_str());
  spectrum_label.DrawLatex(kAnnotationX, 0.66, vertex_label.str().c_str());
  summary_spectrum_plot_pad->cd();
  summary_spectrum_plot_pad->RedrawAxis();
  summary_spectrum_canvas.cd();
  summary_spectrum_canvas.SaveAs((output_base + ".pdf").c_str());

  const auto draw_category_annotation = [&]() {
    TLatex label;
    label.SetNDC();
    label.SetTextAlign(13);
    label.SetTextSize(kTextSize);
    label.DrawLatex(kAnnotationX, 0.96, "#it{#bf{sPHENIX}} Internal");
    label.DrawLatex(kAnnotationX, 0.91, sample_label.c_str());
    label.DrawLatex(kAnnotationX, 0.86, "Denominator: all #pi^{0}-main anchors");
    label.DrawLatex(kAnnotationX, 0.81, anchor_label.str().c_str());
    label.DrawLatex(kAnnotationX, 0.76, min_cluster_energy_label.str().c_str());
    label.DrawLatex(kAnnotationX, 0.71, recovery_label.str().c_str());
    label.DrawLatex(kAnnotationX, 0.66, vertex_label.str().c_str());
  };

  TCanvas detailed_fraction_canvas("c_pythia_pi0_anchor_detailed_category_fractions", "Detailed Pythia pi0 anchor category fractions", kCanvasWidth, kCanvasHeight);
  detailed_fraction_canvas.SetCanvasSize(kCanvasWidth, kCanvasHeight);
  auto detailed_fraction_plot_pad = make_plot_pad("pythia_pi0_anchor_detailed_category_fractions_plot");
  detailed_fractions[0]->SetMinimum(0.0);
  detailed_fractions[0]->SetMaximum(1.05);
  detailed_fractions[0]->Draw("E1");
  for (std::size_t index = 1; index < detailed_fractions.size(); ++index) detailed_fractions[index]->Draw("E1 SAME");
  detailed_fraction_canvas.cd();
  TLegend detailed_fraction_legend(0.50, 0.60, 0.95, 0.97);
  detailed_fraction_legend.SetBorderSize(0);
  detailed_fraction_legend.SetFillStyle(0);
  detailed_fraction_legend.SetTextSize(0.020);
  detailed_fraction_legend.AddEntry(detailed_fractions[0].get(), "Separated", "lep");
  detailed_fraction_legend.AddEntry(detailed_fractions[1].get(), "Merged", "lep");
  detailed_fraction_legend.AddEntry(detailed_fractions[2].get(), "Single contaminated", "lep");
  detailed_fraction_legend.AddEntry(detailed_fractions[3].get(), "Missing: energy threshold", "lep");
  detailed_fraction_legend.AddEntry(detailed_fractions[4].get(), "Missing: displaced partner cluster", "lep");
  detailed_fraction_legend.AddEntry(detailed_fractions[5].get(), "Missing: acceptance", "lep");
  detailed_fraction_legend.AddEntry(detailed_fractions[6].get(), "Missing: no CEMC deposit", "lep");
  detailed_fraction_legend.AddEntry(detailed_fractions[7].get(), "Missing: unclustered deposit", "lep");
  detailed_fraction_legend.AddEntry(detailed_fractions[8].get(), "Missing: match incomplete", "lep");
  detailed_fraction_legend.AddEntry(detailed_fractions[9].get(), "Missing: other", "lep");
  detailed_fraction_legend.AddEntry(detailed_fractions[10].get(), "Other", "lep");
  detailed_fraction_legend.Draw();
  draw_category_annotation();
  detailed_fraction_plot_pad->cd();
  detailed_fraction_plot_pad->RedrawAxis();
  detailed_fraction_canvas.cd();
  detailed_fraction_canvas.SaveAs((output_base + "_category_fractions_detailed.pdf").c_str());

  TCanvas summary_fraction_canvas("c_pythia_pi0_anchor_summary_category_fractions", "Summary Pythia pi0 anchor category fractions", kCanvasWidth, kCanvasHeight);
  summary_fraction_canvas.SetCanvasSize(kCanvasWidth, kCanvasHeight);
  auto summary_fraction_plot_pad = make_plot_pad("pythia_pi0_anchor_summary_category_fractions_plot");
  summary_fractions[0]->SetMinimum(0.0);
  summary_fractions[0]->SetMaximum(1.05);
  summary_fractions[0]->Draw("E1");
  for (std::size_t index = 1; index < summary_fractions.size(); ++index) summary_fractions[index]->Draw("E1 SAME");
  summary_fraction_canvas.cd();
  TLegend summary_fraction_legend(0.55, 0.72, 0.94, 0.95);
  summary_fraction_legend.SetBorderSize(0);
  summary_fraction_legend.SetFillStyle(0);
  summary_fraction_legend.SetTextSize(kTextSize);
  summary_fraction_legend.AddEntry(summary_fractions[0].get(), "Separated", "lep");
  summary_fraction_legend.AddEntry(summary_fractions[1].get(), "Merged", "lep");
  summary_fraction_legend.AddEntry(summary_fractions[2].get(), "Single contaminated", "lep");
  summary_fraction_legend.AddEntry(summary_fractions[3].get(), "Missing", "lep");
  summary_fraction_legend.AddEntry(summary_fractions[4].get(), "Other", "lep");
  summary_fraction_legend.Draw();
  draw_category_annotation();
  summary_fraction_plot_pad->cd();
  summary_fraction_plot_pad->RedrawAxis();
  summary_fraction_canvas.cd();
  summary_fraction_canvas.SaveAs((output_base + "_category_fractions.pdf").c_str());

  TCanvas detailed_fraction_stack_canvas("c_pythia_pi0_anchor_detailed_category_fraction_stack", "Detailed Pythia pi0 anchor category fraction stack", kCanvasWidth, kCanvasHeight);
  detailed_fraction_stack_canvas.SetCanvasSize(kCanvasWidth, kCanvasHeight);
  auto detailed_fraction_stack_plot_pad = make_plot_pad("pythia_pi0_anchor_detailed_category_fraction_stack_plot");
  THStack detailed_fraction_stack("stack_pythia_pi0_anchor_detailed_category_fractions", "");
  for (auto& histogram : detailed_stacked_fractions) detailed_fraction_stack.Add(histogram.get());
  detailed_fraction_stack.SetMinimum(0.0);
  detailed_fraction_stack.SetMaximum(1.05);
  detailed_fraction_stack.Draw("HIST");
  detailed_fraction_stack.GetXaxis()->SetTitle("Anchor Cluster E_{T} [GeV]");
  detailed_fraction_stack.GetYaxis()->SetTitle("Fraction");
  style_plot_axes(detailed_fraction_stack.GetXaxis(), detailed_fraction_stack.GetYaxis());
  detailed_fraction_stack_canvas.cd();
  TLegend detailed_fraction_stack_legend(0.50, 0.60, 0.95, 0.97);
  detailed_fraction_stack_legend.SetBorderSize(0);
  detailed_fraction_stack_legend.SetFillStyle(0);
  detailed_fraction_stack_legend.SetTextSize(0.020);
  detailed_fraction_stack_legend.AddEntry(detailed_stacked_fractions[0].get(), "Separated", "f");
  detailed_fraction_stack_legend.AddEntry(detailed_stacked_fractions[1].get(), "Merged", "f");
  detailed_fraction_stack_legend.AddEntry(detailed_stacked_fractions[2].get(), "Single contaminated", "f");
  detailed_fraction_stack_legend.AddEntry(detailed_stacked_fractions[3].get(), "Missing: energy threshold", "f");
  detailed_fraction_stack_legend.AddEntry(detailed_stacked_fractions[4].get(), "Missing: displaced partner cluster", "f");
  detailed_fraction_stack_legend.AddEntry(detailed_stacked_fractions[5].get(), "Missing: acceptance", "f");
  detailed_fraction_stack_legend.AddEntry(detailed_stacked_fractions[6].get(), "Missing: no CEMC deposit", "f");
  detailed_fraction_stack_legend.AddEntry(detailed_stacked_fractions[7].get(), "Missing: unclustered deposit", "f");
  detailed_fraction_stack_legend.AddEntry(detailed_stacked_fractions[8].get(), "Missing: match incomplete", "f");
  detailed_fraction_stack_legend.AddEntry(detailed_stacked_fractions[9].get(), "Missing: other", "f");
  detailed_fraction_stack_legend.AddEntry(detailed_stacked_fractions[10].get(), "Other", "f");
  detailed_fraction_stack_legend.Draw();
  draw_category_annotation();
  detailed_fraction_stack_plot_pad->cd();
  detailed_fraction_stack_plot_pad->RedrawAxis();
  detailed_fraction_stack_canvas.cd();
  detailed_fraction_stack_canvas.SaveAs((output_base + "_category_fraction_stack_detailed.pdf").c_str());

  TCanvas summary_fraction_stack_canvas("c_pythia_pi0_anchor_summary_category_fraction_stack", "Summary Pythia pi0 anchor category fraction stack", kCanvasWidth, kCanvasHeight);
  summary_fraction_stack_canvas.SetCanvasSize(kCanvasWidth, kCanvasHeight);
  auto summary_fraction_stack_plot_pad = make_plot_pad("pythia_pi0_anchor_summary_category_fraction_stack_plot");
  THStack summary_fraction_stack("stack_pythia_pi0_anchor_summary_category_fractions", "");
  for (auto& histogram : summary_stacked_fractions) summary_fraction_stack.Add(histogram.get());
  summary_fraction_stack.SetMinimum(0.0);
  summary_fraction_stack.SetMaximum(1.05);
  summary_fraction_stack.Draw("HIST");
  summary_fraction_stack.GetXaxis()->SetTitle("Anchor Cluster E_{T} [GeV]");
  summary_fraction_stack.GetYaxis()->SetTitle("Fraction");
  style_plot_axes(summary_fraction_stack.GetXaxis(), summary_fraction_stack.GetYaxis());
  summary_fraction_stack_canvas.cd();
  TLegend summary_fraction_stack_legend(0.55, 0.72, 0.94, 0.95);
  summary_fraction_stack_legend.SetBorderSize(0);
  summary_fraction_stack_legend.SetFillStyle(0);
  summary_fraction_stack_legend.SetTextSize(kTextSize);
  summary_fraction_stack_legend.AddEntry(summary_stacked_fractions[0].get(), "Separated", "f");
  summary_fraction_stack_legend.AddEntry(summary_stacked_fractions[1].get(), "Merged", "f");
  summary_fraction_stack_legend.AddEntry(summary_stacked_fractions[2].get(), "Single contaminated", "f");
  summary_fraction_stack_legend.AddEntry(summary_stacked_fractions[3].get(), "Missing", "f");
  summary_fraction_stack_legend.AddEntry(summary_stacked_fractions[4].get(), "Other", "f");
  summary_fraction_stack_legend.Draw();
  draw_category_annotation();
  summary_fraction_stack_plot_pad->cd();
  summary_fraction_stack_plot_pad->RedrawAxis();
  summary_fraction_stack_canvas.cd();
  summary_fraction_stack_canvas.SaveAs((output_base + "_category_fraction_stack.pdf").c_str());

  TFile output((output_base + ".root").c_str(), "RECREATE");
  if (output.IsZombie()) {
    return 6;
  }
  for (std::size_t index = 0; index < raw.size(); ++index) {
    raw[index]->Write();
    density[index]->Write();
  }
  for (auto& histogram : detailed_fractions) histogram->Write();
  for (auto& histogram : summary_fractions) histogram->Write();

  int output_schema_version = 9;
  long long manifest_begin = partials.front().manifest_begin;
  long long manifest_end = partials.back().manifest_end;
  long long partial_file_count = static_cast<long long>(partials.size());
  long long input_file_count = manifest_end - manifest_begin;
  unsigned char contains_raw_histograms = 1U;
  unsigned char contains_bin_width_normalized_histograms = 1U;
  unsigned char contains_category_fractions = 1U;
  unsigned char contains_summary_category_fractions = 1U;
  TTree metadata("metadata", "Final Pythia pi0 anchor-cluster metadata");
  metadata.Branch("schema_version", &output_schema_version);
  metadata.Branch("manifest_path", &total.manifest_path);
  metadata.Branch("manifest_begin", &manifest_begin);
  metadata.Branch("manifest_end", &manifest_end);
  metadata.Branch("partial_file_count", &partial_file_count);
  metadata.Branch("input_file_count", &input_file_count);
  metadata.Branch("cluster_collection", &total.cluster_collection);
  metadata.Branch("tower_geom_node", &total.tower_geom_node);
  metadata.Branch("classification_unit", &total.classification_unit);
  metadata.Branch("pi0_selection", &total.pi0_selection);
  metadata.Branch("partner_selection", &total.partner_selection);
  metadata.Branch("topology_definition", &total.topology_definition);
  metadata.Branch("topology_priority", &total.topology_priority);
  metadata.Branch("missing_category_priority", &total.missing_category_priority);
  metadata.Branch("response_policy", &total.response_policy);
  metadata.Branch("photon_recovery_policy", &total.photon_recovery_policy);
  metadata.Branch("vertex_selection", &total.vertex_selection);
  metadata.Branch("signal_embedding_id", &total.signal_embedding_id);
  metadata.Branch("n_bins", &total.n_bins);
  metadata.Branch("et_max", &total.et_max);
  metadata.Branch("anchor_cluster_eta_max", &total.anchor_cluster_eta_max);
  metadata.Branch("partner_cluster_eta_max", &total.partner_cluster_eta_max);
  metadata.Branch("cemc_acceptance_eta_max", &total.cemc_acceptance_eta_max);
  metadata.Branch("pre_cemc_interaction_radius", &total.pre_cemc_interaction_radius);
  metadata.Branch("min_cluster_energy", &total.min_cluster_energy);
  metadata.Branch("dominant_fraction_min", &total.dominant_fraction_min);
  metadata.Branch("anchor_pi0_fraction_min", &total.anchor_pi0_fraction_min);
  metadata.Branch("min_energy_contribution_fraction", &total.min_energy_contribution_fraction);
  metadata.Branch("min_photon_energy_recovery", &total.min_photon_energy_recovery);
  metadata.Branch("min_direct_match_cluster_energy_coverage", &total.min_direct_match_cluster_energy_coverage);
  metadata.Branch("missing_diagnostic_max_delta_r", &total.missing_diagnostic_max_delta_r);
  metadata.Branch("enable_missing_diagnostics", &total.enable_missing_diagnostics);
  metadata.Branch("max_abs_vertex_z", &total.max_abs_vertex_z);
  metadata.Branch("pi0_truth_matching_algorithm_version", &total.matcher_version);
  metadata.Branch("pi0_topology_algorithm_version", &total.topology_version);
  metadata.Branch("contains_raw_histograms", &contains_raw_histograms);
  metadata.Branch("contains_bin_width_normalized_histograms", &contains_bin_width_normalized_histograms);
  metadata.Branch("contains_category_fractions", &contains_category_fractions);
  metadata.Branch("contains_summary_category_fractions", &contains_summary_category_fractions);
  metadata.Branch("events_processed", &total.events_processed);
  metadata.Branch("events_written", &total.events_written);
  metadata.Branch("events_invalid", &total.events_invalid);
  metadata.Branch("events_vertex_rejected", &total.events_vertex_rejected);
  metadata.Branch("cluster_considered_count", &total.cluster_considered);
  metadata.Branch("cluster_invalid_truth_count", &total.cluster_invalid_truth);
  metadata.Branch("prompt_cluster_count", &total.prompt_count);
  metadata.Branch("pi0_candidate_g4_decay_count", &total.candidate_g4);
  metadata.Branch("pi0_candidate_generator_decay_count", &total.candidate_generator);
  metadata.Branch("pi0_malformed_daughters_count", &total.malformed_daughters);
  metadata.Branch("anchor_cluster_count", &total.anchor_count);
  metadata.Branch("anchor_g4_decay_count", &total.anchor_g4);
  metadata.Branch("anchor_generator_decay_count", &total.anchor_generator);
  metadata.Branch("anchor_ambiguous_main_count", &total.ambiguous_main);
  metadata.Branch("energy_match_invalid_count", &total.energy_match_invalid);
  metadata.Branch("separated_count", &total.separated_count);
  metadata.Branch("merged_count", &total.merged_count);
  metadata.Branch("single_contaminated_count", &total.single_contaminated_count);
  metadata.Branch("missing_count", &total.missing_count);
  metadata.Branch("missing_energy_threshold_count", &total.missing_energy_threshold_count);
  metadata.Branch("missing_displaced_partner_cluster_count", &total.missing_displaced_partner_cluster_count);
  metadata.Branch("missing_acceptance_count", &total.missing_acceptance_count);
  metadata.Branch("missing_no_cemc_deposit_count", &total.missing_no_cemc_deposit_count);
  metadata.Branch("missing_unclustered_deposit_count", &total.missing_unclustered_deposit_count);
  metadata.Branch("missing_match_incomplete_count", &total.missing_match_incomplete_count);
  metadata.Branch("missing_other_count", &total.missing_other_count);
  metadata.Branch("other_count", &total.other_count);
  metadata.Fill();
  metadata.Write();
  output.Close();
  if (output.TestBit(TFile::kWriteError))
  {
    return 6;
  }

  std::cout
      << "FinalizePythiaPi0AnchorClusterSpectra - partials/files/events"
      << "/vertex-rejected/anchor/separated/merged/single-contaminated/missing(energy/displaced/acceptance/no-CEMC/unclustered/match-incomplete/other)/other = "
      << partial_file_count << "/" << input_file_count << "/"
      << total.events_processed << "/" << total.events_vertex_rejected << "/"
      << total.anchor_count << "/"
      << total.separated_count << "/" << total.merged_count << "/" << total.single_contaminated_count << "/"
      << total.missing_count << "(" << total.missing_energy_threshold_count << "/" << total.missing_displaced_partner_cluster_count << "/"
      << total.missing_acceptance_count << "/" << total.missing_no_cemc_deposit_count << "/" << total.missing_unclustered_deposit_count << "/"
      << total.missing_match_incomplete_count << "/" << total.missing_other_count << ")/"
      << total.other_count << std::endl;
  return 0;
}
