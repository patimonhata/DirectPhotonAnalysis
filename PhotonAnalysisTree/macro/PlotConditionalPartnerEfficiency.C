#include "Utilities/sPhenixStyle.C"

#include <TAttMarker.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TH1D.h>
#include <THStack.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace {
constexpr double pi = 3.14159265358979323846;
const std::vector<double> truth_pt_edges = {3.0, 5.0, 7.0, 9.0,
                                            11.0, 13.0, 15.0};
constexpr std::size_t n_truth_pt_bins = 6U;
const std::array<std::string, n_truth_pt_bins> truth_pt_labels = {
    "3 #leq p_{T}^{truth} < 5 GeV", "5 #leq p_{T}^{truth} < 7 GeV",
    "7 #leq p_{T}^{truth} < 9 GeV",
    "9 #leq p_{T}^{truth} < 11 GeV", "11 #leq p_{T}^{truth} < 13 GeV",
    "13 #leq p_{T}^{truth} #leq 15 GeV"};
const std::array<std::string, n_truth_pt_bins> reco_et_labels = {
    "3 #leq E_{T}^{cluster} < 5 GeV",
    "5 #leq E_{T}^{cluster} < 7 GeV",
    "7 #leq E_{T}^{cluster} < 9 GeV",
    "9 #leq E_{T}^{cluster} < 11 GeV",
    "11 #leq E_{T}^{cluster} < 13 GeV",
    "13 #leq E_{T}^{cluster} #leq 15 GeV"};
constexpr int canvas_columns = 3;
constexpr int canvas_rows = 3;
constexpr int canvas_width = 1350;
constexpr int canvas_height = 1350;
constexpr int conditions_pad = static_cast<int>(n_truth_pt_bins) + 1;
constexpr int legend_pad = conditions_pad + 1;

struct CollectionBranches {
  std::vector<double> *cluster_e = nullptr;
  std::vector<double> *cluster_et = nullptr;
  std::vector<double> *cluster_eta = nullptr;
  std::vector<double> *cluster_phi = nullptr;
  std::vector<unsigned int> *pair_i = nullptr;
  std::vector<unsigned int> *pair_j = nullptr;
  std::vector<double> *pair_mass = nullptr;
};

struct MatchResult {
  bool valid = false;
  std::size_t cluster0 = 0;
  std::size_t cluster1 = 0;
  double delta_r0 = std::numeric_limits<double>::infinity();
  double delta_r1 = std::numeric_limits<double>::infinity();
};

struct StageHistograms {
  std::unique_ptr<TH1D> reference;
  std::unique_ptr<TH1D> denominator;
  std::unique_ptr<TH1D> matched;
  std::unique_ptr<TH1D> removed;
  std::unique_ptr<TH1D> efficiency_matched;
  std::unique_ptr<TH1D> efficiency_removed;
  Long64_t n_reference = 0;
  Long64_t n_denominator = 0;
  Long64_t n_matched = 0;
  Long64_t n_removed = 0;
  Long64_t n_malformed = 0;
};

enum class EventComponent : std::size_t {
  correct_pair = 0,
  merged_candidate = 1,
  individual_anchor = 2,
  other_anchor = 3,
  count = 4
};

constexpr std::size_t n_event_components =
    static_cast<std::size_t>(EventComponent::count);
const std::array<std::string, n_event_components> event_component_names = {
    "correct_pair", "merged_candidate", "individual_anchor_no_partner",
    "other_anchor"};
const std::array<std::string, n_event_components> event_component_labels = {
    "Correct pair", "Merged candidate (no correct pair)",
    "Individual-#gamma anchor, no partner", "Other selected anchor"};

struct EventComponentHistograms {
  std::unique_ptr<TH1D> reference;
  std::unique_ptr<TH1D> in_acceptance;
  std::unique_ptr<TH1D> selected;
  std::array<std::unique_ptr<TH1D>, n_event_components> component;
  std::array<std::unique_ptr<TH1D>, n_event_components> fraction;
  Long64_t n_reference = 0;
  Long64_t n_in_acceptance = 0;
  Long64_t n_selected = 0;
  std::array<Long64_t, n_event_components> n_component = {};
  Long64_t n_malformed = 0;
};

// Reconstruction topologies for every generated pi0 in the requested truth
// acceptance. Unlike EventComponentHistograms, these histograms do not first
// require a selected reconstructed anchor cluster.
struct TruthPi0ComponentHistograms {
  std::unique_ptr<TH1D> reference;
  std::unique_ptr<TH1D> in_truth_eta;
  std::array<std::unique_ptr<TH1D>, n_event_components> component;
  std::array<std::unique_ptr<TH1D>, n_event_components> fraction;
  Long64_t n_reference = 0;
  Long64_t n_in_truth_eta = 0;
  std::array<Long64_t, n_event_components> n_component = {};
  Long64_t n_malformed = 0;
};

const std::array<std::string, n_event_components> truth_pi0_component_names = {
    "separated_pair_candidate", "merged_candidate",
    "individual_cluster_no_pair", "no_matched_topology"};
const std::array<std::string, n_event_components> truth_pi0_component_labels = {
    "Two separated clusters", "Merged-cluster candidate",
    "Individual-#gamma cluster, no pair", "No matched topology"};
constexpr std::size_t n_cut_stages = 4U;

struct CutStageHistograms {
  std::array<std::unique_ptr<TH1D>, n_cut_stages> stage;
  std::array<Long64_t, n_cut_stages> n_stage = {};
  Long64_t n_malformed = 0;
};

enum class ClusterSelection : std::size_t {
  all = 0,
  min_energy = 1,
  eta = 2,
  min_energy_eta = 3,
  count = 4
};

constexpr std::size_t n_cluster_selections =
    static_cast<std::size_t>(ClusterSelection::count);
constexpr std::size_t n_cluster_retentions = n_cluster_selections - 1U;
const std::array<std::string, n_cluster_selections> cluster_selection_names = {
    "all", "min_energy", "eta", "min_energy_eta"};
const std::array<std::string, n_cluster_retentions> cluster_retention_names = {
    "min_energy", "eta", "min_energy_eta"};

constexpr std::size_t selection_index(const ClusterSelection selection) {
  return static_cast<std::size_t>(selection);
}

struct RecoEtHistograms {
  std::array<std::unique_ptr<TH1D>, n_truth_pt_bins> truth_pt_contribution;
  std::array<std::unique_ptr<TH1D>, n_cluster_selections> selection;
  std::array<std::unique_ptr<TH1D>, n_cluster_retentions> retention;
  std::array<Long64_t, n_truth_pt_bins> n_truth_pt_contribution = {};
  std::array<Long64_t, n_cluster_selections> n_selection = {};
};

enum class CoreCondition : std::size_t {
  reference = 0,
  selected = 1,
  matched = 2,
  removed = 3,
  count = 4
};

constexpr std::size_t n_core_conditions =
    static_cast<std::size_t>(CoreCondition::count);
constexpr std::size_t condition_index(const CoreCondition condition) {
  return static_cast<std::size_t>(condition);
}

// Color identifies the core selection role throughout the plot family.
const std::array<int, n_core_conditions> core_condition_colors = {
    kBlack, kGreen + 2, kBlue + 1, kRed + 1};

// Marker families distinguish the counted object and whether the quantity is a
// count or a derived ratio. Exact event-level quantities retain the same marker
// between the count, component, and cut-stage views.
const std::array<int, n_core_conditions> event_count_markers = {
    kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown};
const std::array<int, n_core_conditions> cluster_count_markers = {
    kFullDiamond, kFullCross, kFullStar, kFullCrossX};
const std::array<int, n_core_conditions> event_ratio_markers = {
    kOpenCircle, kOpenSquare, kOpenTriangleUp, kOpenTriangleDown};
const std::array<int, n_core_conditions> cluster_ratio_markers = {
    kOpenDiamond, kOpenSquareDiagonal, kOpenStar, kOpenCross};
const std::array<int, n_event_components> event_component_colors = {
    core_condition_colors[condition_index(CoreCondition::matched)],
    kMagenta - 3, kOrange + 7, kGray + 1};
const std::array<int, n_event_components> event_component_count_markers = {
    event_count_markers[condition_index(CoreCondition::matched)],
    kFullDiamond, kFullStar, kFullCross};
const std::array<int, n_event_components> event_component_ratio_markers = {
    event_ratio_markers[condition_index(CoreCondition::matched)],
    kOpenDiamond, kOpenStar, kOpenCross};


void set_point_style(TH1D &histogram, const int color, const int marker,
                     const int line_style = 1,
                     const double marker_size = 0.65) {
  histogram.SetLineColor(color);
  histogram.SetMarkerColor(color);
  histogram.SetMarkerStyle(marker);
  histogram.SetMarkerSize(marker_size);
  histogram.SetLineStyle(line_style);
  histogram.SetLineWidth(2);
}

double wrap_delta_phi(double value) {
  while (value > pi) {
    value -= 2.0 * pi;
  }
  while (value <= -pi) {
    value += 2.0 * pi;
  }
  return value;
}

double delta_r(const double eta0, const double phi0, const double eta1,
               const double phi1) {
  return std::hypot(eta0 - eta1, wrap_delta_phi(phi0 - phi1));
}

bool make_output_directory(const std::string &output_base) {
  const std::size_t slash = output_base.find_last_of('/');
  if (slash == std::string::npos) {
    return true;
  }
  const std::string directory = output_base.substr(0, slash);
  if (directory.empty() || !gSystem->AccessPathName(directory.c_str())) {
    return true;
  }
  if (gSystem->mkdir(directory.c_str(), true) != 0) {
    std::cerr << "Failed to create output directory " << directory << std::endl;
    return false;
  }
  return true;
}

std::string collection_output_base(const std::string &output_base,
                                   const std::string &collection) {
  const std::size_t slash = output_base.find_last_of("/");
  const std::string directory =
      slash == std::string::npos ? "" : output_base.substr(0, slash + 1U);
  const std::string stem =
      slash == std::string::npos ? output_base : output_base.substr(slash + 1U);
  return directory + collection + "/" + stem + "_" + collection;
}

std::size_t find_truth_pt_bin(const double truth_pt) {
  if (!std::isfinite(truth_pt) || truth_pt < truth_pt_edges.front() ||
      truth_pt > truth_pt_edges.back()) {
    return truth_pt_labels.size();
  }
  if (truth_pt == truth_pt_edges.back()) {
    return truth_pt_labels.size() - 1U;
  }
  const auto upper =
      std::upper_bound(truth_pt_edges.begin(), truth_pt_edges.end(), truth_pt);
  return static_cast<std::size_t>(std::distance(truth_pt_edges.begin(), upper) -
                                  1);
}

std::size_t find_reco_et_bin(const double cluster_et) {
  if (!std::isfinite(cluster_et) || cluster_et < truth_pt_edges.front() ||
      cluster_et > truth_pt_edges.back()) {
    return reco_et_labels.size();
  }
  if (cluster_et == truth_pt_edges.back()) {
    return reco_et_labels.size() - 1U;
  }
  const auto upper =
      std::upper_bound(truth_pt_edges.begin(), truth_pt_edges.end(), cluster_et);
  return static_cast<std::size_t>(
      std::distance(truth_pt_edges.begin(), upper) - 1);
}

bool valid_cluster(const CollectionBranches &branches,
                   const std::size_t index) {
  return branches.cluster_e && branches.cluster_et && branches.cluster_eta &&
         branches.cluster_phi && index < branches.cluster_e->size() &&
         index < branches.cluster_et->size() &&
         index < branches.cluster_eta->size() &&
         index < branches.cluster_phi->size() &&
         std::isfinite(branches.cluster_e->at(index)) &&
         std::isfinite(branches.cluster_et->at(index)) &&
         std::isfinite(branches.cluster_eta->at(index)) &&
         std::isfinite(branches.cluster_phi->at(index));
}
bool passes_min_cluster_energy(const CollectionBranches &branches,
                               const std::size_t index,
                               const double min_cluster_energy) {
  return valid_cluster(branches, index) &&
         branches.cluster_e->at(index) >= min_cluster_energy;
}

bool selected_anchor(const CollectionBranches &branches,
                     const std::size_t index, const double eta_max,
                     const double et_min, const double et_max) {
  return valid_cluster(branches, index) &&
         std::abs(branches.cluster_eta->at(index)) < eta_max &&
         branches.cluster_et->at(index) > et_min &&
         branches.cluster_et->at(index) < et_max;
}

MatchResult match_event(const std::vector<double> &truth_eta,
                        const std::vector<double> &truth_phi,
                        const CollectionBranches &branches,
                        const double min_cluster_energy) {
  MatchResult result;
  if (truth_eta.size() != 2U || truth_phi.size() != 2U ||
      !branches.cluster_et) {
    return result;
  }
  double best_cost = std::numeric_limits<double>::infinity();
  for (std::size_t first = 0; first < branches.cluster_et->size(); ++first) {
    if (!passes_min_cluster_energy(branches, first, min_cluster_energy)) {
      continue;
    }
    const double first_delta_r =
        delta_r(truth_eta[0], truth_phi[0], branches.cluster_eta->at(first),
                branches.cluster_phi->at(first));
    for (std::size_t second = 0; second < branches.cluster_et->size();
         ++second) {
      if (first == second ||
          !passes_min_cluster_energy(branches, second, min_cluster_energy)) {
        continue;
      }
      const double second_delta_r =
          delta_r(truth_eta[1], truth_phi[1], branches.cluster_eta->at(second),
                  branches.cluster_phi->at(second));
      const double cost =
          first_delta_r * first_delta_r + second_delta_r * second_delta_r;
      if (cost < best_cost) {
        best_cost = cost;
        result.valid = true;
        result.cluster0 = first;
        result.cluster1 = second;
        result.delta_r0 = first_delta_r;
        result.delta_r1 = second_delta_r;
      }
    }
  }
  return result;
}

MatchResult match_anchor(const std::size_t anchor,
                         const std::vector<double> &truth_eta,
                         const std::vector<double> &truth_phi,
                         const CollectionBranches &branches,
                         const double min_cluster_energy) {
  MatchResult result;
  if (truth_eta.size() != 2U || truth_phi.size() != 2U ||
      !passes_min_cluster_energy(branches, anchor, min_cluster_energy)) {
    return result;
  }
  double best_cost = std::numeric_limits<double>::infinity();
  for (std::size_t partner = 0; partner < branches.cluster_et->size();
       ++partner) {
    if (partner == anchor ||
        !passes_min_cluster_energy(branches, partner, min_cluster_energy)) {
      continue;
    }
    for (std::size_t anchor_gamma = 0; anchor_gamma < 2U; ++anchor_gamma) {
      const std::size_t partner_gamma = 1U - anchor_gamma;
      const double anchor_delta_r = delta_r(
          truth_eta[anchor_gamma], truth_phi[anchor_gamma],
          branches.cluster_eta->at(anchor), branches.cluster_phi->at(anchor));
      const double partner_delta_r = delta_r(
          truth_eta[partner_gamma], truth_phi[partner_gamma],
          branches.cluster_eta->at(partner), branches.cluster_phi->at(partner));
      const double cost =
          anchor_delta_r * anchor_delta_r + partner_delta_r * partner_delta_r;
      if (cost < best_cost) {
        best_cost = cost;
        result.valid = true;
        result.cluster0 = anchor;
        result.cluster1 = partner;
        result.delta_r0 = anchor_delta_r;
        result.delta_r1 = partner_delta_r;
      }
    }
  }
  return result;
}

bool find_pair_mass(const CollectionBranches &branches,
                    const std::size_t first_cluster,
                    const std::size_t second_cluster, double &mass) {
  if (!branches.pair_i || !branches.pair_j || !branches.pair_mass ||
      branches.pair_i->size() != branches.pair_j->size() ||
      branches.pair_i->size() != branches.pair_mass->size()) {
    return false;
  }
  for (std::size_t pair = 0; pair < branches.pair_i->size(); ++pair) {
    const std::size_t i = branches.pair_i->at(pair);
    const std::size_t j = branches.pair_j->at(pair);
    if ((i == first_cluster && j == second_cluster) ||
        (i == second_cluster && j == first_cluster)) {
      mass = branches.pair_mass->at(pair);
      return std::isfinite(mass);
    }
  }
  return false;
}

StageHistograms make_histograms(const std::string &prefix,
                                const int asymmetry_nbins) {
  StageHistograms histograms;
  const auto make = [&](const std::string &suffix) {
    auto histogram = std::make_unique<TH1D>(
        ("h_" + prefix + "_" + suffix).c_str(), "", asymmetry_nbins, 0.0, 1.0);
    histogram->SetDirectory(nullptr);
    histogram->Sumw2();
    histogram->SetStats(false);
    return histogram;
  };
  histograms.reference = make("reference");
  histograms.denominator = make("denominator");
  histograms.matched = make("matched");
  histograms.removed = make("removed");
  return histograms;
}

std::vector<StageHistograms>
make_truth_pt_histograms(const std::string &prefix, const int asymmetry_nbins) {
  std::vector<StageHistograms> histograms;
  histograms.reserve(truth_pt_labels.size());
  for (std::size_t bin = 0; bin < truth_pt_labels.size(); ++bin) {
    histograms.push_back(make_histograms(
        prefix + "_truth_pt_" + std::to_string(bin), asymmetry_nbins));
  }
  return histograms;
}

EventComponentHistograms make_component_histograms(const std::string &prefix,
                                                   const int asymmetry_nbins) {
  EventComponentHistograms histograms;
  const auto make = [&](const std::string &suffix) {
    auto histogram = std::make_unique<TH1D>(
        ("h_" + prefix + "_" + suffix).c_str(), "", asymmetry_nbins, 0.0, 1.0);
    histogram->SetDirectory(nullptr);
    histogram->Sumw2();
    histogram->SetStats(false);
    return histogram;
  };
  histograms.reference = make("reference");
  histograms.in_acceptance = make("in_acceptance");
  histograms.selected = make("selected");
  for (std::size_t component = 0; component < n_event_components; ++component) {
    histograms.component[component] = make(event_component_names[component]);
  }
  return histograms;
}

std::vector<EventComponentHistograms>
make_truth_pt_component_histograms(const std::string &prefix,
                                   const int asymmetry_nbins) {
  std::vector<EventComponentHistograms> histograms;
  histograms.reserve(truth_pt_labels.size());
  for (std::size_t bin = 0; bin < truth_pt_labels.size(); ++bin) {
    histograms.push_back(make_component_histograms(
        prefix + "_truth_pt_" + std::to_string(bin), asymmetry_nbins));
  }
  return histograms;
}

TruthPi0ComponentHistograms
make_truth_pi0_component_histograms(const std::string &prefix,
                                    const int asymmetry_nbins) {
  TruthPi0ComponentHistograms histograms;
  const auto make = [&](const std::string &suffix) {
    auto histogram = std::make_unique<TH1D>(
        ("h_" + prefix + "_" + suffix).c_str(), "", asymmetry_nbins, 0.0, 1.0);
    histogram->SetDirectory(nullptr);
    histogram->Sumw2();
    histogram->SetStats(false);
    return histogram;
  };
  histograms.reference = make("reference");
  histograms.in_truth_eta = make("in_truth_eta");
  for (std::size_t component = 0; component < n_event_components; ++component) {
    histograms.component[component] =
        make(truth_pi0_component_names[component]);
  }
  return histograms;
}

std::vector<TruthPi0ComponentHistograms>
make_truth_pt_truth_pi0_component_histograms(const std::string &prefix,
                                             const int asymmetry_nbins) {
  std::vector<TruthPi0ComponentHistograms> histograms;
  histograms.reserve(truth_pt_labels.size());
  for (std::size_t bin = 0; bin < truth_pt_labels.size(); ++bin) {
    histograms.push_back(make_truth_pi0_component_histograms(
        prefix + "_truth_pt_" + std::to_string(bin), asymmetry_nbins));
  }
  return histograms;
}
CutStageHistograms
make_cut_stage_histograms(const std::string &prefix,
                          const std::array<std::string, n_cut_stages> &suffixes,
                          const int asymmetry_nbins) {
  CutStageHistograms histograms;
  for (std::size_t stage = 0; stage < n_cut_stages; ++stage) {
    histograms.stage[stage] =
        std::make_unique<TH1D>(("h_" + prefix + "_" + suffixes[stage]).c_str(),
                               "", asymmetry_nbins, 0.0, 1.0);
    histograms.stage[stage]->SetDirectory(nullptr);
    histograms.stage[stage]->Sumw2();
    histograms.stage[stage]->SetStats(false);
  }
  return histograms;
}

std::vector<CutStageHistograms> make_truth_pt_cut_stage_histograms(
    const std::string &prefix,
    const std::array<std::string, n_cut_stages> &suffixes,
    const int asymmetry_nbins) {
  std::vector<CutStageHistograms> histograms;
  histograms.reserve(truth_pt_labels.size());
  for (std::size_t bin = 0; bin < truth_pt_labels.size(); ++bin) {
    histograms.push_back(
        make_cut_stage_histograms(prefix + "_truth_pt_" + std::to_string(bin),
                                  suffixes, asymmetry_nbins));
  }
  return histograms;
}

bool valid_collection_shape(const CollectionBranches &branches) {
  return branches.cluster_e && branches.cluster_et && branches.cluster_eta &&
         branches.cluster_phi &&
         branches.cluster_e->size() == branches.cluster_et->size() &&
         branches.cluster_e->size() == branches.cluster_eta->size() &&
         branches.cluster_e->size() == branches.cluster_phi->size();
}

RecoEtHistograms make_reco_et_histograms(const std::string &prefix,
                                         const int asymmetry_nbins) {
  RecoEtHistograms histograms;
  const auto make = [&](const std::string &suffix) {
    auto histogram = std::make_unique<TH1D>(
        ("h_" + prefix + "_" + suffix).c_str(), "", asymmetry_nbins, 0.0,
        1.0);
    histogram->SetDirectory(nullptr);
    histogram->Sumw2();
    histogram->SetStats(false);
    return histogram;
  };
  for (std::size_t truth_bin = 0; truth_bin < n_truth_pt_bins; ++truth_bin) {
    histograms.truth_pt_contribution[truth_bin] =
        make("truth_pt_" + std::to_string(truth_bin) + "_selected");
  }
  for (std::size_t selection = 0; selection < n_cluster_selections;
       ++selection) {
    histograms.selection[selection] =
        make("selection_" + cluster_selection_names[selection]);
  }
  for (std::size_t retention = 0; retention < n_cluster_retentions;
       ++retention) {
    histograms.retention[retention] =
        make("retention_" + cluster_retention_names[retention]);
  }
  return histograms;
}

std::vector<RecoEtHistograms>
make_reco_et_histogram_set(const std::string &prefix,
                           const int asymmetry_nbins) {
  std::vector<RecoEtHistograms> histograms;
  histograms.reserve(reco_et_labels.size());
  for (std::size_t reco_bin = 0; reco_bin < reco_et_labels.size(); ++reco_bin) {
    histograms.push_back(make_reco_et_histograms(
        prefix + "_reco_et_" + std::to_string(reco_bin), asymmetry_nbins));
  }
  return histograms;
}

bool fill_reco_et_histograms(std::vector<RecoEtHistograms> &histograms,
                             const std::size_t truth_pt_bin,
                             const double truth_alpha,
                             const CollectionBranches &branches,
                             const double min_cluster_energy,
                             const double eta_max) {
  if (!valid_collection_shape(branches)) {
    return false;
  }
  for (std::size_t cluster = 0; cluster < branches.cluster_et->size();
       ++cluster) {
    if (!valid_cluster(branches, cluster)) {
      continue;
    }
    const std::size_t reco_bin =
        find_reco_et_bin(branches.cluster_et->at(cluster));
    if (reco_bin >= histograms.size()) {
      continue;
    }
    RecoEtHistograms &current = histograms[reco_bin];
    const bool passes_energy =
        branches.cluster_e->at(cluster) >= min_cluster_energy;
    const bool passes_eta =
        std::abs(branches.cluster_eta->at(cluster)) < eta_max;
    const std::array<bool, n_cluster_selections> passed = {
        true, passes_energy, passes_eta, passes_energy && passes_eta};
    for (std::size_t selection = 0; selection < n_cluster_selections;
         ++selection) {
      if (passed[selection]) {
        current.selection[selection]->Fill(truth_alpha);
        ++current.n_selection[selection];
      }
    }
    if (passes_energy && passes_eta) {
      current.truth_pt_contribution[truth_pt_bin]->Fill(truth_alpha);
      ++current.n_truth_pt_contribution[truth_pt_bin];
    }
  }
  return true;
}

void make_reco_et_retentions(std::vector<RecoEtHistograms> &histograms) {
  for (RecoEtHistograms &current : histograms) {
    for (std::size_t retention = 0; retention < n_cluster_retentions;
         ++retention) {
      current.retention[retention]->Divide(
          current.selection[retention + 1U].get(),
          current.selection[selection_index(ClusterSelection::all)].get(), 1.0,
          1.0, "B");
    }
  }
}

void fill_event_selection_stages(CutStageHistograms &histograms,
                                 const double truth_alpha,
                                 const CollectionBranches &branches,
                                 const double min_cluster_energy,
                                 const double anchor_eta_max,
                                 const double anchor_et_min,
                                 const double anchor_et_max) {
  ++histograms.n_stage[0];
  histograms.stage[0]->Fill(truth_alpha);
  if (!valid_collection_shape(branches)) {
    ++histograms.n_malformed;
    return;
  }

  bool has_baseline_cluster = false;
  bool has_central_cluster = false;
  bool has_selected_cluster = false;
  for (std::size_t cluster = 0; cluster < branches.cluster_e->size();
       ++cluster) {
    if (!passes_min_cluster_energy(branches, cluster, min_cluster_energy)) {
      continue;
    }
    has_baseline_cluster = true;
    if (!(std::abs(branches.cluster_eta->at(cluster)) < anchor_eta_max)) {
      continue;
    }
    has_central_cluster = true;
    has_selected_cluster |= selected_anchor(branches, cluster, anchor_eta_max,
                                            anchor_et_min, anchor_et_max);
  }
  const std::array<bool, n_cut_stages> passed = {
      true, has_baseline_cluster, has_central_cluster, has_selected_cluster};
  for (std::size_t stage = 1; stage < n_cut_stages; ++stage) {
    if (passed[stage]) {
      ++histograms.n_stage[stage];
      histograms.stage[stage]->Fill(truth_alpha);
    }
  }
}

void fill_pair_selection_stages(
    CutStageHistograms &histograms, const double truth_alpha,
    const std::vector<double> &truth_eta, const std::vector<double> &truth_phi,
    const CollectionBranches &branches, const double min_cluster_energy,
    const double anchor_eta_max, const double anchor_et_min,
    const double anchor_et_max, const double delta_r_cut) {
  ++histograms.n_stage[0];
  histograms.stage[0]->Fill(truth_alpha);
  if (!valid_collection_shape(branches)) {
    ++histograms.n_malformed;
    return;
  }

  const MatchResult match =
      match_event(truth_eta, truth_phi, branches, min_cluster_energy);
  if (!match.valid ||
      !(std::max(match.delta_r0, match.delta_r1) < delta_r_cut)) {
    return;
  }
  ++histograms.n_stage[1];
  histograms.stage[1]->Fill(truth_alpha);

  const bool has_central_endpoint =
      std::abs(branches.cluster_eta->at(match.cluster0)) < anchor_eta_max ||
      std::abs(branches.cluster_eta->at(match.cluster1)) < anchor_eta_max;
  if (!has_central_endpoint) {
    return;
  }
  ++histograms.n_stage[2];
  histograms.stage[2]->Fill(truth_alpha);

  const bool has_selected_endpoint =
      selected_anchor(branches, match.cluster0, anchor_eta_max, anchor_et_min,
                      anchor_et_max) ||
      selected_anchor(branches, match.cluster1, anchor_eta_max, anchor_et_min,
                      anchor_et_max);
  if (has_selected_endpoint) {
    ++histograms.n_stage[3];
    histograms.stage[3]->Fill(truth_alpha);
  }
}

void fill_event_component_reference(EventComponentHistograms &histograms,
                                    const double truth_alpha) {
  ++histograms.n_reference;
  histograms.reference->Fill(truth_alpha);
}

void fill_event_components(
    EventComponentHistograms &histograms, const double truth_alpha,
    const double truth_pt, const std::vector<double> &truth_daughter_pt,
    const std::vector<double> &truth_eta, const std::vector<double> &truth_phi,
    const CollectionBranches &branches, const double anchor_eta_max,
    const double anchor_et_min, const double anchor_et_max,
    const double min_cluster_energy, const double delta_r_cut,
    const double merged_delta_r_cut, const double merged_response_min,
    const double merged_response_max, const double individual_response_min,
    const double individual_response_max) {
  ++histograms.n_in_acceptance;
  histograms.in_acceptance->Fill(truth_alpha);
  if (!valid_collection_shape(branches) || truth_daughter_pt.size() != 2U ||
      truth_eta.size() != 2U || truth_phi.size() != 2U || !(truth_pt > 0.0) ||
      !(truth_daughter_pt[0] > 0.0) || !(truth_daughter_pt[1] > 0.0)) {
    ++histograms.n_malformed;
    return;
  }

  std::vector<std::size_t> selected;
  for (std::size_t cluster = 0; cluster < branches.cluster_et->size();
       ++cluster) {
    if (selected_anchor(branches, cluster, anchor_eta_max, anchor_et_min,
                        anchor_et_max)) {
      selected.push_back(cluster);
    }
  }
  if (selected.empty()) {
    return;
  }
  ++histograms.n_selected;
  histograms.selected->Fill(truth_alpha);

  EventComponent classification = EventComponent::other_anchor;
  const MatchResult pair_match =
      match_event(truth_eta, truth_phi, branches, min_cluster_energy);
  const bool correct_pair =
      pair_match.valid &&
      std::max(pair_match.delta_r0, pair_match.delta_r1) < delta_r_cut &&
      (selected_anchor(branches, pair_match.cluster0, anchor_eta_max,
                       anchor_et_min, anchor_et_max) ||
       selected_anchor(branches, pair_match.cluster1, anchor_eta_max,
                       anchor_et_min, anchor_et_max));

  if (correct_pair) {
    classification = EventComponent::correct_pair;
  } else {
    bool has_merged_candidate = false;
    for (const std::size_t anchor : selected) {
      const double delta_r0 =
          delta_r(truth_eta[0], truth_phi[0], branches.cluster_eta->at(anchor),
                  branches.cluster_phi->at(anchor));
      const double delta_r1 =
          delta_r(truth_eta[1], truth_phi[1], branches.cluster_eta->at(anchor),
                  branches.cluster_phi->at(anchor));
      const double response = branches.cluster_et->at(anchor) / truth_pt;
      has_merged_candidate |=
          std::max(delta_r0, delta_r1) < merged_delta_r_cut &&
          response >= merged_response_min && response <= merged_response_max;
    }

    if (has_merged_candidate) {
      classification = EventComponent::merged_candidate;
    } else {
      bool has_individual_anchor = false;
      for (const std::size_t anchor : selected) {
        for (std::size_t gamma = 0; gamma < 2U; ++gamma) {
          const double anchor_delta_r =
              delta_r(truth_eta[gamma], truth_phi[gamma],
                      branches.cluster_eta->at(anchor),
                      branches.cluster_phi->at(anchor));
          const double response =
              branches.cluster_et->at(anchor) / truth_daughter_pt[gamma];
          has_individual_anchor |= anchor_delta_r < delta_r_cut &&
                                   response >= individual_response_min &&
                                   response <= individual_response_max;
        }
      }
      if (has_individual_anchor) {
        classification = EventComponent::individual_anchor;
      }
    }
  }

  const std::size_t component = static_cast<std::size_t>(classification);
  ++histograms.n_component[component];
  histograms.component[component]->Fill(truth_alpha);
}

void fill_truth_pi0_components(
    TruthPi0ComponentHistograms &histograms, const double truth_alpha,
    const bool passes_truth_eta, const double truth_pt,
    const std::vector<double> &truth_daughter_pt,
    const std::vector<double> &truth_eta, const std::vector<double> &truth_phi,
    const CollectionBranches &branches, const double min_cluster_energy,
    const double delta_r_cut, const double merged_delta_r_cut,
    const double merged_response_min, const double merged_response_max,
    const double individual_response_min,
    const double individual_response_max) {
  ++histograms.n_reference;
  histograms.reference->Fill(truth_alpha);

  if (!passes_truth_eta) {
    return;
  }
  ++histograms.n_in_truth_eta;
  histograms.in_truth_eta->Fill(truth_alpha);

  EventComponent classification = EventComponent::other_anchor;
  if (!valid_collection_shape(branches) || truth_daughter_pt.size() != 2U ||
      truth_eta.size() != 2U || truth_phi.size() != 2U || !(truth_pt > 0.0) ||
      !(truth_daughter_pt[0] > 0.0) || !(truth_daughter_pt[1] > 0.0)) {
    ++histograms.n_malformed;
  } else {
    const MatchResult pair_match =
        match_event(truth_eta, truth_phi, branches, min_cluster_energy);
    const bool separated_pair =
        pair_match.valid &&
        std::max(pair_match.delta_r0, pair_match.delta_r1) < delta_r_cut;

    if (separated_pair) {
      classification = EventComponent::correct_pair;
    } else {
      bool has_merged_candidate = false;
      for (std::size_t cluster = 0; cluster < branches.cluster_et->size();
           ++cluster) {
        if (!valid_cluster(branches, cluster)) {
          continue;
        }
        const double delta_r0 = delta_r(truth_eta[0], truth_phi[0],
                                        branches.cluster_eta->at(cluster),
                                        branches.cluster_phi->at(cluster));
        const double delta_r1 = delta_r(truth_eta[1], truth_phi[1],
                                        branches.cluster_eta->at(cluster),
                                        branches.cluster_phi->at(cluster));
        const double response = branches.cluster_et->at(cluster) / truth_pt;
        has_merged_candidate |=
            std::max(delta_r0, delta_r1) < merged_delta_r_cut &&
            response >= merged_response_min && response <= merged_response_max;
      }

      if (has_merged_candidate) {
        classification = EventComponent::merged_candidate;
      } else {
        bool has_individual_cluster = false;
        for (std::size_t cluster = 0; cluster < branches.cluster_et->size();
             ++cluster) {
          if (!valid_cluster(branches, cluster)) {
            continue;
          }
          for (std::size_t gamma = 0; gamma < 2U; ++gamma) {
            const double cluster_delta_r =
                delta_r(truth_eta[gamma], truth_phi[gamma],
                        branches.cluster_eta->at(cluster),
                        branches.cluster_phi->at(cluster));
            const double response =
                branches.cluster_et->at(cluster) / truth_daughter_pt[gamma];
            has_individual_cluster |= cluster_delta_r < delta_r_cut &&
                                      response >= individual_response_min &&
                                      response <= individual_response_max;
          }
        }
        if (has_individual_cluster) {
          classification = EventComponent::individual_anchor;
        }
      }
    }
  }

  const std::size_t component = static_cast<std::size_t>(classification);
  ++histograms.n_component[component];
  histograms.component[component]->Fill(truth_alpha);
}

void make_truth_pi0_component_fractions(TruthPi0ComponentHistograms &histograms,
                                        const std::string &prefix) {
  for (std::size_t component = 0; component < n_event_components; ++component) {
    histograms.fraction[component].reset(
        static_cast<TH1D *>(histograms.component[component]->Clone(
            ("h_" + prefix + "_fraction_" +
             truth_pi0_component_names[component])
                .c_str())));
    histograms.fraction[component]->SetDirectory(nullptr);
    histograms.fraction[component]->Divide(
        histograms.component[component].get(), histograms.in_truth_eta.get(),
        1.0, 1.0, "B");
  }
}
void make_component_fractions(EventComponentHistograms &histograms,
                              const std::string &prefix) {
  for (std::size_t component = 0; component < n_event_components; ++component) {
    histograms.fraction[component].reset(
        static_cast<TH1D *>(histograms.component[component]->Clone(
            ("h_" + prefix + "_fraction_" + event_component_names[component])
                .c_str())));
    histograms.fraction[component]->SetDirectory(nullptr);
    histograms.fraction[component]->Divide(
        histograms.component[component].get(), histograms.selected.get(), 1.0,
        1.0, "B");
  }
}

void fill_event_histograms(
    StageHistograms &histograms, const double truth_alpha,
    const std::vector<double> &truth_eta, const std::vector<double> &truth_phi,
    const CollectionBranches &branches, const double anchor_eta_max,
    const double anchor_et_min, const double anchor_et_max,
    const double min_cluster_energy, const double delta_r_cut,
    const double mass_min, const double mass_max) {
  ++histograms.n_reference;
  histograms.reference->Fill(truth_alpha);
  if (!valid_collection_shape(branches)) {
    ++histograms.n_malformed;
    return;
  }
  bool has_selected_anchor = false;
  for (std::size_t cluster = 0; cluster < branches.cluster_et->size();
       ++cluster) {
    has_selected_anchor |= selected_anchor(branches, cluster, anchor_eta_max,
                                           anchor_et_min, anchor_et_max);
  }
  if (!has_selected_anchor) {
    return;
  }
  ++histograms.n_denominator;
  histograms.denominator->Fill(truth_alpha);

  const MatchResult match =
      match_event(truth_eta, truth_phi, branches, min_cluster_energy);
  if (!match.valid ||
      !(std::max(match.delta_r0, match.delta_r1) < delta_r_cut) ||
      (!selected_anchor(branches, match.cluster0, anchor_eta_max, anchor_et_min,
                        anchor_et_max) &&
       !selected_anchor(branches, match.cluster1, anchor_eta_max, anchor_et_min,
                        anchor_et_max))) {
    return;
  }
  ++histograms.n_matched;
  histograms.matched->Fill(truth_alpha);
  double mass = -999.0;
  if (!find_pair_mass(branches, match.cluster0, match.cluster1, mass)) {
    ++histograms.n_malformed;
    return;
  }
  if (mass >= mass_min && mass <= mass_max) {
    ++histograms.n_removed;
    histograms.removed->Fill(truth_alpha);
  }
}

void fill_cluster_histograms(
    StageHistograms &histograms, const double truth_alpha,
    const std::vector<double> &truth_eta, const std::vector<double> &truth_phi,
    const CollectionBranches &branches, const double anchor_eta_max,
    const double anchor_et_min, const double anchor_et_max,
    const double min_cluster_energy, const double delta_r_cut,
    const double mass_min, const double mass_max) {
  if (!valid_collection_shape(branches)) {
    ++histograms.n_malformed;
    return;
  }
  for (std::size_t anchor = 0; anchor < branches.cluster_et->size(); ++anchor) {
    if (!valid_cluster(branches, anchor) ||
        !(std::abs(branches.cluster_eta->at(anchor)) < anchor_eta_max)) {
      continue;
    }
    ++histograms.n_reference;
    histograms.reference->Fill(truth_alpha);
    if (!selected_anchor(branches, anchor, anchor_eta_max, anchor_et_min,
                         anchor_et_max)) {
      continue;
    }
    ++histograms.n_denominator;
    histograms.denominator->Fill(truth_alpha);
    const MatchResult match = match_anchor(anchor, truth_eta, truth_phi,
                                           branches, min_cluster_energy);
    if (!match.valid ||
        !(std::max(match.delta_r0, match.delta_r1) < delta_r_cut)) {
      continue;
    }
    ++histograms.n_matched;
    histograms.matched->Fill(truth_alpha);
    double mass = -999.0;
    if (!find_pair_mass(branches, match.cluster0, match.cluster1, mass)) {
      ++histograms.n_malformed;
      continue;
    }
    if (mass >= mass_min && mass <= mass_max) {
      ++histograms.n_removed;
      histograms.removed->Fill(truth_alpha);
    }
  }
}

void make_efficiencies(StageHistograms &histograms, const std::string &prefix) {
  histograms.efficiency_matched.reset(
      static_cast<TH1D *>(histograms.matched->Clone(
          ("h_" + prefix + "_efficiency_matched").c_str())));
  histograms.efficiency_removed.reset(
      static_cast<TH1D *>(histograms.removed->Clone(
          ("h_" + prefix + "_efficiency_removed").c_str())));
  histograms.efficiency_matched->SetDirectory(nullptr);
  histograms.efficiency_removed->SetDirectory(nullptr);
  histograms.efficiency_matched->Divide(
      histograms.matched.get(), histograms.denominator.get(), 1.0, 1.0, "B");
  histograms.efficiency_removed->Divide(
      histograms.removed.get(), histograms.denominator.get(), 1.0, 1.0, "B");
}

void style_counts(StageHistograms &histograms, const bool event_family) {
  const auto &markers =
      event_family ? event_count_markers : cluster_count_markers;
  const std::array<TH1D *, n_core_conditions> styled_histograms = {
      histograms.reference.get(), histograms.denominator.get(),
      histograms.matched.get(), histograms.removed.get()};
  for (std::size_t condition = 0; condition < n_core_conditions; ++condition) {
    TH1D *histogram = styled_histograms[condition];
    set_point_style(*histogram, core_condition_colors[condition],
                    markers[condition]);
    histogram->GetXaxis()->SetTitle("Truth #alpha");
    histogram->GetYaxis()->SetTitle(event_family ? "Events" : "Clusters");
  }
}

void draw_information(StageHistograms &histograms,
                      const std::string &collection_label,
                      const bool event_family, const bool efficiency,
                      const double anchor_eta_max, const double anchor_et_min,
                      const double anchor_et_max, const double delta_r_cut,
                      const double mass_min, const double mass_max,
                      const bool conditions) {
  if (conditions) {
    TLatex label;
    label.SetNDC();
    label.SetTextAlign(13);
    label.SetTextSize(0.044);
    label.DrawLatex(0.14, 0.92, "#it{#bf{sPHENIX}} Internal");
    label.DrawLatex(0.14, 0.84,
                    ("Single #pi^{0} gun, " + collection_label).c_str());
    label.DrawLatex(0.14, 0.76,
                    event_family ? "Event-conditional efficiency"
                                 : "Anchor-cluster conditional efficiency");
    label.DrawLatex(
        0.14, 0.68,
        Form("|#eta_{anchor}| < %.1f, %.0f < E_{T}^{anchor} < %.0f GeV",
             anchor_eta_max, anchor_et_min, anchor_et_max));
    label.DrawLatex(0.14, 0.60,
                    Form("max(#DeltaR_{0},#DeltaR_{1}) < %.3f", delta_r_cut));
    label.DrawLatex(
        0.14, 0.52,
        Form("%.2f #leq m_{#gamma#gamma} #leq %.2f GeV", mass_min, mass_max));
    return;
  }

  TLegend legend(0.10, 0.25, 0.95, 0.75);
  legend.SetTextSize(0.040);
  if (efficiency) {
    legend.AddEntry(histograms.efficiency_matched.get(),
                    event_family ? "Correct pair / selected events"
                                 : "Correct partner / selected clusters",
                    "lep");
    legend.AddEntry(histograms.efficiency_removed.get(),
                    event_family ? "Mass-window pair / selected events"
                                 : "Removed / selected clusters",
                    "lep");
  } else if (event_family) {
    legend.AddEntry(histograms.reference.get(),
                    "Truth #pi^{0} #rightarrow #gamma#gamma in acceptance",
                    "lep");
    legend.AddEntry(histograms.denominator.get(),
                    "Events with selected cluster", "lep");
    legend.AddEntry(histograms.matched.get(),
                    "Events with correct matched pair", "lep");
    legend.AddEntry(histograms.removed.get(), "+ pair mass window", "lep");
  } else {
    legend.AddEntry(histograms.reference.get(), "Central clusters", "lep");
    legend.AddEntry(histograms.denominator.get(), "Selected anchor clusters",
                    "lep");
    legend.AddEntry(histograms.matched.get(), "Anchors with correct partner",
                    "lep");
    legend.AddEntry(histograms.removed.get(), "Anchors removed by mass window",
                    "lep");
  }
  legend.DrawClone();
}

void draw_counts(std::vector<StageHistograms> &histograms,
                 const std::string &collection_label, const bool event_family,
                 const std::string &output_path, const double anchor_eta_max,
                 const double anchor_et_min, const double anchor_et_max,
                 const double delta_r_cut, const double mass_min,
                 const double mass_max) {
  TCanvas canvas(
      ("c_counts_" + collection_label + (event_family ? "_event" : "_cluster"))
          .c_str(),
      "Conditional counts by truth pT", canvas_width, canvas_height);
  canvas.Divide(canvas_columns, canvas_rows);
  for (std::size_t bin = 0; bin < histograms.size(); ++bin) {
    canvas.cd(static_cast<int>(bin + 1U));
    StageHistograms &current = histograms[bin];
    style_counts(current, event_family);
    current.reference->SetMinimum(0.0);
    current.reference->SetMaximum(current.reference->GetMaximum() > 0.0
                                      ? 1.30 * current.reference->GetMaximum()
                                      : 1.0);
    current.reference->Draw("E1");
    current.denominator->Draw("E1 SAME");
    current.matched->Draw("E1 SAME");
    current.removed->Draw("E1 SAME");
    TLatex panel;
    panel.SetNDC();
    panel.SetTextAlign(13);
    panel.SetTextSize(0.045);
    panel.DrawLatex(0.18, 0.92, truth_pt_labels[bin].c_str());
    gPad->RedrawAxis();
  }
  canvas.cd(conditions_pad);
  draw_information(histograms.front(), collection_label, event_family, false,
                   anchor_eta_max, anchor_et_min, anchor_et_max, delta_r_cut,
                   mass_min, mass_max, true);
  canvas.cd(legend_pad);
  draw_information(histograms.front(), collection_label, event_family, false,
                   anchor_eta_max, anchor_et_min, anchor_et_max, delta_r_cut,
                   mass_min, mass_max, false);
  canvas.SaveAs(output_path.c_str());
}

void draw_efficiencies(std::vector<StageHistograms> &histograms,
                       const std::string &collection_label,
                       const bool event_family, const std::string &output_path,
                       const double anchor_eta_max, const double anchor_et_min,
                       const double anchor_et_max, const double delta_r_cut,
                       const double mass_min, const double mass_max) {
  TCanvas canvas(("c_efficiency_" + collection_label +
                  (event_family ? "_event" : "_cluster"))
                     .c_str(),
                 "Conditional efficiencies by truth pT", canvas_width,
                 canvas_height);
  canvas.Divide(canvas_columns, canvas_rows);
  for (std::size_t bin = 0; bin < histograms.size(); ++bin) {
    canvas.cd(static_cast<int>(bin + 1U));
    TH1D *matched = histograms[bin].efficiency_matched.get();
    TH1D *removed = histograms[bin].efficiency_removed.get();
    const auto &markers =
        event_family ? event_ratio_markers : cluster_ratio_markers;
    set_point_style(
        *matched,
        core_condition_colors[condition_index(CoreCondition::matched)],
        markers[condition_index(CoreCondition::matched)], 2, 0.9);
    set_point_style(
        *removed,
        core_condition_colors[condition_index(CoreCondition::removed)],
        markers[condition_index(CoreCondition::removed)], 2, 0.9);
    for (TH1D *histogram : {matched, removed}) {
      histogram->SetStats(false);
      histogram->GetXaxis()->SetTitle("Truth #alpha");
      histogram->GetYaxis()->SetTitle("Conditional efficiency");
      histogram->SetMinimum(0.0);
      histogram->SetMaximum(1.05);
    }
    matched->Draw("E1");
    removed->Draw("E1 SAME");
    TLatex panel;
    panel.SetNDC();
    panel.SetTextAlign(13);
    panel.SetTextSize(0.045);
    panel.DrawLatex(0.18, 0.92, truth_pt_labels[bin].c_str());
    gPad->RedrawAxis();
  }
  canvas.cd(conditions_pad);
  draw_information(histograms.front(), collection_label, event_family, true,
                   anchor_eta_max, anchor_et_min, anchor_et_max, delta_r_cut,
                   mass_min, mass_max, true);
  canvas.cd(legend_pad);
  draw_information(histograms.front(), collection_label, event_family, true,
                   anchor_eta_max, anchor_et_min, anchor_et_max, delta_r_cut,
                   mass_min, mass_max, false);
  canvas.SaveAs(output_path.c_str());
}

void style_component_histogram(TH1D &histogram, const int color,
                               const int marker, const int line_style = 1,
                               const double marker_size = 0.65) {
  set_point_style(histogram, color, marker, line_style, marker_size);
  histogram.SetFillColorAlpha(color, 0.45);
  histogram.GetXaxis()->SetTitle("Truth #alpha");
}

void draw_component_information(
    EventComponentHistograms &histograms, const std::string &collection_label,
    const bool fractions, const double anchor_eta_max,
    const double anchor_et_min, const double anchor_et_max,
    const double delta_r_cut, const double merged_delta_r_cut,
    const double merged_response_min, const double merged_response_max,
    const double individual_response_min, const double individual_response_max,
    const bool conditions) {
  if (conditions) {
    TLatex label;
    label.SetNDC();
    label.SetTextAlign(13);
    label.SetTextSize(0.039);
    label.DrawLatex(0.10, 0.94, "#it{#bf{sPHENIX}} Internal");
    label.DrawLatex(
        0.10, 0.87,
        ("Single #pi^{0} gun, " + collection_label +
         (fractions ? ", selected-event fractions" : ", event components"))
            .c_str());
    label.DrawLatex(
        0.10, 0.80,
        Form("|#eta_{anchor}| < %.1f, %.0f < E_{T}^{anchor} < %.0f GeV",
             anchor_eta_max, anchor_et_min, anchor_et_max));
    label.DrawLatex(0.10, 0.73,
                    Form("Pair/individual match: #DeltaR < %.3f", delta_r_cut));
    label.DrawLatex(0.10, 0.66,
                    Form("Merged: max #DeltaR < %.3f, %.1f < "
                         "E_{T}^{anchor}/p_{T}^{#pi^{0}} < %.1f",
                         merged_delta_r_cut, merged_response_min,
                         merged_response_max));
    label.DrawLatex(
        0.10, 0.59,
        Form("Individual: %.1f < E_{T}^{anchor}/p_{T}^{#gamma} < %.1f",
             individual_response_min, individual_response_max));
    label.DrawLatex(
        0.10, 0.52,
        "Exclusive priority: correct > merged > individual > other");

    return;
  }

  TLegend legend(0.08, 0.16, 0.96, 0.84);
  legend.SetTextSize(0.040);
  if (!fractions) {
    legend.AddEntry(histograms.reference.get(),
                    "Truth #pi^{0} #rightarrow #gamma#gamma", "lep");
    legend.AddEntry(histograms.in_acceptance.get(),
                    "Truth #pi^{0} #rightarrow #gamma#gamma in acceptance",
                    "lep");
    legend.AddEntry(histograms.selected.get(), "Events with selected cluster",
                    "lep");
  }
  for (std::size_t component = 0; component < n_event_components; ++component) {
    legend.AddEntry(fractions ? histograms.fraction[component].get()
                              : histograms.component[component].get(),
                    event_component_labels[component].c_str(),
                    fractions ? "lep" : "f");
  }
  legend.DrawClone();
}

void draw_component_counts(
    std::vector<EventComponentHistograms> &histograms,
    const std::string &collection_label, const std::string &output_path,
    const double anchor_eta_max, const double anchor_et_min,
    const double anchor_et_max, const double delta_r_cut,
    const double merged_delta_r_cut, const double merged_response_min,
    const double merged_response_max, const double individual_response_min,
    const double individual_response_max) {
  const std::array<std::size_t, n_event_components> stack_order = {0, 1, 2, 3};
  TCanvas canvas(("c_event_components_" + collection_label).c_str(),
                 "Selected-event components by truth pT", canvas_width,
                 canvas_height);
  canvas.Divide(canvas_columns, canvas_rows);
  std::vector<std::unique_ptr<THStack>> stacks;
  stacks.reserve(histograms.size());
  for (std::size_t bin = 0; bin < histograms.size(); ++bin) {
    canvas.cd(static_cast<int>(bin + 1U));
    EventComponentHistograms &current = histograms[bin];
    set_point_style(
        *current.reference,
        core_condition_colors[condition_index(CoreCondition::reference)],
        event_count_markers[condition_index(CoreCondition::reference)]);
    current.reference->GetXaxis()->SetTitle("Truth #alpha");
    current.reference->GetYaxis()->SetTitle("Events");
    current.reference->SetMinimum(0.0);
    current.reference->SetMaximum(current.reference->GetMaximum() > 0.0
                                      ? 1.30 * current.reference->GetMaximum()
                                      : 1.0);
    set_point_style(*current.in_acceptance, TColor::GetColor("#17BECF"),
                    kOpenSquare);
    set_point_style(
        *current.selected,
        core_condition_colors[condition_index(CoreCondition::selected)],
        event_count_markers[condition_index(CoreCondition::selected)]);
    for (std::size_t component = 0; component < n_event_components;
         ++component) {
      style_component_histogram(*current.component[component],
                                event_component_colors[component],
                                event_component_count_markers[component]);
    }
    stacks.push_back(std::make_unique<THStack>(
        ("stack_" + collection_label + "_" + std::to_string(bin)).c_str(), ""));
    for (const std::size_t component : stack_order) {
      stacks.back()->Add(current.component[component].get());
    }
    current.reference->Draw("E1");
    stacks.back()->Draw("HIST SAME");
    current.selected->Draw("E1 SAME");
    current.in_acceptance->Draw("E1 SAME");
    current.reference->Draw("E1 SAME");
    TLatex panel;
    panel.SetNDC();
    panel.SetTextAlign(13);
    panel.SetTextSize(0.045);
    panel.DrawLatex(0.18, 0.92, truth_pt_labels[bin].c_str());
    gPad->RedrawAxis();
  }
  canvas.cd(conditions_pad);
  draw_component_information(
      histograms.front(), collection_label, false, anchor_eta_max,
      anchor_et_min, anchor_et_max, delta_r_cut, merged_delta_r_cut,
      merged_response_min, merged_response_max, individual_response_min,
      individual_response_max, true);
  canvas.cd(legend_pad);
  draw_component_information(
      histograms.front(), collection_label, false, anchor_eta_max,
      anchor_et_min, anchor_et_max, delta_r_cut, merged_delta_r_cut,
      merged_response_min, merged_response_max, individual_response_min,
      individual_response_max, false);
  canvas.SaveAs(output_path.c_str());
}

void draw_component_fractions(
    std::vector<EventComponentHistograms> &histograms,
    const std::string &collection_label, const std::string &output_path,
    const double anchor_eta_max, const double anchor_et_min,
    const double anchor_et_max, const double delta_r_cut,
    const double merged_delta_r_cut, const double merged_response_min,
    const double merged_response_max, const double individual_response_min,
    const double individual_response_max) {
  TCanvas canvas(("c_event_component_fractions_" + collection_label).c_str(),
                 "Selected-event component fractions by truth pT", canvas_width,
                 canvas_height);
  canvas.Divide(canvas_columns, canvas_rows);
  for (std::size_t bin = 0; bin < histograms.size(); ++bin) {
    canvas.cd(static_cast<int>(bin + 1U));
    EventComponentHistograms &current = histograms[bin];
    for (std::size_t component = 0; component < n_event_components;
         ++component) {
      TH1D &fraction = *current.fraction[component];
      style_component_histogram(fraction, event_component_colors[component],
                                event_component_ratio_markers[component], 2,
                                0.9);
      fraction.SetFillStyle(0);
      fraction.GetYaxis()->SetTitle("Fraction of selected events");
      fraction.SetMinimum(0.0);
      fraction.SetMaximum(1.05);
      fraction.Draw(component == 0U ? "E1" : "E1 SAME");
    }
    TLatex panel;
    panel.SetNDC();
    panel.SetTextAlign(13);
    panel.SetTextSize(0.045);
    panel.DrawLatex(0.18, 0.92, truth_pt_labels[bin].c_str());
    gPad->RedrawAxis();
  }
  canvas.cd(conditions_pad);
  draw_component_information(
      histograms.front(), collection_label, true, anchor_eta_max, anchor_et_min,
      anchor_et_max, delta_r_cut, merged_delta_r_cut, merged_response_min,
      merged_response_max, individual_response_min, individual_response_max,
      true);
  canvas.cd(legend_pad);
  draw_component_information(
      histograms.front(), collection_label, true, anchor_eta_max, anchor_et_min,
      anchor_et_max, delta_r_cut, merged_delta_r_cut, merged_response_min,
      merged_response_max, individual_response_min, individual_response_max,
      false);
  canvas.SaveAs(output_path.c_str());
}

void draw_truth_pi0_component_information(
    TruthPi0ComponentHistograms &histograms,
    const std::string &collection_label, const bool fractions,
    const double truth_eta_max, const double min_cluster_energy,
    const double delta_r_cut, const double merged_delta_r_cut,
    const double merged_response_min, const double merged_response_max,
    const double individual_response_min, const double individual_response_max,
    const bool conditions) {
  if (conditions) {
    TLatex label;
    label.SetNDC();
    label.SetTextAlign(13);
    label.SetTextSize(0.036);
    label.DrawLatex(0.08, 0.95, "#it{#bf{sPHENIX}} Internal");
    label.DrawLatex(0.08, 0.89,
                    ("Single #pi^{0} gun, " + collection_label +
                     (fractions ? ", truth-denominator fractions"
                                : ", truth-denominator components"))
                        .c_str());
    label.DrawLatex(
        0.08, 0.83,
        Form("|#eta_{truth}^{#pi^{0}}| < %.1f; no reco #eta or E_{T} cut",
             truth_eta_max));
    label.DrawLatex(0.08, 0.77,
                    Form("Separated: E_{cluster} #geq %.3g GeV, #DeltaR < %.3f",
                         min_cluster_energy, delta_r_cut));
    label.DrawLatex(0.08, 0.71,
                    Form("Merged: max #DeltaR < %.3f, %.1f < "
                         "E_{T}^{cluster}/p_{T}^{#pi^{0}} < %.1f",
                         merged_delta_r_cut, merged_response_min,
                         merged_response_max));
    label.DrawLatex(0.08, 0.65,
                    Form("Individual: #DeltaR < %.3f, %.1f < "
                         "E_{T}^{cluster}/p_{T}^{#gamma} < %.1f",
                         delta_r_cut, individual_response_min,
                         individual_response_max));
    label.DrawLatex(
        0.08, 0.59,
        "Exclusive priority: separated > merged > individual > none");
    return;
  }

  TLegend legend(0.08, 0.15, 0.96, 0.85);
  legend.SetTextSize(0.038);
  if (!fractions) {
    legend.AddEntry(histograms.reference.get(),
                    "Generated #pi^{0} #rightarrow #gamma#gamma", "lep");
    legend.AddEntry(histograms.in_truth_eta.get(),
                    "Generated #pi^{0} #rightarrow #gamma#gamma with "
                    "|#eta_{truth}^{#pi^{0}}| < cut",
                    "lep");
  }
  for (std::size_t component = 0; component < n_event_components; ++component) {
    legend.AddEntry(fractions ? histograms.fraction[component].get()
                              : histograms.component[component].get(),
                    truth_pi0_component_labels[component].c_str(),
                    fractions ? "lep" : "f");
  }
  legend.DrawClone();
}

void draw_truth_pi0_component_counts(
    std::vector<TruthPi0ComponentHistograms> &histograms,
    const std::string &collection_label, const std::string &output_path,
    const double truth_eta_max, const double min_cluster_energy,
    const double delta_r_cut, const double merged_delta_r_cut,
    const double merged_response_min, const double merged_response_max,
    const double individual_response_min,
    const double individual_response_max) {
  const std::array<std::size_t, n_event_components> stack_order = {0, 1, 2, 3};
  TCanvas canvas(("c_truth_pi0_event_components_" + collection_label).c_str(),
                 "Reconstruction topologies for central truth pi0 by truth pT",
                 canvas_width, canvas_height);
  canvas.Divide(canvas_columns, canvas_rows);
  std::vector<std::unique_ptr<THStack>> stacks;
  stacks.reserve(histograms.size());
  for (std::size_t bin = 0; bin < histograms.size(); ++bin) {
    canvas.cd(static_cast<int>(bin + 1U));
    TruthPi0ComponentHistograms &current = histograms[bin];
    set_point_style(
        *current.reference,
        core_condition_colors[condition_index(CoreCondition::reference)],
        event_count_markers[condition_index(CoreCondition::reference)]);
    current.reference->GetXaxis()->SetTitle("Truth #alpha");
    current.reference->GetYaxis()->SetTitle("Generated #pi^{0}");
    current.reference->SetMinimum(0.0);
    current.reference->SetMaximum(current.reference->GetMaximum() > 0.0
                                      ? 1.30 * current.reference->GetMaximum()
                                      : 1.0);
    set_point_style(
        *current.in_truth_eta,
        core_condition_colors[condition_index(CoreCondition::selected)],
        event_count_markers[condition_index(CoreCondition::selected)]);
    for (std::size_t component = 0; component < n_event_components;
         ++component) {
      style_component_histogram(*current.component[component],
                                event_component_colors[component],
                                event_component_count_markers[component]);
    }
    stacks.push_back(std::make_unique<THStack>(
        ("stack_truth_pi0_" + collection_label + "_" + std::to_string(bin))
            .c_str(),
        ""));
    for (const std::size_t component : stack_order) {
      stacks.back()->Add(current.component[component].get());
    }
    current.reference->Draw("E1");
    stacks.back()->Draw("HIST SAME");
    current.in_truth_eta->Draw("E1 SAME");
    current.reference->Draw("E1 SAME");
    TLatex panel;
    panel.SetNDC();
    panel.SetTextAlign(13);
    panel.SetTextSize(0.045);
    panel.DrawLatex(0.18, 0.92, truth_pt_labels[bin].c_str());
    gPad->RedrawAxis();
  }
  canvas.cd(conditions_pad);
  draw_truth_pi0_component_information(
      histograms.front(), collection_label, false, truth_eta_max,
      min_cluster_energy, delta_r_cut, merged_delta_r_cut, merged_response_min,
      merged_response_max, individual_response_min, individual_response_max,
      true);
  canvas.cd(legend_pad);
  draw_truth_pi0_component_information(
      histograms.front(), collection_label, false, truth_eta_max,
      min_cluster_energy, delta_r_cut, merged_delta_r_cut, merged_response_min,
      merged_response_max, individual_response_min, individual_response_max,
      false);
  canvas.SaveAs(output_path.c_str());
}

void draw_truth_pi0_component_fractions(
    std::vector<TruthPi0ComponentHistograms> &histograms,
    const std::string &collection_label, const std::string &output_path,
    const double truth_eta_max, const double min_cluster_energy,
    const double delta_r_cut, const double merged_delta_r_cut,
    const double merged_response_min, const double merged_response_max,
    const double individual_response_min,
    const double individual_response_max) {
  TCanvas canvas(
      ("c_truth_pi0_event_component_fractions_" + collection_label).c_str(),
      "Reconstruction-topology fractions for central truth pi0 by truth pT",
      canvas_width, canvas_height);
  canvas.Divide(canvas_columns, canvas_rows);
  for (std::size_t bin = 0; bin < histograms.size(); ++bin) {
    canvas.cd(static_cast<int>(bin + 1U));
    TruthPi0ComponentHistograms &current = histograms[bin];
    for (std::size_t component = 0; component < n_event_components;
         ++component) {
      TH1D &fraction = *current.fraction[component];
      style_component_histogram(fraction, event_component_colors[component],
                                event_component_ratio_markers[component], 2,
                                0.9);
      fraction.SetFillStyle(0);
      fraction.GetYaxis()->SetTitle("Fraction of central truth #pi^{0}");
      fraction.SetMinimum(0.0);
      fraction.SetMaximum(1.05);
      fraction.Draw(component == 0U ? "E1" : "E1 SAME");
    }
    TLatex panel;
    panel.SetNDC();
    panel.SetTextAlign(13);
    panel.SetTextSize(0.045);
    panel.DrawLatex(0.18, 0.92, truth_pt_labels[bin].c_str());
    gPad->RedrawAxis();
  }
  canvas.cd(conditions_pad);
  draw_truth_pi0_component_information(
      histograms.front(), collection_label, true, truth_eta_max,
      min_cluster_energy, delta_r_cut, merged_delta_r_cut, merged_response_min,
      merged_response_max, individual_response_min, individual_response_max,
      true);
  canvas.cd(legend_pad);
  draw_truth_pi0_component_information(
      histograms.front(), collection_label, true, truth_eta_max,
      min_cluster_energy, delta_r_cut, merged_delta_r_cut, merged_response_min,
      merged_response_max, individual_response_min, individual_response_max,
      false);
  canvas.SaveAs(output_path.c_str());
}
void draw_cut_stage_information(
    CutStageHistograms &histograms, const std::string &collection_label,
    const bool pair_family, const double min_cluster_energy,
    const double anchor_eta_max, const double anchor_et_min,
    const double anchor_et_max, const double delta_r_cut,
    const bool conditions) {
  if (conditions) {
    TLatex label;
    label.SetNDC();
    label.SetTextAlign(13);
    label.SetTextSize(0.041);
    label.DrawLatex(0.10, 0.93, "#it{#bf{sPHENIX}} Internal");
    label.DrawLatex(0.10, 0.85,
                    ("Single #pi^{0} gun, " + collection_label).c_str());
    label.DrawLatex(0.10, 0.77,
                    pair_family ? "Matched-pair selection stages"
                                : "Event cluster-selection stages");
    label.DrawLatex(0.10, 0.69,
                    Form("E_{cluster} #geq %.3g GeV", min_cluster_energy));
    label.DrawLatex(0.10, 0.61,
                    Form("|#eta| < %.1f, %.0f < E_{T} < %.0f GeV",
                         anchor_eta_max, anchor_et_min, anchor_et_max));
    if (pair_family) {
      label.DrawLatex(0.10, 0.53,
                      Form("max(#DeltaR_{0},#DeltaR_{1}) < %.3f", delta_r_cut));
    }

    return;
  }

  TLegend legend(0.08, 0.20, 0.96, 0.80);
  legend.SetTextSize(0.038);
  legend.AddEntry(histograms.stage[0].get(),
                  "Truth #pi^{0} #rightarrow #gamma#gamma in acceptance",
                  "lep");
  if (pair_family) {
    legend.AddEntry(histograms.stage[1].get(), "#DeltaR-matched pair", "lep");
    legend.AddEntry(histograms.stage[2].get(),
                    Form("+ #geq1 endpoint with |#eta| < %.1f", anchor_eta_max),
                    "lep");
    legend.AddEntry(histograms.stage[3].get(),
                    Form("+ same endpoint with %.0f < E_{T} < %.0f GeV",
                         anchor_et_min, anchor_et_max),
                    "lep");
  } else {
    legend.AddEntry(histograms.stage[1].get(),
                    "Events with #geq1 cluster above E threshold", "lep");
    legend.AddEntry(histograms.stage[2].get(),
                    Form("+ #geq1 cluster with |#eta| < %.1f", anchor_eta_max),
                    "lep");
    legend.AddEntry(histograms.stage[3].get(),
                    Form("+ same cluster with %.0f < E_{T} < %.0f GeV",
                         anchor_et_min, anchor_et_max),
                    "lep");
  }
  legend.DrawClone();
}

void draw_cut_stages(std::vector<CutStageHistograms> &histograms,
                     const std::string &collection_label,
                     const bool pair_family, const std::string &output_path,
                     const double min_cluster_energy,
                     const double anchor_eta_max, const double anchor_et_min,
                     const double anchor_et_max, const double delta_r_cut) {
  // The two cut-stage families have different physics at the same ordinal
  // stage, so give intermediate stages family-specific colors. The reference
  // and final stages reuse the exact styles of their equivalent event counts.
  const std::array<int, n_cut_stages> event_colors = {
      core_condition_colors[condition_index(CoreCondition::reference)],
      TColor::GetColor("#8C564B"), TColor::GetColor("#17BECF"),
      core_condition_colors[condition_index(CoreCondition::selected)]};
  const std::array<int, n_cut_stages> event_markers = {
      event_count_markers[condition_index(CoreCondition::reference)],
      kOpenCircle, kOpenSquare,
      event_count_markers[condition_index(CoreCondition::selected)]};
  const std::array<int, n_cut_stages> pair_colors = {
      core_condition_colors[condition_index(CoreCondition::reference)],
      TColor::GetColor("#9467BD"), TColor::GetColor("#BCBD22"),
      core_condition_colors[condition_index(CoreCondition::matched)]};
  const std::array<int, n_cut_stages> pair_markers = {
      event_count_markers[condition_index(CoreCondition::reference)],
      kOpenDiamond, kOpenCross,
      event_count_markers[condition_index(CoreCondition::matched)]};
  const auto &colors = pair_family ? pair_colors : event_colors;
  const auto &markers = pair_family ? pair_markers : event_markers;
  TCanvas canvas(("c_" + collection_label +
                  (pair_family ? "_pair_stages" : "_event_stages"))
                     .c_str(),
                 "Selection stages by truth pT", canvas_width, canvas_height);
  canvas.Divide(canvas_columns, canvas_rows);
  for (std::size_t bin = 0; bin < histograms.size(); ++bin) {
    canvas.cd(static_cast<int>(bin + 1U));
    for (std::size_t stage = 0; stage < n_cut_stages; ++stage) {
      TH1D &histogram = *histograms[bin].stage[stage];
      set_point_style(histogram, colors[stage], markers[stage]);
      histogram.GetXaxis()->SetTitle("Truth #alpha");
      histogram.GetYaxis()->SetTitle("Events");
      histogram.SetMinimum(0.0);
      histogram.SetMaximum(histograms[bin].stage[0]->GetMaximum() > 0.0
                               ? 1.30 * histograms[bin].stage[0]->GetMaximum()
                               : 1.0);
      histogram.Draw(stage == 0U ? "E1" : "E1 SAME");
    }
    histograms[bin].stage[0]->Draw("E1 SAME");
    TLatex panel;
    panel.SetNDC();
    panel.SetTextAlign(13);
    panel.SetTextSize(0.045);
    panel.DrawLatex(0.18, 0.92, truth_pt_labels[bin].c_str());
    gPad->RedrawAxis();
  }
  canvas.cd(conditions_pad);
  draw_cut_stage_information(histograms.front(), collection_label, pair_family,
                             min_cluster_energy, anchor_eta_max, anchor_et_min,
                             anchor_et_max, delta_r_cut, true);
  canvas.cd(legend_pad);
  draw_cut_stage_information(histograms.front(), collection_label, pair_family,
                             min_cluster_energy, anchor_eta_max, anchor_et_min,
                             anchor_et_max, delta_r_cut, false);
  canvas.SaveAs(output_path.c_str());
}

void write_cut_stage_histograms(TFile &output,
                                std::vector<CutStageHistograms> &histograms) {
  output.cd();
  for (CutStageHistograms &current : histograms) {
    for (std::size_t stage = 0; stage < n_cut_stages; ++stage) {
      current.stage[stage]->Write();
    }
  }
}

bool print_cut_stage_summary(
    const std::string &label,
    const std::vector<CutStageHistograms> &histograms) {
  bool nesting_ok = true;
  for (std::size_t bin = 0; bin < histograms.size(); ++bin) {
    const CutStageHistograms &current = histograms[bin];
    bool bin_ok = true;
    for (std::size_t stage = 1; stage < n_cut_stages; ++stage) {
      bin_ok &= current.n_stage[stage] <= current.n_stage[stage - 1U];
    }
    nesting_ok &= bin_ok;
    std::cout << label << " " << truth_pt_labels[bin]
              << " stage0/stage1/stage2/stage3/malformed = "
              << current.n_stage[0] << "/" << current.n_stage[1] << "/"
              << current.n_stage[2] << "/" << current.n_stage[3] << "/"
              << current.n_malformed << " nesting=" << (bin_ok ? "OK" : "FAIL")
              << std::endl;
  }
  return nesting_ok;
}

void write_component_histograms(
    TFile &output, std::vector<EventComponentHistograms> &histograms) {
  output.cd();
  for (EventComponentHistograms &current : histograms) {
    current.reference->Write();
    current.in_acceptance->Write();
    current.selected->Write();
    for (std::size_t component = 0; component < n_event_components;
         ++component) {
      current.component[component]->Write();
      current.fraction[component]->Write();
    }
  }
}

void write_truth_pi0_component_histograms(
    TFile &output, std::vector<TruthPi0ComponentHistograms> &histograms) {
  output.cd();
  for (TruthPi0ComponentHistograms &current : histograms) {
    current.reference->Write();
    current.in_truth_eta->Write();
    for (std::size_t component = 0; component < n_event_components;
         ++component) {
      current.component[component]->Write();
      current.fraction[component]->Write();
    }
  }
}
void write_histograms(TFile &output, std::vector<StageHistograms> &histograms) {
  output.cd();
  for (StageHistograms &current : histograms) {
    for (TH1D *histogram :
         {current.reference.get(), current.denominator.get(),
          current.matched.get(), current.removed.get(),
          current.efficiency_matched.get(), current.efficiency_removed.get()}) {
      histogram->Write();
    }
  }
}

void print_summary(const std::string &label,
                   const std::vector<StageHistograms> &histograms) {
  for (std::size_t bin = 0; bin < histograms.size(); ++bin) {
    const StageHistograms &current = histograms[bin];
    std::cout << label << " " << truth_pt_labels[bin]
              << " reference/denominator/matched/removed/malformed = "
              << current.n_reference << "/" << current.n_denominator << "/"
              << current.n_matched << "/" << current.n_removed << "/"
              << current.n_malformed << std::endl;
  }
}

bool print_component_summary(
    const std::string &label,
    const std::vector<EventComponentHistograms> &histograms) {
  bool closure_ok = true;
  for (std::size_t bin = 0; bin < histograms.size(); ++bin) {
    const EventComponentHistograms &current = histograms[bin];
    Long64_t component_sum = 0;
    for (const Long64_t count : current.n_component) {
      component_sum += count;
    }
    const bool bin_ok = component_sum == current.n_selected &&
                        current.n_selected <= current.n_in_acceptance &&
                        current.n_in_acceptance <= current.n_reference;
    closure_ok &= bin_ok;
    std::cout
        << label << " " << truth_pt_labels[bin]
        << " reference/acceptance/selected/correct/merged/individual/other/"
           "malformed = "
        << current.n_reference << "/" << current.n_in_acceptance << "/"
        << current.n_selected;
    for (const Long64_t count : current.n_component) {
      std::cout << "/" << count;
    }
    std::cout << "/" << current.n_malformed
              << " nesting/closure=" << (bin_ok ? "OK" : "FAIL") << std::endl;
  }
  return closure_ok;
}

bool print_truth_pi0_component_summary(
    const std::string &label,
    const std::vector<TruthPi0ComponentHistograms> &histograms) {
  bool closure_ok = true;
  for (std::size_t bin = 0; bin < histograms.size(); ++bin) {
    const TruthPi0ComponentHistograms &current = histograms[bin];
    Long64_t component_sum = 0;
    for (const Long64_t count : current.n_component) {
      component_sum += count;
    }
    const bool bin_ok = component_sum == current.n_in_truth_eta &&
                        current.n_in_truth_eta <= current.n_reference;
    closure_ok &= bin_ok;
    std::cout << label << " " << truth_pt_labels[bin]
              << " reference/truth_eta/separated/merged/individual/none/"
                 "malformed = "
              << current.n_reference << "/" << current.n_in_truth_eta;
    for (const Long64_t count : current.n_component) {
      std::cout << "/" << count;
    }
    std::cout << "/" << current.n_malformed
              << " nesting/closure=" << (bin_ok ? "OK" : "FAIL") << std::endl;
  }
  return closure_ok;
}
void draw_reco_et_information(RecoEtHistograms &histograms,
                              const std::string &collection_label,
                              const bool retention,
                              const double min_cluster_energy,
                              const double eta_max, const bool conditions) {
  if (conditions) {
    TLatex label;
    label.SetNDC();
    label.SetTextAlign(13);
    label.SetTextSize(0.041);
    label.DrawLatex(0.09, 0.94, "#it{#bf{sPHENIX}} Internal");
    label.DrawLatex(0.09, 0.86,
                    ("Single #pi^{0} gun, " + collection_label).c_str());
    label.DrawLatex(0.09, 0.78,
                    retention
                        ? "Cluster-selection retention by reconstructed E_{T}"
                        : "Truth-p_{T} contributions by reconstructed E_{T}");
    label.DrawLatex(0.09, 0.70,
                    Form("E_{cluster} #geq %.3g GeV, |#eta_{cluster}| < %.1f",
                         min_cluster_energy, eta_max));

    return;
  }

  if (retention) {
    TLegend legend(0.08, 0.25, 0.96, 0.75);
    legend.SetTextSize(0.040);
    legend.AddEntry(histograms.retention[0].get(), "E cut only / no cut",
                    "lep");
    legend.AddEntry(histograms.retention[1].get(), "#eta cut only / no cut",
                    "lep");
    legend.AddEntry(histograms.retention[2].get(), "Both cuts / no cut", "lep");
    legend.DrawClone();
    return;
  }

  TLegend legend(0.08, 0.12, 0.96, 0.88);
  legend.SetTextSize(0.036);
  for (std::size_t truth_bin = 0; truth_bin < n_truth_pt_bins; ++truth_bin) {
    legend.AddEntry(histograms.truth_pt_contribution[truth_bin].get(),
                    truth_pt_labels[truth_bin].c_str(), "f");
  }
  legend.AddEntry(
      histograms.selection[selection_index(ClusterSelection::min_energy_eta)]
          .get(),
      "Total after both cuts", "lep");
  legend.DrawClone();
}

void draw_reco_et_truth_pt_stack(std::vector<RecoEtHistograms> &histograms,
                                 const std::string &collection_label,
                                 const std::string &output_path,
                                 const double min_cluster_energy,
                                 const double eta_max) {
  // A dedicated ordered palette prevents truth-pT bins from borrowing the
  // categorical colors used for selections and event classifications.
  const std::array<int, n_truth_pt_bins> colors = {
      TColor::GetColor("#440154"), TColor::GetColor("#414487"),
      TColor::GetColor("#2A788E"), TColor::GetColor("#22A884"),
      TColor::GetColor("#7AD151"), TColor::GetColor("#FDE725")};
  TCanvas canvas(("c_reco_et_truth_pt_stack_" + collection_label).c_str(),
                 "Truth-pT contributions in reconstructed-ET panels",
                 canvas_width, canvas_height);
  canvas.Divide(canvas_columns, canvas_rows);
  std::vector<std::unique_ptr<THStack>> stacks;
  stacks.reserve(histograms.size());
  for (std::size_t reco_bin = 0; reco_bin < histograms.size(); ++reco_bin) {
    canvas.cd(static_cast<int>(reco_bin + 1U));
    RecoEtHistograms &current = histograms[reco_bin];
    TH1D &total =
        *current.selection[selection_index(ClusterSelection::min_energy_eta)];
    set_point_style(total, kBlack, kFullDiamond);
    total.GetXaxis()->SetTitle("Truth #alpha");
    total.GetYaxis()->SetTitle("Clusters");
    total.SetMinimum(0.0);
    total.SetMaximum(total.GetMaximum() > 0.0 ? 1.30 * total.GetMaximum()
                                              : 1.0);

    stacks.push_back(std::make_unique<THStack>(
        ("stack_reco_et_" + collection_label + "_" + std::to_string(reco_bin))
            .c_str(),
        ""));
    for (std::size_t truth_bin = 0; truth_bin < n_truth_pt_bins; ++truth_bin) {
      TH1D &contribution = *current.truth_pt_contribution[truth_bin];
      contribution.SetLineColor(colors[truth_bin]);
      contribution.SetFillColorAlpha(colors[truth_bin], 0.55);
      contribution.SetLineWidth(2);
      stacks.back()->Add(&contribution);
    }
    total.Draw("E1");
    stacks.back()->Draw("HIST SAME");
    total.Draw("E1 SAME");
    TLatex panel;
    panel.SetNDC();
    panel.SetTextAlign(13);
    panel.SetTextSize(0.045);
    panel.DrawLatex(0.18, 0.92, reco_et_labels[reco_bin].c_str());
    gPad->RedrawAxis();
  }
  canvas.cd(conditions_pad);
  draw_reco_et_information(histograms.front(), collection_label, false,
                           min_cluster_energy, eta_max, true);
  canvas.cd(legend_pad);
  draw_reco_et_information(histograms.front(), collection_label, false,
                           min_cluster_energy, eta_max, false);
  canvas.SaveAs(output_path.c_str());
}

void draw_reco_et_retentions(std::vector<RecoEtHistograms> &histograms,
                             const std::string &collection_label,
                             const std::string &output_path,
                             const double min_cluster_energy,
                             const double eta_max) {
  const std::array<int, n_cluster_retentions> colors = {kBlue + 1, kOrange + 7,
                                                        kGreen + 2};
  const std::array<int, n_cluster_retentions> markers = {24, 25, 20};
  TCanvas canvas(("c_reco_et_retention_" + collection_label).c_str(),
                 "Cluster-selection retention in reconstructed-ET panels",
                 canvas_width, canvas_height);
  canvas.Divide(canvas_columns, canvas_rows);
  for (std::size_t reco_bin = 0; reco_bin < histograms.size(); ++reco_bin) {
    canvas.cd(static_cast<int>(reco_bin + 1U));
    for (std::size_t retention = 0; retention < n_cluster_retentions;
         ++retention) {
      TH1D &histogram = *histograms[reco_bin].retention[retention];
      histogram.SetLineColor(colors[retention]);
      histogram.SetMarkerColor(colors[retention]);
      histogram.SetMarkerStyle(markers[retention]);
      histogram.SetMarkerSize(0.7);
      histogram.SetLineWidth(2);
      histogram.GetXaxis()->SetTitle("Truth #alpha");
      histogram.GetYaxis()->SetTitle("Fraction retained");
      histogram.SetMinimum(0.0);
      histogram.SetMaximum(1.05);
      histogram.Draw(retention == 0U ? "E1" : "E1 SAME");
    }
    TLatex panel;
    panel.SetNDC();
    panel.SetTextAlign(13);
    panel.SetTextSize(0.045);
    panel.DrawLatex(0.18, 0.92, reco_et_labels[reco_bin].c_str());
    gPad->RedrawAxis();
  }
  canvas.cd(conditions_pad);
  draw_reco_et_information(histograms.front(), collection_label, true,
                           min_cluster_energy, eta_max, true);
  canvas.cd(legend_pad);
  draw_reco_et_information(histograms.front(), collection_label, true,
                           min_cluster_energy, eta_max, false);
  canvas.SaveAs(output_path.c_str());
}

void write_reco_et_histograms(TFile &output,
                              std::vector<RecoEtHistograms> &histograms) {
  output.cd();
  for (RecoEtHistograms &current : histograms) {
    for (auto &histogram : current.truth_pt_contribution) {
      histogram->Write();
    }
    for (auto &histogram : current.selection) {
      histogram->Write();
    }
    for (auto &histogram : current.retention) {
      histogram->Write();
    }
  }
}

bool print_reco_et_summary(const std::string &label,
                           const std::vector<RecoEtHistograms> &histograms) {
  bool selection_ok = true;
  for (std::size_t reco_bin = 0; reco_bin < histograms.size(); ++reco_bin) {
    const RecoEtHistograms &current = histograms[reco_bin];
    Long64_t contribution_sum = 0;
    for (const Long64_t count : current.n_truth_pt_contribution) {
      contribution_sum += count;
    }
    const Long64_t all =
        current.n_selection[selection_index(ClusterSelection::all)];
    const Long64_t energy =
        current.n_selection[selection_index(ClusterSelection::min_energy)];
    const Long64_t eta =
        current.n_selection[selection_index(ClusterSelection::eta)];
    const Long64_t both = current.n_selection[
        selection_index(ClusterSelection::min_energy_eta)];
    const bool bin_ok = contribution_sum == both && energy <= all && eta <= all &&
                        both <= energy && both <= eta;
    selection_ok &= bin_ok;
    std::cout << label << " " << reco_et_labels[reco_bin]
              << " all/energy/eta/both/contribution_sum = " << all << "/"
              << energy << "/" << eta << "/" << both << "/"
              << contribution_sum << " closure=" << (bin_ok ? "OK" : "FAIL")
              << std::endl;
  }
  return selection_ok;
}

} // namespace

int PlotConditionalPartnerEfficiency(
    const std::string input_path =
        "PhotonAnalysisTree/output/merged/"
        "100kevents_pi0_3to15GeV_etapm1_vertexpm60.root",
    const std::string output_base =
        "PhotonAnalysisTree/output/plots/conditional_efficiency/"
        "conditional_partner_et3to5",
    const double anchor_eta_max = 0.7, const double anchor_et_min = 3.0,
    const double anchor_et_max = 5.0, const double delta_r_cut = 0.03,
    const double mass_window_min = 0.10, const double mass_window_max = 0.18,
    const int asymmetry_nbins = 20, const double merged_delta_r_cut = 0.06,
    const double merged_response_min = 0.5,
    const double merged_response_max = 1.5,
    const double individual_response_min = 0.5,
    const double individual_response_max = 1.5,
    const double min_cluster_energy = 0.3, const double truth_eta_max = 0.7) {
  if (input_path.empty() || output_base.empty() || !(anchor_eta_max > 0.0) ||
      !(anchor_et_min >= 0.0 && anchor_et_min < anchor_et_max) ||
      !(delta_r_cut > 0.0) || !(merged_delta_r_cut > 0.0) ||
      !(merged_response_min >= 0.0 &&
        merged_response_min < merged_response_max) ||
      !(individual_response_min >= 0.0 &&
        individual_response_min < individual_response_max) ||
      !(min_cluster_energy >= 0.0) || !(truth_eta_max > 0.0) ||
      !(mass_window_min >= 0.0 && mass_window_min < mass_window_max) ||
      asymmetry_nbins <= 0) {
    std::cerr << "PlotConditionalPartnerEfficiency - invalid argument"
              << std::endl;
    return 1;
  }
  SetsPhenixStyle();
  TH1::AddDirectory(false);
  std::unique_ptr<TFile> input(TFile::Open(input_path.c_str(), "READ"));
  if (!input || input->IsZombie()) {
    std::cerr << "Failed to open " << input_path << std::endl;
    return 2;
  }
  TTree *tree = input->Get<TTree>("event_tree");
  if (!tree) {
    std::cerr << "Missing TTree event_tree in " << input_path << std::endl;
    return 3;
  }

  unsigned char truth_valid = 0;
  unsigned char truth_is_pi0_to_2gamma = 0;
  unsigned char truth_both_gamma_in_acceptance = 0;
  double truth_pt = -999.0;
  double truth_pi0_eta = -999.0;
  double truth_alpha = -999.0;
  std::vector<double> *truth_eta = nullptr;
  std::vector<double> *truth_phi = nullptr;
  std::vector<double> *truth_daughter_pt = nullptr;
  std::vector<unsigned char> *truth_projection_valid = nullptr;
  CollectionBranches split;
  CollectionBranches nosplit;
  tree->SetBranchStatus("*", false);
  bool branches_ok = true;
  const auto bind = [&](const std::string &name, auto *address) {
    if (!tree->GetBranch(name.c_str())) {
      std::cerr << "Missing branch " << name << std::endl;
      branches_ok = false;
      return;
    }
    tree->SetBranchStatus(name.c_str(), true);
    if (tree->SetBranchAddress(name.c_str(), address) < 0) {
      std::cerr << "Failed to bind branch " << name << std::endl;
      branches_ok = false;
    }
  };
  const auto bind_collection = [&](const std::string &prefix,
                                   CollectionBranches &branches) {
    bind(prefix + "_cluster_e", &branches.cluster_e);
    bind(prefix + "_cluster_et", &branches.cluster_et);
    bind(prefix + "_cluster_eta", &branches.cluster_eta);
    bind(prefix + "_cluster_phi", &branches.cluster_phi);
    bind(prefix + "_pair_cluster_i", &branches.pair_i);
    bind(prefix + "_pair_cluster_j", &branches.pair_j);
    bind(prefix + "_pair_m_gg", &branches.pair_mass);
  };
  bind("truth_valid", &truth_valid);
  bind("truth_is_pi0_to_2gamma", &truth_is_pi0_to_2gamma);
  bind("truth_both_gamma_in_acceptance", &truth_both_gamma_in_acceptance);
  bind("truth_pt", &truth_pt);
  bind("truth_eta", &truth_pi0_eta);
  bind("truth_pair_e_asym", &truth_alpha);
  bind("truth_daughter_pt", &truth_daughter_pt);
  bind("truth_daughter_eta", &truth_eta);
  bind("truth_daughter_phi", &truth_phi);
  bind("truth_daughter_projection_valid", &truth_projection_valid);
  bind_collection("split", split);
  bind_collection("nosplit", nosplit);
  if (!branches_ok) {
    return 4;
  }

  std::vector<StageHistograms> split_event =
      make_truth_pt_histograms("split_event", asymmetry_nbins);
  std::vector<StageHistograms> split_cluster =
      make_truth_pt_histograms("split_cluster", asymmetry_nbins);
  std::vector<StageHistograms> nosplit_event =
      make_truth_pt_histograms("nosplit_event", asymmetry_nbins);
  std::vector<StageHistograms> nosplit_cluster =
      make_truth_pt_histograms("nosplit_cluster", asymmetry_nbins);
  std::vector<EventComponentHistograms> split_components =
      make_truth_pt_component_histograms("split_event_component",
                                         asymmetry_nbins);
  std::vector<EventComponentHistograms> nosplit_components =
      make_truth_pt_component_histograms("nosplit_event_component",
                                         asymmetry_nbins);
  std::vector<TruthPi0ComponentHistograms> split_truth_pi0_components =
      make_truth_pt_truth_pi0_component_histograms(
          "split_central_truth_pi0_component", asymmetry_nbins);
  std::vector<TruthPi0ComponentHistograms> nosplit_truth_pi0_components =
      make_truth_pt_truth_pi0_component_histograms(
          "nosplit_central_truth_pi0_component", asymmetry_nbins);
  const std::array<std::string, n_cut_stages> event_stage_suffixes = {
      "truth", "cluster_energy", "central_cluster", "selected_cluster"};
  const std::array<std::string, n_cut_stages> pair_stage_suffixes = {
      "truth", "matched_pair", "matched_pair_central", "matched_pair_selected"};
  std::vector<CutStageHistograms> split_event_stages =
      make_truth_pt_cut_stage_histograms("split_event_stage",
                                         event_stage_suffixes, asymmetry_nbins);
  std::vector<CutStageHistograms> split_pair_stages =
      make_truth_pt_cut_stage_histograms("split_pair_stage",
                                         pair_stage_suffixes, asymmetry_nbins);
  std::vector<CutStageHistograms> nosplit_event_stages =
      make_truth_pt_cut_stage_histograms("nosplit_event_stage",
                                         event_stage_suffixes, asymmetry_nbins);
  std::vector<CutStageHistograms> nosplit_pair_stages =
      make_truth_pt_cut_stage_histograms("nosplit_pair_stage",
                                         pair_stage_suffixes, asymmetry_nbins);
  std::vector<RecoEtHistograms> split_reco_et =
      make_reco_et_histogram_set("split_cluster", asymmetry_nbins);
  std::vector<RecoEtHistograms> nosplit_reco_et =
      make_reco_et_histogram_set("nosplit_cluster", asymmetry_nbins);

  Long64_t invalid_truth_shape = 0;
  Long64_t split_reco_et_malformed = 0;
  Long64_t nosplit_reco_et_malformed = 0;
  Long64_t truth_pt_outside_bins = 0;
  for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry) {
    if (entry % 5000 == 0) {
      std::cout << "PlotConditionalPartnerEfficiency - entry " << entry << " / "
                << tree->GetEntries() << std::endl;
    }
    tree->GetEntry(entry);
    if (!truth_valid || !truth_is_pi0_to_2gamma) {
      continue;
    }
    if (!std::isfinite(truth_alpha) || truth_alpha < 0.0 || truth_alpha > 1.0 ||
        !std::isfinite(truth_pi0_eta) || !truth_eta || !truth_phi ||
        !truth_daughter_pt || truth_eta->size() != 2U ||
        truth_phi->size() != 2U || truth_daughter_pt->size() != 2U ||
        !std::isfinite(truth_daughter_pt->at(0)) ||
        !std::isfinite(truth_daughter_pt->at(1))) {
      ++invalid_truth_shape;
      continue;
    }
    const std::size_t bin = find_truth_pt_bin(truth_pt);
    if (bin >= truth_pt_labels.size()) {
      ++truth_pt_outside_bins;
      continue;
    }

    fill_event_component_reference(split_components[bin], truth_alpha);
    fill_event_component_reference(nosplit_components[bin], truth_alpha);

    const bool passes_truth_eta = std::abs(truth_pi0_eta) < truth_eta_max;
    fill_truth_pi0_components(split_truth_pi0_components[bin], truth_alpha,
                              passes_truth_eta, truth_pt, *truth_daughter_pt,
                              *truth_eta, *truth_phi, split, min_cluster_energy,
                              delta_r_cut, merged_delta_r_cut,
                              merged_response_min, merged_response_max,
                              individual_response_min, individual_response_max);
    fill_truth_pi0_components(
        nosplit_truth_pi0_components[bin], truth_alpha, passes_truth_eta,
        truth_pt, *truth_daughter_pt, *truth_eta, *truth_phi, nosplit,
        min_cluster_energy, delta_r_cut, merged_delta_r_cut,
        merged_response_min, merged_response_max, individual_response_min,
        individual_response_max);

    // Preserve the original acceptance requirements for the existing plot
    // families. They are deliberately not part of the new truth-pi0
    // denominator above.
    if (!truth_both_gamma_in_acceptance || !truth_projection_valid ||
        truth_projection_valid->size() != 2U ||
        !truth_projection_valid->at(0) || !truth_projection_valid->at(1)) {
      continue;
    }
    if (!fill_reco_et_histograms(split_reco_et, bin, truth_alpha, split,
                                  min_cluster_energy, anchor_eta_max)) {
      ++split_reco_et_malformed;
    }
    if (!fill_reco_et_histograms(nosplit_reco_et, bin, truth_alpha, nosplit,
                                  min_cluster_energy, anchor_eta_max)) {
      ++nosplit_reco_et_malformed;
    }
    fill_event_histograms(split_event[bin], truth_alpha, *truth_eta, *truth_phi,
                          split, anchor_eta_max, anchor_et_min, anchor_et_max,
                          min_cluster_energy, delta_r_cut, mass_window_min,
                          mass_window_max);
    fill_cluster_histograms(split_cluster[bin], truth_alpha, *truth_eta,
                            *truth_phi, split, anchor_eta_max, anchor_et_min,
                            anchor_et_max, min_cluster_energy, delta_r_cut,
                            mass_window_min, mass_window_max);
    fill_event_histograms(nosplit_event[bin], truth_alpha, *truth_eta,
                          *truth_phi, nosplit, anchor_eta_max, anchor_et_min,
                          anchor_et_max, min_cluster_energy, delta_r_cut,
                          mass_window_min, mass_window_max);
    fill_cluster_histograms(nosplit_cluster[bin], truth_alpha, *truth_eta,
                            *truth_phi, nosplit, anchor_eta_max, anchor_et_min,
                            anchor_et_max, min_cluster_energy, delta_r_cut,
                            mass_window_min, mass_window_max);
    fill_event_components(split_components[bin], truth_alpha, truth_pt,
                          *truth_daughter_pt, *truth_eta, *truth_phi, split,
                          anchor_eta_max, anchor_et_min, anchor_et_max,
                          min_cluster_energy, delta_r_cut, merged_delta_r_cut,
                          merged_response_min, merged_response_max,
                          individual_response_min, individual_response_max);
    fill_event_components(nosplit_components[bin], truth_alpha, truth_pt,
                          *truth_daughter_pt, *truth_eta, *truth_phi, nosplit,
                          anchor_eta_max, anchor_et_min, anchor_et_max,
                          min_cluster_energy, delta_r_cut, merged_delta_r_cut,
                          merged_response_min, merged_response_max,
                          individual_response_min, individual_response_max);
    fill_event_selection_stages(split_event_stages[bin], truth_alpha, split,
                                min_cluster_energy, anchor_eta_max,
                                anchor_et_min, anchor_et_max);
    fill_pair_selection_stages(split_pair_stages[bin], truth_alpha, *truth_eta,
                               *truth_phi, split, min_cluster_energy,
                               anchor_eta_max, anchor_et_min, anchor_et_max,
                               delta_r_cut);
    fill_event_selection_stages(nosplit_event_stages[bin], truth_alpha, nosplit,
                                min_cluster_energy, anchor_eta_max,
                                anchor_et_min, anchor_et_max);
    fill_pair_selection_stages(nosplit_pair_stages[bin], truth_alpha,
                               *truth_eta, *truth_phi, nosplit,
                               min_cluster_energy, anchor_eta_max,
                               anchor_et_min, anchor_et_max, delta_r_cut);
  }

  make_reco_et_retentions(split_reco_et);
  make_reco_et_retentions(nosplit_reco_et);
  for (std::size_t bin = 0; bin < truth_pt_labels.size(); ++bin) {
    make_efficiencies(split_event[bin],
                      "split_event_truth_pt_" + std::to_string(bin));
    make_efficiencies(split_cluster[bin],
                      "split_cluster_truth_pt_" + std::to_string(bin));
    make_efficiencies(nosplit_event[bin],
                      "nosplit_event_truth_pt_" + std::to_string(bin));
    make_efficiencies(nosplit_cluster[bin],
                      "nosplit_cluster_truth_pt_" + std::to_string(bin));
    make_component_fractions(split_components[bin],
                             "split_event_component_truth_pt_" +
                                 std::to_string(bin));
    make_component_fractions(nosplit_components[bin],
                             "nosplit_event_component_truth_pt_" +
                                 std::to_string(bin));
    make_truth_pi0_component_fractions(
        split_truth_pi0_components[bin],
        "split_central_truth_pi0_component_truth_pt_" + std::to_string(bin));
    make_truth_pi0_component_fractions(
        nosplit_truth_pi0_components[bin],
        "nosplit_central_truth_pi0_component_truth_pt_" + std::to_string(bin));
  }
  const std::string split_output_base =
      collection_output_base(output_base, "split");
  const std::string nosplit_output_base =
      collection_output_base(output_base, "nosplit");
  if (!make_output_directory(output_base) ||
      !make_output_directory(split_output_base) ||
      !make_output_directory(nosplit_output_base)) {
    return 5;
  }
  const auto draw_collection = [&](std::vector<StageHistograms> &event,
                                   std::vector<StageHistograms> &cluster,
                                   const std::string &collection_base,
                                   const std::string &label) {
    draw_counts(event, label, true,
                collection_base + "_event_counts_truth_pt.pdf",
                anchor_eta_max, anchor_et_min, anchor_et_max, delta_r_cut,
                mass_window_min, mass_window_max);
    draw_efficiencies(event, label, true,
                      collection_base +
                          "_event_efficiency_truth_pt.pdf",
                      anchor_eta_max, anchor_et_min, anchor_et_max, delta_r_cut,
                      mass_window_min, mass_window_max);
    draw_counts(cluster, label, false,
                collection_base + "_cluster_counts_truth_pt.pdf",
                anchor_eta_max, anchor_et_min, anchor_et_max, delta_r_cut,
                mass_window_min, mass_window_max);
    draw_efficiencies(cluster, label, false,
                      collection_base +
                          "_cluster_efficiency_truth_pt.pdf",
                      anchor_eta_max, anchor_et_min, anchor_et_max, delta_r_cut,
                      mass_window_min, mass_window_max);
  };
  draw_collection(split_event, split_cluster, split_output_base, "SPLIT");
  draw_collection(nosplit_event, nosplit_cluster, nosplit_output_base,
                  "NO_SPLIT");
  const auto draw_components =
      [anchor_eta_max, anchor_et_min, anchor_et_max, delta_r_cut,
       merged_delta_r_cut, merged_response_min, merged_response_max,
       individual_response_min, individual_response_max](
          std::vector<EventComponentHistograms> &components,
          const std::string &collection_base, const std::string &label) {
        draw_component_counts(
            components, label,
            collection_base + "_event_components_truth_pt.pdf",
            anchor_eta_max, anchor_et_min, anchor_et_max, delta_r_cut,
            merged_delta_r_cut, merged_response_min, merged_response_max,
            individual_response_min, individual_response_max);
        draw_component_fractions(
            components, label,
            collection_base +
                "_event_component_fractions_truth_pt.pdf",
            anchor_eta_max, anchor_et_min, anchor_et_max, delta_r_cut,
            merged_delta_r_cut, merged_response_min, merged_response_max,
            individual_response_min, individual_response_max);
      };
  draw_components(split_components, split_output_base, "SPLIT");
  draw_components(nosplit_components, nosplit_output_base, "NO_SPLIT");
  const auto draw_truth_pi0_components =
      [truth_eta_max, min_cluster_energy, delta_r_cut, merged_delta_r_cut,
       merged_response_min, merged_response_max, individual_response_min,
       individual_response_max](
          std::vector<TruthPi0ComponentHistograms> &components,
          const std::string &collection_base, const std::string &label) {
        draw_truth_pi0_component_counts(
            components, label,
            collection_base +
                "_central_truth_pi0_event_components_truth_pt.pdf",
            truth_eta_max, min_cluster_energy, delta_r_cut, merged_delta_r_cut,
            merged_response_min, merged_response_max, individual_response_min,
            individual_response_max);
        draw_truth_pi0_component_fractions(
            components, label,
            collection_base +
                "_central_truth_pi0_event_component_fractions_truth_pt.pdf",
            truth_eta_max, min_cluster_energy, delta_r_cut, merged_delta_r_cut,
            merged_response_min, merged_response_max, individual_response_min,
            individual_response_max);
      };
  draw_truth_pi0_components(split_truth_pi0_components, split_output_base,
                            "SPLIT");
  draw_truth_pi0_components(nosplit_truth_pi0_components, nosplit_output_base,
                            "NO_SPLIT");
  draw_cut_stages(split_event_stages, "SPLIT", false,
                  split_output_base + "_event_selection_stages_truth_pt.pdf",
                  min_cluster_energy, anchor_eta_max, anchor_et_min,
                  anchor_et_max, delta_r_cut);
  draw_cut_stages(split_pair_stages, "SPLIT", true,
                  split_output_base +
                      "_matched_pair_selection_stages_truth_pt.pdf",
                  min_cluster_energy, anchor_eta_max, anchor_et_min,
                  anchor_et_max, delta_r_cut);
  draw_cut_stages(nosplit_event_stages, "NO_SPLIT", false,
                  nosplit_output_base +
                      "_event_selection_stages_truth_pt.pdf",
                  min_cluster_energy, anchor_eta_max, anchor_et_min,
                  anchor_et_max, delta_r_cut);
  draw_cut_stages(nosplit_pair_stages, "NO_SPLIT", true,
                  nosplit_output_base +
                      "_matched_pair_selection_stages_truth_pt.pdf",
                  min_cluster_energy, anchor_eta_max, anchor_et_min,
                  anchor_et_max, delta_r_cut);
  draw_reco_et_truth_pt_stack(
      split_reco_et, "SPLIT",
      split_output_base + "_cluster_truth_pt_stack_reco_et.pdf",
      min_cluster_energy, anchor_eta_max);
  draw_reco_et_retentions(
      split_reco_et, "SPLIT",
      split_output_base + "_cluster_selection_retention_reco_et.pdf",
      min_cluster_energy, anchor_eta_max);
  draw_reco_et_truth_pt_stack(
      nosplit_reco_et, "NO_SPLIT",
      nosplit_output_base + "_cluster_truth_pt_stack_reco_et.pdf",
      min_cluster_energy, anchor_eta_max);
  draw_reco_et_retentions(
      nosplit_reco_et, "NO_SPLIT",
      nosplit_output_base + "_cluster_selection_retention_reco_et.pdf",
      min_cluster_energy, anchor_eta_max);

  TFile output((output_base + ".root").c_str(), "RECREATE");
  if (output.IsZombie()) {
    std::cerr << "Failed to create " << output_base << ".root" << std::endl;
    return 6;
  }
  write_histograms(output, split_event);
  write_histograms(output, split_cluster);
  write_histograms(output, nosplit_event);
  write_histograms(output, nosplit_cluster);
  TNamed source("input_path", input_path.c_str());
  write_component_histograms(output, split_components);
  write_component_histograms(output, nosplit_components);
  write_truth_pi0_component_histograms(output, split_truth_pi0_components);
  write_truth_pi0_component_histograms(output, nosplit_truth_pi0_components);
  write_cut_stage_histograms(output, split_event_stages);
  write_cut_stage_histograms(output, split_pair_stages);
  write_cut_stage_histograms(output, nosplit_event_stages);
  write_cut_stage_histograms(output, nosplit_pair_stages);
  write_reco_et_histograms(output, split_reco_et);
  write_reco_et_histograms(output, nosplit_reco_et);
  TParameter<double> stored_eta_max("anchor_eta_max", anchor_eta_max);
  TParameter<double> stored_truth_eta_max("truth_pi0_eta_max", truth_eta_max);
  TParameter<double> stored_et_min("anchor_et_min", anchor_et_min);
  TParameter<double> stored_et_max("anchor_et_max", anchor_et_max);
  TParameter<double> stored_delta_r("delta_r_cut", delta_r_cut);
  TParameter<double> stored_mass_min("mass_window_min", mass_window_min);
  TParameter<double> stored_mass_max("mass_window_max", mass_window_max);
  source.Write();
  TParameter<double> stored_merged_delta_r("merged_delta_r_cut",
                                           merged_delta_r_cut);
  TParameter<double> stored_merged_response_min("merged_response_min",
                                                merged_response_min);
  TParameter<double> stored_merged_response_max("merged_response_max",
                                                merged_response_max);
  TParameter<double> stored_individual_response_min("individual_response_min",
                                                    individual_response_min);
  TParameter<double> stored_individual_response_max("individual_response_max",
                                                    individual_response_max);
  stored_eta_max.Write();
  stored_truth_eta_max.Write();
  TParameter<double> stored_min_cluster_energy("analysis_min_cluster_energy",
                                               min_cluster_energy);
  stored_et_min.Write();
  stored_et_max.Write();
  stored_delta_r.Write();
  stored_mass_min.Write();
  stored_mass_max.Write();
  stored_merged_delta_r.Write();
  stored_merged_response_min.Write();
  stored_merged_response_max.Write();
  stored_individual_response_min.Write();
  stored_individual_response_max.Write();
  stored_min_cluster_energy.Write();
  output.Close();

  print_summary("SPLIT event", split_event);
  print_summary("SPLIT cluster", split_cluster);
  print_summary("NO_SPLIT event", nosplit_event);
  print_summary("NO_SPLIT cluster", nosplit_cluster);
  const bool component_closure_ok =
      print_component_summary("SPLIT component", split_components) &&
      print_component_summary("NO_SPLIT component", nosplit_components);
  const bool truth_pi0_component_closure_ok =
      print_truth_pi0_component_summary("SPLIT central truth pi0 component",
                                        split_truth_pi0_components) &&
      print_truth_pi0_component_summary("NO_SPLIT central truth pi0 component",
                                        nosplit_truth_pi0_components);
  const bool reco_et_selection_ok =
      print_reco_et_summary("SPLIT reco-ET", split_reco_et) &&
      print_reco_et_summary("NO_SPLIT reco-ET", nosplit_reco_et);
  bool cut_stage_nesting_ok = true;
  cut_stage_nesting_ok &=
      print_cut_stage_summary("SPLIT event stages", split_event_stages);
  cut_stage_nesting_ok &=
      print_cut_stage_summary("SPLIT pair stages", split_pair_stages);
  cut_stage_nesting_ok &=
      print_cut_stage_summary("NO_SPLIT event stages", nosplit_event_stages);
  cut_stage_nesting_ok &=
      print_cut_stage_summary("NO_SPLIT pair stages", nosplit_pair_stages);
  std::cout << "PlotConditionalPartnerEfficiency - invalid truth shape / truth "
               "pT outside bins = "
            << invalid_truth_shape << " / " << truth_pt_outside_bins
            << std::endl;
  std::cout << "Reco-ET malformed SPLIT / NO_SPLIT collections = "
            << split_reco_et_malformed << " / " << nosplit_reco_et_malformed
            << std::endl;
  std::cout << "Wrote " << output_base << ".root and twenty-four PDF plots"
            << std::endl;
  Long64_t malformed = invalid_truth_shape + split_reco_et_malformed +
                       nosplit_reco_et_malformed;
  for (const auto *family :
       {&split_event, &split_cluster, &nosplit_event, &nosplit_cluster}) {
    for (const StageHistograms &histograms : *family) {
      malformed += histograms.n_malformed;
    }
  }
  for (const auto *family : {&split_components, &nosplit_components}) {
    for (const EventComponentHistograms &histograms : *family) {
      malformed += histograms.n_malformed;
    }
  }
  for (const auto *family :
       {&split_truth_pi0_components, &nosplit_truth_pi0_components}) {
    for (const TruthPi0ComponentHistograms &histograms : *family) {
      malformed += histograms.n_malformed;
    }
  }
  for (const auto *family : {&split_event_stages, &split_pair_stages,
                             &nosplit_event_stages, &nosplit_pair_stages}) {
    for (const CutStageHistograms &histograms : *family) {
      malformed += histograms.n_malformed;
    }
  }
  return malformed == 0 && component_closure_ok && cut_stage_nesting_ok &&
                 reco_et_selection_ok && truth_pi0_component_closure_ok
             ? 0
             : 7;
}
