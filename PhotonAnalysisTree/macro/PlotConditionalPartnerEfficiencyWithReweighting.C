#include "Utilities/sPhenixStyle.C"

#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TH1D.h>
#include <THStack.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TNamed.h>
#include <TPad.h>
#include <TParameter.h>
#include <TString.h>
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

namespace conditional_partner_reweighting {

// const std::array<double, 13> truth_pt_edges = {
//     3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0,
//     10.0, 11.0, 12.0, 13.0, 14.0, 15.0};
const std::array<double, 39> truth_pt_edges = {
    2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
    11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0,
    19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0,
    27.0, 28.0, 29.0, 30.0, 31.0, 32.0, 33.0, 34.0,
    35.0, 36.0, 37.0, 38.0, 39.0, 40.0};
constexpr std::size_t n_truth_pt_bins = truth_pt_edges.size() - 1U;
constexpr double anchor_et_bin_width = 0.2;
constexpr double anchor_et_histogram_min = 0.0;
constexpr double anchor_et_histogram_max = 40.0;
constexpr int n_anchor_et_bins = 200;

enum class Topology : std::size_t {
  split = 0,
  merged = 1,
  missing = 2,
  other = 3,
  count = 4
};

constexpr std::size_t n_topologies = static_cast<std::size_t>(Topology::count);
const std::array<std::string, n_topologies> topology_names = {
    "separated_pair_candidate", "merged_candidate",
    "individual_cluster_no_pair", "no_matched_topology"};
const std::array<std::string, n_topologies> topology_labels = {
    "Separated", "Merged", "Missing partner", "Other"};
const std::array<int, n_topologies> topology_colors = {
    kAzure + 7, kMagenta + 1, kGreen + 2, kGray + 2};

constexpr int plot_canvas_width = 1100;
constexpr int plot_canvas_height = 900;
constexpr double plot_area_top = 0.52;
constexpr double plot_left_margin = 0.13;
constexpr double plot_right_margin = 0.04;
constexpr double plot_bottom_margin = 0.16;
constexpr double plot_top_margin = 0.04;
constexpr double plot_annotation_x = 0.06;
constexpr double plot_legend_x = 0.55;
constexpr double plot_text_size = 0.026;

std::unique_ptr<TPad> make_plot_pad(const std::string &name, const bool log_y = false) {
  auto pad = std::make_unique<TPad>(name.c_str(), "", 0.0, 0.0, 1.0, plot_area_top);
  pad->SetLeftMargin(plot_left_margin);
  pad->SetRightMargin(plot_right_margin);
  pad->SetBottomMargin(plot_bottom_margin);
  pad->SetTopMargin(plot_top_margin);
  pad->SetTicks(1, 1);
  pad->SetLogy(log_y);
  pad->Draw();
  pad->cd();
  return pad;
}

struct CollectionBranches {
  std::vector<double> *cluster_e = nullptr;
  std::vector<double> *cluster_et = nullptr;
  std::vector<double> *cluster_eta = nullptr;
  std::vector<double> *cluster_phi = nullptr;
  std::vector<unsigned char> *truth_match_valid = nullptr;
  std::vector<float> *truth_gamma0_edep = nullptr;
  std::vector<float> *truth_gamma1_edep = nullptr;
  std::vector<float> *truth_gamma0_fraction = nullptr;
  std::vector<float> *truth_gamma1_fraction = nullptr;
  std::vector<float> *truth_gamma0_recovery = nullptr;
  std::vector<float> *truth_gamma1_recovery = nullptr;
};

struct WeightedHistograms {
  std::unique_ptr<TH1D> total;
  std::array<std::unique_ptr<TH1D>, n_topologies> component;
  std::array<std::unique_ptr<TH1D>, n_topologies> fraction;
  Long64_t raw_total = 0;
  std::array<Long64_t, n_topologies> raw_component = {};
  Long64_t malformed = 0;
  double sumw = 0.0;
  double sumw2 = 0.0;
};

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
  const std::size_t slash = output_base.find_last_of('/');
  const std::string directory =
      slash == std::string::npos ? "" : output_base.substr(0, slash + 1U);
  const std::string stem =
      slash == std::string::npos ? output_base : output_base.substr(slash + 1U);
  return directory + collection + "/" + stem + "_" + collection;
}

bool valid_collection_shape(const CollectionBranches &branches) {
  return branches.cluster_e && branches.cluster_et && branches.cluster_eta &&
         branches.cluster_phi &&
         branches.cluster_e->size() == branches.cluster_et->size() &&
         branches.cluster_e->size() == branches.cluster_eta->size() &&
         branches.cluster_e->size() == branches.cluster_phi->size();
}

bool valid_contribution_shape(const CollectionBranches &branches) {
  if (!valid_collection_shape(branches) || !branches.truth_match_valid ||
      !branches.truth_gamma0_edep || !branches.truth_gamma1_edep ||
      !branches.truth_gamma0_fraction || !branches.truth_gamma1_fraction ||
      !branches.truth_gamma0_recovery || !branches.truth_gamma1_recovery) {
    return false;
  }
  const std::size_t size = branches.cluster_e->size();
  return branches.truth_match_valid->size() == size &&
         branches.truth_gamma0_edep->size() == size &&
         branches.truth_gamma1_edep->size() == size &&
         branches.truth_gamma0_fraction->size() == size &&
         branches.truth_gamma1_fraction->size() == size &&
         branches.truth_gamma0_recovery->size() == size &&
         branches.truth_gamma1_recovery->size() == size;
}

bool valid_cluster(const CollectionBranches &branches,
                   const std::size_t cluster) {
  return valid_collection_shape(branches) &&
         cluster < branches.cluster_e->size() &&
         std::isfinite(branches.cluster_e->at(cluster)) &&
         std::isfinite(branches.cluster_et->at(cluster)) &&
         std::isfinite(branches.cluster_eta->at(cluster)) &&
         std::isfinite(branches.cluster_phi->at(cluster));
}

bool has_gamma_contribution(const CollectionBranches &branches,
                            const std::size_t cluster,
                            const std::size_t gamma,
                            const double min_fraction) {
  if (!valid_contribution_shape(branches) ||
      cluster >= branches.cluster_e->size() ||
      !branches.truth_match_valid->at(cluster)) {
    return false;
  }
  const float edep = gamma == 0U ? branches.truth_gamma0_edep->at(cluster)
                                 : branches.truth_gamma1_edep->at(cluster);
  const float fraction =
      gamma == 0U ? branches.truth_gamma0_fraction->at(cluster)
                  : branches.truth_gamma1_fraction->at(cluster);
  return std::isfinite(edep) && std::isfinite(fraction) && edep > 0.0F &&
         fraction > min_fraction;
}

Topology classify_energy_contribution(
    const double truth_pt, const std::vector<double> &truth_daughter_pt,
    const CollectionBranches &branches, const double min_cluster_energy,
    const double min_contribution_fraction,
    const double merged_response_min, const double merged_response_max,
    const double individual_response_min,
    const double individual_response_max, bool &malformed) {
  malformed = false;
  if (!valid_contribution_shape(branches) ||
      truth_daughter_pt.size() != 2U || !(truth_pt > 0.0) ||
      !(truth_daughter_pt[0] > 0.0) || !(truth_daughter_pt[1] > 0.0)) {
    malformed = true;
    return Topology::other;
  }

  for (std::size_t cluster = 0; cluster < branches.cluster_e->size();
       ++cluster) {
    const std::array<float, 6> values = {
        branches.truth_gamma0_edep->at(cluster),
        branches.truth_gamma1_edep->at(cluster),
        branches.truth_gamma0_fraction->at(cluster),
        branches.truth_gamma1_fraction->at(cluster),
        branches.truth_gamma0_recovery->at(cluster),
        branches.truth_gamma1_recovery->at(cluster)};
    for (const float value : values) {
      if (!std::isfinite(value) || value < 0.0F) {
        malformed = true;
        return Topology::other;
      }
    }
  }

  constexpr std::size_t invalid_cluster = std::numeric_limits<std::size_t>::max();
  std::array<std::size_t, 2> maximum_deposit_cluster = {
      invalid_cluster, invalid_cluster};
  std::array<float, 2> maximum_deposit = {-1.0F, -1.0F};
  for (std::size_t gamma = 0; gamma < 2U; ++gamma) {
    for (std::size_t cluster = 0; cluster < branches.cluster_e->size();
         ++cluster) {
      if (!valid_cluster(branches, cluster) ||
          !has_gamma_contribution(branches, cluster, gamma,
                                  min_contribution_fraction)) {
        continue;
      }
      const float deposit =
          gamma == 0U ? branches.truth_gamma0_edep->at(cluster)
                      : branches.truth_gamma1_edep->at(cluster);
      if (deposit > maximum_deposit[gamma]) {
        maximum_deposit[gamma] = deposit;
        maximum_deposit_cluster[gamma] = cluster;
      }
    }
  }

  const bool gamma0_matched = maximum_deposit_cluster[0] != invalid_cluster;
  const bool gamma1_matched = maximum_deposit_cluster[1] != invalid_cluster;
  if (gamma0_matched && gamma1_matched &&
      maximum_deposit_cluster[0] != maximum_deposit_cluster[1] &&
      branches.cluster_e->at(maximum_deposit_cluster[0]) >=
          min_cluster_energy &&
      branches.cluster_e->at(maximum_deposit_cluster[1]) >=
          min_cluster_energy) {
    return Topology::split;
  }

  if (gamma0_matched && gamma1_matched &&
      maximum_deposit_cluster[0] == maximum_deposit_cluster[1]) {
    const std::size_t cluster = maximum_deposit_cluster[0];
    const double response = branches.cluster_et->at(cluster) / truth_pt;
    if (response >= merged_response_min &&
        response <= merged_response_max) {
      return Topology::merged;
    }
  }

  for (std::size_t gamma = 0; gamma < 2U; ++gamma) {
    const std::size_t cluster = maximum_deposit_cluster[gamma];
    if (cluster == invalid_cluster) {
      continue;
    }
    const double response = branches.cluster_et->at(cluster) / truth_daughter_pt[gamma];
    if (response >= individual_response_min &&
        response <= individual_response_max) {
      return Topology::missing;
    }
  }
  return Topology::other;
}

std::unique_ptr<TH1D> make_histogram(const std::string &name) {
  auto histogram = std::make_unique<TH1D>(name.c_str(), "", static_cast<int>(n_truth_pt_bins), truth_pt_edges.data());
  histogram->SetDirectory(nullptr);
  histogram->Sumw2();
  histogram->SetStats(false);
  return histogram;
}

std::unique_ptr<TH1D> make_anchor_et_histogram(
    const std::string &name) {
  auto histogram = std::make_unique<TH1D>(
      name.c_str(), "", n_anchor_et_bins, anchor_et_histogram_min,
      anchor_et_histogram_max);
  histogram->SetDirectory(nullptr);
  histogram->Sumw2();
  histogram->SetStats(false);
  return histogram;
}

WeightedHistograms make_histograms(const std::string &collection) {
  WeightedHistograms histograms;
  const std::string prefix =
      "h_" + collection + "_energy_contribution_central_truth_pi0_";
  histograms.total =
      make_histogram(prefix + "component_total_vs_truth_pt");
  for (std::size_t component = 0; component < n_topologies; ++component) {
    histograms.component[component] = make_histogram(
        prefix + topology_names[component] + "_vs_truth_pt");
  }
  return histograms;
}

WeightedHistograms make_anchor_et_histograms(
    const std::string &collection) {
  WeightedHistograms histograms;
  const std::string prefix =
      "h_" + collection +
      "_energy_contribution_central_truth_pi0_anchor_";
  histograms.total = make_anchor_et_histogram(
      prefix + "component_total_vs_anchor_cluster_et");
  for (std::size_t component = 0; component < n_topologies; ++component) {
    histograms.component[component] = make_anchor_et_histogram(
        prefix + topology_names[component] + "_vs_anchor_cluster_et");
  }
  return histograms;
}

double histogram_fill_pt(const double truth_pt) {
  if (truth_pt == truth_pt_edges.back()) {
    return std::nextafter(truth_pt_edges.back(), truth_pt_edges.front());
  }
  return truth_pt;
}

void fill_histograms(WeightedHistograms &histograms,
                     const double truth_pt,
                     const std::vector<double> &truth_daughter_pt,
                     const CollectionBranches &branches,
                     const double event_weight,
                     const double min_cluster_energy,
                     const double min_contribution_fraction,
                     const double merged_response_min,
                     const double merged_response_max,
                     const double individual_response_min,
                     const double individual_response_max) {
  bool malformed = false;
  const Topology topology = classify_energy_contribution(
      truth_pt, truth_daughter_pt, branches, min_cluster_energy,
      min_contribution_fraction, merged_response_min, merged_response_max,
      individual_response_min, individual_response_max, malformed);
  const std::size_t component = static_cast<std::size_t>(topology);
  const double fill_pt = histogram_fill_pt(truth_pt);
  histograms.total->Fill(fill_pt, event_weight);
  histograms.component[component]->Fill(fill_pt, event_weight);
  ++histograms.raw_total;
  ++histograms.raw_component[component];
  histograms.malformed += malformed;
  histograms.sumw += event_weight;
  histograms.sumw2 += event_weight * event_weight;
}

void fill_anchor_et_histograms(
    WeightedHistograms &histograms, const double truth_pt,
    const std::vector<double> &truth_daughter_pt,
    const CollectionBranches &branches, const double event_weight,
    const double anchor_eta_max, const double min_cluster_energy,
    const double min_contribution_fraction,
    const double merged_response_min, const double merged_response_max,
    const double individual_response_min,
    const double individual_response_max) {
  bool malformed = false;
  const Topology topology = classify_energy_contribution(
      truth_pt, truth_daughter_pt, branches, min_cluster_energy,
      min_contribution_fraction, merged_response_min, merged_response_max,
      individual_response_min, individual_response_max, malformed);
  const std::size_t component = static_cast<std::size_t>(topology);
  if (!valid_collection_shape(branches)) {
    return;
  }
  for (std::size_t cluster = 0; cluster < branches.cluster_et->size();
       ++cluster) {
    if (!valid_cluster(branches, cluster) ||
        !(std::abs(branches.cluster_eta->at(cluster)) < anchor_eta_max)) {
      continue;
    }
    const double anchor_et = branches.cluster_et->at(cluster);
    histograms.total->Fill(anchor_et, event_weight);
    histograms.component[component]->Fill(anchor_et, event_weight);
    ++histograms.raw_total;
    ++histograms.raw_component[component];
    histograms.malformed += malformed;
    histograms.sumw += event_weight;
    histograms.sumw2 += event_weight * event_weight;
  }
}

void make_fractions(WeightedHistograms &histograms,
                    const std::string &collection) {
  const std::string prefix =
      "h_" + collection + "_energy_contribution_central_truth_pi0_";
  for (std::size_t component = 0; component < n_topologies; ++component) {
    histograms.fraction[component].reset(
        static_cast<TH1D *>(histograms.component[component]->Clone(
            (prefix + topology_names[component] + "_fraction_vs_truth_pt")
                .c_str())));
    TH1D &fraction = *histograms.fraction[component];
    fraction.SetDirectory(nullptr);
    fraction.Divide(histograms.component[component].get(),
                    histograms.total.get(), 1.0, 1.0, "B");
    fraction.SetStats(false);
    fraction.SetLineColor(topology_colors[component]);
    fraction.SetMarkerColor(topology_colors[component]);
    fraction.SetMarkerStyle(20 + static_cast<int>(component));
    fraction.SetMarkerSize(0.9);
    fraction.SetLineWidth(2);
  }
}

void make_anchor_et_fractions(WeightedHistograms &histograms,
                              const std::string &collection) {
  const std::string prefix =
      "h_" + collection +
      "_energy_contribution_central_truth_pi0_anchor_";
  for (std::size_t component = 0; component < n_topologies; ++component) {
    histograms.fraction[component].reset(
        static_cast<TH1D *>(histograms.component[component]->Clone(
            (prefix + topology_names[component] +
             "_fraction_vs_anchor_cluster_et")
                .c_str())));
    TH1D &fraction = *histograms.fraction[component];
    fraction.SetDirectory(nullptr);
    fraction.Divide(histograms.component[component].get(),
                    histograms.total.get(), 1.0, 1.0, "B");
    fraction.SetStats(false);
    fraction.SetLineColor(topology_colors[component]);
    fraction.SetMarkerColor(topology_colors[component]);
    fraction.SetMarkerStyle(20 + static_cast<int>(component));
    fraction.SetMarkerSize(0.9);
    fraction.SetLineWidth(2);
  }
}

void style_histograms(WeightedHistograms &histograms) {
  histograms.total->SetLineColor(kBlue + 1);
  histograms.total->SetMarkerColor(kBlue + 1);
  histograms.total->SetMarkerStyle(kFullSquare);
  histograms.total->SetMarkerSize(0.9);
  histograms.total->SetLineWidth(2);
  for (std::size_t component = 0; component < n_topologies; ++component) {
    histograms.component[component]->SetLineColor(topology_colors[component]);
    histograms.component[component]->SetFillColor(topology_colors[component]);
    histograms.component[component]->SetLineWidth(1);
  }
}

double minimum_positive_bin(const TH1D &histogram) {
  double minimum = std::numeric_limits<double>::infinity();
  for (int bin = 1; bin <= histogram.GetNbinsX(); ++bin) {
    const double content = histogram.GetBinContent(bin);
    if (content > 0.0 && content < minimum) {
      minimum = content;
    }
  }
  return minimum;
}

void draw_histograms(WeightedHistograms &histograms,
                     const std::string &collection_label,
                     const std::string &output_path,
                     const double truth_eta_max,
                     const double min_cluster_energy,
                     const double min_contribution_fraction,
                     const double merged_response_min,
                     const double merged_response_max,
                     const double individual_response_min,
                     const double individual_response_max,
                     const double weight_exponent,
                     const double weight_reference_pt,
                     const bool log_y,
                     const bool anchor_et_axis = false,
                     const double anchor_eta_max = 0.0) {
  style_histograms(histograms);
  const std::string axis_suffix =
      anchor_et_axis ? "_anchor_et" : "_truth_pt";
  auto total = std::unique_ptr<TH1D>(
      static_cast<TH1D *>(histograms.total->Clone(
          (std::string(histograms.total->GetName()) +
           (log_y ? "_draw_log" : "_draw_linear"))
              .c_str())));
  total->SetDirectory(nullptr);
  total->GetXaxis()->SetTitle(
      anchor_et_axis ? "Anchor cluster E_{T} [GeV]"
                     : "Truth p_{T}^{#pi^{0}} [GeV]");
  total->GetYaxis()->SetTitle(
      anchor_et_axis
          ? Form("Weighted anchor clusters / %.2g GeV (a.u.)",
                 total->GetXaxis()->GetBinWidth(1))
          : "Weighted generated #pi^{0} / 1 GeV (a.u.)");
  total->GetXaxis()->CenterTitle();
  total->GetYaxis()->CenterTitle();
  const double maximum = total->GetMaximum();
  if (log_y) {
    const double minimum = minimum_positive_bin(*total);
    total->SetMinimum(std::isfinite(minimum) ? 0.5 * minimum : 1.0e-12);
    total->SetMaximum(maximum > 0.0 ? 20.0 * maximum : 1.0);
  } else {
    total->SetMinimum(0.0);
    total->SetMaximum(maximum > 0.0 ? 2.5 * maximum : 1.0);
  }

  THStack stack(
      ("stack_" + collection_label + axis_suffix +
       (log_y ? "_reweighted_log" : "_reweighted_linear"))
          .c_str(),
      "");
  for (auto &component : histograms.component) {
    stack.Add(component.get());
  }

  TCanvas canvas(
      ("c_" + collection_label + axis_suffix +
       (log_y ? "_reweighted_log" : "_reweighted_linear"))
          .c_str(),
      anchor_et_axis
          ? "Weighted anchor-cluster ET by event topology"
          : "Weighted central truth pi0 reconstruction topologies",
      plot_canvas_width, plot_canvas_height);
  canvas.SetCanvasSize(plot_canvas_width, plot_canvas_height);
  auto plot_pad = make_plot_pad(std::string(canvas.GetName()) + "_plot", log_y);
  total->Draw("E1");
  stack.Draw("HIST SAME");
  total->Draw("E1 SAME");

  canvas.cd();
  TLatex label;
  label.SetNDC();
  label.SetTextAlign(13);
  label.SetTextSize(plot_text_size);
  label.DrawLatex(plot_annotation_x, 0.96, "#it{#bf{sPHENIX}} Internal");
  label.DrawLatex(
      plot_annotation_x, 0.92, ("Single #pi^{0} gun, " + collection_label).c_str());
  if (anchor_et_axis) {
    label.DrawLatex(plot_annotation_x, 0.88, "Event category per selected anchor");
    label.DrawLatex(
        plot_annotation_x, 0.84,
        Form("Truth/anchor |#eta| < %.1f / %.1f",
             truth_eta_max, anchor_eta_max));
    label.DrawLatex(
        plot_annotation_x, 0.80,
        Form("Displayed anchor E_{T}: [%.3g, %.3g) GeV",
             anchor_et_histogram_min, anchor_et_histogram_max));
    label.DrawLatex(
        plot_annotation_x, 0.76,
        Form("Anchor weight: (p_{T}^{truth}/%.3g GeV)^{%.3f}",
             weight_reference_pt, weight_exponent));
  } else {
    label.DrawLatex(plot_annotation_x, 0.88, "Energy-deposit match; truth #alpha integrated");
    label.DrawLatex(plot_annotation_x, 0.84, Form("|#eta_{truth}^{#pi^{0}}| < %.1f; no reco cuts", truth_eta_max));
    label.DrawLatex(
        plot_annotation_x, 0.80,
        Form("w = (p_{T}^{truth}/%.3g GeV)^{%.3f}",
             weight_reference_pt, weight_exponent));
  }
  const double match_y = anchor_et_axis ? 0.72 : 0.76;
  const double split_y = anchor_et_axis ? 0.68 : 0.72;
  const double merged_y = anchor_et_axis ? 0.64 : 0.68;
  const double missing_y = anchor_et_axis ? 0.60 : 0.64;
  label.DrawLatex(plot_annotation_x, match_y, Form("Max-E_{dep} match; f_{#gamma} > %.3g", min_contribution_fraction));
  label.DrawLatex(plot_annotation_x, split_y, Form("Separated: two clusters; E_{clus} #geq %.3g GeV", min_cluster_energy));
  label.DrawLatex(
      plot_annotation_x, merged_y,
      Form("Merged: same match; R_{#pi^{0}} [%.1f, %.1f]",
           merged_response_min, merged_response_max));
  label.DrawLatex(
      plot_annotation_x, missing_y,
      Form("Missing: one match; R_{#gamma} [%.1f, %.1f]",
           individual_response_min, individual_response_max));

  TLegend legend(plot_legend_x, 0.67, 0.94, 0.95);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetTextSize(plot_text_size);
  legend.AddEntry(
      total.get(),
      anchor_et_axis ? "Central anchor clusters"
                     : "Truth #pi^{0} in acceptance", "lep");
  for (std::size_t component = 0; component < n_topologies; ++component) {
    legend.AddEntry(histograms.component[component].get(),
                    topology_labels[component].c_str(), "f");
  }
  legend.DrawClone();
  plot_pad->cd();
  plot_pad->RedrawAxis();
  canvas.cd();
  canvas.SaveAs(output_path.c_str());
}

void draw_fraction_information(
    const std::string &collection_label, const double truth_eta_max,
    const double min_cluster_energy,
    const double min_contribution_fraction,
    const double merged_response_min, const double merged_response_max,
    const double individual_response_min,
    const double individual_response_max,
    const double weight_exponent, const double weight_reference_pt,
    const bool anchor_et_axis = false,
    const double anchor_eta_max = 0.0) {
  TLatex label;
  label.SetNDC();
  label.SetTextAlign(13);
  label.SetTextSize(plot_text_size);
  label.DrawLatex(plot_annotation_x, 0.96, "#it{#bf{sPHENIX}} Internal");
  label.DrawLatex(
      plot_annotation_x, 0.92, ("Single #pi^{0} gun, " + collection_label).c_str());
  if (anchor_et_axis) {
    label.DrawLatex(plot_annotation_x, 0.88, "Event category per selected anchor");
    label.DrawLatex(plot_annotation_x, 0.84, "Denominator: weighted selected anchors");
    label.DrawLatex(
        plot_annotation_x, 0.80,
        Form("Truth/anchor |#eta| < %.1f / %.1f",
             truth_eta_max, anchor_eta_max));
    label.DrawLatex(
        plot_annotation_x, 0.76,
        Form("Displayed anchor E_{T}: [%.3g, %.3g) GeV",
             anchor_et_histogram_min, anchor_et_histogram_max));
    label.DrawLatex(
        plot_annotation_x, 0.72,
        Form("Anchor weight: (p_{T}^{truth}/%.3g GeV)^{%.3f}",
             weight_reference_pt, weight_exponent));
  } else {
    label.DrawLatex(plot_annotation_x, 0.88, "Energy-deposit match; truth #alpha integrated");
    label.DrawLatex(plot_annotation_x, 0.84, "Denominator: weighted accepted truth #pi^{0}");
    label.DrawLatex(plot_annotation_x, 0.80, Form("|#eta_{truth}^{#pi^{0}}| < %.1f; no reco cuts", truth_eta_max));
    label.DrawLatex(
        plot_annotation_x, 0.76,
        Form("w = (p_{T}^{truth}/%.3g GeV)^{%.3f}",
             weight_reference_pt, weight_exponent));
  }
  const double match_y = anchor_et_axis ? 0.68 : 0.72;
  const double split_y = anchor_et_axis ? 0.64 : 0.68;
  const double merged_y = anchor_et_axis ? 0.60 : 0.64;
  const double missing_y = anchor_et_axis ? 0.56 : 0.60;
  label.DrawLatex(plot_annotation_x, match_y, Form("Max-E_{dep} match; f_{#gamma} > %.3g", min_contribution_fraction));
  label.DrawLatex(plot_annotation_x, split_y, Form("Separated: two clusters; E_{clus} #geq %.3g GeV", min_cluster_energy));
  label.DrawLatex(
      plot_annotation_x, merged_y,
      Form("Merged: same match; R_{#pi^{0}} [%.1f, %.1f]",
           merged_response_min, merged_response_max));
  label.DrawLatex(
      plot_annotation_x, missing_y,
      Form("Missing: one match; R_{#gamma} [%.1f, %.1f]",
           individual_response_min, individual_response_max));
}

void draw_fraction_overlay(
    WeightedHistograms &histograms, const std::string &collection_label,
    const std::string &output_path, const double truth_eta_max,
    const double min_cluster_energy,
    const double min_contribution_fraction,
    const double merged_response_min, const double merged_response_max,
    const double individual_response_min,
    const double individual_response_max,
    const double weight_exponent, const double weight_reference_pt,
    const bool anchor_et_axis = false,
    const double anchor_eta_max = 0.0) {
  TH1D &frame = *histograms.fraction.front();
  frame.SetMinimum(0.0);
  frame.SetMaximum(1.05);
  frame.GetXaxis()->SetTitle(
      anchor_et_axis ? "Anchor cluster E_{T} [GeV]"
                     : "Truth p_{T}^{#pi^{0}} [GeV]");
  frame.GetYaxis()->SetTitle(
      anchor_et_axis ? "Category / selected anchor clusters"
                     : "Category / truth #pi^{0} in acceptance");
  frame.GetXaxis()->CenterTitle();
  frame.GetYaxis()->CenterTitle();

  TCanvas canvas(
      ("c_" + collection_label +
       (anchor_et_axis ? "_anchor_et" : "_truth_pt") +
       "_reweighted_category_fractions")
          .c_str(),
      anchor_et_axis ? "Weighted anchor category fractions"
                     : "Weighted truth-pi0 category fractions",
      plot_canvas_width, plot_canvas_height);
  canvas.SetCanvasSize(plot_canvas_width, plot_canvas_height);
  auto plot_pad = make_plot_pad(std::string(canvas.GetName()) + "_plot");
  frame.Draw("E1");
  for (std::size_t component = 1; component < n_topologies; ++component) {
    histograms.fraction[component]->Draw("E1 SAME");
  }

  canvas.cd();
  TLegend legend(plot_legend_x, 0.72, 0.94, 0.95);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetTextSize(plot_text_size);
  for (std::size_t component = 0; component < n_topologies; ++component) {
    legend.AddEntry(histograms.fraction[component].get(),
                    topology_labels[component].c_str(), "lep");
  }
  legend.DrawClone();
  draw_fraction_information(
      collection_label, truth_eta_max, min_cluster_energy,
      min_contribution_fraction, merged_response_min, merged_response_max,
      individual_response_min, individual_response_max, weight_exponent,
      weight_reference_pt, anchor_et_axis, anchor_eta_max);
  plot_pad->cd();
  plot_pad->RedrawAxis();
  canvas.cd();
  canvas.SaveAs(output_path.c_str());
}

void draw_fraction_stack(
    const WeightedHistograms &histograms,
    const std::string &collection_label, const std::string &output_path,
    const double truth_eta_max, const double min_cluster_energy,
    const double min_contribution_fraction,
    const double merged_response_min, const double merged_response_max,
    const double individual_response_min,
    const double individual_response_max,
    const double weight_exponent, const double weight_reference_pt,
    const bool anchor_et_axis = false,
    const double anchor_eta_max = 0.0) {
  std::array<std::unique_ptr<TH1D>, n_topologies> stacked;
  THStack stack(
      ("stack_" + collection_label +
       (anchor_et_axis ? "_anchor_et" : "_truth_pt") +
       "_reweighted_category_fractions")
          .c_str(),
      "");
  for (std::size_t component = 0; component < n_topologies; ++component) {
    stacked[component].reset(
        static_cast<TH1D *>(histograms.fraction[component]->Clone(
            (std::string(histograms.fraction[component]->GetName()) +
             "_stack_component_draw")
                .c_str())));
    stacked[component]->SetDirectory(nullptr);
    stacked[component]->SetFillColor(topology_colors[component]);
    stacked[component]->SetLineColor(kBlack);
    stacked[component]->SetLineWidth(1);
    stacked[component]->SetMarkerStyle(0);
    stack.Add(stacked[component].get());
  }

  TCanvas canvas(
      ("c_" + collection_label +
       (anchor_et_axis ? "_anchor_et" : "_truth_pt") +
       "_reweighted_category_fraction_stack")
          .c_str(),
      anchor_et_axis ? "Weighted anchor category fraction stack"
                     : "Weighted truth-pi0 category fraction stack",
      plot_canvas_width, plot_canvas_height);
  canvas.SetCanvasSize(plot_canvas_width, plot_canvas_height);
  auto plot_pad = make_plot_pad(std::string(canvas.GetName()) + "_plot");
  stack.SetMinimum(0.0);
  stack.SetMaximum(1.05);
  stack.Draw("HIST");
  stack.GetXaxis()->SetTitle(
      anchor_et_axis ? "Anchor cluster E_{T} [GeV]"
                     : "Truth p_{T}^{#pi^{0}} [GeV]");
  stack.GetYaxis()->SetTitle(
      anchor_et_axis ? "Category / selected anchor clusters"
                     : "Category / truth #pi^{0} in acceptance");
  stack.GetXaxis()->CenterTitle();
  stack.GetYaxis()->CenterTitle();

  canvas.cd();
  TLegend legend(plot_legend_x, 0.72, 0.94, 0.95);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetTextSize(plot_text_size);
  for (std::size_t component = 0; component < n_topologies; ++component) {
    legend.AddEntry(stacked[component].get(),
                    topology_labels[component].c_str(), "f");
  }
  legend.DrawClone();
  draw_fraction_information(
      collection_label, truth_eta_max, min_cluster_energy,
      min_contribution_fraction, merged_response_min, merged_response_max,
      individual_response_min, individual_response_max, weight_exponent,
      weight_reference_pt, anchor_et_axis, anchor_eta_max);
  plot_pad->cd();
  plot_pad->RedrawAxis();
  canvas.cd();
  canvas.SaveAs(output_path.c_str());
}


bool check_closure(const WeightedHistograms &histograms,
                   const std::string &label) {
  Long64_t raw_sum = 0;
  for (const Long64_t count : histograms.raw_component) {
    raw_sum += count;
  }
  bool ok = raw_sum == histograms.raw_total;
  for (int bin = 0; bin <= histograms.total->GetNbinsX() + 1;
       ++bin) {
    double component_sum = 0.0;
    for (const auto &component : histograms.component) {
      component_sum += component->GetBinContent(bin);
    }
    const double total = histograms.total->GetBinContent(bin);
    const double tolerance =
        128.0 * std::numeric_limits<double>::epsilon() *
        std::max(1.0, std::abs(total));
    if (std::abs(component_sum - total) > tolerance) {
      std::cerr << label << " weighted closure failure in bin " << bin
                << ": components/total = " << component_sum << "/"
                << total << std::endl;
      ok = false;
    }
    if (total > 0.0) {
      double fraction_sum = 0.0;
      for (const auto &fraction : histograms.fraction) {
        if (!fraction || !std::isfinite(fraction->GetBinContent(bin)) ||
            !std::isfinite(fraction->GetBinError(bin))) {
          ok = false;
          continue;
        }
        fraction_sum += fraction->GetBinContent(bin);
      }
      if (std::abs(fraction_sum - 1.0) > 1.0e-12) {
        std::cerr << label << " fraction closure failure in bin " << bin
                  << ": sum = " << fraction_sum << std::endl;
        ok = false;
      }
    }
  }
  std::cout << label << " raw/weighted/sumw2/malformed = "
            << histograms.raw_total << "/" << histograms.sumw << "/"
            << histograms.sumw2 << "/" << histograms.malformed << std::endl;
  return ok;
}

void write_histograms(TFile &output,
                      const WeightedHistograms &histograms) {
  output.cd();
  histograms.total->Write();
  for (const auto &component : histograms.component) {
    component->Write();
  }
  for (const auto &fraction : histograms.fraction) {
    fraction->Write();
  }
}

}  // namespace conditional_partner_reweighting

int PlotConditionalPartnerEfficiencyWithReweighting(
    const std::string input_path =
        // "PhotonAnalysisTree/output/merged/newscheme_100kevents_pi0_3to15GeV_etapm1_vertexpm60.root",
        "PhotonAnalysisTree/output/merged/newscheme_496p5kevents_pi0_2to40GeV_etapm0p7_vertexpm60.root",
    const std::string output_base =
        "PhotonAnalysisTree/output/plots/conditional_efficiency_reweighted/ET3to5GeV_threshold100MeV_widepT/conditional_partner_et3to5_reweighted",
    const double anchor_eta_max = 0.7, const double anchor_et_min = 3.0,
    const double anchor_et_max = 5.0,
    const double min_cluster_energy = 0.1,
    const double truth_eta_max = 0.7,
    const double min_contribution_fraction = 0.0,
    const double merged_response_min = 0.5,
    const double merged_response_max = 1.5,
    const double individual_response_min = 0.5,
    const double individual_response_max = 1.5,
    const double pt_weight_exponent = -8.122,
    const double weight_reference_pt = 3.0) {
  using namespace conditional_partner_reweighting;

  if (input_path.empty() || output_base.empty() ||
      !(anchor_eta_max > 0.0) ||
      !(min_cluster_energy >= 0.0) || !(truth_eta_max > 0.0) ||
      !(min_contribution_fraction >= 0.0 &&
        min_contribution_fraction < 1.0) ||
      !(merged_response_min >= 0.0 &&
        merged_response_min < merged_response_max) ||
      !(individual_response_min >= 0.0 &&
        individual_response_min < individual_response_max) ||
      !std::isfinite(pt_weight_exponent) ||
      !(weight_reference_pt > 0.0)) {
    std::cerr
        << "PlotConditionalPartnerEfficiencyWithReweighting - invalid argument"
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
  double truth_pt = -999.0;
  double truth_pi0_eta = -999.0;
  double truth_alpha = -999.0;
  std::vector<double> *truth_daughter_pt = nullptr;
  std::vector<double> *truth_daughter_eta = nullptr;
  std::vector<double> *truth_daughter_phi = nullptr;
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
    bind(prefix + "_cluster_truth_match_valid", &branches.truth_match_valid);
    bind(prefix + "_cluster_truth_gamma0_edep",
         &branches.truth_gamma0_edep);
    bind(prefix + "_cluster_truth_gamma1_edep",
         &branches.truth_gamma1_edep);
    bind(prefix + "_cluster_truth_gamma0_fraction",
         &branches.truth_gamma0_fraction);
    bind(prefix + "_cluster_truth_gamma1_fraction",
         &branches.truth_gamma1_fraction);
    bind(prefix + "_cluster_truth_gamma0_recovery",
         &branches.truth_gamma0_recovery);
    bind(prefix + "_cluster_truth_gamma1_recovery",
         &branches.truth_gamma1_recovery);
  };

  bind("truth_valid", &truth_valid);
  bind("truth_is_pi0_to_2gamma", &truth_is_pi0_to_2gamma);
  bind("truth_pt", &truth_pt);
  bind("truth_eta", &truth_pi0_eta);
  bind("truth_pair_e_asym", &truth_alpha);
  bind("truth_daughter_pt", &truth_daughter_pt);
  bind("truth_daughter_eta", &truth_daughter_eta);
  bind("truth_daughter_phi", &truth_daughter_phi);
  bind_collection("split", split);
  bind_collection("nosplit", nosplit);
  if (!branches_ok) {
    return 4;
  }

  WeightedHistograms split_histograms = make_histograms("split");
  WeightedHistograms nosplit_histograms = make_histograms("nosplit");
  WeightedHistograms split_anchor_et_histograms = make_anchor_et_histograms("split");
  WeightedHistograms nosplit_anchor_et_histograms = make_anchor_et_histograms("nosplit");
  Long64_t invalid_truth_shape = 0;
  Long64_t truth_pt_outside_range = 0;
  Long64_t invalid_weight = 0;

  const Long64_t entries = tree->GetEntries();
  for (Long64_t entry = 0; entry < entries; ++entry) {
    if (entry % 5000 == 0) {
      std::cout
          << "PlotConditionalPartnerEfficiencyWithReweighting - entry "
          << entry << " / " << entries << std::endl;
    }
    tree->GetEntry(entry);
    if (!truth_valid || !truth_is_pi0_to_2gamma) {
      continue;
    }
    if (!std::isfinite(truth_pt) || !std::isfinite(truth_pi0_eta) ||
        !std::isfinite(truth_alpha) || truth_alpha < 0.0 ||
        truth_alpha > 1.0 || !truth_daughter_pt || !truth_daughter_eta ||
        !truth_daughter_phi || truth_daughter_pt->size() != 2U ||
        truth_daughter_eta->size() != 2U ||
        truth_daughter_phi->size() != 2U ||
        !std::isfinite(truth_daughter_pt->at(0)) ||
        !std::isfinite(truth_daughter_pt->at(1))) {
      ++invalid_truth_shape;
      continue;
    }
    if (truth_pt < truth_pt_edges.front() ||
        truth_pt > truth_pt_edges.back()) {
      ++truth_pt_outside_range;
      continue;
    }
    if (!(std::abs(truth_pi0_eta) < truth_eta_max)) {
      continue;
    }

    const double event_weight = std::pow(truth_pt / weight_reference_pt, pt_weight_exponent);
    if (!std::isfinite(event_weight) || !(event_weight > 0.0)) {
      ++invalid_weight;
      continue;
    }
    fill_histograms(
        split_histograms, truth_pt, *truth_daughter_pt, split, event_weight,
        min_cluster_energy, min_contribution_fraction, merged_response_min,
        merged_response_max, individual_response_min, individual_response_max);
    fill_histograms(
        nosplit_histograms, truth_pt, *truth_daughter_pt, nosplit,
        event_weight, min_cluster_energy, min_contribution_fraction,
        merged_response_min, merged_response_max, individual_response_min,
        individual_response_max);
    fill_anchor_et_histograms(
        split_anchor_et_histograms, truth_pt, *truth_daughter_pt, split,
        event_weight, anchor_eta_max, min_cluster_energy,
        min_contribution_fraction, merged_response_min, merged_response_max,
        individual_response_min,
        individual_response_max);
    fill_anchor_et_histograms(
        nosplit_anchor_et_histograms, truth_pt, *truth_daughter_pt, nosplit,
        event_weight, anchor_eta_max, min_cluster_energy,
        min_contribution_fraction, merged_response_min, merged_response_max,
        individual_response_min,
        individual_response_max);
  }

  make_fractions(split_histograms, "split");
  make_fractions(nosplit_histograms, "nosplit");
  make_anchor_et_fractions(split_anchor_et_histograms, "split");
  make_anchor_et_fractions(nosplit_anchor_et_histograms, "nosplit");
  const bool split_truth_closure = check_closure(split_histograms, "SPLIT truth-pT");
  const bool nosplit_truth_closure = check_closure(nosplit_histograms, "NO-SPLIT truth-pT");
  const bool split_anchor_closure = check_closure(split_anchor_et_histograms, "SPLIT anchor-ET");
  const bool nosplit_anchor_closure = check_closure(nosplit_anchor_et_histograms, "NO-SPLIT anchor-ET");
  const bool closure_ok = split_truth_closure && nosplit_truth_closure &&
                          split_anchor_closure && nosplit_anchor_closure;
  const std::string split_output_base = collection_output_base(output_base, "split");
  const std::string nosplit_output_base = collection_output_base(output_base, "nosplit");
  if (!make_output_directory(output_base) ||
      !make_output_directory(split_output_base) ||
      !make_output_directory(nosplit_output_base)) {
    return 5;
  }

  const auto draw_collection =
      [&](WeightedHistograms &histograms,
          WeightedHistograms &anchor_et_histograms,
          const std::string &collection_label,
          const std::string &collection_base) {
        const std::string plot_base =
            collection_base +
            "_central_truth_pi0_energy_contribution_event_components_vs_truth_pt";
        draw_histograms(
            histograms, collection_label, plot_base + ".pdf", truth_eta_max,
            min_cluster_energy, min_contribution_fraction,
            merged_response_min, merged_response_max,
            individual_response_min, individual_response_max,
            pt_weight_exponent, weight_reference_pt, false);
        draw_histograms(
            histograms, collection_label, plot_base + "_logy.pdf",
            truth_eta_max, min_cluster_energy, min_contribution_fraction,
            merged_response_min, merged_response_max,
            individual_response_min, individual_response_max,
            pt_weight_exponent, weight_reference_pt, true);
        const std::string fraction_base =
            collection_base +
            "_central_truth_pi0_energy_contribution_event_component";
        draw_fraction_overlay(
            histograms, collection_label,
            fraction_base + "_fractions_vs_truth_pt.pdf", truth_eta_max,
            min_cluster_energy, min_contribution_fraction,
            merged_response_min, merged_response_max,
            individual_response_min, individual_response_max,
            pt_weight_exponent, weight_reference_pt);
        draw_fraction_stack(
            histograms, collection_label,
            fraction_base + "_fraction_stack_vs_truth_pt.pdf", truth_eta_max,
            min_cluster_energy, min_contribution_fraction,
            merged_response_min, merged_response_max,
            individual_response_min, individual_response_max,
            pt_weight_exponent, weight_reference_pt);
        const std::string anchor_plot_base =
            collection_base +
            "_central_truth_pi0_energy_contribution_event_components_vs_anchor_cluster_et";
        draw_histograms(
            anchor_et_histograms, collection_label,
            anchor_plot_base + ".pdf", truth_eta_max, min_cluster_energy,
            min_contribution_fraction, merged_response_min,
            merged_response_max, individual_response_min,
            individual_response_max, pt_weight_exponent,
            weight_reference_pt, false, true, anchor_eta_max);
        draw_histograms(
            anchor_et_histograms, collection_label,
            anchor_plot_base + "_logy.pdf", truth_eta_max,
            min_cluster_energy, min_contribution_fraction,
            merged_response_min, merged_response_max,
            individual_response_min, individual_response_max,
            pt_weight_exponent, weight_reference_pt, true, true, anchor_eta_max);
        const std::string anchor_fraction_base =
            collection_base +
            "_central_truth_pi0_energy_contribution_event_component";
        draw_fraction_overlay(
            anchor_et_histograms, collection_label,
            anchor_fraction_base + "_fractions_vs_anchor_cluster_et.pdf",
            truth_eta_max, min_cluster_energy, min_contribution_fraction,
            merged_response_min, merged_response_max,
            individual_response_min, individual_response_max,
            pt_weight_exponent, weight_reference_pt, true, anchor_eta_max);
        draw_fraction_stack(
            anchor_et_histograms, collection_label,
            anchor_fraction_base + "_fraction_stack_vs_anchor_cluster_et.pdf",
            truth_eta_max, min_cluster_energy, min_contribution_fraction,
            merged_response_min, merged_response_max,
            individual_response_min, individual_response_max,
            pt_weight_exponent, weight_reference_pt, true, anchor_eta_max);
      };
  draw_collection(split_histograms, split_anchor_et_histograms, "SPLIT",
                  split_output_base);
  draw_collection(nosplit_histograms, nosplit_anchor_et_histograms,
                  "NO_SPLIT", nosplit_output_base);

  TFile output((output_base + ".root").c_str(), "RECREATE");
  if (output.IsZombie()) {
    std::cerr << "Failed to create " << output_base << ".root" << std::endl;
    return 6;
  }
  write_histograms(output, split_histograms);
  write_histograms(output, nosplit_histograms);
  write_histograms(output, split_anchor_et_histograms);
  write_histograms(output, nosplit_anchor_et_histograms);

  TNamed("input_path", input_path.c_str()).Write();
  TNamed("weight_formula",
         "pow(truth_pt / weight_reference_pt, pt_weight_exponent)")
      .Write();
  TNamed("topology_definition",
         "maximum_energy_deposit_cluster_per_truth_pi0_daughter")
      .Write();
  TNamed("topology_mapping",
         "split=separated_pair_candidate;merged=merged_candidate;"
         "missing=individual_cluster_no_pair;other=no_matched_topology")
      .Write();
  TNamed("anchor_selection_note",
         "every valid cluster satisfying only the anchor eta cut is an anchor; "
         "there is no anchor ET selection; "
         "the collection event-level topology is assigned to every anchor; "
         "no anchor truth-contributor or recovery requirement")
      .Write();
  TParameter<double>("anchor_eta_max", anchor_eta_max).Write();
  TParameter<double>("anchor_et_min_unused", anchor_et_min).Write();
  TParameter<double>("anchor_et_max_unused", anchor_et_max).Write();
  TParameter<double>("anchor_et_histogram_min",
                     anchor_et_histogram_min)
      .Write();
  TParameter<double>("anchor_et_histogram_max",
                     anchor_et_histogram_max)
      .Write();
  TParameter<double>("anchor_et_histogram_bin_width",
                     split_anchor_et_histograms.total->GetBinWidth(1))
      .Write();
  TParameter<double>("truth_pi0_eta_max", truth_eta_max).Write();
  TParameter<double>("analysis_min_cluster_energy", min_cluster_energy).Write();
  TParameter<double>("min_energy_contribution_fraction",
                     min_contribution_fraction)
      .Write();
  TParameter<double>("merged_response_min", merged_response_min).Write();
  TParameter<double>("merged_response_max", merged_response_max).Write();
  TParameter<double>("individual_response_min", individual_response_min)
      .Write();
  TParameter<double>("individual_response_max", individual_response_max)
      .Write();
  TParameter<int>("contains_category_fractions", 1).Write();
  TParameter<int>("contains_anchor_et_spectra", 1).Write();
  TParameter<double>("pt_weight_exponent", pt_weight_exponent).Write();
  TParameter<double>("weight_reference_pt", weight_reference_pt).Write();
  TParameter<double>("split_sum_weights", split_histograms.sumw).Write();
  TParameter<double>("split_sum_weights_squared", split_histograms.sumw2)
      .Write();
  TParameter<double>("nosplit_sum_weights", nosplit_histograms.sumw).Write();
  TParameter<double>("nosplit_sum_weights_squared", nosplit_histograms.sumw2)
      .Write();
  TParameter<Long64_t>("split_raw_count", split_histograms.raw_total).Write();
  TParameter<Long64_t>("nosplit_raw_count", nosplit_histograms.raw_total)
      .Write();
  TParameter<Long64_t>("split_malformed_count", split_histograms.malformed)
      .Write();
  TParameter<Long64_t>("nosplit_malformed_count",
                       nosplit_histograms.malformed)
      .Write();
  TParameter<double>("split_anchor_sum_weights",
                     split_anchor_et_histograms.sumw)
      .Write();
  TParameter<double>("split_anchor_sum_weights_squared",
                     split_anchor_et_histograms.sumw2)
      .Write();
  TParameter<double>("nosplit_anchor_sum_weights",
                     nosplit_anchor_et_histograms.sumw)
      .Write();
  TParameter<double>("nosplit_anchor_sum_weights_squared",
                     nosplit_anchor_et_histograms.sumw2)
      .Write();
  TParameter<Long64_t>("split_anchor_raw_count",
                       split_anchor_et_histograms.raw_total)
      .Write();
  TParameter<Long64_t>("nosplit_anchor_raw_count",
                       nosplit_anchor_et_histograms.raw_total)
      .Write();
  TParameter<Long64_t>("split_anchor_malformed_count",
                       split_anchor_et_histograms.malformed)
      .Write();
  TParameter<Long64_t>("nosplit_anchor_malformed_count",
                       nosplit_anchor_et_histograms.malformed)
      .Write();
  TParameter<Long64_t>("invalid_truth_shape_count", invalid_truth_shape)
      .Write();
  TParameter<Long64_t>("truth_pt_outside_range_count",
                       truth_pt_outside_range)
      .Write();
  TParameter<Long64_t>("invalid_weight_count", invalid_weight).Write();
  output.Close();

  std::cout
      << "PlotConditionalPartnerEfficiencyWithReweighting - invalid truth / "
         "outside pT / invalid weight = "
      << invalid_truth_shape << " / " << truth_pt_outside_range << " / "
      << invalid_weight << std::endl;
  std::cout << "Wrote " << output_base
            << ".root and sixteen weighted PDF plots" << std::endl;
  return closure_ok && invalid_truth_shape == 0 &&
                 truth_pt_outside_range == 0 && invalid_weight == 0 &&
                 split_histograms.malformed == 0 &&
                 nosplit_histograms.malformed == 0 &&
                 split_anchor_et_histograms.malformed == 0 &&
                 nosplit_anchor_et_histograms.malformed == 0
             ? 0
             : 7;
}
