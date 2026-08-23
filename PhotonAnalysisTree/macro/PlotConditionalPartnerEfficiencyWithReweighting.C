#include "Utilities/sPhenixStyle.C"

#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TH1D.h>
#include <THStack.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace conditional_partner_reweighting {

const std::array<double, 13> truth_pt_edges = {
    3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0,
    10.0, 11.0, 12.0, 13.0, 14.0, 15.0};
constexpr std::size_t n_truth_pt_bins = truth_pt_edges.size() - 1U;

enum class Topology : std::size_t {
  split = 0,
  merged = 1,
  missing = 2,
  other = 3,
  count = 4
};

constexpr std::size_t n_topologies =
    static_cast<std::size_t>(Topology::count);
const std::array<std::string, n_topologies> topology_names = {
    "separated_pair_candidate", "merged_candidate",
    "individual_cluster_no_pair", "no_matched_topology"};
const std::array<std::string, n_topologies> topology_labels = {
    "Split (two separated clusters)", "Merged",
    "Missing (individual cluster only)", "Other"};
const std::array<int, n_topologies> topology_colors = {
    kBlue + 1, kMagenta - 3, kOrange + 7, kGray + 1};

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

  constexpr std::size_t invalid_cluster =
      std::numeric_limits<std::size_t>::max();
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
    const double response =
        branches.cluster_et->at(cluster) / truth_daughter_pt[gamma];
    if (response >= individual_response_min &&
        response <= individual_response_max) {
      return Topology::missing;
    }
  }
  return Topology::other;
}

std::unique_ptr<TH1D> make_histogram(const std::string &name) {
  auto histogram = std::make_unique<TH1D>(
      name.c_str(), "", static_cast<int>(n_truth_pt_bins),
      truth_pt_edges.data());
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

void style_histograms(WeightedHistograms &histograms) {
  histograms.total->SetLineColor(kGreen + 2);
  histograms.total->SetMarkerColor(kGreen + 2);
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
                     const bool log_y) {
  style_histograms(histograms);
  auto total = std::unique_ptr<TH1D>(
      static_cast<TH1D *>(histograms.total->Clone(
          (std::string(histograms.total->GetName()) +
           (log_y ? "_draw_log" : "_draw_linear"))
              .c_str())));
  total->SetDirectory(nullptr);
  total->GetXaxis()->SetTitle("Truth p_{T}^{#pi^{0}} [GeV]");
  total->GetYaxis()->SetTitle(
      "Weighted generated #pi^{0} / 1 GeV (a.u.)");
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
      ("stack_" + collection_label +
       (log_y ? "_reweighted_log" : "_reweighted_linear"))
          .c_str(),
      "");
  for (auto &component : histograms.component) {
    stack.Add(component.get());
  }

  TCanvas canvas(
      ("c_" + collection_label +
       (log_y ? "_reweighted_log" : "_reweighted_linear"))
          .c_str(),
      "Weighted central truth pi0 reconstruction topologies", 1400, 800);
  canvas.SetLeftMargin(0.14);
  canvas.SetRightMargin(0.36);
  canvas.SetBottomMargin(0.12);
  canvas.SetTopMargin(0.04);
  canvas.SetLogy(log_y);
  total->Draw("E1");
  stack.Draw("HIST SAME");
  total->Draw("E1 SAME");

  TLatex label;
  label.SetNDC();
  label.SetTextAlign(13);
  label.SetTextSize(0.021);
  label.DrawLatex(0.67, 0.93, "#it{#bf{sPHENIX}} Internal");
  label.DrawLatex(
      0.67, 0.88, ("Single #pi^{0} gun, " + collection_label).c_str());
  label.DrawLatex(
      0.67, 0.83, "Energy-deposit matching; truth #alpha integrated");
  label.DrawLatex(
      0.67, 0.78,
      Form("|#eta_{truth}^{#pi^{0}}| < %.1f; no reco #eta or E_{T} cut",
           truth_eta_max));
  label.DrawLatex(
      0.67, 0.73,
      Form("Weight: (p_{T}^{truth}/%.3g GeV)^{%.3f}",
           weight_reference_pt, weight_exponent));
  label.DrawLatex(
      0.67, 0.65,
      Form("Maximum-E_{dep} match; f_{#gamma} > %.3g",
           min_contribution_fraction));
  label.DrawLatex(
      0.67, 0.61,
      Form("Split: distinct matches; E_{cluster} #geq %.3g GeV",
           min_cluster_energy));
  label.DrawLatex(
      0.67, 0.57,
      Form("Merged: same match; response [%.1f, %.1f]",
           merged_response_min, merged_response_max));
  label.DrawLatex(
      0.67, 0.53,
      Form("Missing: one match; response [%.1f, %.1f]",
           individual_response_min, individual_response_max));

  TLegend legend(0.67, 0.14, 0.98, 0.43);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetTextSize(0.022);
  legend.AddEntry(total.get(), "Truth #pi^{0} in acceptance", "lep");
  for (std::size_t component = 0; component < n_topologies; ++component) {
    legend.AddEntry(histograms.component[component].get(),
                    topology_labels[component].c_str(), "f");
  }
  legend.DrawClone();
  gPad->RedrawAxis();
  canvas.SaveAs(output_path.c_str());
}

void draw_fraction_information(
    const std::string &collection_label, const double truth_eta_max,
    const double min_cluster_energy,
    const double min_contribution_fraction,
    const double merged_response_min, const double merged_response_max,
    const double individual_response_min,
    const double individual_response_max,
    const double weight_exponent, const double weight_reference_pt) {
  TLatex label;
  label.SetNDC();
  label.SetTextAlign(13);
  label.SetTextSize(0.021);
  label.DrawLatex(0.67, 0.93, "#it{#bf{sPHENIX}} Internal");
  label.DrawLatex(
      0.67, 0.88, ("Single #pi^{0} gun, " + collection_label).c_str());
  label.DrawLatex(
      0.67, 0.83, "Energy-deposit matching; truth #alpha integrated");
  label.DrawLatex(
      0.67, 0.78, "Denominator: weighted truth #pi^{0} in acceptance");
  label.DrawLatex(
      0.67, 0.73,
      Form("|#eta_{truth}^{#pi^{0}}| < %.1f; no reco #eta or E_{T} cut",
           truth_eta_max));
  label.DrawLatex(
      0.67, 0.68,
      Form("Weight: (p_{T}^{truth}/%.3g GeV)^{%.3f}",
           weight_reference_pt, weight_exponent));
  label.DrawLatex(
      0.67, 0.61,
      Form("Maximum-E_{dep} match; f_{#gamma} > %.3g",
           min_contribution_fraction));
  label.DrawLatex(
      0.67, 0.56,
      Form("Split: distinct matches; E_{cluster} #geq %.3g GeV",
           min_cluster_energy));
  label.DrawLatex(
      0.67, 0.51,
      Form("Merged: same match; response [%.1f, %.1f]",
           merged_response_min, merged_response_max));
  label.DrawLatex(
      0.67, 0.46,
      Form("Missing: one match; response [%.1f, %.1f]",
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
    const double weight_exponent, const double weight_reference_pt) {
  TH1D &frame = *histograms.fraction.front();
  frame.SetMinimum(0.0);
  frame.SetMaximum(1.05);
  frame.GetXaxis()->SetTitle("Truth p_{T}^{#pi^{0}} [GeV]");
  frame.GetYaxis()->SetTitle("Category / truth #pi^{0} in acceptance");
  frame.GetXaxis()->CenterTitle();
  frame.GetYaxis()->CenterTitle();

  TCanvas canvas(
      ("c_" + collection_label + "_reweighted_category_fractions").c_str(),
      "Weighted truth-pi0 category fractions", 1400, 800);
  canvas.SetLeftMargin(0.14);
  canvas.SetRightMargin(0.36);
  canvas.SetBottomMargin(0.12);
  canvas.SetTopMargin(0.04);
  frame.Draw("E1");
  for (std::size_t component = 1; component < n_topologies; ++component) {
    histograms.fraction[component]->Draw("E1 SAME");
  }

  TLegend legend(0.67, 0.14, 0.98, 0.40);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetTextSize(0.022);
  for (std::size_t component = 0; component < n_topologies; ++component) {
    legend.AddEntry(histograms.fraction[component].get(),
                    topology_labels[component].c_str(), "lep");
  }
  legend.DrawClone();
  draw_fraction_information(
      collection_label, truth_eta_max, min_cluster_energy,
      min_contribution_fraction, merged_response_min, merged_response_max,
      individual_response_min, individual_response_max, weight_exponent,
      weight_reference_pt);
  gPad->RedrawAxis();
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
    const double weight_exponent, const double weight_reference_pt) {
  std::array<std::unique_ptr<TH1D>, n_topologies> stacked;
  THStack stack(
      ("stack_" + collection_label + "_reweighted_category_fractions").c_str(),
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
      ("c_" + collection_label + "_reweighted_category_fraction_stack")
          .c_str(),
      "Weighted truth-pi0 category fraction stack", 1400, 800);
  canvas.SetLeftMargin(0.14);
  canvas.SetRightMargin(0.36);
  canvas.SetBottomMargin(0.12);
  canvas.SetTopMargin(0.04);
  stack.SetMinimum(0.0);
  stack.SetMaximum(1.05);
  stack.Draw("HIST");
  stack.GetXaxis()->SetTitle("Truth p_{T}^{#pi^{0}} [GeV]");
  stack.GetYaxis()->SetTitle("Category / truth #pi^{0} in acceptance");
  stack.GetXaxis()->CenterTitle();
  stack.GetYaxis()->CenterTitle();

  TLegend legend(0.67, 0.14, 0.98, 0.40);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetTextSize(0.022);
  for (std::size_t component = 0; component < n_topologies; ++component) {
    legend.AddEntry(stacked[component].get(),
                    topology_labels[component].c_str(), "f");
  }
  legend.DrawClone();
  draw_fraction_information(
      collection_label, truth_eta_max, min_cluster_energy,
      min_contribution_fraction, merged_response_min, merged_response_max,
      individual_response_min, individual_response_max, weight_exponent,
      weight_reference_pt);
  gPad->RedrawAxis();
  canvas.SaveAs(output_path.c_str());
}


bool check_closure(const WeightedHistograms &histograms,
                   const std::string &label) {
  Long64_t raw_sum = 0;
  for (const Long64_t count : histograms.raw_component) {
    raw_sum += count;
  }
  bool ok = raw_sum == histograms.raw_total;
  for (int bin = 1; bin <= histograms.total->GetNbinsX(); ++bin) {
    double component_sum = 0.0;
    for (const auto &component : histograms.component) {
      component_sum += component->GetBinContent(bin);
    }
    const double total = histograms.total->GetBinContent(bin);
    const double tolerance =
        32.0 * std::numeric_limits<double>::epsilon() *
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
        "PhotonAnalysisTree/output/merged/"
        "newscheme_100kevents_pi0_3to15GeV_etapm1_vertexpm60.root",
    const std::string output_base =
        "PhotonAnalysisTree/output/plots/conditional_efficiency_reweighted/"
        "ET3to5GeV_threshold100MeV/conditional_partner_et3to5_reweighted",
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
      !(anchor_et_min >= 0.0 && anchor_et_min < anchor_et_max) ||
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

    const double event_weight =
        std::pow(truth_pt / weight_reference_pt, pt_weight_exponent);
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
  }

  make_fractions(split_histograms, "split");
  make_fractions(nosplit_histograms, "nosplit");
  const bool closure_ok =
      check_closure(split_histograms, "SPLIT") &&
      check_closure(nosplit_histograms, "NO_SPLIT");
  const std::string split_output_base =
      collection_output_base(output_base, "split");
  const std::string nosplit_output_base =
      collection_output_base(output_base, "nosplit");
  if (!make_output_directory(output_base) ||
      !make_output_directory(split_output_base) ||
      !make_output_directory(nosplit_output_base)) {
    return 5;
  }

  const auto draw_collection =
      [&](WeightedHistograms &histograms,
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
      };
  draw_collection(split_histograms, "SPLIT", split_output_base);
  draw_collection(nosplit_histograms, "NO_SPLIT", nosplit_output_base);

  TFile output((output_base + ".root").c_str(), "RECREATE");
  if (output.IsZombie()) {
    std::cerr << "Failed to create " << output_base << ".root" << std::endl;
    return 6;
  }
  write_histograms(output, split_histograms);
  write_histograms(output, nosplit_histograms);

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
         "anchor eta and ET arguments are retained for interface compatibility "
         "but are not applied to this truth-pi0-denominator plot")
      .Write();
  TParameter<double>("anchor_eta_max_unused", anchor_eta_max).Write();
  TParameter<double>("anchor_et_min_unused", anchor_et_min).Write();
  TParameter<double>("anchor_et_max_unused", anchor_et_max).Write();
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
            << ".root and eight weighted PDF plots" << std::endl;
  return closure_ok && invalid_truth_shape == 0 &&
                 truth_pt_outside_range == 0 && invalid_weight == 0 &&
                 split_histograms.malformed == 0 &&
                 nosplit_histograms.malformed == 0
             ? 0
             : 7;
}
