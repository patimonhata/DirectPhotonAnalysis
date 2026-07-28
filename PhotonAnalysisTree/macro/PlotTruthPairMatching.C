#include "Utilities/sPhenixStyle.C"

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
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

struct CollectionBranches {
  std::vector<double> *cluster_e = nullptr;
  std::vector<double> *cluster_eta = nullptr;
  std::vector<double> *cluster_phi = nullptr;
  std::vector<unsigned int> *pair_i = nullptr;
  std::vector<unsigned int> *pair_j = nullptr;
  std::vector<double> *pair_mass = nullptr;
};

struct MatchResult {
  bool valid = false;
  std::size_t cluster_for_gamma0 = 0;
  std::size_t cluster_for_gamma1 = 0;
  double delta_r0 = std::numeric_limits<double>::infinity();
  double delta_r1 = std::numeric_limits<double>::infinity();
};

struct CollectionHistograms {
  std::unique_ptr<TH1D> truth;
  std::unique_ptr<TH1D> matched;
  std::unique_ptr<TH1D> matched_mass_window;
  std::unique_ptr<TH1D> efficiency_matched;
  std::unique_ptr<TH1D> efficiency_mass_window;
  std::unique_ptr<TH1D> delta_r_each;
  std::unique_ptr<TH1D> delta_r_max;
  std::unique_ptr<TH1D> matched_mass;

  Long64_t truth_events = 0;
  Long64_t events_with_two_clusters = 0;
  Long64_t matched_events = 0;
  Long64_t mass_window_events = 0;
  Long64_t malformed_events = 0;
  Long64_t delta_r_overflow_events = 0;
};

const std::vector<double> truth_pt_edges = {5.0, 7.0, 9.0, 11.0, 13.0, 15.0};
const std::array<std::string, 5> truth_pt_labels = {
    "5 #leq p_{T}^{truth} < 7 GeV", "7 #leq p_{T}^{truth} < 9 GeV",
    "9 #leq p_{T}^{truth} < 11 GeV", "11 #leq p_{T}^{truth} < 13 GeV",
    "13 #leq p_{T}^{truth} #leq 15 GeV"};

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

CollectionHistograms make_histograms(const std::string &prefix,
                                     const int asymmetry_nbins) {
  CollectionHistograms histograms;
  const auto make_asymmetry_histogram = [&](const std::string &name) {
    auto histogram = std::make_unique<TH1D>(
        ("h_" + prefix + "_" + name).c_str(), "", asymmetry_nbins, 0.0, 1.0);
    histogram->SetDirectory(nullptr);
    histogram->Sumw2();
    histogram->SetStats(false);
    return histogram;
  };

  histograms.truth = make_asymmetry_histogram("truth");
  histograms.matched = make_asymmetry_histogram("matched");
  histograms.matched_mass_window =
      make_asymmetry_histogram("matched_mass_window");

  histograms.delta_r_each = std::make_unique<TH1D>(
      ("h_" + prefix + "_delta_r_each").c_str(), "", 200, 0.0, 1.0);
  histograms.delta_r_max = std::make_unique<TH1D>(
      ("h_" + prefix + "_delta_r_max").c_str(), "", 200, 0.0, 1.0);
  histograms.matched_mass = std::make_unique<TH1D>(
      ("h_" + prefix + "_matched_mass").c_str(), "", 150, 0.0, 0.30);
  for (TH1D *histogram :
       {histograms.delta_r_each.get(), histograms.delta_r_max.get(),
        histograms.matched_mass.get()}) {
    histogram->SetDirectory(nullptr);
    histogram->Sumw2();
    histogram->SetStats(false);
  }
  return histograms;
}

std::vector<CollectionHistograms>
make_truth_pt_histograms(const std::string &prefix, const int asymmetry_nbins) {
  std::vector<CollectionHistograms> histograms;
  histograms.reserve(truth_pt_labels.size());
  for (std::size_t bin = 0; bin < truth_pt_labels.size(); ++bin) {
    histograms.push_back(make_histograms(
        prefix + "_truth_pt_" + std::to_string(bin), asymmetry_nbins));
  }
  return histograms;
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

MatchResult match_truth_to_clusters(const std::vector<double> &truth_eta,
                                    const std::vector<double> &truth_phi,
                                    const CollectionBranches &branches,
                                    const double min_cluster_energy) {
  MatchResult result;
  if (truth_eta.size() != 2U || truth_phi.size() != 2U || !branches.cluster_e ||
      !branches.cluster_eta || !branches.cluster_phi) {
    return result;
  }

  const std::size_t ncluster = branches.cluster_e->size();
  if (branches.cluster_eta->size() != ncluster ||
      branches.cluster_phi->size() != ncluster) {
    return result;
  }

  std::vector<std::size_t> selected_clusters;
  selected_clusters.reserve(ncluster);
  for (std::size_t cluster = 0; cluster < ncluster; ++cluster) {
    const double energy = branches.cluster_e->at(cluster);
    const double eta = branches.cluster_eta->at(cluster);
    const double phi = branches.cluster_phi->at(cluster);
    if (std::isfinite(energy) && energy >= min_cluster_energy &&
        std::isfinite(eta) && std::isfinite(phi)) {
      selected_clusters.push_back(cluster);
    }
  }

  if (selected_clusters.size() < 2U) {
    return result;
  }

  double best_cost = std::numeric_limits<double>::infinity();
  for (const std::size_t cluster0 : selected_clusters) {
    const double delta_r0 =
        delta_r(truth_eta[0], truth_phi[0], branches.cluster_eta->at(cluster0),
                branches.cluster_phi->at(cluster0));
    for (const std::size_t cluster1 : selected_clusters) {
      if (cluster0 == cluster1) {
        continue;
      }

      const double delta_r1 = delta_r(truth_eta[1], truth_phi[1],
                                      branches.cluster_eta->at(cluster1),
                                      branches.cluster_phi->at(cluster1));
      const double cost = delta_r0 * delta_r0 + delta_r1 * delta_r1;
      if (cost < best_cost) {
        best_cost = cost;
        result.valid = true;
        result.cluster_for_gamma0 = cluster0;
        result.cluster_for_gamma1 = cluster1;
        result.delta_r0 = delta_r0;
        result.delta_r1 = delta_r1;
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

void make_efficiencies(CollectionHistograms &histograms,
                       const std::string &prefix) {
  histograms.efficiency_matched.reset(
      static_cast<TH1D *>(histograms.matched->Clone(
          ("h_" + prefix + "_efficiency_matched").c_str())));
  histograms.efficiency_mass_window.reset(
      static_cast<TH1D *>(histograms.matched_mass_window->Clone(
          ("h_" + prefix + "_efficiency_mass_window").c_str())));
  histograms.efficiency_matched->SetDirectory(nullptr);
  histograms.efficiency_mass_window->SetDirectory(nullptr);
  histograms.efficiency_matched->Divide(histograms.matched.get(),
                                        histograms.truth.get(), 1.0, 1.0, "B");
  histograms.efficiency_mass_window->Divide(
      histograms.matched_mass_window.get(), histograms.truth.get(), 1.0, 1.0,
      "B");
}

void style_count_histograms(CollectionHistograms &histograms) {
  histograms.truth->SetLineColor(kBlack);
  histograms.truth->SetMarkerColor(kBlack);
  histograms.truth->SetMarkerStyle(20);
  histograms.matched->SetLineColor(kBlue + 1);
  histograms.matched->SetMarkerColor(kBlue + 1);
  histograms.matched->SetMarkerStyle(21);
  histograms.matched_mass_window->SetLineColor(kRed + 1);
  histograms.matched_mass_window->SetMarkerColor(kRed + 1);
  histograms.matched_mass_window->SetMarkerStyle(22);
  for (TH1D *histogram : {histograms.truth.get(), histograms.matched.get(),
                          histograms.matched_mass_window.get()}) {
    histogram->SetLineWidth(2);
    histogram->SetMarkerSize(0.8);
    histogram->GetXaxis()->SetTitle(
        "Truth energy asymmetry |E_{1}-E_{2}|/(E_{1}+E_{2})");
    histogram->GetYaxis()->SetTitle("Events (truth pairs)");
  }
}

void draw_labels(const std::string &collection_label, const double delta_r_cut,
                 const double mass_min, const double mass_max,
                 const double min_cluster_energy) {
  TLatex label;
  label.SetNDC();
  label.SetTextAlign(13);
  label.SetTextSize(0.035);
  label.DrawLatex(0.18, 0.92, "#it{#bf{sPHENIX}} Internal");
  label.DrawLatex(0.18, 0.87,
                  ("Single #pi^{0} gun, " + collection_label).c_str());
  label.DrawLatex(0.18, 0.82,
                  Form("max(#DeltaR_{0},#DeltaR_{1}) < %.3f", delta_r_cut));
  label.DrawLatex(
      0.18, 0.77,
      Form("%.2f #leq m_{#gamma#gamma} #leq %.2f GeV", mass_min, mass_max));
  label.DrawLatex(0.18, 0.72,
                  Form("E_{cluster} #geq %.3g GeV", min_cluster_energy));
}

void draw_count_plot(CollectionHistograms &histograms,
                     const std::string &collection_label,
                     const std::string &output_path, const double delta_r_cut,
                     const double mass_min, const double mass_max,
                     const double min_cluster_energy) {
  style_count_histograms(histograms);
  const double maximum = histograms.truth->GetMaximum();
  histograms.truth->SetMinimum(0.0);
  histograms.truth->SetMaximum(maximum > 0.0 ? 1.35 * maximum : 1.0);

  TCanvas canvas(("c_" + collection_label + "_counts").c_str(),
                 "Truth pair matching counts", 1000, 800);
  histograms.truth->Draw("E1");
  histograms.matched->Draw("E1 SAME");
  histograms.matched_mass_window->Draw("E1 SAME");

  TLegend legend(0.53, 0.72, 0.91, 0.89);
  legend.AddEntry(histograms.truth.get(),
                  "Truth #pi^{0}#rightarrow#gamma#gamma in acceptance", "lep");
  legend.AddEntry(histograms.matched.get(), "#DeltaR-matched pair", "lep");
  legend.AddEntry(histograms.matched_mass_window.get(),
                  "#DeltaR-matched + mass window", "lep");
  legend.Draw();
  draw_labels(collection_label, delta_r_cut, mass_min, mass_max,
              min_cluster_energy);
  canvas.RedrawAxis();
  canvas.SaveAs(output_path.c_str());
}

void draw_efficiency_plot(CollectionHistograms &histograms,
                          const std::string &collection_label,
                          const std::string &output_path,
                          const double delta_r_cut, const double mass_min,
                          const double mass_max,
                          const double min_cluster_energy) {
  TH1D *matched = histograms.efficiency_matched.get();
  TH1D *mass_window = histograms.efficiency_mass_window.get();
  matched->SetStats(false);
  matched->SetLineColor(kBlue + 1);
  matched->SetMarkerColor(kBlue + 1);
  matched->SetMarkerStyle(21);
  mass_window->SetStats(false);
  mass_window->SetLineColor(kRed + 1);
  mass_window->SetMarkerColor(kRed + 1);
  mass_window->SetMarkerStyle(22);
  for (TH1D *histogram : {matched, mass_window}) {
    histogram->SetLineWidth(2);
    histogram->SetMarkerSize(0.8);
    histogram->GetXaxis()->SetTitle(
        "Truth energy asymmetry |E_{1}-E_{2}|/(E_{1}+E_{2})");
    histogram->GetYaxis()->SetTitle("Fraction of accepted truth pairs");
    histogram->SetMinimum(0.0);
    histogram->SetMaximum(1.05);
  }

  TCanvas canvas(("c_" + collection_label + "_efficiency").c_str(),
                 "Truth pair matching efficiency", 1000, 800);
  matched->Draw("E1");
  mass_window->Draw("E1 SAME");
  TLegend legend(0.56, 0.76, 0.91, 0.89);
  legend.AddEntry(matched, "#DeltaR-matched / truth", "lep");
  legend.AddEntry(mass_window, "(#DeltaR-matched + mass window) / truth",
                  "lep");
  legend.Draw();
  draw_labels(collection_label, delta_r_cut, mass_min, mass_max,
              min_cluster_energy);
  canvas.RedrawAxis();
  canvas.SaveAs(output_path.c_str());
}

void draw_delta_r_plot(CollectionHistograms &histograms,
                       const std::string &collection_label,
                       const std::string &output_path, const double delta_r_cut,
                       const double min_cluster_energy) {
  TH1D *each = histograms.delta_r_each.get();
  TH1D *maximum = histograms.delta_r_max.get();
  each->SetLineColor(kGray + 2);
  each->SetLineStyle(2);
  each->SetLineWidth(2);
  maximum->SetLineColor(kBlue + 1);
  maximum->SetLineWidth(2);
  maximum->GetXaxis()->SetTitle("#DeltaR(truth #gamma, matched cluster)");
  maximum->GetYaxis()->SetTitle("Entries / 0.005");
  maximum->GetXaxis()->SetRangeUser(0.0, 0.20);
  maximum->SetMinimum(0.5);
  const double plot_maximum =
      std::max(each->GetMaximum(), maximum->GetMaximum());
  maximum->SetMaximum(plot_maximum > 0.0 ? 3.0 * plot_maximum : 1.0);

  TCanvas canvas(("c_" + collection_label + "_delta_r").c_str(),
                 "Truth-cluster delta R", 1000, 800);
  canvas.SetLogy();
  maximum->Draw("HIST");
  each->Draw("HIST SAME");
  TLine cut_line(delta_r_cut, 0.5, delta_r_cut, maximum->GetMaximum());
  cut_line.SetLineColor(kRed + 1);
  cut_line.SetLineStyle(7);
  cut_line.SetLineWidth(2);
  cut_line.Draw();

  TLegend legend(0.54, 0.72, 0.91, 0.89);
  legend.AddEntry(maximum, "max(#DeltaR_{0}, #DeltaR_{1}); one/event", "l");
  legend.AddEntry(each, "Individual assignments; two/event", "l");
  legend.AddEntry(&cut_line, Form("#DeltaR cut = %.3f", delta_r_cut), "l");
  legend.Draw();

  TLatex label;
  label.SetNDC();
  label.SetTextAlign(13);
  label.SetTextSize(0.035);
  label.DrawLatex(0.18, 0.92, "#it{#bf{sPHENIX}} Internal");
  label.DrawLatex(0.18, 0.87,
                  ("Single #pi^{0} gun, " + collection_label).c_str());
  label.DrawLatex(0.18, 0.82,
                  Form("E_{cluster} #geq %.3g GeV, before #DeltaR cut",
                       min_cluster_energy));
  canvas.RedrawAxis();
  canvas.SaveAs(output_path.c_str());
}

void draw_truth_pt_information_panel(CollectionHistograms &histograms,
                                     const std::string &collection_label,
                                     const std::string &plot_type,
                                     const double delta_r_cut,
                                     const double mass_min,
                                     const double mass_max,
                                     const double min_cluster_energy) {
  TLatex label;
  label.SetNDC();
  label.SetTextAlign(13);
  label.SetTextSize(0.045);
  label.DrawLatex(0.15, 0.91, "#it{#bf{sPHENIX}} Internal");
  label.DrawLatex(0.15, 0.83,
                  ("Single #pi^{0} gun, " + collection_label).c_str());
  label.DrawLatex(0.15, 0.75, "Binned in truth #pi^{0} p_{T}");
  label.DrawLatex(0.15, 0.67,
                  Form("max(#DeltaR_{0},#DeltaR_{1}) < %.3f", delta_r_cut));
  label.DrawLatex(
      0.15, 0.59,
      Form("%.2f #leq m_{#gamma#gamma} #leq %.2f GeV", mass_min, mass_max));
  label.DrawLatex(0.15, 0.51,
                  Form("E_{cluster} #geq %.3g GeV", min_cluster_energy));

  TLegend legend(0.15, 0.18, 0.91, 0.43);
  if (plot_type == "counts") {
    legend.AddEntry(histograms.truth.get(),
                    "Truth #pi^{0}#rightarrow#gamma#gamma in acceptance",
                    "lep");
    legend.AddEntry(histograms.matched.get(), "#DeltaR-matched pair", "lep");
    legend.AddEntry(histograms.matched_mass_window.get(),
                    "#DeltaR-matched + mass window", "lep");
  } else if (plot_type == "efficiency") {
    legend.AddEntry(histograms.efficiency_matched.get(),
                    "#DeltaR-matched / truth", "lep");
    legend.AddEntry(histograms.efficiency_mass_window.get(),
                    "(#DeltaR-matched + mass window) / truth", "lep");
  } else {
    legend.AddEntry(histograms.delta_r_max.get(),
                    "max(#DeltaR_{0}, #DeltaR_{1}); one/event", "l");
    legend.AddEntry(histograms.delta_r_each.get(),
                    "Individual assignments; two/event", "l");
  }
  legend.DrawClone();
}

void draw_truth_pt_counts(std::vector<CollectionHistograms> &histograms,
                          const std::string &collection_label,
                          const std::string &output_path,
                          const double delta_r_cut, const double mass_min,
                          const double mass_max,
                          const double min_cluster_energy) {
  TCanvas canvas(("c_" + collection_label + "_counts_truth_pt").c_str(),
                 "Truth pair matching counts by truth pT", 1500, 900);
  canvas.Divide(3, 2);
  for (std::size_t bin = 0; bin < histograms.size(); ++bin) {
    canvas.cd(static_cast<int>(bin + 1U));
    CollectionHistograms &current = histograms[bin];
    style_count_histograms(current);
    const double maximum = current.truth->GetMaximum();
    current.truth->SetMinimum(0.0);
    current.truth->SetMaximum(maximum > 0.0 ? 1.30 * maximum : 1.0);
    current.truth->Draw("E1");
    current.matched->Draw("E1 SAME");
    current.matched_mass_window->Draw("E1 SAME");
    TLatex label;
    label.SetNDC();
    label.SetTextAlign(13);
    label.SetTextSize(0.045);
    label.DrawLatex(0.18, 0.92, truth_pt_labels[bin].c_str());
    gPad->RedrawAxis();
  }
  canvas.cd(6);
  draw_truth_pt_information_panel(histograms.front(), collection_label,
                                  "counts", delta_r_cut, mass_min, mass_max,
                                  min_cluster_energy);
  canvas.SaveAs(output_path.c_str());
}

void draw_truth_pt_efficiencies(std::vector<CollectionHistograms> &histograms,
                                const std::string &collection_label,
                                const std::string &output_path,
                                const double delta_r_cut, const double mass_min,
                                const double mass_max,
                                const double min_cluster_energy) {
  TCanvas canvas(("c_" + collection_label + "_efficiency_truth_pt").c_str(),
                 "Truth pair matching efficiency by truth pT", 1500, 900);
  canvas.Divide(3, 2);
  for (std::size_t bin = 0; bin < histograms.size(); ++bin) {
    canvas.cd(static_cast<int>(bin + 1U));
    TH1D *matched = histograms[bin].efficiency_matched.get();
    TH1D *mass_window = histograms[bin].efficiency_mass_window.get();
    matched->SetStats(false);
    matched->SetLineColor(kBlue + 1);
    matched->SetMarkerColor(kBlue + 1);
    matched->SetMarkerStyle(21);
    mass_window->SetStats(false);
    mass_window->SetLineColor(kRed + 1);
    mass_window->SetMarkerColor(kRed + 1);
    mass_window->SetMarkerStyle(22);
    for (TH1D *histogram : {matched, mass_window}) {
      histogram->SetLineWidth(2);
      histogram->SetMarkerSize(0.7);
      histogram->GetXaxis()->SetTitle(
          "Truth energy asymmetry |E_{1}-E_{2}|/(E_{1}+E_{2})");
      histogram->GetYaxis()->SetTitle("Fraction of accepted truth pairs");
      histogram->SetMinimum(0.0);
      histogram->SetMaximum(1.05);
    }
    matched->Draw("E1");
    mass_window->Draw("E1 SAME");
    TLatex label;
    label.SetNDC();
    label.SetTextAlign(13);
    label.SetTextSize(0.045);
    label.DrawLatex(0.18, 0.92, truth_pt_labels[bin].c_str());
    gPad->RedrawAxis();
  }
  canvas.cd(6);
  draw_truth_pt_information_panel(histograms.front(), collection_label,
                                  "efficiency", delta_r_cut, mass_min, mass_max,
                                  min_cluster_energy);
  canvas.SaveAs(output_path.c_str());
}

void draw_truth_pt_delta_r(std::vector<CollectionHistograms> &histograms,
                           const std::string &collection_label,
                           const std::string &output_path,
                           const double delta_r_cut, const double mass_min,
                           const double mass_max,
                           const double min_cluster_energy) {
  TCanvas canvas(("c_" + collection_label + "_delta_r_truth_pt").c_str(),
                 "Truth-cluster delta R by truth pT", 1500, 900);
  canvas.Divide(3, 2);
  std::vector<std::unique_ptr<TLine>> cut_lines;
  cut_lines.reserve(histograms.size());
  for (std::size_t bin = 0; bin < histograms.size(); ++bin) {
    canvas.cd(static_cast<int>(bin + 1U));
    gPad->SetLogy();
    TH1D *each = histograms[bin].delta_r_each.get();
    TH1D *maximum = histograms[bin].delta_r_max.get();
    each->SetLineColor(kGray + 2);
    each->SetLineStyle(2);
    each->SetLineWidth(2);
    maximum->SetLineColor(kBlue + 1);
    maximum->SetLineWidth(2);
    maximum->GetXaxis()->SetTitle("#DeltaR(truth #gamma, matched cluster)");
    maximum->GetYaxis()->SetTitle("Entries / 0.005");
    maximum->GetXaxis()->SetRangeUser(0.0, 0.20);
    maximum->SetMinimum(0.5);
    const double plot_maximum =
        std::max(each->GetMaximum(), maximum->GetMaximum());
    maximum->SetMaximum(plot_maximum > 0.0 ? 3.0 * plot_maximum : 1.0);
    maximum->Draw("HIST");
    each->Draw("HIST SAME");
    cut_lines.push_back(std::make_unique<TLine>(delta_r_cut, 0.5, delta_r_cut,
                                                maximum->GetMaximum()));
    cut_lines.back()->SetLineColor(kRed + 1);
    cut_lines.back()->SetLineStyle(7);
    cut_lines.back()->SetLineWidth(2);
    cut_lines.back()->Draw();
    TLatex label;
    label.SetNDC();
    label.SetTextAlign(13);
    label.SetTextSize(0.045);
    label.DrawLatex(0.18, 0.92, truth_pt_labels[bin].c_str());
    gPad->RedrawAxis();
  }
  canvas.cd(6);
  draw_truth_pt_information_panel(histograms.front(), collection_label,
                                  "delta_r", delta_r_cut, mass_min, mass_max,
                                  min_cluster_energy);
  canvas.SaveAs(output_path.c_str());
}

void write_histograms(TFile &output, CollectionHistograms &histograms) {
  output.cd();
  for (TH1D *histogram :
       {histograms.truth.get(), histograms.matched.get(),
        histograms.matched_mass_window.get(),
        histograms.efficiency_matched.get(),
        histograms.efficiency_mass_window.get(), histograms.delta_r_each.get(),
        histograms.delta_r_max.get(), histograms.matched_mass.get()}) {
    histogram->Write();
  }
}

void write_histograms(TFile &output,
                      std::vector<CollectionHistograms> &histograms) {
  for (CollectionHistograms &current : histograms) {
    write_histograms(output, current);
  }
}

void print_summary(const std::string &collection,
                   const CollectionHistograms &histograms) {
  std::cout << collection
            << " truth/two-cluster/deltaR-matched/mass-window/malformed/"
               "deltaR-overflow = "
            << histograms.truth_events << "/"
            << histograms.events_with_two_clusters << "/"
            << histograms.matched_events << "/" << histograms.mass_window_events
            << "/" << histograms.malformed_events << "/"
            << histograms.delta_r_overflow_events << std::endl;
}
} // namespace

int PlotTruthPairMatching(
    const std::string input_path = "PhotonAnalysisTree/output/merged/100kevents_eta_5to15GeV_etapm1.root",
    const std::string output_base =
        "PhotonAnalysisTree/output/plots/pair_matching_efficiency/"
        "truth_pair_matching",
    const double delta_r_cut = 0.03, const double mass_window_min = 0.10,
    const double mass_window_max = 0.18, const double min_cluster_energy = 0.1,
    const int asymmetry_nbins = 20) {
  if (input_path.empty() || output_base.empty() || !(delta_r_cut > 0.0) ||
      !(mass_window_min >= 0.0 && mass_window_min < mass_window_max) ||
      !(min_cluster_energy >= 0.0) || asymmetry_nbins <= 0) {
    std::cerr << "PlotTruthPairMatching - invalid argument" << std::endl;
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
  double truth_pair_e_asym = -999.0;
  std::vector<double> *truth_eta = nullptr;
  std::vector<double> *truth_phi = nullptr;
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
  bind("truth_pair_e_asym", &truth_pair_e_asym);
  bind("truth_daughter_eta", &truth_eta);
  bind("truth_daughter_phi", &truth_phi);
  bind("truth_daughter_projection_valid", &truth_projection_valid);
  bind_collection("split", split);
  bind_collection("nosplit", nosplit);
  if (!branches_ok) {
    return 4;
  }

  CollectionHistograms split_histograms =
      make_histograms("split", asymmetry_nbins);
  CollectionHistograms nosplit_histograms =
      make_histograms("nosplit", asymmetry_nbins);
  std::vector<CollectionHistograms> split_truth_pt_histograms =
      make_truth_pt_histograms("split", asymmetry_nbins);
  std::vector<CollectionHistograms> nosplit_truth_pt_histograms =
      make_truth_pt_histograms("nosplit", asymmetry_nbins);

  const auto process_collection = [&](const CollectionBranches &branches,
                                      CollectionHistograms &histograms) {
    ++histograms.truth_events;
    histograms.truth->Fill(truth_pair_e_asym);

    const MatchResult match = match_truth_to_clusters(
        *truth_eta, *truth_phi, branches, min_cluster_energy);
    if (!match.valid) {
      const bool branch_shapes_valid =
          branches.cluster_e && branches.cluster_eta && branches.cluster_phi &&
          branches.cluster_e->size() == branches.cluster_eta->size() &&
          branches.cluster_e->size() == branches.cluster_phi->size();
      if (!branch_shapes_valid) {
        ++histograms.malformed_events;
      }
      return;
    }
    ++histograms.events_with_two_clusters;

    histograms.delta_r_each->Fill(match.delta_r0);
    histograms.delta_r_each->Fill(match.delta_r1);
    const double maximum_delta_r = std::max(match.delta_r0, match.delta_r1);
    histograms.delta_r_max->Fill(maximum_delta_r);
    if (maximum_delta_r >= histograms.delta_r_max->GetXaxis()->GetXmax()) {
      ++histograms.delta_r_overflow_events;
    }
    if (!(maximum_delta_r < delta_r_cut)) {
      return;
    }

    ++histograms.matched_events;
    histograms.matched->Fill(truth_pair_e_asym);
    double mass = -999.0;
    if (!find_pair_mass(branches, match.cluster_for_gamma0,
                        match.cluster_for_gamma1, mass)) {
      ++histograms.malformed_events;
      return;
    }
    histograms.matched_mass->Fill(mass);
    if (mass >= mass_window_min && mass <= mass_window_max) {
      ++histograms.mass_window_events;
      histograms.matched_mass_window->Fill(truth_pair_e_asym);
    }
  };

  const Long64_t entries = tree->GetEntries();
  Long64_t invalid_truth_shape = 0;
  Long64_t truth_pt_outside_bins = 0;
  for (Long64_t entry = 0; entry < entries; ++entry) {
    if (entry % 5000 == 0) {
      std::cout << "PlotTruthPairMatching - entry " << entry << " / " << entries
                << std::endl;
    }
    tree->GetEntry(entry);
    if (!truth_valid || !truth_is_pi0_to_2gamma ||
        !truth_both_gamma_in_acceptance) {
      continue;
    }
    if (!std::isfinite(truth_pair_e_asym) || truth_pair_e_asym < 0.0 ||
        truth_pair_e_asym > 1.0 || !truth_eta || !truth_phi ||
        !truth_projection_valid || truth_eta->size() != 2U ||
        truth_phi->size() != 2U ||
        truth_projection_valid->size() != 2U ||
        !truth_projection_valid->at(0) || !truth_projection_valid->at(1)) {
      ++invalid_truth_shape;
      continue;
    }

    process_collection(split, split_histograms);
    process_collection(nosplit, nosplit_histograms);
    const std::size_t truth_pt_bin = find_truth_pt_bin(truth_pt);
    if (truth_pt_bin < truth_pt_labels.size()) {
      process_collection(split, split_truth_pt_histograms[truth_pt_bin]);
      process_collection(nosplit, nosplit_truth_pt_histograms[truth_pt_bin]);
    } else {
      ++truth_pt_outside_bins;
    }
  }

  make_efficiencies(split_histograms, "split");
  make_efficiencies(nosplit_histograms, "nosplit");
  for (std::size_t bin = 0; bin < truth_pt_labels.size(); ++bin) {
    make_efficiencies(split_truth_pt_histograms[bin],
                      "split_truth_pt_" + std::to_string(bin));
    make_efficiencies(nosplit_truth_pt_histograms[bin],
                      "nosplit_truth_pt_" + std::to_string(bin));
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
  draw_count_plot(split_histograms, "SPLIT", split_output_base + "_counts.pdf",
                  delta_r_cut, mass_window_min, mass_window_max,
                  min_cluster_energy);
  draw_efficiency_plot(split_histograms, "SPLIT",
                       split_output_base + "_efficiency.pdf", delta_r_cut,
                       mass_window_min, mass_window_max, min_cluster_energy);
  draw_delta_r_plot(split_histograms, "SPLIT",
                    split_output_base + "_delta_r.pdf", delta_r_cut,
                    min_cluster_energy);

  draw_count_plot(nosplit_histograms, "NO_SPLIT",
                  nosplit_output_base + "_counts.pdf", delta_r_cut,
                  mass_window_min, mass_window_max, min_cluster_energy);
  draw_efficiency_plot(nosplit_histograms, "NO_SPLIT",
                       nosplit_output_base + "_efficiency.pdf", delta_r_cut,
                       mass_window_min, mass_window_max, min_cluster_energy);
  draw_delta_r_plot(nosplit_histograms, "NO_SPLIT",
                    nosplit_output_base + "_delta_r.pdf", delta_r_cut,
                    min_cluster_energy);
  draw_truth_pt_counts(split_truth_pt_histograms, "SPLIT",
                       split_output_base + "_counts_truth_pt.pdf", delta_r_cut,
                       mass_window_min, mass_window_max, min_cluster_energy);
  draw_truth_pt_efficiencies(split_truth_pt_histograms, "SPLIT",
                             split_output_base + "_efficiency_truth_pt.pdf",
                             delta_r_cut, mass_window_min, mass_window_max,
                             min_cluster_energy);
  draw_truth_pt_delta_r(split_truth_pt_histograms, "SPLIT",
                        split_output_base + "_delta_r_truth_pt.pdf",
                        delta_r_cut, mass_window_min, mass_window_max,
                        min_cluster_energy);
  draw_truth_pt_counts(nosplit_truth_pt_histograms, "NO_SPLIT",
                       nosplit_output_base + "_counts_truth_pt.pdf",
                       delta_r_cut, mass_window_min, mass_window_max,
                       min_cluster_energy);
  draw_truth_pt_efficiencies(nosplit_truth_pt_histograms, "NO_SPLIT",
                             nosplit_output_base + "_efficiency_truth_pt.pdf",
                             delta_r_cut, mass_window_min, mass_window_max,
                             min_cluster_energy);
  draw_truth_pt_delta_r(nosplit_truth_pt_histograms, "NO_SPLIT",
                        nosplit_output_base + "_delta_r_truth_pt.pdf",
                        delta_r_cut, mass_window_min, mass_window_max,
                        min_cluster_energy);

  TFile output((output_base + ".root").c_str(), "RECREATE");
  if (output.IsZombie()) {
    std::cerr << "Failed to create " << output_base << ".root" << std::endl;
    return 6;
  }
  write_histograms(output, split_histograms);
  write_histograms(output, nosplit_histograms);
  write_histograms(output, split_truth_pt_histograms);
  write_histograms(output, nosplit_truth_pt_histograms);
  output.Close();

  std::cout << "PlotTruthPairMatching - input events/invalid truth shape = "
            << entries << "/" << invalid_truth_shape << std::endl;
  print_summary("SPLIT", split_histograms);
  print_summary("NO_SPLIT", nosplit_histograms);
  for (std::size_t bin = 0; bin < truth_pt_labels.size(); ++bin) {
    print_summary("SPLIT " + truth_pt_labels[bin],
                  split_truth_pt_histograms[bin]);
    print_summary("NO_SPLIT " + truth_pt_labels[bin],
                  nosplit_truth_pt_histograms[bin]);
  }
  std::cout << "PlotTruthPairMatching - truth pT outside configured bins = "
            << truth_pt_outside_bins << std::endl;
  std::cout << "Wrote " << output_base << ".root and twelve PDF plots"
            << std::endl;

  const Long64_t malformed = split_histograms.malformed_events +
                             nosplit_histograms.malformed_events +
                             invalid_truth_shape;
  return malformed == 0 ? 0 : 7;
}
