#include "Utilities/sPhenixStyle.C"

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
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
const std::vector<double> truth_pt_edges = {5.0, 7.0, 9.0, 11.0, 13.0, 15.0};
const std::array<std::string, 5> truth_pt_labels = {
    "5 #leq p_{T}^{truth} < 7 GeV", "7 #leq p_{T}^{truth} < 9 GeV",
    "9 #leq p_{T}^{truth} < 11 GeV", "11 #leq p_{T}^{truth} < 13 GeV",
    "13 #leq p_{T}^{truth} #leq 15 GeV"};

struct CollectionBranches {
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

bool valid_cluster(const CollectionBranches &branches,
                   const std::size_t index) {
  return branches.cluster_et && branches.cluster_eta && branches.cluster_phi &&
         index < branches.cluster_et->size() &&
         index < branches.cluster_eta->size() &&
         index < branches.cluster_phi->size() &&
         std::isfinite(branches.cluster_et->at(index)) &&
         std::isfinite(branches.cluster_eta->at(index)) &&
         std::isfinite(branches.cluster_phi->at(index));
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
                        const CollectionBranches &branches) {
  MatchResult result;
  if (truth_eta.size() != 2U || truth_phi.size() != 2U ||
      !branches.cluster_et) {
    return result;
  }
  double best_cost = std::numeric_limits<double>::infinity();
  for (std::size_t first = 0; first < branches.cluster_et->size(); ++first) {
    if (!valid_cluster(branches, first)) {
      continue;
    }
    const double first_delta_r =
        delta_r(truth_eta[0], truth_phi[0], branches.cluster_eta->at(first),
                branches.cluster_phi->at(first));
    for (std::size_t second = 0; second < branches.cluster_et->size();
         ++second) {
      if (first == second || !valid_cluster(branches, second)) {
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
                         const CollectionBranches &branches) {
  MatchResult result;
  if (truth_eta.size() != 2U || truth_phi.size() != 2U ||
      !valid_cluster(branches, anchor)) {
    return result;
  }
  double best_cost = std::numeric_limits<double>::infinity();
  for (std::size_t partner = 0; partner < branches.cluster_et->size();
       ++partner) {
    if (partner == anchor || !valid_cluster(branches, partner)) {
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

void fill_event_histograms(
    StageHistograms &histograms, const double truth_alpha,
    const std::vector<double> &truth_eta, const std::vector<double> &truth_phi,
    const CollectionBranches &branches, const double anchor_eta_max,
    const double anchor_et_min, const double anchor_et_max,
    const double delta_r_cut, const double mass_min, const double mass_max) {
  ++histograms.n_reference;
  histograms.reference->Fill(truth_alpha);
  if (!branches.cluster_et || !branches.cluster_eta || !branches.cluster_phi ||
      branches.cluster_et->size() != branches.cluster_eta->size() ||
      branches.cluster_et->size() != branches.cluster_phi->size()) {
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

  const MatchResult match = match_event(truth_eta, truth_phi, branches);
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
    const double delta_r_cut, const double mass_min, const double mass_max) {
  if (!branches.cluster_et || !branches.cluster_eta || !branches.cluster_phi ||
      branches.cluster_et->size() != branches.cluster_eta->size() ||
      branches.cluster_et->size() != branches.cluster_phi->size()) {
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
    const MatchResult match =
        match_anchor(anchor, truth_eta, truth_phi, branches);
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
  histograms.reference->SetLineColor(kBlack);
  histograms.reference->SetMarkerColor(kBlack);
  histograms.denominator->SetLineColor(kGreen + 2);
  histograms.denominator->SetMarkerColor(kGreen + 2);
  histograms.matched->SetLineColor(kBlue + 1);
  histograms.matched->SetMarkerColor(kBlue + 1);
  histograms.removed->SetLineColor(kRed + 1);
  histograms.removed->SetMarkerColor(kRed + 1);
  int marker = 20;
  for (TH1D *histogram :
       {histograms.reference.get(), histograms.denominator.get(),
        histograms.matched.get(), histograms.removed.get()}) {
    histogram->SetMarkerStyle(marker++);
    histogram->SetMarkerSize(0.65);
    histogram->SetLineWidth(2);
    histogram->GetXaxis()->SetTitle("Truth #alpha");
    histogram->GetYaxis()->SetTitle(event_family ? "Events" : "Clusters");
  }
}

void draw_information(StageHistograms &histograms,
                      const std::string &collection_label,
                      const bool event_family, const bool efficiency,
                      const double anchor_eta_max, const double anchor_et_min,
                      const double anchor_et_max, const double delta_r_cut,
                      const double mass_min, const double mass_max) {
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
  TLegend legend(0.14, 0.12, 0.93, 0.43);
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
    legend.AddEntry(histograms.reference.get(), "Accepted truth events", "lep");
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
      "Conditional counts by truth pT", 1500, 900);
  canvas.Divide(3, 2);
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
  canvas.cd(6);
  draw_information(histograms.front(), collection_label, event_family, false,
                   anchor_eta_max, anchor_et_min, anchor_et_max, delta_r_cut,
                   mass_min, mass_max);
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
                 "Conditional efficiencies by truth pT", 1500, 900);
  canvas.Divide(3, 2);
  for (std::size_t bin = 0; bin < histograms.size(); ++bin) {
    canvas.cd(static_cast<int>(bin + 1U));
    TH1D *matched = histograms[bin].efficiency_matched.get();
    TH1D *removed = histograms[bin].efficiency_removed.get();
    matched->SetLineColor(kBlue + 1);
    matched->SetMarkerColor(kBlue + 1);
    matched->SetMarkerStyle(22);
    removed->SetLineColor(kRed + 1);
    removed->SetMarkerColor(kRed + 1);
    removed->SetMarkerStyle(23);
    for (TH1D *histogram : {matched, removed}) {
      histogram->SetStats(false);
      histogram->SetLineWidth(2);
      histogram->SetMarkerSize(0.7);
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
  canvas.cd(6);
  draw_information(histograms.front(), collection_label, event_family, true,
                   anchor_eta_max, anchor_et_min, anchor_et_max, delta_r_cut,
                   mass_min, mass_max);
  canvas.SaveAs(output_path.c_str());
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
} // namespace

int PlotConditionalPartnerEfficiency(
    const std::string input_path = "PhotonAnalysisTree/output/merged/100kevents_pi0_5to15GeV_etapm1.root",
    const std::string output_base =
        "PhotonAnalysisTree/output/plots/conditional_partner_et5to7",
    const double anchor_eta_max = 0.7, const double anchor_et_min = 5.0,
    const double anchor_et_max = 7.0, const double delta_r_cut = 0.03,
    const double mass_window_min = 0.10, const double mass_window_max = 0.18,
    const int asymmetry_nbins = 20) {
  if (input_path.empty() || output_base.empty() || !(anchor_eta_max > 0.0) ||
      !(anchor_et_min >= 0.0 && anchor_et_min < anchor_et_max) ||
      !(delta_r_cut > 0.0) ||
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
  double truth_alpha = -999.0;
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
  bind("truth_pair_e_asym", &truth_alpha);
  bind("truth_daughter_projection_eta", &truth_eta);
  bind("truth_daughter_projection_phi", &truth_phi);
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

  Long64_t invalid_truth_shape = 0;
  Long64_t truth_pt_outside_bins = 0;
  for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry) {
    if (entry % 5000 == 0) {
      std::cout << "PlotConditionalPartnerEfficiency - entry " << entry << " / "
                << tree->GetEntries() << std::endl;
    }
    tree->GetEntry(entry);
    if (!truth_valid || !truth_is_pi0_to_2gamma ||
        !truth_both_gamma_in_acceptance) {
      continue;
    }
    if (!std::isfinite(truth_alpha) || truth_alpha < 0.0 || truth_alpha > 1.0 ||
        !truth_eta || !truth_phi || !truth_projection_valid ||
        truth_eta->size() != 2U || truth_phi->size() != 2U ||
        truth_projection_valid->size() != 2U ||
        !truth_projection_valid->at(0) || !truth_projection_valid->at(1)) {
      ++invalid_truth_shape;
      continue;
    }
    const std::size_t bin = find_truth_pt_bin(truth_pt);
    if (bin >= truth_pt_labels.size()) {
      ++truth_pt_outside_bins;
      continue;
    }
    fill_event_histograms(split_event[bin], truth_alpha, *truth_eta, *truth_phi,
                          split, anchor_eta_max, anchor_et_min, anchor_et_max,
                          delta_r_cut, mass_window_min, mass_window_max);
    fill_cluster_histograms(split_cluster[bin], truth_alpha, *truth_eta,
                            *truth_phi, split, anchor_eta_max, anchor_et_min,
                            anchor_et_max, delta_r_cut, mass_window_min,
                            mass_window_max);
    fill_event_histograms(nosplit_event[bin], truth_alpha, *truth_eta,
                          *truth_phi, nosplit, anchor_eta_max, anchor_et_min,
                          anchor_et_max, delta_r_cut, mass_window_min,
                          mass_window_max);
    fill_cluster_histograms(nosplit_cluster[bin], truth_alpha, *truth_eta,
                            *truth_phi, nosplit, anchor_eta_max, anchor_et_min,
                            anchor_et_max, delta_r_cut, mass_window_min,
                            mass_window_max);
  }

  for (std::size_t bin = 0; bin < truth_pt_labels.size(); ++bin) {
    make_efficiencies(split_event[bin],
                      "split_event_truth_pt_" + std::to_string(bin));
    make_efficiencies(split_cluster[bin],
                      "split_cluster_truth_pt_" + std::to_string(bin));
    make_efficiencies(nosplit_event[bin],
                      "nosplit_event_truth_pt_" + std::to_string(bin));
    make_efficiencies(nosplit_cluster[bin],
                      "nosplit_cluster_truth_pt_" + std::to_string(bin));
  }
  if (!make_output_directory(output_base)) {
    return 5;
  }
  const auto draw_collection = [&](std::vector<StageHistograms> &event,
                                   std::vector<StageHistograms> &cluster,
                                   const std::string &prefix,
                                   const std::string &label) {
    draw_counts(event, label, true,
                output_base + "_" + prefix + "_event_counts_truth_pt.pdf",
                anchor_eta_max, anchor_et_min, anchor_et_max, delta_r_cut,
                mass_window_min, mass_window_max);
    draw_efficiencies(event, label, true,
                      output_base + "_" + prefix +
                          "_event_efficiency_truth_pt.pdf",
                      anchor_eta_max, anchor_et_min, anchor_et_max, delta_r_cut,
                      mass_window_min, mass_window_max);
    draw_counts(cluster, label, false,
                output_base + "_" + prefix + "_cluster_counts_truth_pt.pdf",
                anchor_eta_max, anchor_et_min, anchor_et_max, delta_r_cut,
                mass_window_min, mass_window_max);
    draw_efficiencies(cluster, label, false,
                      output_base + "_" + prefix +
                          "_cluster_efficiency_truth_pt.pdf",
                      anchor_eta_max, anchor_et_min, anchor_et_max, delta_r_cut,
                      mass_window_min, mass_window_max);
  };
  draw_collection(split_event, split_cluster, "split", "SPLIT");
  draw_collection(nosplit_event, nosplit_cluster, "nosplit", "NO_SPLIT");

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
  TParameter<double> stored_eta_max("anchor_eta_max", anchor_eta_max);
  TParameter<double> stored_et_min("anchor_et_min", anchor_et_min);
  TParameter<double> stored_et_max("anchor_et_max", anchor_et_max);
  TParameter<double> stored_delta_r("delta_r_cut", delta_r_cut);
  TParameter<double> stored_mass_min("mass_window_min", mass_window_min);
  TParameter<double> stored_mass_max("mass_window_max", mass_window_max);
  source.Write();
  stored_eta_max.Write();
  stored_et_min.Write();
  stored_et_max.Write();
  stored_delta_r.Write();
  stored_mass_min.Write();
  stored_mass_max.Write();
  output.Close();

  print_summary("SPLIT event", split_event);
  print_summary("SPLIT cluster", split_cluster);
  print_summary("NO_SPLIT event", nosplit_event);
  print_summary("NO_SPLIT cluster", nosplit_cluster);
  std::cout << "PlotConditionalPartnerEfficiency - invalid truth shape / truth "
               "pT outside bins = "
            << invalid_truth_shape << " / " << truth_pt_outside_bins
            << std::endl;
  std::cout << "Wrote " << output_base << ".root and eight PDF plots"
            << std::endl;
  Long64_t malformed = invalid_truth_shape;
  for (const auto *family :
       {&split_event, &split_cluster, &nosplit_event, &nosplit_cluster}) {
    for (const StageHistograms &histograms : *family) {
      malformed += histograms.n_malformed;
    }
  }
  return malformed == 0 ? 0 : 7;
}
