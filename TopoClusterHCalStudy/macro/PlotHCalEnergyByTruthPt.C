#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TStyle.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace
{
constexpr int n_truth_pt_bins = 10;
constexpr double truth_pt_bin_width = 1.0;
constexpr double emcal_energy_threshold = 10.0;

constexpr int n_hcal_energy_bins = 30;
constexpr double hcal_energy_min = 0.0;
constexpr double hcal_energy_max = 3.0;

using HistogramArray = std::array<std::unique_ptr<TH1D>, n_truth_pt_bins>;

bool load_hcal_energy_tree(TFile &file, TTree *&tree, const std::string &input_file)
{
  if (file.IsZombie())
  {
    std::cerr << "PlotHCalEnergyByTruthPt - could not open input file: "
              << input_file << std::endl;
    return false;
  }

  file.GetObject("topocluster_tree", tree);
  if (!tree || tree->GetEntries() <= 0 || !tree->GetBranch("truth_pt") ||
      !tree->GetBranch("truth_eta") ||
      !tree->GetBranch("truth_energy_asymmetry") ||
      !tree->GetBranch("emcal_energy") ||
      !tree->GetBranch("hcal_total_energy"))
  {
    std::cerr << "PlotHCalEnergyByTruthPt - missing tree, entries, or required branches in: "
              << input_file << std::endl;
    return false;
  }
  return true;
}

bool determine_truth_pt_minimum(TTree &gamma_tree, TTree &pi0_tree, int &truth_pt_minimum) {
  const double observed_minimum = std::min(gamma_tree.GetMinimum("truth_pt"), pi0_tree.GetMinimum("truth_pt"));
  const double observed_maximum = std::max(gamma_tree.GetMaximum("truth_pt"), pi0_tree.GetMaximum("truth_pt"));

  if (observed_minimum >= 25.0 && observed_maximum <= 35.0) {
    truth_pt_minimum = 25;
  } else if (observed_minimum >= 35.0 && observed_maximum <= 45.0) {
    truth_pt_minimum = 35;
  } else {
    std::cerr << "PlotHCalEnergyByTruthPt - expected truth pT range 25-35 or 35-45 GeV, observed "
              << observed_minimum << "-" << observed_maximum << " GeV"
              << std::endl;
    return false;
  }

  std::cout << "PlotHCalEnergyByTruthPt - observed truth pT range: " << observed_minimum << "-" << observed_maximum << " GeV" << std::endl;
  return true;
}

void create_histograms(HistogramArray &histograms, const std::string &sample, const int truth_pt_minimum) {
  for (int bin = 0; bin < n_truth_pt_bins; ++bin) {
    const int pt_low = truth_pt_minimum + bin;
    const std::string name = sample + "_hcal_energy_pt_" + std::to_string(pt_low) + "_" + std::to_string(pt_low + 1);
    histograms[bin] = std::make_unique<TH1D>(
        name.c_str(),
        "",
        n_hcal_energy_bins,
        hcal_energy_min,
        hcal_energy_max);
    histograms[bin]->SetDirectory(nullptr);
    histograms[bin]->Sumw2();
  }
}

bool fill_histograms(
    TTree &tree,
    HistogramArray &histograms,
    const int truth_pt_minimum,
    const double truth_abs_eta_minimum,
    const double truth_abs_eta_maximum,
    const bool apply_truth_energy_asymmetry_cut,
    const double truth_energy_asymmetry_minimum,
    const double truth_energy_asymmetry_maximum)
{
  float truth_pt = 0.0F;
  float truth_eta = 0.0F;
  float truth_energy_asymmetry = -1.0F;
  std::vector<float> *emcal_energy = nullptr;
  std::vector<float> *hcal_total_energy = nullptr;
  tree.SetBranchAddress("truth_pt", &truth_pt);
  tree.SetBranchAddress("truth_eta", &truth_eta);
  tree.SetBranchAddress("truth_energy_asymmetry", &truth_energy_asymmetry);
  tree.SetBranchAddress("emcal_energy", &emcal_energy);
  tree.SetBranchAddress("hcal_total_energy", &hcal_total_energy);

  const Long64_t n_entries = tree.GetEntries();
  for (Long64_t entry = 0; entry < n_entries; ++entry) {
    tree.GetEntry(entry);
    const double truth_abs_eta = std::abs(truth_eta);
    if (truth_abs_eta < truth_abs_eta_minimum || truth_abs_eta >= truth_abs_eta_maximum) {
      continue;
    }
    if (apply_truth_energy_asymmetry_cut &&
        (truth_energy_asymmetry < truth_energy_asymmetry_minimum || truth_energy_asymmetry >= truth_energy_asymmetry_maximum)) {
      continue;
    }
    const int bin = static_cast<int>(std::floor((truth_pt - truth_pt_minimum) / truth_pt_bin_width));
    if (bin < 0 || bin >= n_truth_pt_bins) {
      continue;
    }
    if (!emcal_energy || !hcal_total_energy || emcal_energy->size() != hcal_total_energy->size()) {
      std::cerr << "PlotHCalEnergyByTruthPt - inconsistent cluster vectors at entry " << entry << std::endl;
      tree.ResetBranchAddresses();
      return false;
    }

    for (std::size_t cluster = 0; cluster < emcal_energy->size(); ++cluster) {
      if (emcal_energy->at(cluster) <= emcal_energy_threshold) {
        continue;
      }
      histograms[bin]->Fill(hcal_total_energy->at(cluster));
    }
  }
  tree.ResetBranchAddresses();

  for (int bin = 0; bin < n_truth_pt_bins; ++bin) {
    const double count = histograms[bin]->GetEntries();
    const double integral = histograms[bin]->Integral();
    if (count <= 0.0 || integral <= 0.0) {
      const double pt_low = truth_pt_minimum + bin * truth_pt_bin_width;
      const double pt_high = pt_low + truth_pt_bin_width;
      std::cerr << "PlotHCalEnergyByTruthPt - no selected TopoClusters for truth pT bin " << pt_low << "-" << pt_high << " GeV" << std::endl;
      return false;
    }
    // Normalize to every selected cluster, including the overflow above the
    // displayed 3 GeV range, so the visible shapes retain their tail fraction.
    histograms[bin]->Scale(1.0 / count);
  }
  return true;
}

void draw_panel(
    TH1D &gamma_histogram,
    TH1D &pi0_histogram,
    const int pt_low,
    const double truth_abs_eta_minimum,
    const double truth_abs_eta_maximum,
    const double truth_energy_asymmetry_minimum,
    const double truth_energy_asymmetry_maximum)
{
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.04);
  gPad->SetBottomMargin(0.14);
  gPad->SetTopMargin(0.10);

  gamma_histogram.SetTitle(Form(
      "%d #leq p_{T}^{truth} < %d GeV;E_{HCalIn}^{cluster} + E_{HCalOut}^{cluster} [GeV];Normalized TopoClusters", pt_low, pt_low + 1));
  gamma_histogram.SetLineColor(kAzure + 2);
  gamma_histogram.SetLineWidth(3);
  pi0_histogram.SetLineColor(kOrange + 7);
  pi0_histogram.SetLineWidth(3);

  const double common_maximum = std::max(gamma_histogram.GetMaximum(), pi0_histogram.GetMaximum());
  gamma_histogram.SetMinimum(1.0e-4);
  gamma_histogram.SetMaximum(2.5 * common_maximum);
  gamma_histogram.GetXaxis()->CenterTitle();
  gamma_histogram.GetYaxis()->CenterTitle();
  gamma_histogram.GetXaxis()->SetTitleSize(0.045);
  gamma_histogram.GetYaxis()->SetTitleSize(0.045);
  gamma_histogram.GetXaxis()->SetLabelSize(0.040);
  gamma_histogram.GetYaxis()->SetLabelSize(0.040);
  gamma_histogram.GetXaxis()->SetTitleOffset(1.25);
  gamma_histogram.GetYaxis()->SetTitleOffset(1.45);
  gamma_histogram.Draw("HIST");
  pi0_histogram.Draw("HIST SAME");

  TLegend legend(0.38, 0.57, 0.94, 0.89);
  legend.SetHeader(Form(
      "#splitline{%.2g #leq |#eta^{truth}| < %.2g, E_{EMCal}^{cluster} > %.0f GeV}{#pi^{0}: %.2g #leq A_{E}^{truth} < %.2g}",
      truth_abs_eta_minimum,
      truth_abs_eta_maximum,
      emcal_energy_threshold,
      truth_energy_asymmetry_minimum,
      truth_energy_asymmetry_maximum));
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetTextSize(0.038);
  legend.AddEntry(&gamma_histogram, Form("Single #gamma (N=%.0f)", gamma_histogram.GetEntries()), "l");
  legend.AddEntry(&pi0_histogram, Form("Single #pi^{0} (N=%.0f)", pi0_histogram.GetEntries()), "l");
  legend.DrawClone();
}
}  // namespace

int PlotHCalEnergyByTruthPt(
    const std::string gamma_input_file = "/sphenix/user/ryotaro/DirectPhotonAnalysis/TopoClusterHCalStudy/output/merge/35to45GeV/topocluster_hcal_gamma_merged.root",
    const std::string pi0_input_file = "/sphenix/user/ryotaro/DirectPhotonAnalysis/TopoClusterHCalStudy/output/merge/35to45GeV/topocluster_hcal_pi0_merged.root",
    const std::string output_directory = "/sphenix/user/ryotaro/DirectPhotonAnalysis/TopoClusterHCalStudy/output/plot/35to45GeV",
    const double truth_abs_eta_minimum = 0.0,
    const double truth_abs_eta_maximum = 0.1,
    const double truth_energy_asymmetry_minimum = 0.0,
    const double truth_energy_asymmetry_maximum = 1.)
{
  if (!std::isfinite(truth_abs_eta_minimum) ||
      !std::isfinite(truth_abs_eta_maximum) ||
      truth_abs_eta_minimum < 0.0 ||
      truth_abs_eta_maximum <= truth_abs_eta_minimum ||
      !std::isfinite(truth_energy_asymmetry_minimum) ||
      !std::isfinite(truth_energy_asymmetry_maximum) ||
      truth_energy_asymmetry_minimum < 0.0 ||
      truth_energy_asymmetry_maximum > 1.0 ||
      truth_energy_asymmetry_maximum <= truth_energy_asymmetry_minimum)
  {
    std::cerr << "PlotHCalEnergyByTruthPt - invalid truth eta or energy asymmetry range" << std::endl;
    return EXIT_FAILURE;
  }

  const std::filesystem::path output_path(output_directory);
  std::error_code directory_error;
  std::filesystem::create_directories(output_path, directory_error);
  if (directory_error) {
    std::cerr << "PlotHCalEnergyByTruthPt - could not create output directory: " << directory_error.message() << std::endl;
    return EXIT_FAILURE;
  }

  TFile gamma_file(gamma_input_file.c_str(), "READ");
  TFile pi0_file(pi0_input_file.c_str(), "READ");
  TTree *gamma_tree = nullptr;
  TTree *pi0_tree = nullptr;
  if (!load_hcal_energy_tree(gamma_file, gamma_tree, gamma_input_file) || !load_hcal_energy_tree(pi0_file, pi0_tree, pi0_input_file)) {
    return EXIT_FAILURE;
  }

  int truth_pt_minimum = 0;
  if (!determine_truth_pt_minimum(*gamma_tree, *pi0_tree, truth_pt_minimum)) {
    return EXIT_FAILURE;
  }

  HistogramArray gamma_histograms;
  HistogramArray pi0_histograms;
  create_histograms(gamma_histograms, "gamma", truth_pt_minimum);
  create_histograms(pi0_histograms, "pi0", truth_pt_minimum);
  if (!fill_histograms(
          *gamma_tree,
          gamma_histograms,
          truth_pt_minimum,
          truth_abs_eta_minimum,
          truth_abs_eta_maximum,
          false,
          truth_energy_asymmetry_minimum,
          truth_energy_asymmetry_maximum) ||
      !fill_histograms(
          *pi0_tree,
          pi0_histograms,
          truth_pt_minimum,
          truth_abs_eta_minimum,
          truth_abs_eta_maximum,
          true,
          truth_energy_asymmetry_minimum,
          truth_energy_asymmetry_maximum))
  {
    return EXIT_FAILURE;
  }

  gStyle->SetOptStat(0);
  TCanvas canvas("hcal_energy_by_truth_pt", "HCAL energy by truth pT", 2000, 1400);
  canvas.Divide(5, 2, 0.001, 0.001);

  for (int bin = 0; bin < n_truth_pt_bins; ++bin) {
    canvas.cd(bin + 1);
    draw_panel(
        *gamma_histograms[bin],
        *pi0_histograms[bin],
        truth_pt_minimum + bin,
        truth_abs_eta_minimum,
        truth_abs_eta_maximum,
        truth_energy_asymmetry_minimum,
        truth_energy_asymmetry_maximum);
  }

  const int truth_pt_maximum = truth_pt_minimum + n_truth_pt_bins;
  std::string eta_tag = Form("%.3gto%.3g", truth_abs_eta_minimum, truth_abs_eta_maximum);
  std::replace(eta_tag.begin(), eta_tag.end(), '.', 'p');
  std::string asymmetry_tag = Form("%.3gto%.3g", truth_energy_asymmetry_minimum, truth_energy_asymmetry_maximum);
  std::replace(asymmetry_tag.begin(), asymmetry_tag.end(), '.', 'p');
  const std::string output_name = "hcal_total_energy_truth_pt_" + std::to_string(truth_pt_minimum) + "to" + std::to_string(truth_pt_maximum) +
                                  "_abseta_" + eta_tag + "_pi0ae_" + asymmetry_tag + ".pdf";
  const std::string output_file = (output_path / output_name).string();
  canvas.SaveAs(output_file.c_str());

  for (int bin = 0; bin < n_truth_pt_bins; ++bin) {
    std::cout << "PlotHCalEnergyByTruthPt - "
              << truth_pt_minimum + bin << "-"
              << truth_pt_minimum + bin + 1 << " GeV: gamma="
              << gamma_histograms[bin]->GetEntries()
              << ", pi0=" << pi0_histograms[bin]->GetEntries()
              << std::endl;
  }
  std::cout << "PlotHCalEnergyByTruthPt - output: "
            << output_file << std::endl;
  return EXIT_SUCCESS;
}
