#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TColor.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TROOT.h>

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>

namespace
{
constexpr int n_emcal_bins = 120;
constexpr double emcal_min = 0.0;
constexpr double emcal_max = 60.0;
constexpr int n_hcal_bins = 160;
constexpr double hcal_min = 0.0;
constexpr double hcal_max = 80.0;

constexpr int n_ratio_bins = 100;
constexpr double ratio_min = 0.0;
constexpr double ratio_max = 2.0;
constexpr double ratio_emcal_threshold = 10.0;

bool load_correlation_tree(TFile &file, TTree *&tree, const std::string &input_file) {
  if (file.IsZombie()) {
    std::cerr << "PlotEMCalHCalCorrelation - could not open input file: " << input_file << std::endl;
    return false;
  }

  file.GetObject("topocluster_tree", tree);
  if (!tree || tree->GetEntries() <= 0 || !tree->GetBranch("emcal_energy") || !tree->GetBranch("hcal_total_energy")) {
    std::cerr << "PlotEMCalHCalCorrelation - missing tree, entries, or required branches in: " << input_file << std::endl;
    return false;
  }
  return true;
}

bool fill_histograms(TTree &tree, TH2D &correlation_histogram, TH1D &ratio_histogram) {
  const std::string correlation_expression = "hcal_total_energy:emcal_energy>>" + std::string(correlation_histogram.GetName());
  if (tree.Draw(correlation_expression.c_str(), "", "goff") < 0) {
    std::cerr << "PlotEMCalHCalCorrelation - correlation draw failed" << std::endl;
    return false;
  }

  const std::string ratio_expression = "hcal_total_energy/emcal_energy>>" + std::string(ratio_histogram.GetName());
  const std::string ratio_selection =  "emcal_energy>" + std::to_string(ratio_emcal_threshold);
  if (tree.Draw(ratio_expression.c_str(), ratio_selection.c_str(), "goff") <= 0 || ratio_histogram.Integral() <= 0.0) {
    std::cerr << "PlotEMCalHCalCorrelation - ratio draw failed" << std::endl;
    return false;
  }

  ratio_histogram.Scale(1.0 / ratio_histogram.Integral());
  return true;
}

void draw_correlation_panel(
    TH2D &histogram,
    const char *sample_label,
    const Long64_t n_events,
    const double common_maximum)
{
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.16);
  gPad->SetBottomMargin(0.12);
  gPad->SetTopMargin(0.10);
  gPad->SetLogz();

  histogram.SetTitle(Form("%s;E_{EMCal}^{cluster} [GeV];E_{HCalIn}^{cluster} + E_{HCalOut}^{cluster} [GeV]", sample_label));
  histogram.SetMinimum(1.0);
  histogram.SetMaximum(common_maximum);
  histogram.GetXaxis()->CenterTitle();
  histogram.GetYaxis()->CenterTitle();
  histogram.GetZaxis()->SetTitle("TopoClusters / bin");
  histogram.GetZaxis()->CenterTitle();
  histogram.GetXaxis()->SetTitleOffset(1.15);
  histogram.GetYaxis()->SetTitleOffset(1.25);
  histogram.GetZaxis()->SetTitleOffset(1.20);
  histogram.Draw("COLZ");

  TLatex label;
  label.SetNDC();
  label.SetTextSize(0.032);
  label.SetTextAlign(13);
  label.DrawLatex(0.14, 0.87, Form("25 < E_{particle} < 35 GeV, events = %lld, clusters = %lld", n_events, static_cast<Long64_t>(histogram.GetEntries())));
}

void draw_ratio_comparison(TH1D &gamma_histogram, TH1D &pi0_histogram) {
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.04);
  gPad->SetBottomMargin(0.12);
  gPad->SetTopMargin(0.08);
  gPad->SetLogy();

  gamma_histogram.SetTitle("; (E_{HCalIn}^{cluster} + E_{HCalOut}^{cluster}) / E_{EMCal}^{cluster}; Normalized TopoClusters");
  gamma_histogram.SetLineColor(kAzure + 2);
  gamma_histogram.SetLineWidth(3);
  pi0_histogram.SetLineColor(kOrange + 7);
  pi0_histogram.SetLineWidth(3);

  const double common_maximum = std::max(gamma_histogram.GetMaximum(), pi0_histogram.GetMaximum());
  gamma_histogram.SetMinimum(1.0e-5);
  gamma_histogram.SetMaximum(2.0 * common_maximum);
  gamma_histogram.GetXaxis()->CenterTitle();
  gamma_histogram.GetYaxis()->CenterTitle();
  gamma_histogram.GetXaxis()->SetTitleOffset(1.15);
  gamma_histogram.GetYaxis()->SetTitleOffset(1.35);
  gamma_histogram.Draw("HIST");
  pi0_histogram.Draw("HIST SAME");

  TLegend legend(0.56, 0.70, 0.92, 0.88);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetTextSize(0.035);
  legend.AddEntry(
      &gamma_histogram,
      Form("Single #gamma (%.0f clusters)", gamma_histogram.GetEntries()),
      "l");
  legend.AddEntry(
      &pi0_histogram,
      Form("Single #pi^{0} (%.0f clusters)", pi0_histogram.GetEntries()),
      "l");
  legend.DrawClone();

  TLatex label;
  label.SetNDC();
  label.SetTextSize(0.036);
  label.SetTextAlign(13);
  label.DrawLatex(0.15, 0.88, "25 < E_{particle} < 35 GeV");
  label.DrawLatex(0.15, 0.83, Form("E_{EMCal}^{cluster} > %.0f GeV", ratio_emcal_threshold));
}
}  // namespace

int PlotEMCalHCalCorrelation(
    const std::string gamma_input_file = "/sphenix/user/ryotaro/DirectPhotonAnalysis/TopoClusterHCalStudy/output/merge/25to35GeV/topocluster_hcal_gamma_merged.root",
    const std::string pi0_input_file =   "/sphenix/user/ryotaro/DirectPhotonAnalysis/TopoClusterHCalStudy/output/merge/25to35GeV/topocluster_hcal_pi0_merged.root",
    const std::string output_directory = "/sphenix/user/ryotaro/DirectPhotonAnalysis/TopoClusterHCalStudy/output/plot/25to35GeV")
{
  const std::filesystem::path output_path(output_directory);
  std::error_code directory_error;
  std::filesystem::create_directories(output_path, directory_error);
  if (directory_error) {
    std::cerr << "PlotEMCalHCalCorrelation - could not create output directory: " << directory_error.message() << std::endl;
    return EXIT_FAILURE;
  }

  TFile gamma_file(gamma_input_file.c_str(), "READ");
  TFile pi0_file(pi0_input_file.c_str(), "READ");
  TTree *gamma_tree = nullptr;
  TTree *pi0_tree = nullptr;
  if (!load_correlation_tree(gamma_file, gamma_tree, gamma_input_file) || !load_correlation_tree(pi0_file, pi0_tree, pi0_input_file)) {
    return EXIT_FAILURE;
  }

  gROOT->cd();
  TH2D gamma_correlation("gamma_emcal_hcal_correlation", "", n_emcal_bins, emcal_min, emcal_max, n_hcal_bins, hcal_min, hcal_max);
  TH2D pi0_correlation("pi0_emcal_hcal_correlation", "", n_emcal_bins, emcal_min, emcal_max, n_hcal_bins, hcal_min, hcal_max);
  TH1D gamma_ratio("gamma_hcal_emcal_ratio", "", n_ratio_bins, ratio_min, ratio_max);
  TH1D pi0_ratio("pi0_hcal_emcal_ratio", "", n_ratio_bins, ratio_min, ratio_max);

  if (!fill_histograms(*gamma_tree, gamma_correlation, gamma_ratio) || !fill_histograms(*pi0_tree, pi0_correlation, pi0_ratio)) {
    return EXIT_FAILURE;
  }

  const Long64_t gamma_events = gamma_tree->GetEntries();
  const Long64_t pi0_events = pi0_tree->GetEntries();

  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBird);
  gStyle->SetNumberContours(255);

  const double common_correlation_maximum = std::max(gamma_correlation.GetMaximum(), pi0_correlation.GetMaximum());

  TCanvas correlation_canvas("emcal_hcal_correlation", "EMCal-HCal correlation", 1500, 680);
  correlation_canvas.Divide(2, 1, 0.001, 0.001);
  correlation_canvas.cd(1);
  draw_correlation_panel(gamma_correlation, "Single #gamma", gamma_events, common_correlation_maximum);
  correlation_canvas.cd(2);
  draw_correlation_panel(pi0_correlation, "Single #pi^{0}", pi0_events, common_correlation_maximum);

  TCanvas ratio_canvas("hcal_emcal_ratio", "HCAL-EMCal energy ratio", 900, 700);
  draw_ratio_comparison(gamma_ratio, pi0_ratio);

  const std::string correlation_output = (output_path / "emcal_hcal_correlation.pdf").string();
  const std::string ratio_output = (output_path / "hcal_emcal_ratio.pdf").string();
  correlation_canvas.SaveAs(correlation_output.c_str());
  ratio_canvas.SaveAs(ratio_output.c_str());

  std::cout << "PlotEMCalHCalCorrelation - gamma: " << gamma_events
            << " events, " << gamma_correlation.GetEntries() << " clusters, "
            << gamma_ratio.GetEntries() << " ratio entries\n"
            << "PlotEMCalHCalCorrelation - pi0: " << pi0_events
            << " events, " << pi0_correlation.GetEntries() << " clusters, "
            << pi0_ratio.GetEntries() << " ratio entries\n"
            << "PlotEMCalHCalCorrelation - output: " << correlation_output
            << '\n'
            << "PlotEMCalHCalCorrelation - output: " << ratio_output
            << std::endl;
  return EXIT_SUCCESS;
}
