#include <TCanvas.h>
#include <TChain.h>
#include <TH1D.h>
#include <TLatex.h>
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

namespace
{
constexpr int n_truth_pt_bins = 10;
constexpr double truth_pt_bin_width = 1.0;
constexpr double emcal_energy_threshold = 10.0;

constexpr int n_hcal_energy_bins = 120;
constexpr double hcal_energy_min = 0.0;
constexpr double hcal_energy_max = 3.0;

using HistogramArray = std::array<std::unique_ptr<TH1D>, n_truth_pt_bins>;
using CountArray = std::array<Long64_t, n_truth_pt_bins>;

bool load_chain(TChain &chain, const std::string &input_pattern)
{
  if (chain.Add(input_pattern.c_str()) == 0)
  {
    std::cerr << "PlotHCalEnergyByTruthPt - no input files matched: "
              << input_pattern << std::endl;
    return false;
  }

  if (chain.GetEntries() <= 0 || !chain.GetBranch("truth_pt") ||
      !chain.GetBranch("emcal_energy") ||
      !chain.GetBranch("hcal_total_energy"))
  {
    std::cerr << "PlotHCalEnergyByTruthPt - missing entries or required branches in: "
              << input_pattern << std::endl;
    return false;
  }
  return true;
}

bool determine_truth_pt_minimum(
    TChain &gamma_chain,
    TChain &pi0_chain,
    int &truth_pt_minimum)
{
  const double observed_minimum = std::min(
      gamma_chain.GetMinimum("truth_pt"),
      pi0_chain.GetMinimum("truth_pt"));
  const double observed_maximum = std::max(
      gamma_chain.GetMaximum("truth_pt"),
      pi0_chain.GetMaximum("truth_pt"));

  if (observed_minimum >= 25.0 && observed_maximum <= 35.0)
  {
    truth_pt_minimum = 25;
  }
  else if (observed_minimum >= 35.0 && observed_maximum <= 45.0)
  {
    truth_pt_minimum = 35;
  }
  else
  {
    std::cerr << "PlotHCalEnergyByTruthPt - expected truth pT range 25-35 or 35-45 GeV, observed "
              << observed_minimum << "-" << observed_maximum << " GeV"
              << std::endl;
    return false;
  }

  std::cout << "PlotHCalEnergyByTruthPt - observed truth pT range: "
            << observed_minimum << "-" << observed_maximum << " GeV"
            << std::endl;
  return true;
}

void create_histograms(
    HistogramArray &histograms,
    const std::string &sample,
    const int truth_pt_minimum)
{
  for (int bin = 0; bin < n_truth_pt_bins; ++bin)
  {
    const int pt_low = truth_pt_minimum + bin;
    const std::string name = sample + "_hcal_energy_pt_" +
        std::to_string(pt_low) + "_" + std::to_string(pt_low + 1);
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
    TChain &chain,
    HistogramArray &histograms,
    CountArray &counts,
    const int truth_pt_minimum)
{
  float truth_pt = 0.0F;
  std::vector<float> *emcal_energy = nullptr;
  std::vector<float> *hcal_total_energy = nullptr;
  chain.SetBranchAddress("truth_pt", &truth_pt);
  chain.SetBranchAddress("emcal_energy", &emcal_energy);
  chain.SetBranchAddress("hcal_total_energy", &hcal_total_energy);

  const Long64_t n_entries = chain.GetEntries();
  for (Long64_t entry = 0; entry < n_entries; ++entry)
  {
    chain.GetEntry(entry);
    const int bin = static_cast<int>(
        std::floor((truth_pt - truth_pt_minimum) / truth_pt_bin_width));
    if (bin < 0 || bin >= n_truth_pt_bins)
    {
      continue;
    }
    if (!emcal_energy || !hcal_total_energy ||
        emcal_energy->size() != hcal_total_energy->size())
    {
      std::cerr << "PlotHCalEnergyByTruthPt - inconsistent cluster vectors at entry "
                << entry << std::endl;
      chain.ResetBranchAddresses();
      return false;
    }

    for (std::size_t cluster = 0; cluster < emcal_energy->size(); ++cluster)
    {
      if (emcal_energy->at(cluster) <= emcal_energy_threshold)
      {
        continue;
      }
      histograms[bin]->Fill(hcal_total_energy->at(cluster));
      ++counts[bin];
    }
  }
  chain.ResetBranchAddresses();

  for (int bin = 0; bin < n_truth_pt_bins; ++bin)
  {
    const double integral = histograms[bin]->Integral();
    if (counts[bin] <= 0 || integral <= 0.0)
    {
      const double pt_low = truth_pt_minimum + bin * truth_pt_bin_width;
      const double pt_high = pt_low + truth_pt_bin_width;
      std::cerr << "PlotHCalEnergyByTruthPt - no selected TopoClusters for truth pT bin "
                << pt_low << "-" << pt_high << " GeV" << std::endl;
      return false;
    }
    histograms[bin]->Scale(1.0 / integral);
  }
  return true;
}

void draw_panel(
    TH1D &gamma_histogram,
    TH1D &pi0_histogram,
    const Long64_t gamma_count,
    const Long64_t pi0_count,
    const int pt_low)
{
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.04);
  gPad->SetBottomMargin(0.14);
  gPad->SetTopMargin(0.10);
  gPad->SetLogy();

  gamma_histogram.SetTitle(Form(
      "%d #leq p_{T}^{truth} < %d GeV;E_{HCalIn}^{cluster} + E_{HCalOut}^{cluster} [GeV];Normalized TopoClusters",
      pt_low,
      pt_low + 1));
  gamma_histogram.SetLineColor(kAzure + 2);
  gamma_histogram.SetLineWidth(3);
  pi0_histogram.SetLineColor(kOrange + 7);
  pi0_histogram.SetLineWidth(3);

  const double common_maximum = std::max(
      gamma_histogram.GetMaximum(),
      pi0_histogram.GetMaximum());
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

  TLegend legend(0.48, 0.67, 0.94, 0.89);
  legend.SetHeader(Form("E_{EMCal}^{cluster} > %.0f GeV", emcal_energy_threshold));
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetTextSize(0.038);
  legend.AddEntry(
      &gamma_histogram,
      Form("Single #gamma (N=%lld)", gamma_count),
      "l");
  legend.AddEntry(
      &pi0_histogram,
      Form("Single #pi^{0} (N=%lld)", pi0_count),
      "l");
  legend.DrawClone();
}
}  // namespace

int PlotHCalEnergyByTruthPt(
    const std::string gamma_input_pattern = "/sphenix/user/ryotaro/DirectPhotonAnalysis/TopoClusterHCalStudy/output/merge/25to35GeV/topocluster_hcal_gamma_merged.root",
    const std::string pi0_input_pattern = "/sphenix/user/ryotaro/DirectPhotonAnalysis/TopoClusterHCalStudy/output/merge/25to35GeV/topocluster_hcal_pi0_merged.root",
    const std::string output_directory = "/sphenix/user/ryotaro/DirectPhotonAnalysis/TopoClusterHCalStudy/output/plot/25to35GeV")
{
  const std::filesystem::path output_path(output_directory);
  std::error_code directory_error;
  std::filesystem::create_directories(output_path, directory_error);
  if (directory_error)
  {
    std::cerr << "PlotHCalEnergyByTruthPt - could not create output directory: "
              << directory_error.message() << std::endl;
    return EXIT_FAILURE;
  }

  TChain gamma_chain("topocluster_tree");
  TChain pi0_chain("topocluster_tree");
  if (!load_chain(gamma_chain, gamma_input_pattern) ||
      !load_chain(pi0_chain, pi0_input_pattern))
  {
    return EXIT_FAILURE;
  }

  int truth_pt_minimum = 0;
  if (!determine_truth_pt_minimum(
          gamma_chain,
          pi0_chain,
          truth_pt_minimum))
  {
    return EXIT_FAILURE;
  }

  HistogramArray gamma_histograms;
  HistogramArray pi0_histograms;
  CountArray gamma_counts{};
  CountArray pi0_counts{};
  create_histograms(gamma_histograms, "gamma", truth_pt_minimum);
  create_histograms(pi0_histograms, "pi0", truth_pt_minimum);
  if (!fill_histograms(
          gamma_chain,
          gamma_histograms,
          gamma_counts,
          truth_pt_minimum) ||
      !fill_histograms(
          pi0_chain,
          pi0_histograms,
          pi0_counts,
          truth_pt_minimum))
  {
    return EXIT_FAILURE;
  }

  gStyle->SetOptStat(0);
  TCanvas canvas(
      "hcal_energy_by_truth_pt",
      "HCAL energy by truth pT",
      2000,
      1400);
  canvas.Divide(5, 2, 0.001, 0.001);

  for (int bin = 0; bin < n_truth_pt_bins; ++bin)
  {
    canvas.cd(bin + 1);
    draw_panel(
        *gamma_histograms[bin],
        *pi0_histograms[bin],
        gamma_counts[bin],
        pi0_counts[bin],
        truth_pt_minimum + bin);
  }

  const int truth_pt_maximum = truth_pt_minimum + n_truth_pt_bins;
  const std::string output_file =
      (output_path /
       ("hcal_total_energy_truth_pt_" +
        std::to_string(truth_pt_minimum) + "to" +
        std::to_string(truth_pt_maximum) + ".pdf"))
          .string();
  canvas.SaveAs(output_file.c_str());

  for (int bin = 0; bin < n_truth_pt_bins; ++bin)
  {
    std::cout << "PlotHCalEnergyByTruthPt - "
              << truth_pt_minimum + bin << "-"
              << truth_pt_minimum + bin + 1 << " GeV: gamma="
              << gamma_counts[bin] << ", pi0=" << pi0_counts[bin]
              << std::endl;
  }
  std::cout << "PlotHCalEnergyByTruthPt - output: "
            << output_file << std::endl;
  return EXIT_SUCCESS;
}
