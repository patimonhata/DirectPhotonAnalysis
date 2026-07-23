#include "Utility/sPhenixStyle.C"

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TString.h>

#include <array>
#include <iostream>
#include <memory>
#include <string>

int PlotPairEnergyAsymmetry(
    const std::string input_directory = "Pi0Reconstruction/output",
    const std::string output_base =
        "Pi0Reconstruction/output/5GeV_pi0_eta0_pair_energy_asymmetry_SPLIT")
{
  SetsPhenixStyle();
  TH1::AddDirectory(false);

  const std::array<std::string, 3> input_names = {
      "5GeV_pi0_eta0_over100MeV_rehist_SPLIT.root",
      "5GeV_pi0_eta0_over300MeV_rehist_SPLIT.root",
      "5GeV_pi0_eta0_over500MeV_rehist_SPLIT.root"};
  const std::array<int, 3> energy_thresholds_mev = {100, 300, 500};
  const std::array<Color_t, 3> colors = {
      kRed + 1,
      kBlue + 1,
      kBlack};

  std::array<std::unique_ptr<TH1D>, 3> histograms;
  for (std::size_t i = 0; i < input_names.size(); ++i)
  {
    const std::string input_path = input_directory + "/" + input_names[i];
    std::unique_ptr<TFile> input(TFile::Open(input_path.c_str(), "READ"));
    if (!input || input->IsZombie())
    {
      std::cerr << "Failed to open " << input_path << std::endl;
      return 1;
    }

    TH1D *source = dynamic_cast<TH1D *>(input->Get("h_pair_e_asym"));
    if (!source)
    {
      std::cerr << "Missing TH1D h_pair_e_asym in " << input_path << std::endl;
      return 1;
    }

    histograms[i].reset(
        static_cast<TH1D *>(source->Clone(("h_pair_e_asym_" + std::to_string(i)).c_str())));
    histograms[i]->SetDirectory(nullptr);
    histograms[i]->SetStats(false);
    histograms[i]->SetLineColor(colors[i]);
    histograms[i]->SetFillStyle(0);
  }

  TCanvas canvas("c_pair_energy_asymmetry", "Pair energy asymmetry", 1000, 800);

  histograms[0]->SetTitle("");
  histograms[0]->GetXaxis()->SetTitle("|E_{1}-E_{2}|/(E_{1}+E_{2})");
  histograms[0]->GetYaxis()->SetTitle("Pairs");
  histograms[0]->GetXaxis()->SetRangeUser(0.0, 1.0);
  histograms[0]->SetMinimum(0.0);
  histograms[0]->SetMaximum(4000.0);
  histograms[0]->Draw("HIST");
  histograms[1]->Draw("HIST SAME");
  histograms[2]->Draw("HIST SAME");

  TLine truth_line(0.0, 2000.0, 1.0, 2000.0);
  truth_line.SetLineColor(kGray + 2);
  truth_line.SetLineStyle(2);
  truth_line.SetLineWidth(2);
  truth_line.Draw("SAME");

  TLegend legend(0.20, 0.57, 0.95, 0.82);
  legend.SetFillColor(kWhite);
  legend.SetFillStyle(1001);
  legend.SetBorderSize(0);
  legend.SetHeader("#varepsilon #approx entries / 10^{5}", "L");
  const std::array<std::size_t, 3> legend_order = {2, 1, 0};
  for (const std::size_t i : legend_order)
  {
    const double efficiency_percent = histograms[i]->GetEntries() / 1000.0;
    const TString legend_label = TString::Format(
        "E_{cluster}>%d MeV, #varepsilon=%.1f%%",
        energy_thresholds_mev[i],
        efficiency_percent);
    legend.AddEntry(histograms[i].get(), legend_label.Data(), "l");
  }
  legend.AddEntry(&truth_line, "Expected truth distribution", "l");
  legend.Draw();

  TLatex label;
  label.SetNDC();
  label.SetTextAlign(13);
  label.DrawLatex(0.20, 0.90, "#it{#bf{sPHENIX}} Internal");
  label.DrawLatex(0.52, 0.90, "Single #pi^{0}, 5 GeV, #eta = 0");

  canvas.RedrawAxis();
  canvas.SaveAs((output_base + ".pdf").c_str());

  std::cout << "Saved " << output_base << ".pdf and " << output_base << ".png"
            << std::endl;
  return 0;
}
