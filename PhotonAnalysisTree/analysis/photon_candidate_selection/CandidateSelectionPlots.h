#ifndef PHOTON_CANDIDATE_CANDIDATESELECTIONPLOTS_H
#define PHOTON_CANDIDATE_CANDIDATESELECTIONPLOTS_H

#include "CandidateSelectionCommon.h"

namespace
{
std::unique_ptr<TPad> make_plot_pad(const std::string& name, bool log_y = false, double height = kPlotPadHeight, double width = 1.0)
{
  auto pad = std::make_unique<TPad>(name.c_str(), "", 0.0, 0.0, width, height);
  pad->SetLogy(log_y);
  pad->Draw();
  pad->cd();
  return pad;
}

void style_axes(TH1* histogram, const char* y_title)
{
  histogram->GetXaxis()->SetTitle("Candidate Cluster E_{T} [GeV]");
  histogram->GetYaxis()->SetTitle(y_title);
  for (TAxis* axis : {histogram->GetXaxis(), histogram->GetYaxis()})
  {
    axis->CenterTitle();
  }
}

void scale_axes_to_canvas(TAxis* x_axis, TAxis* y_axis, double pad_height)
{
  const double scale = 1.0 / pad_height;
  x_axis->SetLabelSize(gStyle->GetLabelSize("X") * scale);
  x_axis->SetTitleSize(gStyle->GetTitleSize("X") * scale);
  y_axis->SetLabelSize(gStyle->GetLabelSize("Y") * scale);
  y_axis->SetTitleSize(gStyle->GetTitleSize("Y") * scale);
}

void scale_axes_to_canvas(TH1* histogram, double pad_height)
{
  scale_axes_to_canvas(histogram->GetXaxis(), histogram->GetYaxis(), pad_height);
}

void scale_legend_to_canvas(TLegend& legend, double pad_height)
{
  legend.SetTextSize(gStyle->GetLegendTextSize() / pad_height);
}

std::unique_ptr<TH1D> make_plot_frame(const TH1D& source, const std::string& name, const char* y_title, double minimum, double maximum)
{
  auto frame = std::unique_ptr<TH1D>(static_cast<TH1D*>(source.Clone(name.c_str())));
  frame->SetDirectory(nullptr);
  frame->Reset("ICES");
  frame->SetBins(source.GetNbinsX(), 0.0, kPlotEtMaximum);
  frame->SetMinimum(minimum);
  frame->SetMaximum(maximum);
  style_axes(frame.get(), y_title);
  return frame;
}

struct PlotCaption
{
  std::string family;
  std::string selection;
  double min_cluster_energy = 0.0;
  double pi0_partner_min_energy = 0.0;
  double eta_partner_min_energy = 0.0;
  double pi0_mass_min = 0.0;
  double pi0_mass_max = 0.0;
  double eta_mass_min = 0.0;
  double eta_mass_max = 0.0;
};
PlotCaption plot_caption;

std::vector<std::string> caption_lines(const PlotCaption& value)
{
  const auto number = [](double x) { std::ostringstream out; out << std::fixed << std::setprecision(2) << x; return out.str(); };
  const std::string sample = value.family == "jet" ? "Pythia 8 p+p Jet samples" : "Pythia 8 p+p PhotonJet samples";
  const std::string pi0_tagging = "#pi^{0} tagging: E_{partner}>" + number(value.pi0_partner_min_energy) + " GeV; " + number(value.pi0_mass_min) + "< m <" + number(value.pi0_mass_max) + " GeV";
  const std::string eta_tagging = "#eta  tagging: E_{partner}>" + number(value.eta_partner_min_energy) + " GeV; " + number(value.eta_mass_min) + "< m <" + number(value.eta_mass_max) + " GeV";
  const std::string candidate = value.selection == kSelectionLabels.front()
      ? "Candidate Cluster: 5 < E_{T} < 35 GeV, |#eta| < 0.7"
      : "Candidate Cluster: after " + value.selection;
  return {
      "#it{#bf{sPHENIX}} Internal",
      // "#it{#bf{sPHENIX}} Simulation Work-in-Progress",
      sample + ", |z_{vertex}^{truth}| < 60 cm",
      candidate,
      pi0_tagging,
      eta_tagging};
}

void draw_annotations()
{
  TLatex label;
  label.SetNDC();
  label.SetTextAlign(13);
  label.SetTextColor(kBlack);
  const auto lines = caption_lines(plot_caption);
  constexpr double line_spacing = 0.06;
  for (std::size_t i = 0; i < lines.size(); ++i)
  {
    const double y = 0.975 - line_spacing * i - (i >= 2 ? 0.01 : 0.0);
    label.DrawLatex(0.055, y, lines[i].c_str());
  }
}

// Presentation-only labels: keep stored category names and work-in-progress plots unchanged.
std::string compact_legend_label(std::string label)
{
  for (const std::string word : {"Single contaminated", "single contaminated"})
  {
    const auto pos = label.find(word);
    if (pos != std::string::npos) label.replace(pos, word.size(), word.front() == 'S' ? "Contaminated" : "contaminated");
  }
  return label;
}

void add_readable_legend_entry(TLegend& legend, TObject* object, const std::string& source, const char* option)
{
  const std::string label = compact_legend_label(source);
  auto split = label.find(", mass ");
  if (split != std::string::npos)
  {
    legend.AddEntry(object, label.substr(0, split).c_str(), option);
    legend.AddEntry(static_cast<TObject*>(nullptr), ("  " + label.substr(split + 2)).c_str(), "");
  }
  else if ((split = label.find(" or no CEMC deposit")) != std::string::npos)
  {
    legend.AddEntry(object, label.substr(0, split).c_str(), option);
    legend.AddEntry(static_cast<TObject*>(nullptr), "  or no CEMC deposit", "");
  }
  else legend.AddEntry(object, label.c_str(), option);
}

constexpr int kDetailedCanvasWidth = 2000;
constexpr double kDetailedPlotWidth = 1600.0 * kCompositionPlotWidth / kDetailedCanvasWidth;
constexpr double kDetailedLegendLeft = kDetailedPlotWidth + 0.02;

void draw_spectrum(const std::array<std::unique_ptr<TH1D>, kSpectrumCount>& density, const std::vector<std::size_t>& indices,
                   const std::string& output_path, bool detailed)
{
  SetsPhenixStyle();
  const bool side_legend = detailed;
  TCanvas canvas(("c_" + std::string(detailed ? "detailed" : "summary") + "_region_a_topology_spectrum").c_str(), "", side_legend ? kDetailedCanvasWidth : kCanvasWidth, kCanvasHeight);
  auto plot_pad = make_plot_pad(std::string(detailed ? "detailed" : "summary") + "_region_a_topology_spectrum_pad", true, kPlotPadHeight,
                                side_legend ? kDetailedPlotWidth : 1.0);
  scale_axes_to_canvas(density[indices.front()].get(), kPlotPadHeight);
  density[indices.front()]->GetXaxis()->SetTitleOffset(1.0);
  density[indices.front()]->GetYaxis()->SetTitleOffset(1.0);
  density[indices.front()]->SetMinimum(kSpectrumYMinimum);
  density[indices.front()]->SetMaximum(kSpectrumYMaximum);
  density[indices.front()]->Draw("HIST");
  for (std::size_t position = 1; position < indices.size(); ++position) density[indices[position]]->Draw("HIST SAME");
  if (side_legend) canvas.cd();
  TLegend legend(side_legend ? kDetailedLegendLeft : 0.46, side_legend ? 0.05 : 0.38, side_legend ? 0.99 : 0.98, side_legend ? 0.95 : 0.94);
  if (!side_legend) scale_legend_to_canvas(legend, kPlotPadHeight);
  legend.SetFillStyle(0);
  for (std::size_t index : indices) add_readable_legend_entry(legend, density[index].get(), kLabels[index], "l");
  legend.Draw();
  plot_pad->RedrawAxis();
  canvas.cd();
  draw_annotations();
  canvas.SaveAs(output_path.c_str());
}

std::vector<std::unique_ptr<TH1D>> make_fractions(const std::array<std::unique_ptr<TH1D>, kSpectrumCount>& weighted_pb,
                                                  const std::vector<std::size_t>& indices, const std::string& prefix)
{
  std::vector<std::unique_ptr<TH1D>> result;
  for (std::size_t index : indices)
  {
    const std::string name = prefix + kKeys[index] + "_fraction";
    auto histogram = std::unique_ptr<TH1D>(static_cast<TH1D*>(weighted_pb[index]->Clone(name.c_str())));
    histogram->SetDirectory(nullptr);
    histogram->Reset("ICES");
    for (int bin = 0; bin <= histogram->GetNbinsX() + 1; ++bin)
    {
      const double denominator = weighted_pb[1]->GetBinContent(bin);
      if (denominator == 0.0) continue;
      const double numerator = weighted_pb[index]->GetBinContent(bin);
      const double value = numerator / denominator;
      const double numerator_w2 = std::pow(weighted_pb[index]->GetBinError(bin), 2);
      const double denominator_w2 = std::pow(weighted_pb[1]->GetBinError(bin), 2);
      const double complement_w2 = std::max(0.0, denominator_w2 - numerator_w2);
      const double variance = (std::pow(1.0 - value, 2) * numerator_w2 + std::pow(value, 2) * complement_w2) / std::pow(denominator, 2);
      histogram->SetBinContent(bin, value);
      histogram->SetBinError(bin, std::sqrt(std::max(0.0, variance)));
    }
    histogram->SetLineColor(kColors[index]);
    histogram->SetMarkerColor(kColors[index]);
    style_axes(histogram.get(), "Fraction");
    result.push_back(std::move(histogram));
  }
  return result;
}

void draw_fraction_lines(const std::vector<std::unique_ptr<TH1D>>& fractions, const std::vector<std::size_t>& indices,
                         const std::string& output_path, bool detailed)
{
  SetsPhenixStyle();
  const bool side_legend = detailed;
  TCanvas canvas(("c_" + std::string(detailed ? "detailed" : "summary") + "_region_a_topology_fraction").c_str(), "", side_legend ? kDetailedCanvasWidth : kCanvasWidth, kCanvasHeight);
  auto plot_pad = make_plot_pad(std::string(detailed ? "detailed" : "summary") + "_region_a_topology_fraction_pad", false, kPlotPadHeight,
                                side_legend ? kDetailedPlotWidth : 1.0);
  scale_axes_to_canvas(fractions.front().get(), kPlotPadHeight);
  fractions.front()->GetXaxis()->SetTitleOffset(1.0);
  fractions.front()->GetYaxis()->SetTitleOffset(1.0);
  fractions.front()->SetMinimum(0.0);
  fractions.front()->SetMaximum(1.05);
  fractions.front()->Draw("E1");
  for (std::size_t index = 1; index < fractions.size(); ++index) fractions[index]->Draw("E1 SAME");
  if (side_legend) canvas.cd();
  TLegend legend(side_legend ? kDetailedLegendLeft : 0.46, side_legend ? 0.10 : 0.50, side_legend ? 0.99 : 0.98, side_legend ? 0.90 : 0.94);
  if (!side_legend) scale_legend_to_canvas(legend, kPlotPadHeight);
  legend.SetFillStyle(0);
  for (std::size_t index = 0; index < fractions.size(); ++index) add_readable_legend_entry(legend, fractions[index].get(), kLabels[indices[index]], "lep");
  legend.Draw();
  plot_pad->RedrawAxis();
  canvas.cd();
  draw_annotations();
  canvas.SaveAs(output_path.c_str());
}

void draw_fraction_stack(std::vector<std::unique_ptr<TH1D>>& fractions, const std::vector<std::size_t>& indices,
                         const std::string& output_path, bool detailed)
{
  SetsPhenixStyle();
  TCanvas canvas(("c_" + std::string(detailed ? "detailed" : "summary") + "_region_a_topology_fraction_stack").c_str(), "", detailed ? kDetailedCanvasWidth : 1600, kCanvasHeight);
  auto plot_pad = make_plot_pad(std::string(detailed ? "detailed" : "summary") + "_region_a_topology_fraction_stack_pad", false, kPlotPadHeight,
                                detailed ? kDetailedPlotWidth : kCompositionPlotWidth);
  THStack stack(("stack_" + std::string(detailed ? "detailed" : "summary") + "_region_a_topology_fraction").c_str(), "");
  for (std::size_t index = 0; index < fractions.size(); ++index)
  {
    fractions[index]->SetFillColor(kColors[indices[index]]);
    fractions[index]->SetLineColor(kBlack);
    stack.Add(fractions[index].get());
  }
  stack.SetMinimum(0.0);
  stack.SetMaximum(1.05);
  stack.Draw("HIST");
  scale_axes_to_canvas(stack.GetXaxis(), stack.GetYaxis(), kPlotPadHeight);
  stack.GetXaxis()->SetTitleOffset(1.0);
  stack.GetYaxis()->SetTitleOffset(1.0);
  stack.GetXaxis()->SetTitle("Candidate Cluster E_{T} [GeV]");
  stack.GetYaxis()->SetTitle("Fraction");
  stack.GetXaxis()->CenterTitle();
  stack.GetYaxis()->CenterTitle();
  canvas.cd();
  TLegend legend(detailed ? kDetailedLegendLeft : 0.72, detailed ? 0.05 : 0.30, 0.99, detailed ? 0.95 : 0.70);
  legend.SetFillStyle(0);
  for (std::size_t index = fractions.size(); index-- > 0;) add_readable_legend_entry(legend, fractions[index].get(), kLabels[indices[index]], "f");
  legend.Draw();
  plot_pad->RedrawAxis();
  canvas.cd();
  draw_annotations();
  canvas.SaveAs(output_path.c_str());
}
}

#endif
