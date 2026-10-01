#ifndef PI0_BDT_HISTOGRAMS_H
#define PI0_BDT_HISTOGRAMS_H
#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TNamed.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <memory>
#include <string>
#include <vector>

namespace pi0_bdt
{
constexpr std::array<const char*, 5> keys = {"separated", "merged", "single_contaminated", "missing", "other"};
constexpr std::array<const char*, 5> labels = {"Separated", "Merged", "Single contaminated", "Missing", "Other"};
constexpr std::array<const char*, 4> isolation = {"all", "isolated", "nonisolated", "gap"};
constexpr std::array<int, 5> colors = {kAzure + 7, kMagenta + 1, kCyan + 2, kGreen + 2, kGray + 1};
struct Histograms
{
  struct Cell
  {
    TH1D count, weighted, summary;
    TH2D score_et;
    Cell(const std::string& name): count((name + "_count").c_str(), ";BDT score;Clusters / bin", 100, 0, 1),
      weighted((name + "_pb").c_str(), ";BDT score;Weighted yield [pb/bin]", 100, 0, 1),
      summary((name + "_summary").c_str(), "", 4, 0, 4),
      score_et((name + "_score_et_pb").c_str(), ";BDT score;Anchor E_{T} [GeV]", 100, 0, 1, 70, 0, 35)
    {
      for (TH1* h : {static_cast<TH1*>(&count), static_cast<TH1*>(&weighted), static_cast<TH1*>(&summary), static_cast<TH1*>(&score_et)})
      { h->SetDirectory(nullptr); h->Sumw2(); }
      summary.GetXaxis()->SetBinLabel(1, "selected_count"); summary.GetXaxis()->SetBinLabel(2, "invalid_bdt_count");
      summary.GetXaxis()->SetBinLabel(3, "selected_pb"); summary.GetXaxis()->SetBinLabel(4, "invalid_bdt_pb");
    }
  };
  std::array<std::array<std::unique_ptr<Cell>, 5>, 4> cells;
  Histograms()
  {
    for (int i = 0; i < 4; ++i) for (int k = 0; k < 5; ++k)
      cells[i][k] = std::make_unique<Cell>(std::string("h_") + isolation[i] + "_" + keys[k]);
  }
  bool fill(int category, double score, bool valid, bool isolated, bool nonisolated, double et, double weight)
  {
    if (category < 0 || category >= 5 || (isolated && nonisolated) || !std::isfinite(weight)) return false;
    const int group = isolated ? 1 : nonisolated ? 2 : 3;
    for (int i : {0, group})
    {
      auto& c = *cells[i][category]; c.summary.Fill(0.5); c.summary.Fill(2.5, weight);
      if (!valid || !std::isfinite(score)) { c.summary.Fill(1.5); c.summary.Fill(3.5, weight); continue; }
      // Include the probability endpoint 1 in the final visible bin; preserve true overflow.
      const double x = score == 1.0 ? std::nextafter(1.0, 0.0) : score;
      c.count.Fill(x); c.weighted.Fill(x, weight); c.score_et.Fill(x, et, weight);
    }
    return true;
  }
  bool write(TDirectory& directory) const
  {
    TDirectory::TContext context(&directory);
    for (const auto& row : cells) for (const auto& c : row)
      if (c->count.Write() <= 0 || c->weighted.Write() <= 0 || c->summary.Write() <= 0 || c->score_et.Write() <= 0) return false;
    return true;
  }
  bool add(TDirectory& directory)
  {
    for (auto& row : cells) for (auto& c : row)
      for (TH1* target : {static_cast<TH1*>(&c->count), static_cast<TH1*>(&c->weighted), static_cast<TH1*>(&c->summary), static_cast<TH1*>(&c->score_et)})
      {
        auto* source = directory.Get<TH1>(target->GetName());
        if (!source || source->GetDimension() != target->GetDimension() || source->GetNcells() != target->GetNcells()) return false;
        if (source->GetXaxis()->GetXmin() != target->GetXaxis()->GetXmin() || source->GetXaxis()->GetXmax() != target->GetXaxis()->GetXmax() ||
            source->GetYaxis()->GetXmin() != target->GetYaxis()->GetXmin() || source->GetYaxis()->GetXmax() != target->GetYaxis()->GetXmax()) return false;
        if (!target->Add(source)) return false;
      }
    return true;
  }
  bool valid() const
  {
    for (int k = 0; k < 5; ++k)
      for (int kind = 0; kind < 3; ++kind)
      {
        auto hist = [&](int i) -> const TH1D& { return kind == 0 ? cells[i][k]->count : kind == 1 ? cells[i][k]->weighted : cells[i][k]->summary; };
        for (int b = 0; b <= hist(0).GetNbinsX() + 1; ++b)
        {
          double sum = 0; for (int i = 1; i < 4; ++i) sum += hist(i).GetBinContent(b);
          if (std::abs(hist(0).GetBinContent(b) - sum) > 1e-8 * std::max(1.0, std::abs(sum))) return false;
        }
      }
    return true;
  }
  void panel(int group, bool normalized, const std::string& selection, bool requires_isolation) const
  {
    std::vector<std::unique_ptr<TH1D>> plots;
    double maximum = 0; double selected = 0; double invalid = 0;
    for (int k = 0; k < 5; ++k)
    {
      const auto& c = *cells[group][k]; selected += c.summary.GetBinContent(1); invalid += c.summary.GetBinContent(2);
      auto h = std::unique_ptr<TH1D>(static_cast<TH1D*>(c.weighted.Clone())); h->SetDirectory(nullptr); h->SetStats(false);
      if (normalized && h->Integral() > 0) h->Scale(1.0 / h->Integral());
      h->SetLineColor(colors[k]); h->SetLineWidth(2); maximum = std::max(maximum, h->GetMaximum()); plots.push_back(std::move(h));
    }
    gPad->SetLeftMargin(0.16); gPad->SetRightMargin(0.04); gPad->SetBottomMargin(0.14); gPad->SetTopMargin(0.20);
    gPad->SetLogy(!normalized && maximum > 0);
    auto& first = *plots.front(); first.SetTitle("");
    for (TAxis* axis : {first.GetXaxis(), first.GetYaxis()})
    { axis->SetTitleFont(42); axis->SetLabelFont(42); axis->SetTitleSize(0.045); axis->SetLabelSize(0.04); }
    first.GetYaxis()->SetTitleOffset(1.55); first.GetXaxis()->SetTitleOffset(1.1);
    first.GetYaxis()->SetTitle(normalized ? "Unit-area fraction / bin" : "Weighted yield [pb/bin]");
    first.SetMaximum(maximum > 0 ? maximum * (normalized ? 1.6 : 30) : 1);
    first.SetMinimum(!normalized && maximum > 0 ? maximum * 1e-6 : 0);
    first.DrawCopy("HIST");
    TLegend legend(0.59, 0.55, 0.95, 0.78); legend.SetBorderSize(0); legend.SetTextSize(0.032); legend.SetTextFont(42); legend.SetFillStyle(0);
    for (int k = 0; k < 5; ++k) { auto* drawn = plots[k]->DrawCopy("HIST SAME"); legend.AddEntry(drawn, labels[k], "l"); }
    legend.DrawClone();
    TLatex text; text.SetNDC(); text.SetTextFont(42); text.SetTextSize(0.034);
    const char* group_labels[] = {"All", "Isolated", "Non-isolated", "Gap"};
    text.DrawLatex(0.16, 0.95, (selection + " | " + group_labels[group]).c_str());
    text.SetTextSize(0.03);
    std::string counts = "Selected N=" + std::to_string(static_cast<long long>(selected)) + ", invalid BDT=" + std::to_string(static_cast<long long>(invalid));
    if (group == 0)
    {
      double gap = 0; for (const auto& cell : cells[3]) gap += cell->summary.GetBinContent(1);
      counts += ", gap=" + std::to_string(static_cast<long long>(gap));
    }
    text.DrawLatex(0.16, 0.88, counts.c_str());
    if (selected == 0)
    {
      text.DrawLatex(0.24, 0.46, "0 entries");
      if (requires_isolation && group == 2) text.DrawLatex(0.24, 0.40, "Isolated required by selection");
    }
    gPad->RedrawAxis();
  }
  bool export_files(const std::string& directory, const std::string& selection, bool requires_isolation) const
  {
    if (!valid()) return false;
    TFile output((directory + "/pi0_bdt_score.root").c_str(), "RECREATE");
    if (output.IsZombie() || !write(output)) return false;
    TNamed definition("population", "composition pi0: valid anchor, main fraction > 0.5, no multiple-origin majority"); definition.Write();
    TNamed cut("selection", selection.c_str()); cut.Write();
    TNamed normalization("normalization", "cross-section and generator weighted; unit area uses visible score bins; invalid BDT counted separately"); normalization.Write();
    output.Close();
    std::ofstream csv(directory + "/pi0_bdt_score_summary.csv");
    csv << "isolation,topology,selected_count,invalid_bdt_count,selected_pb,invalid_bdt_pb,valid_bdt_pb,underflow_count,overflow_count\n" << std::setprecision(17);
    for (int i = 0; i < 4; ++i) for (int k = 0; k < 5; ++k)
    {
      const auto& c = *cells[i][k]; csv << isolation[i] << ',' << keys[k];
      for (int b = 1; b <= 4; ++b) csv << ',' << c.summary.GetBinContent(b);
      csv << ',' << c.weighted.Integral(0, 101) << ',' << c.count.GetBinContent(0) << ',' << c.count.GetBinContent(101) << '\n';
    }
    if (!csv.good()) return false;
    TCanvas standard("pi0_bdt_standard", "", 1600, 700); standard.Divide(2, 1);
    for (int j = 0; j < 2; ++j) { standard.cd(j + 1); panel(0, j == 1, selection, requires_isolation); }
    standard.SaveAs((directory + "/pi0_bdt_score_by_topology.pdf").c_str());
    TCanvas split("pi0_bdt_isolation", "", 2100, 1300); split.Divide(3, 2);
    for (int j = 0; j < 2; ++j) for (int i = 0; i < 3; ++i) { split.cd(j * 3 + i + 1); panel(i, j == 1, selection, requires_isolation); }
    split.SaveAs((directory + "/pi0_bdt_score_by_topology_isolation.pdf").c_str());
    return true;
  }
};
}
#endif
