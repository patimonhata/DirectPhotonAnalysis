#ifndef PI0_BDT_HISTOGRAMS_H
#define PI0_BDT_HISTOGRAMS_H
#include "../../macro/Utilities/sPhenixStyle.C"
#include <TCanvas.h>
#include <TPad.h>
#include <TTree.h>
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
#include <sstream>
#include <string>
#include <vector>

namespace pi0_bdt
{
constexpr std::array<const char*, 5> keys = {"separated", "merged", "single_contaminated", "missing", "other"};
constexpr std::array<const char*, 5> labels = {"Separated", "Merged", "Single contaminated", "Missing", "Other"};
constexpr std::array<const char*, 4> isolation = {"all", "isolated", "nonisolated", "gap"};
constexpr std::array<double, 6> et_edges = {5, 6, 8, 10, 18, 30};
constexpr std::array<int, 5> colors = {kAzure + 7, kMagenta + 1, kCyan + 2, kGreen + 2, kGray + 1};
// Reconstruct the page caption from existing merge metadata for standalone redraws.
inline std::vector<std::string> read_plot_caption(TDirectory& input, const std::string& selection)
{
  std::vector<std::string> result = {"#it{#bf{sPHENIX}} Internal", "Candidate Cluster: after " + selection};
  auto* tree = input.Get<TTree>("metadata");
  if (!tree || tree->GetEntries() != 1) return result;
  std::string family_value;
  std::string* family = &family_value;
  const std::array<const char*, 6> names = {"pi0_partner_min_energy", "eta_partner_min_energy", "pi0_mass_min", "pi0_mass_max", "eta_mass_min", "eta_mass_max"};
  std::array<double, 6> values{};
  if (!tree->GetBranch("family") || tree->SetBranchAddress("family", &family) < 0) return result;
  bool valid = true;
  for (std::size_t i = 0; i < names.size(); ++i)
    valid = (tree->GetBranch(names[i]) && tree->SetBranchAddress(names[i], &values[i]) >= 0) && valid;
  valid = valid && tree->GetEntry(0) > 0 && family;
  if (valid)
  {
    const auto number = [](double x) { std::ostringstream out; out << std::fixed << std::setprecision(2) << x; return out.str(); };
    result = {
        "#it{#bf{sPHENIX}} Internal",
        std::string(*family == "jet" ? "Pythia 8 p+p Jet samples" : "Pythia 8 p+p PhotonJet samples") + ", |z_{vertex}^{truth}| < 60 cm",
        selection == "Kinematic" ? "Candidate Cluster: 5 < E_{T} < 35 GeV, |#eta| < 0.7" : "Candidate Cluster: after " + selection,
        "#pi^{0} tagging: E_{partner}>" + number(values[0]) + " GeV; " + number(values[2]) + "< m <" + number(values[3]) + " GeV",
        "#eta  tagging: E_{partner}>" + number(values[1]) + " GeV; " + number(values[4]) + "< m <" + number(values[5]) + " GeV"};
  }
  tree->ResetBranchAddresses();
  return result;
}

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
  // slice=-1 preserves the original full-ET histogram, including 30--35 GeV.
  std::unique_ptr<TH1D> score_histogram(int group, int category, int slice = -1) const
  {
    const auto& cell = *cells.at(group).at(category);
    std::unique_ptr<TH1D> result;
    if (slice < 0) result.reset(static_cast<TH1D*>(cell.weighted.Clone()));
    else
    {
      const double low = et_edges.at(slice), high = et_edges.at(slice + 1);
      const auto* axis = cell.score_et.GetYaxis();
      const int first = axis->FindFixBin(low), last = axis->FindFixBin(high) - 1;
      const std::string name = std::string(cell.weighted.GetName()) + "_et_slice_" + std::to_string(slice);
      result.reset(cell.score_et.ProjectionX(name.c_str(), first, last, "e"));
    }
    result->SetDirectory(nullptr);
    return result;
  }
  void panel(int group, bool normalized, bool requires_isolation, int slice = -1) const
  {
    std::vector<std::unique_ptr<TH1D>> plots;
    double maximum = 0; double selected = 0; double invalid = 0; bool has_scores = false;
    for (int k = 0; k < 5; ++k)
    {
      const auto& c = *cells[group][k]; selected += c.summary.GetBinContent(1); invalid += c.summary.GetBinContent(2);
      auto h = score_histogram(group, k, slice); h->SetStats(false);
      for (int b = 0; b <= h->GetNbinsX() + 1; ++b)
        has_scores = has_scores || h->GetBinContent(b) != 0.0 || h->GetBinError(b) != 0.0;
      if (normalized && h->Integral() > 0) h->Scale(1.0 / h->Integral());
      h->SetLineColor(colors[k]); h->SetLineWidth(2); maximum = std::max(maximum, h->GetMaximum()); plots.push_back(std::move(h));
    }
    gPad->SetLeftMargin(0.18); gPad->SetRightMargin(0.04); gPad->SetBottomMargin(0.16); gPad->SetTopMargin(0.23);
    gPad->SetLogy(!normalized && maximum > 0);
    auto& first = *plots.front(); first.SetTitle("");
    for (TAxis* axis : {first.GetXaxis(), first.GetYaxis()})
    { axis->SetTitleFont(42); axis->SetLabelFont(42); axis->SetTitleSize(0.05); axis->SetLabelSize(0.05); axis->CenterTitle(); }
    first.GetYaxis()->SetTitleOffset(1.5); first.GetXaxis()->SetTitleOffset(1.1);
    first.GetYaxis()->SetTitle(normalized ? "Unit-area fraction / bin" : "Weighted yield [pb/bin]");
    first.SetMaximum(maximum > 0 ? maximum * (normalized ? 1.6 : 30) : 1);
    first.SetMinimum(!normalized && maximum > 0 ? maximum * 1e-6 : 0);
    first.DrawCopy("HIST");
    TLegend legend(0.22, 0.43, 0.65, 0.74); legend.SetBorderSize(0); legend.SetTextSize(0.05); legend.SetTextFont(42); legend.SetFillStyle(0);
    for (int k = 0; k < 5; ++k) { auto* drawn = plots[k]->DrawCopy("HIST SAME"); legend.AddEntry(drawn, k == 2 ? "Contaminated" : labels[k], "l"); }
    legend.DrawClone();
    TLatex text; text.SetNDC(); text.SetTextFont(42); text.SetTextSize(0.05);
    const char* group_labels[] = {"All", "Isolated", "Non-isolated", "Gap"};
    text.DrawLatex(0.18, 0.94, group_labels[group]);
    if (slice < 0)
    {
      text.SetTextSize(0.04);
      text.DrawLatex(0.18, 0.87, ("Selected N=" + std::to_string(static_cast<long long>(selected))).c_str());
      std::string counts = "Invalid BDT=" + std::to_string(static_cast<long long>(invalid));
      if (group == 0)
      {
        double gap = 0; for (const auto& cell : cells[3]) gap += cell->summary.GetBinContent(1);
        counts += ", gap=" + std::to_string(static_cast<long long>(gap));
      }
      text.DrawLatex(0.18, 0.81, counts.c_str());
    }
    text.SetTextSize(0.045);
    if (slice < 0 ? selected == 0 : !has_scores)
    {
      text.DrawLatex(0.24, 0.34, slice < 0 || (requires_isolation && group == 2) ? "0 entries" : "No valid BDT entries");
      if (requires_isolation && group == 2) text.DrawLatex(0.24, 0.27, "Isolated required by selection");
    }
    gPad->RedrawAxis();
  }
  bool export_files(const std::string& directory, const std::string& selection, bool requires_isolation, const std::vector<std::string>& caption = {}) const
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
    return export_plots(directory, selection, requires_isolation, caption);
  }
  bool export_plots(const std::string& directory, const std::string& selection, bool requires_isolation, const std::vector<std::string>& caption = {}) const
  {
    if (!valid()) return false;
    SetsPhenixStyle();
    const std::string standard_path = directory + "/pi0_bdt_score_by_topology.pdf";
    const std::string isolation_path = directory + "/pi0_bdt_score_by_topology_isolation.pdf";
    for (int book = 0; book < 2; ++book)
    {
      const bool split = book == 1;
      const std::string& path = split ? isolation_path : standard_path;
      TCanvas canvas(split ? "pi0_bdt_isolation" : "pi0_bdt_standard", "", split ? 2100 : 1600, split ? 1600 : 1000);
      canvas.Print((path + "[").c_str());
      for (int page = 0; page < 6; ++page)
      {
        const int slice = page < 5 ? page : -1;
        canvas.Clear();
        canvas.cd();
        TLatex heading; heading.SetNDC(); heading.SetTextFont(42); heading.SetTextAlign(13);
        heading.SetTextSize(split ? 0.028 : 0.035);
        const std::vector<std::string> lines = caption.empty()
            ? std::vector<std::string>{"#it{#bf{sPHENIX}} Internal", "Candidate Cluster: after " + selection} : caption;
        const double spacing = split ? 0.04 : 0.05;
        for (std::size_t i = 0; i < lines.size(); ++i) heading.DrawLatex(0.055, 0.98 - spacing * i, lines[i].c_str());
        std::ostringstream range;
        if (slice < 0) range << "All E_{T}: 5 < E_{T} < 35 GeV";
        else range << et_edges[slice] << (slice == 0 ? " < E_{T} < " : " #leq E_{T} < ") << et_edges[slice + 1] << " GeV";
        heading.DrawLatex(0.66, 0.98, range.str().c_str());
        TPad panels("bdt_panels", "", 0, 0, 1, split ? 0.76 : 0.68);
        panels.Draw(); panels.Divide(split ? 3 : 2, split ? 2 : 1, 0.006, 0.008);
        if (split)
        {
          for (int j = 0; j < 2; ++j) for (int i = 0; i < 3; ++i)
          { panels.cd(j * 3 + i + 1); panel(i, j == 1, requires_isolation, slice); }
        }
        else
          for (int j = 0; j < 2; ++j) { panels.cd(j + 1); panel(0, j == 1, requires_isolation, slice); }
        canvas.cd();
        canvas.Print(path.c_str());
      }
      canvas.Print((path + "]").c_str());
    }
    return true;
  }
};
}
#endif
