#ifndef TAGGING_WORKING_POINT_SCAN_H
#define TAGGING_WORKING_POINT_SCAN_H
#include "../tagging_mass_diagnostic/TaggingMassDiagnostic.h"
#include <TNamed.h>
#include <fstream>

namespace working_point_scan
{
using tagging_mass::Metadata;
constexpr int npoints = 6, nmasks = 1 << npoints, reference_bit = 4;
constexpr std::array<const char*, npoints> keys = {"e0p15_m0p08_0p18", "e0p15_m0p08_0p20", "e0p12_m0p08_0p18", "e0p12_m0p08_0p20", "e0p5_m0p10_0p20", "e0p15_m0p10_0p20"};
constexpr std::array<double, npoints> cuts = {.15, .15, .12, .12, .5, .15};
constexpr std::array<double, npoints> lower = {.08, .08, .08, .08, .10, .10};
constexpr std::array<double, npoints> upper = {.18, .20, .18, .20, .20, .20};
constexpr const char* definition = "working_point_scan_v2; pi0: E>.15/.15/.12/.12/.5/.15,m>.08/.08/.08/.08/.10/.10,m<.18/.20/.18/.20/.20/.20; reference=4; eta: E>.5,.45<m<.65; saved Region A; exclusive majority composition; pi0 counts are clusters";

struct Histograms
{
  TH2D prompt_count{"prompt_survival_mask_count", "", 40, 0, 40, nmasks, -.5, nmasks - .5};
  TH2D prompt_pb{"prompt_survival_mask_pb", "", 40, 0, 40, nmasks, -.5, nmasks - .5};
  TH2D background_count{"background_survival_mask_count", "", 40, 0, 40, nmasks, -.5, nmasks - .5};
  TH2D background_pb{"background_survival_mask_pb", "", 40, 0, 40, nmasks, -.5, nmasks - .5};
  TH2D pi0_count{"pi0_survival_mask_count", "", 40, 0, 40, nmasks, -.5, nmasks - .5};
  TH2D pi0_pb{"pi0_survival_mask_pb", "", 40, 0, 40, nmasks, -.5, nmasks - .5};
  std::array<TH2D*, 6> objects() { return {&prompt_count, &prompt_pb, &background_count, &background_pb, &pi0_count, &pi0_pb}; }
  Histograms() { for (auto* h : objects()) { h->SetDirectory(nullptr); h->Sumw2(); } }
  void fill(bool prompt, bool pi0, double et, unsigned mask, double weight)
  { (prompt ? prompt_count : background_count).Fill(et, mask); (prompt ? prompt_pb : background_pb).Fill(et, mask, weight); if (pi0) { pi0_count.Fill(et, mask); pi0_pb.Fill(et, mask, weight); } }
  void write() { for (auto* h : objects()) h->Write(); TNamed("analysis_definition", definition).Write(); }
  bool add(TFile& f)
  {
    auto* tag = f.Get<TNamed>("analysis_definition");
    if (!tag || std::string(tag->GetTitle()) != definition) return false;
    for (auto* h : objects()) if (!tagging_mass::Histograms::add_checked(h, f.Get<TH1>(h->GetName()))) return false;
    return true;
  }
};

struct Yield { double value = 0, variance = 0; };
inline Yield yield(const TH2D& h, double lo, double hi, int bit = -1, int removed_from = -1)
{
  Yield out;
  for (int x = h.GetXaxis()->FindFixBin(lo); x <= h.GetXaxis()->FindFixBin(std::nextafter(hi, lo)); ++x)
    for (int mask = 0; mask < nmasks; ++mask)
    {
      const bool include = removed_from >= 0 ? ((mask & (1 << removed_from)) && !(mask & (1 << bit))) : bit < 0 || (mask & (1 << bit));
      if (include) { out.value += h.GetBinContent(x, mask + 1); out.variance += std::pow(h.GetBinError(x, mask + 1), 2); }
    }
  return out;
}
inline double fraction_error(Yield numerator, Yield denominator)
{
  if (denominator.value == 0) return 0;
  const double f = numerator.value / denominator.value;
  return std::sqrt(std::max(0., (1 - 2 * f) * numerator.variance + f * f * denominator.variance)) / std::abs(denominator.value);
}
inline void summarize(Histograms& h, const std::string& path)
{
  std::ofstream out(path);
  out << "et_low\tet_high\tworking_point\tprompt_before_pb\tbackground_before_pb\tprompt_after_pb\tbackground_after_pb\tpurity\tpurity_error\tprompt_survival\tprompt_survival_error\tbackground_survival\textra_prompt_removed_pb\textra_background_removed_pb\textra_prompt_removed_count\textra_background_removed_count\tprompt_recovered_pb\tprompt_net_loss_fraction\tprompt_gross_loss_fraction\tprompt_recovery_fraction\tpi0_after_count\tpi0_after_pb\tpi0_fraction\tselected_after_count\tselected_after_pb\n";
  out << std::setprecision(12);
  for (auto range : {std::make_pair(5.,6.), {6.,8.}, {8.,10.}, {10.,18.}, {18.,30.}, {5.,30.}})
    for (int w = 0; w < npoints; ++w)
    {
      const auto p0 = yield(h.prompt_pb, range.first, range.second), b0 = yield(h.background_pb, range.first, range.second);
      const auto p = yield(h.prompt_pb, range.first, range.second, w), b = yield(h.background_pb, range.first, range.second, w);
      const Yield all{p.value + b.value, p.variance + b.variance};
      out << range.first << '\t' << range.second << '\t' << keys[w] << '\t' << p0.value << '\t' << b0.value << '\t' << p.value << '\t' << b.value << '\t'
          << (all.value ? p.value / all.value : 0) << '\t' << fraction_error(p, all) << '\t' << (p0.value ? p.value / p0.value : 0) << '\t'
          << fraction_error(p, p0) << '\t' << (b0.value ? b.value / b0.value : 0) << '\t'
          << yield(h.prompt_pb, range.first, range.second, w, reference_bit).value << '\t' << yield(h.background_pb, range.first, range.second, w, reference_bit).value << '\t'
          << yield(h.prompt_count, range.first, range.second, w, reference_bit).value << '\t' << yield(h.background_count, range.first, range.second, w, reference_bit).value << '\t';
      const auto ref = yield(h.prompt_pb, range.first, range.second, reference_bit);
      const auto lost = yield(h.prompt_pb, range.first, range.second, w, reference_bit);
      const auto recovered = yield(h.prompt_pb, range.first, range.second, reference_bit, w);
      const auto pi0 = yield(h.pi0_pb, range.first, range.second, w);
      out << recovered.value << '\t' << (ref.value ? 1 - p.value / ref.value : 0) << '\t'
          << (ref.value ? lost.value / ref.value : 0) << '\t' << (ref.value ? recovered.value / ref.value : 0) << '\t'
          << yield(h.pi0_count, range.first, range.second, w).value << '\t' << pi0.value << '\t' << (all.value ? pi0.value / all.value : 0) << '\t'
          << yield(h.prompt_count, range.first, range.second, w).value + yield(h.background_count, range.first, range.second, w).value << '\t' << all.value << '\n';
    }
}
inline void plot(Histograms& h, const std::string& output)
{
  SetsPhenixStyle();
  for (bool purity : {true, false})
  {
    TCanvas canvas(purity ? "scan_purity" : "scan_prompt", "", 1000, 800);
    std::array<std::unique_ptr<TH1D>, npoints> curves;
    TLegend legend(.46, .16, .91, .44); legend.SetFillStyle(0); legend.SetTextSize(.025);
    for (int w = 0; w < npoints; ++w)
    {
      curves[w] = std::make_unique<TH1D>((std::string("scan_") + keys[w]).c_str(), "", 40, 0, 40);
      curves[w]->SetDirectory(nullptr);
      for (int bin = 1; bin <= 40; ++bin)
      {
        const auto p = yield(h.prompt_pb, bin-1, bin, w);
        const auto b = yield(h.background_pb, bin-1, bin, w);
        const auto den = purity ? Yield{p.value+b.value,p.variance+b.variance} : yield(h.prompt_pb, bin-1, bin);
        curves[w]->SetBinContent(bin, den.value ? p.value/den.value : 0); curves[w]->SetBinError(bin, fraction_error(p,den));
      }
      curves[w]->SetLineColor(tagging_mass::kColors[w]); curves[w]->SetLineWidth(2); curves[w]->SetStats(false);
      curves[w]->SetMinimum(0); curves[w]->SetMaximum(1.05); curves[w]->GetXaxis()->SetTitle("Candidate cluster E_{T} [GeV]");
      curves[w]->GetXaxis()->SetRangeUser(5, 18); curves[w]->GetYaxis()->SetTitle(purity ? "Prompt purity after veto" : "Prompt survival relative to Region A");
      curves[w]->Draw(w ? "HIST SAME" : "HIST");
      std::ostringstream label; label << "E > " << cuts[w] << "; " << lower[w] << " < m < " << upper[w] << " GeV";
      legend.AddEntry(curves[w].get(), label.str().c_str(), "l");
    }
    legend.Draw(); canvas.SaveAs((output + (purity ? "/purity.pdf" : "/prompt_survival.pdf")).c_str());
  }
}
}
#endif
