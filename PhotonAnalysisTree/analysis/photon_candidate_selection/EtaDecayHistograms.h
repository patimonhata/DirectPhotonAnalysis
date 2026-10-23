#ifndef ETA_DECAY_HISTOGRAMS_H
#define ETA_DECAY_HISTOGRAMS_H
#include "CandidateSelectionCommon.h"
namespace eta_decay
{
constexpr std::array<const char*, 3> keys = {"eta_three_pi0", "eta_charged_pions_pi0", "eta_unresolved"};
constexpr std::array<const char*, 3> labels = {"#eta#rightarrow3#pi^{0} (no #pi^{0} majority)", "#eta#rightarrow#pi^{+}#pi^{-}#pi^{0} (no #pi^{0} majority)", "#eta: unresolved daughter shares"};
constexpr std::array<int, 3> colors = {kOrange + 1, kOrange + 3, kYellow + 2};
struct Histograms
{
  std::array<std::unique_ptr<TH1D>, 3> counts, weighted;
  Histograms(int bins, double max)
  {
    for (std::size_t i = 0; i < keys.size(); ++i)
    {
      counts[i] = std::make_unique<TH1D>((std::string("h_candidate_") + keys[i] + "_et_count").c_str(), "", bins, 0, max);
      weighted[i] = std::make_unique<TH1D>((std::string("h_candidate_") + keys[i] + "_et_pb").c_str(), "", bins, 0, max);
      counts[i]->SetDirectory(nullptr); weighted[i]->SetDirectory(nullptr); counts[i]->Sumw2(); weighted[i]->Sumw2();
    }
  }
  void fill(int mode, bool unresolved, double et, double weight)
  {
    const int index = unresolved ? 2 : mode == 2 ? 0 : mode == 3 ? 1 : -1;
    if (index >= 0) { counts[index]->Fill(et); weighted[index]->Fill(et, weight); }
  }
  bool write(TDirectory& directory) const
  {
    directory.cd();
    for (std::size_t i = 0; i < keys.size(); ++i) { counts[i]->Write(); weighted[i]->Write(); }
    return true;
  }
  bool add(TDirectory& directory)
  {
    for (std::size_t i = 0; i < keys.size(); ++i)
    {
      auto* c = directory.Get<TH1D>(counts[i]->GetName()); auto* w = directory.Get<TH1D>(weighted[i]->GetName());
      if (!c || !w || c->GetNbinsX() != counts[i]->GetNbinsX() || w->GetNbinsX() != weighted[i]->GetNbinsX() ||
          c->GetXaxis()->GetXmax() != counts[i]->GetXaxis()->GetXmax() || w->GetXaxis()->GetXmax() != weighted[i]->GetXaxis()->GetXmax() ||
          !counts[i]->Add(c) || !weighted[i]->Add(w)) return false;
    }
    return true;
  }
  bool valid(const TH1D& eta_count, const TH1D& eta_weight, const TH1D& direct_count, const TH1D& direct_weight) const
  {
    for (int bin = 0; bin <= eta_count.GetNbinsX() + 1; ++bin)
    {
      double n = direct_count.GetBinContent(bin), w2 = std::pow(direct_weight.GetBinError(bin), 2);
      for (std::size_t i = 0; i < keys.size(); ++i) { n += counts[i]->GetBinContent(bin); w2 += std::pow(weighted[i]->GetBinError(bin), 2); }
      if (n > eta_count.GetBinContent(bin) || w2 > std::pow(eta_weight.GetBinError(bin), 2) + 1e-8 * std::max(1., w2)) return false;
    }
    return true;
  }
};
using Selections = std::array<std::unique_ptr<Histograms>, kSelectionCount>;
}
#endif
