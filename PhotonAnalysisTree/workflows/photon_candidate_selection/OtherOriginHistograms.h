#ifndef OTHER_ORIGIN_HISTOGRAMS_H
#define OTHER_ORIGIN_HISTOGRAMS_H
#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include <TH2D.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <vector>
namespace other_origin
{
constexpr int kAlgorithmVersion = 1;
enum Category : std::size_t { charged_hadron, electron, muon, decay_photon, other_particle, mixed, unresolved, category_count };
constexpr std::array<const char*, category_count> keys = {"charged_hadron", "electron", "muon", "decay_photon", "other_particle", "mixed", "unresolved"};
constexpr std::array<const char*, category_count> labels = {
    "Other: #pi^{#pm}, K^{#pm}, p/#bar{p}", "Other: e^{#pm}", "Other: #mu^{#pm}", "Other: hadron-decay #gamma (not #pi^{0}/#eta)",
    "Other: remaining particle majority", "Other: no single-particle majority", "Other: unresolved truth"};
constexpr std::array<int, category_count> colors = {kRed - 4, kPink + 1, kViolet - 3, kYellow + 2, kGray + 1, kAzure - 4, kGray + 3};
constexpr std::array<const char*, 3> summary_keys = {"stable_charged", "decay_photon", "residual"};
constexpr std::array<const char*, 3> summary_labels = {"Other: long-lived charged", "Other: hadron-decay #gamma (not #pi^{0}/#eta)", "Other: remaining"};
constexpr std::array<int, 3> summary_colors = {kRed - 4, kYellow + 2, kGray + 1};
inline int summary_category(std::size_t category) { return category <= muon ? 0 : category == decay_photon ? 1 : 2; }
inline bool hadron_parent(int pdg)
{
  const auto* particle = TDatabasePDG::Instance()->GetParticle(pdg);
  if (!particle) return false;
  const std::string type = particle->ParticleClass();
  return type.find("Meson") != std::string::npos || type.find("Baryon") != std::string::npos;
}
struct Result
{
  Category category = unresolved;
  double main_fraction = 0, other_gamma_fraction = 0;
  int pdg = 0, parent = 0;
};
inline Result classify(bool truth_valid, unsigned int begin, unsigned int end, const std::vector<int>& pdg,
                       const std::vector<int>& embedding, const std::vector<float>& fraction, const std::vector<int>& parent,
                       const std::vector<int>& photon_category, bool origin_overlap = false)
{
  Result result;
  if (!truth_valid || begin >= end || end > pdg.size() || pdg.size() != fraction.size() || pdg.size() != embedding.size() ||
      pdg.size() != parent.size() || pdg.size() != photon_category.size() || origin_overlap) return result;
  unsigned int leading = end;
  for (unsigned int i = begin; i < end; ++i)
  {
    if (!std::isfinite(fraction[i]) || fraction[i] < 0 || fraction[i] > 1.0 + 1e-5) return result;
    if (fraction[i] > result.main_fraction) { result.main_fraction = fraction[i]; leading = i; }
  }
  if (leading == end) return result;
  result.pdg = pdg[leading]; result.parent = parent[leading];
  for (unsigned int i = begin; i < end; ++i) if (i != leading && pdg[i] == 22) result.other_gamma_fraction += fraction[i];
  if (result.main_fraction <= .5) { result.category = mixed; return result; }
  if (embedding[leading] != 1) return result;
  const int code = std::abs(result.pdg);
  if (code == 211 || code == 321 || code == 2212) result.category = charged_hadron;
  else if (code == 11) result.category = electron;
  else if (code == 13) result.category = muon;
  else if (code == 22 && photon_category[leading] == 3 && std::abs(result.parent) != 111 && std::abs(result.parent) != 221 && hadron_parent(result.parent))
    result.category = decay_photon;
  else result.category = other_particle;
  return result;
}
struct Histograms
{
  std::array<std::unique_ptr<TH1D>, category_count> counts, weighted, pdg;
  std::array<std::unique_ptr<TH2D>, category_count> purity, other_gamma, score;
  std::unique_ptr<TH1D> decay_parent;
  Histograms(int bins, double max)
  {
    for (std::size_t i = 0; i < category_count; ++i)
    {
      const std::string prefix = std::string("h_other_") + keys[i];
      counts[i] = std::make_unique<TH1D>((prefix + "_et_count").c_str(), "", bins, 0, max);
      weighted[i] = std::make_unique<TH1D>((prefix + "_et_pb").c_str(), "", bins, 0, max);
      pdg[i] = std::make_unique<TH1D>((prefix + "_leading_pdg_pb").c_str(), "", 12001, -6000.5, 6000.5);
      purity[i] = std::make_unique<TH2D>((prefix + "_main_fraction_vs_et_pb").c_str(), "", bins, 0, max, 50, 0, 1.00001);
      other_gamma[i] = std::make_unique<TH2D>((prefix + "_other_primary_gamma_fraction_vs_et_pb").c_str(), "", bins, 0, max, 50, 0, 1.00001);
      score[i] = std::make_unique<TH2D>((prefix + "_bdt_score_vs_et_pb").c_str(), "", bins, 0, max, 100, 0, 1.00001);
      const std::array<TH1*, 6> histograms = {counts[i].get(), weighted[i].get(), pdg[i].get(), purity[i].get(), other_gamma[i].get(), score[i].get()};
      for (auto* h : histograms)
      { h->SetDirectory(nullptr); h->Sumw2(); }
    }
    decay_parent = std::make_unique<TH1D>("h_other_decay_photon_parent_pdg_pb", "", 12001, -6000.5, 6000.5);
    decay_parent->SetDirectory(nullptr); decay_parent->Sumw2();
  }
  void fill(const Result& r, double et, double weight, double bdt, bool bdt_valid)
  {
    const auto i = r.category;
    counts[i]->Fill(et); weighted[i]->Fill(et, weight); pdg[i]->Fill(r.pdg, weight);
    purity[i]->Fill(et, r.main_fraction, weight); other_gamma[i]->Fill(et, r.other_gamma_fraction, weight);
    if (bdt_valid && std::isfinite(bdt) && bdt >= 0 && bdt <= 1) score[i]->Fill(et, bdt, weight);
    if (i == decay_photon) decay_parent->Fill(r.parent, weight);
  }
  std::vector<TH1*> objects() const
  {
    std::vector<TH1*> result;
    for (std::size_t i = 0; i < category_count; ++i)
    { result.insert(result.end(), {counts[i].get(), weighted[i].get(), pdg[i].get(), purity[i].get(), other_gamma[i].get(), score[i].get()}); }
    result.push_back(decay_parent.get()); return result;
  }
  bool write(TDirectory& directory) const { directory.cd(); for (auto* h : objects()) if (h->Write() <= 0) return false; return true; }
  bool add(TDirectory& directory)
  {
    for (auto* h : objects())
    {
      auto* source = directory.Get<TH1>(h->GetName());
      if (!source || source->IsA() != h->IsA() || source->GetNbinsX() != h->GetNbinsX() || source->GetNbinsY() != h->GetNbinsY() ||
          source->GetXaxis()->GetXmin() != h->GetXaxis()->GetXmin() || source->GetXaxis()->GetXmax() != h->GetXaxis()->GetXmax() ||
          source->GetYaxis()->GetXmin() != h->GetYaxis()->GetXmin() || source->GetYaxis()->GetXmax() != h->GetYaxis()->GetXmax() || !h->Add(source)) return false;
    }
    return true;
  }
  bool valid(const TH1D& total_n, const TH1D& total_w) const
  {
    for (int bin = 0; bin <= total_n.GetNbinsX() + 1; ++bin)
    {
      double n = 0, w = 0, w2 = 0;
      for (std::size_t i = 0; i < category_count; ++i)
      { n += counts[i]->GetBinContent(bin); w += weighted[i]->GetBinContent(bin); w2 += std::pow(weighted[i]->GetBinError(bin), 2); }
      const auto close = [](double a, double b) { return std::abs(a - b) <= 1e-8 * std::max({1., std::abs(a), std::abs(b)}); };
      if (n != total_n.GetBinContent(bin) || !close(w, total_w.GetBinContent(bin)) || !close(w2, std::pow(total_w.GetBinError(bin), 2))) return false;
    }
    return true;
  }
};
using Selections = std::array<std::unique_ptr<Histograms>, kSelectionCount>;
}
#endif
