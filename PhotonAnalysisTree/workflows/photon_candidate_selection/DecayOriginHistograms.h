#ifndef DECAY_ORIGIN_HISTOGRAMS_H
#define DECAY_ORIGIN_HISTOGRAMS_H
#include <set>
#include <tuple>
namespace decay_origin
{
constexpr int kAlgorithmVersion = 1;
constexpr double kMaxProperLengthCm = 1e-4, kMaxFlightCm = 1e-2;
enum Category : std::size_t { pi0_gg, eta_gg, decay_gamma, charged, nonprompt, remaining, mixed, unresolved, category_count };
constexpr std::array<const char*, category_count> keys = {"pi0_gg", "eta_gg", "decay_gamma", "charged", "nonprompt", "remaining", "mixed", "unresolved"};
constexpr std::array<const char*, category_count> labels = {
    "Other: #pi^{0}#rightarrow#gamma#gamma majority", "Other: #eta#rightarrow#gamma#gamma majority", "Other: hadron-decay #gamma majority",
    "Other: long-lived charged", "Other: long-lived/displaced parent", "Other: remaining particle majority", "Other: no particle majority", "Other: unresolved"};
constexpr std::array<int, category_count> colors = {kBlue - 7, kOrange - 3, kYellow + 2, kRed - 4, kViolet - 3, kGray + 1, kAzure - 4, kGray + 3};
struct Reader
{
  std::vector<unsigned char>* valid = nullptr;
  std::vector<unsigned int>* offset = nullptr;
#define DECAY_INTS(X) X(pdg) X(embedding) X(source) X(id) X(root_track_id) X(parent_pdg) X(parent_source) X(parent_id) X(status) X(pathway) X(diphoton) X(long_lived) X(daughter_count) X(daughter_pi0_count) X(daughter_gamma_count) X(daughter_charged_count)
#define DECAY_FLOATS(X) X(fraction) X(gamma0_fraction) X(gamma1_fraction) X(flight_cm) X(parent_flight_cm) X(proper_length_cm)
#define DECLARE_INT(member) std::vector<int>* member = nullptr;
#define DECLARE_FLOAT(member) std::vector<float>* member = nullptr;
  DECAY_INTS(DECLARE_INT) DECAY_FLOATS(DECLARE_FLOAT)
#undef DECLARE_INT
#undef DECLARE_FLOAT
  bool bind(TTree* tree)
  {
    bool ok = true;
#define BIND(member) ok &= bind_active(tree, "split_cluster_decay_origin_" #member, &member);
    BIND(valid) BIND(offset) DECAY_INTS(BIND) DECAY_FLOATS(BIND)
#undef BIND
    return ok;
  }
  bool consistent(std::size_t ncluster, const std::vector<unsigned int>& original_offset, const std::vector<float>& original_fraction) const
  {
    if (!valid || valid->size() != ncluster || !offset || offset->size() != ncluster + 1 || offset->front() != 0 || !pdg) return false;
    const auto size = pdg->size(); if (offset->back() != size) return false;
#define SIZE(member) if (!member || member->size() != size) return false;
    DECAY_INTS(SIZE) DECAY_FLOATS(SIZE)
#undef SIZE
    for (std::size_t c = 0; c < ncluster; ++c)
    {
      if ((*offset)[c] > (*offset)[c + 1] || (*offset)[c + 1] > size) return false;
      double before = 0, after = 0; std::set<std::tuple<int, int, int>> identities;
      for (auto i = original_offset[c]; i < original_offset[c + 1]; ++i) before += original_fraction[i];
      for (auto i = (*offset)[c]; i < (*offset)[c + 1]; ++i)
      {
        const double f = (*fraction)[i], a = (*gamma0_fraction)[i], b = (*gamma1_fraction)[i];
        if (!std::isfinite(f) || f < 0 || f > 1.00001 || !identities.emplace((*embedding)[i], (*source)[i], (*id)[i]).second) return false;
        if (!std::isfinite(a) || !std::isfinite(b) || !((a == -1 && b == -1) || (a >= 0 && b >= 0 && std::abs(a + b - f) <= 1e-5))) return false;
        if ((*status)[i] < -3 || (*status)[i] > 1 || ((*source)[i] != 1 && (*source)[i] != 2)) return false;
        for (double value : {(*flight_cm)[i], (*parent_flight_cm)[i], (*proper_length_cm)[i]}) if (!std::isfinite(value) || value < -1) return false;
        after += f;
      }
      if (std::abs(before - after) > 1e-5 * std::max(1., std::abs(before))) return false;
    }
    return true;
  }
};
#undef DECAY_INTS
#undef DECAY_FLOATS
struct Result
{
  Category category = unresolved;
  int pdg = 0, parent_pdg = 0, status = -1, pathway = 0;
  double fraction = 0, flight_cm = -1, parent_flight_cm = -1, second_gamma_fraction = -1;
};
inline Result classify(const Reader& r, std::size_t cluster)
{
  Result result; if (!(*r.valid)[cluster]) return result;
  unsigned int leading = (*r.offset)[cluster + 1]; double unknown = 0;
  for (auto i = (*r.offset)[cluster]; i < (*r.offset)[cluster + 1]; ++i)
  {
    if ((*r.status)[i] < 0) unknown += (*r.fraction)[i];
    if ((*r.fraction)[i] > result.fraction) { leading = i; result.fraction = (*r.fraction)[i]; }
  }
  if (leading == (*r.offset)[cluster + 1]) return result;
  result.pdg = (*r.pdg)[leading]; result.parent_pdg = (*r.parent_pdg)[leading]; result.status = (*r.status)[leading];
  result.pathway = (*r.pathway)[leading]; result.flight_cm = (*r.flight_cm)[leading]; result.parent_flight_cm = (*r.parent_flight_cm)[leading];
  result.second_gamma_fraction = std::min((*r.gamma0_fraction)[leading], (*r.gamma1_fraction)[leading]);
  if (result.fraction <= .5) { result.category = unknown > 0 ? unresolved : mixed; return result; }
  if ((*r.embedding)[leading] != 1 || result.status < 0) return result;
  const int pdg = std::abs(result.pdg);
  if (pdg == 11 || pdg == 13 || pdg == 211 || pdg == 321 || pdg == 2212) result.category = charged;
  else if ((*r.long_lived)[leading]) result.category = nonprompt;
  else if (pdg == 111 && (*r.diphoton)[leading]) result.category = pi0_gg;
  else if (pdg == 221 && (*r.diphoton)[leading]) result.category = eta_gg;
  else if (pdg == 22 && std::abs(result.parent_pdg) != 111 && std::abs(result.parent_pdg) != 221 && other_origin::hadron_parent(result.parent_pdg)) result.category = decay_gamma;
  else result.category = remaining;
  return result;
}
struct Histograms
{
  std::array<std::unique_ptr<TH1D>, category_count> counts, weighted, pdg, parent, status, pathway;
  std::array<std::unique_ptr<TH2D>, category_count> purity, score, flight, parent_flight, second_gamma;
  Histograms(int bins, double max)
  {
    for (std::size_t i = 0; i < category_count; ++i)
    {
      const std::string p = std::string("h_decay_origin_") + keys[i];
      counts[i] = std::make_unique<TH1D>((p + "_et_count").c_str(), "", bins, 0, max);
      weighted[i] = std::make_unique<TH1D>((p + "_et_pb").c_str(), "", bins, 0, max);
      pdg[i] = std::make_unique<TH1D>((p + "_pdg_pb").c_str(), "", 12001, -6000.5, 6000.5);
      parent[i] = std::make_unique<TH1D>((p + "_parent_pdg_pb").c_str(), "", 12001, -6000.5, 6000.5);
      status[i] = std::make_unique<TH1D>((p + "_status_pb").c_str(), "", 5, -3.5, 1.5);
      pathway[i] = std::make_unique<TH1D>((p + "_pathway_pb").c_str(), "", 4, -.5, 3.5);
      purity[i] = std::make_unique<TH2D>((p + "_fraction_vs_et_pb").c_str(), "", bins, 0, max, 50, 0, 1.00001);
      score[i] = std::make_unique<TH2D>((p + "_bdt_score_vs_et_pb").c_str(), "", bins, 0, max, 100, 0, 1.00001);
      flight[i] = std::make_unique<TH2D>((p + "_log10_flight_cm_vs_et_pb").c_str(), "", bins, 0, max, 76, -16, 3);
      parent_flight[i] = std::make_unique<TH2D>((p + "_log10_parent_flight_cm_vs_et_pb").c_str(), "", bins, 0, max, 76, -16, 3);
      second_gamma[i] = std::make_unique<TH2D>((p + "_second_gamma_fraction_vs_et_pb").c_str(), "", bins, 0, max, 50, 0, 1.00001);
    }
    for (auto* h : objects()) { h->SetDirectory(nullptr); h->Sumw2(); }
  }
  std::vector<TH1*> objects() const
  {
    std::vector<TH1*> result;
    for (std::size_t i = 0; i < category_count; ++i)
      result.insert(result.end(), {counts[i].get(), weighted[i].get(), pdg[i].get(), parent[i].get(), status[i].get(), pathway[i].get(),
                                 purity[i].get(), score[i].get(), flight[i].get(), parent_flight[i].get(), second_gamma[i].get()});
    return result;
  }
  void fill(const Result& r, double et, double weight, double bdt, bool bdt_valid)
  {
    const auto i = r.category;
    counts[i]->Fill(et); weighted[i]->Fill(et, weight); pdg[i]->Fill(r.pdg, weight); parent[i]->Fill(r.parent_pdg, weight);
    status[i]->Fill(r.status, weight); pathway[i]->Fill(r.pathway, weight); purity[i]->Fill(et, r.fraction, weight);
    if (bdt_valid && std::isfinite(bdt) && bdt >= 0 && bdt <= 1) score[i]->Fill(et, bdt, weight);
    const auto log_flight = [](double f) { return f < 0 ? -17. : f == 0 ? -15.999 : std::log10(f); };
    flight[i]->Fill(et, log_flight(r.flight_cm), weight); parent_flight[i]->Fill(et, log_flight(r.parent_flight_cm), weight);
    second_gamma[i]->Fill(et, r.second_gamma_fraction, weight); // Unknown photon shares deliberately go to underflow.
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
  bool valid(const TH1D& n, const TH1D& w) const
  {
    for (int bin = 0; bin <= n.GetNbinsX() + 1; ++bin)
    {
      double count = 0, weight = 0, error = 0;
      for (std::size_t i = 0; i < category_count; ++i)
      { count += counts[i]->GetBinContent(bin); weight += weighted[i]->GetBinContent(bin); error += std::pow(weighted[i]->GetBinError(bin), 2); }
      const auto close = [](double a, double b) { return std::abs(a - b) <= 1e-8 * std::max({1., std::abs(a), std::abs(b)}); };
      if (count != n.GetBinContent(bin) || !close(weight, w.GetBinContent(bin)) || !close(error, std::pow(w.GetBinError(bin), 2))) return false;
    }
    return true;
  }
};
using Selections = std::array<std::unique_ptr<Histograms>, kSelectionCount>;
}
#endif
