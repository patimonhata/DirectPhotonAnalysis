#ifndef COMMON_COMPOSITION_HISTOGRAMS_H
#define COMMON_COMPOSITION_HISTOGRAMS_H
#include "CandidateComposition.h"
#include "../../src/CommonTruthComposition.h"
#include <TH3D.h>
#include <THnSparse.h>
namespace common_composition
{
using photon_tree::composition::Main;
using photon_tree::composition::Applicability;
constexpr std::array<const char*, 27> keys = {"all", "prompt", "pi0_separated", "pi0_merged", "pi0_single_contaminated", "pi0_missing", "pi0_other",
  "pi0_not_applicable", "pi0_unknown", "eta_separated", "eta_merged", "eta_single_contaminated", "eta_missing", "eta_other", "eta_gg_unknown",
  "eta_parent", "eta_daughters_unresolved", "decay_photon", "charged", "other_parent", "other_parent_unresolved", "residual_particle",
  "mixed", "unknown", "non_signal", "photon_origin_unknown", "eta_gg_not_applicable"};
constexpr std::array<const char*, 27> labels = {"All", "Prompt #gamma", "#pi^{0}: separated", "#pi^{0}: merged", "#pi^{0}: single contaminated",
  "#pi^{0}: missing", "#pi^{0}: other", "#pi^{0}: topology N/A", "#pi^{0}: topology unknown", "#eta#rightarrow#gamma#gamma: separated",
  "#eta#rightarrow#gamma#gamma: merged", "#eta#rightarrow#gamma#gamma: single contaminated", "#eta#rightarrow#gamma#gamma: missing",
  "#eta#rightarrow#gamma#gamma: other", "#eta#rightarrow#gamma#gamma: topology unknown", "#eta: parent majority", "#eta: daughter shares unknown",
  "Hadron decay #gamma", "Charged particle", "Other parent majority", "Other parent: daughter shares unknown", "Residual particle",
  "Mixed", "Unknown", "Non-signal", "#gamma: origin unknown", "#eta#rightarrow#gamma#gamma: topology N/A"};
inline int topology_index(int topology) { return topology == 1 ? 0 : topology == 2 ? 1 : topology == 4 ? 2 : topology == 3 ? 3 : topology == 0 ? 4 : -1; }
inline int detail(int main, int pdg, int applicability, int topology)
{
  const auto m = static_cast<Main>(main);
  if (m == Main::prompt) return 1;
  if (m == Main::pi0) return applicability == 1 ? topology_index(topology) < 0 ? -1 : 2 + topology_index(topology) : applicability == 2 ? 7 : 8;
  if (m == Main::eta_gg) return applicability == 1 ? topology_index(topology) < 0 ? -1 : 9 + topology_index(topology) : applicability == 2 ? 26 : 14;
  if (m == Main::parent) return pdg == 221 ? 15 : 19;
  if (m == Main::parent_unresolved) return pdg == 221 ? 16 : 20;
  if (m == Main::decay_photon) return 17;
  if (m == Main::charged) return 18;
  if (m == Main::particle) return 21;
  if (m == Main::mixed) return 22;
  if (m == Main::non_signal) return 24;
  if (m == Main::photon_unknown) return 25;
  return 23;
}
struct Reader
{
  std::vector<int> *main = nullptr, *selected_pdg = nullptr, *applicability = nullptr, *topology_reason = nullptr, *legacy_origin = nullptr, *matching_status = nullptr, *anchor_matching_status = nullptr;
  std::vector<unsigned char> *long_lived_ancestor = nullptr, *anchor_matching_usable = nullptr;
  std::vector<unsigned int> *quality = nullptr, *offset = nullptr, *ancestor_offset = nullptr;
  std::vector<int> *table = nullptr, *pdg = nullptr, *embedding = nullptr, *source = nullptr, *id = nullptr, *ancestor_pdg = nullptr;
  std::vector<double> *total = nullptr, *attributed = nullptr, *unresolved = nullptr, *max_photon = nullptr, *other_photons = nullptr, *charged = nullptr, *minimum_coverage = nullptr, *anchor_matching_coverage = nullptr;
  std::vector<double> *deposit = nullptr, *flight_cm = nullptr, *proper_length_cm = nullptr, *ancestor_flight_cm = nullptr, *ancestor_proper_length_cm = nullptr;
  bool bind(TTree* tree)
  {
    bool ok = true;
#define B(member) ok &= bind_active(tree, "split_cluster_common_" #member, &member)
    B(long_lived_ancestor); B(main); B(selected_pdg); B(applicability); B(topology_reason); B(legacy_origin); B(quality); B(offset); B(ancestor_offset);
    B(table); B(pdg); B(embedding); B(source); B(id); B(ancestor_pdg); B(total); B(attributed); B(unresolved); B(max_photon); B(other_photons); B(charged); B(matching_status); B(minimum_coverage);
    B(anchor_matching_status); B(anchor_matching_coverage); B(anchor_matching_usable);
    B(deposit); B(flight_cm); B(proper_length_cm); B(ancestor_flight_cm); B(ancestor_proper_length_cm);
#undef B
    return ok;
  }
  bool valid(std::size_t n) const
  {
    if (!main || !offset || !id || !ancestor_offset || !ancestor_pdg || offset->size() != n + 1 || offset->front() || offset->back() != id->size() ||
        ancestor_offset->size() != id->size() + 1 || ancestor_offset->front() || ancestor_offset->back() != ancestor_pdg->size() ||
        !std::is_sorted(offset->begin(), offset->end()) || !std::is_sorted(ancestor_offset->begin(), ancestor_offset->end())) return false;
#define N(member) if (!member || member->size() != n) return false
    N(long_lived_ancestor); N(main); N(selected_pdg); N(applicability); N(topology_reason); N(legacy_origin); N(quality); N(total); N(attributed); N(unresolved); N(max_photon); N(other_photons); N(charged); N(matching_status); N(minimum_coverage);
    N(anchor_matching_status); N(anchor_matching_coverage); N(anchor_matching_usable);
#undef N
#define N(member) if (!member || member->size() != id->size()) return false
    N(table); N(pdg); N(embedding); N(source); N(deposit); N(flight_cm); N(proper_length_cm);
#undef N
    for (std::size_t c = 0; c < n; ++c)
    {
      const int status = (*anchor_matching_status)[c];
      const double coverage = (*anchor_matching_coverage)[c];
      const auto usable = (*anchor_matching_usable)[c];
      if (status < -1 || status > 2 || !std::isfinite(coverage) || usable > 1 ||
          (status == -1 ? coverage != -1 || usable : coverage < 0) ||
          ((*applicability)[c] == 1 && (status < 0 || !usable))) return false;
      double sum = 0;
      for (auto i = (*offset)[c]; i < (*offset)[c + 1]; ++i) if ((*table)[i] == 0) sum += (*deposit)[i];
      if (!std::isfinite((*total)[c]) || (*total)[c] < 0 || !same_double(sum, (*attributed)[c]) ||
          (std::abs((*attributed)[c] - (*total)[c]) > 1e-5 * std::max(1e-30, (*total)[c]) && (*quality)[c] == 0) ||
          (*main)[c] < 0 || (*main)[c] > 11 || (*legacy_origin)[c] < 0 || (*legacy_origin)[c] > 3 || (*applicability)[c] < 0 || (*applicability)[c] > 2) return false;
    }
    return true;
  }
};
struct Histograms
{
  std::array<std::unique_ptr<TH1D>, keys.size()> counts, weighted;
  std::array<std::unique_ptr<TH1D>, 4> legacy_counts, legacy_weighted;
  TH3D migration_count, migration_weighted;
  std::unique_ptr<THnSparseD> migration_reasons_count, migration_reasons_weighted;
  std::unique_ptr<THnSparseD> ancestry_count, ancestry_weighted, matching_count, matching_weighted, anchor_matching_count, anchor_matching_weighted;
  std::unique_ptr<TH1D> unresolved_count, unresolved_weighted;
  Histograms(int bins, double max) :
    migration_count("migration_count", "All origins;Reconstructed E_{T};Legacy origin;Common detail", bins, 0, max, 4, 0, 4, keys.size()-1, 1, keys.size()),
    migration_weighted("migration_weighted", "All origins;Reconstructed E_{T};Legacy origin;Common detail", bins, 0, max, 4, 0, 4, keys.size()-1, 1, keys.size()),
    ancestry_count(make_sparse("ancestry_count", bins, max, true)), ancestry_weighted(make_sparse("ancestry_weighted", bins, max, true)),
    matching_count(make_sparse("matching_count", bins, max, false)), matching_weighted(make_sparse("matching_weighted", bins, max, false)),
    anchor_matching_count(make_sparse("anchor_matching_count", bins, max, false)),
    anchor_matching_weighted(make_sparse("anchor_matching_weighted", bins, max, false))
  {
    migration_count.SetDirectory(nullptr); migration_weighted.SetDirectory(nullptr); migration_count.Sumw2(); migration_weighted.Sumw2();
    for (std::size_t i = 0; i < keys.size(); ++i)
    { counts[i] = make(std::string(keys[i]) + "_count", bins, max); weighted[i] = make(std::string(keys[i]) + "_weighted", bins, max); }
    for (int i = 0; i < 4; ++i)
    { legacy_counts[i] = make("legacy_" + std::to_string(i) + "_count", bins, max); legacy_weighted[i] = make("legacy_" + std::to_string(i) + "_weighted", bins, max); }
    const int reason_bins[] = {bins, 4, int(keys.size()-1), 32};
    const double reason_low[] = {0,0,1,-.5}, reason_high[] = {max,4,double(keys.size()),31.5};
    migration_reasons_count = std::make_unique<THnSparseD>("migration_reasons_count", "ET,legacy,common_detail,reason_bits(1 boundary/2 unit_priority/4 long_lived/8 topology/16 missing_info)",4,reason_bins,reason_low,reason_high);
    migration_reasons_weighted = std::make_unique<THnSparseD>("migration_reasons_weighted", migration_reasons_count->GetTitle(),4,reason_bins,reason_low,reason_high);
    migration_reasons_count->Sumw2(); migration_reasons_weighted->Sumw2();
    unresolved_count = make("unresolved_deposit_fraction_sum_count", bins, max); unresolved_weighted = make("unresolved_deposit_fraction_sum_weighted", bins, max);
  }
  static std::unique_ptr<TH1D> make(const std::string& name, int bins, double max)
  { auto h = std::make_unique<TH1D>(name.c_str(), "", bins, 0, max); h->SetDirectory(nullptr); h->Sumw2(); return h; }
  static std::unique_ptr<THnSparseD> make_sparse(const char* name, int bins, double max, bool ancestry)
  {
    const int n[] = {bins, ancestry ? 20001 : 3, ancestry ? 3 : 101, ancestry ? 256 : 8};
    const double low[] = {0, ancestry ? -10000.5 : -0.5, ancestry ? -0.5 : 0, -0.5};
    const double high[] = {max, ancestry ? 10000.5 : 2.5, ancestry ? 2.5 : 1.01, ancestry ? 255.5 : 7.5};
    auto h = std::make_unique<THnSparseD>(name, ancestry ? "ET,PDG,relation(0 unit/1 parent/2 ancestor),quality" : "ET,status,coverage,topology_reason", 4, n, low, high); h->Sumw2(); return h;
  }
  void fill(const Reader& r, std::size_t c, int category, double et, double weight, int status = -1, double coverage = -1)
  {
    counts[0]->Fill(et); weighted[0]->Fill(et, weight); counts[category]->Fill(et); weighted[category]->Fill(et, weight);
    const int old = (*r.legacy_origin)[c]; legacy_counts[old]->Fill(et); legacy_weighted[old]->Fill(et, weight);
    migration_count.Fill(et, old + .5, category + .5); migration_weighted.Fill(et, old + .5, category + .5, weight);
    const auto main = static_cast<Main>((*r.main)[c]);
    const int reason = (old != photon_tree::composition::overview(main, (*r.selected_pdg)[c]) ? 1 : 0) |
        (main == Main::charged || main == Main::decay_photon ? 2 : 0) | ((*r.long_lived_ancestor)[c] ? 4 : 0) |
        ((*r.applicability)[c] != 1 && (main == Main::pi0 || main == Main::eta_gg) ? 8 : 0) |
        ((*r.quality)[c] || main == Main::parent_unresolved || main == Main::unknown ? 16 : 0);
    const double migration_reason[] = {et, old+.5, category+.5, double(reason)};
    migration_reasons_count->Fill(migration_reason); migration_reasons_weighted->Fill(migration_reason, weight);
    const double fraction = (*r.total)[c] > 0 ? (*r.unresolved)[c] / (*r.total)[c] : 0;
    unresolved_count->Fill(et, fraction); unresolved_weighted->Fill(et, fraction * weight);
    if (status < 0) { status = (*r.matching_status)[c]; coverage = (*r.minimum_coverage)[c]; }
    const double match[] = {et, double(status), coverage, double((*r.topology_reason)[c])}; matching_count->Fill(match); matching_weighted->Fill(match, weight);
    if ((*r.anchor_matching_status)[c] >= 0)
    {
      const double local[] = {et, double((*r.anchor_matching_status)[c]), (*r.anchor_matching_coverage)[c], double((*r.topology_reason)[c])};
      anchor_matching_count->Fill(local); anchor_matching_weighted->Fill(local, weight);
    }
    std::set<int> ancestors;
    for (auto i = (*r.offset)[c]; i < (*r.offset)[c + 1]; ++i)
    {
      if ((*r.table)[i] == 1) continue;
      const double value[] = {et, double((*r.pdg)[i]), (*r.table)[i] == 0 ? 0. : 1., double((*r.quality)[c])};
      ancestry_count->Fill(value); ancestry_weighted->Fill(value, weight);
      for (auto j = (*r.ancestor_offset)[i]; j < (*r.ancestor_offset)[i + 1]; ++j) ancestors.insert((*r.ancestor_pdg)[j]);
    }
    for (int pdg : ancestors) { const double v[] = {et, double(pdg), 2., double((*r.quality)[c])}; ancestry_count->Fill(v); ancestry_weighted->Fill(v, weight); }
  }
  bool valid() const
  {
    for (int b = 0; b <= counts[0]->GetNbinsX() + 1; ++b)
    {
      for (bool w : {false, true})
      {
        const auto& hs = w ? weighted : counts; double sum = 0, sum2 = 0, old = 0, old2 = 0;
        for (std::size_t k = 1; k < keys.size(); ++k)
        {
          double transition = 0, transition2 = 0;
          for (int origin = 1; origin <= 4; ++origin)
          {
            const auto& m = w ? migration_weighted : migration_count; transition += m.GetBinContent(b, origin, k);
            transition2 += std::pow(m.GetBinError(m.GetBin(b, origin, k)), 2);
          }
          if (!same_weighted(transition, hs[k]->GetBinContent(b)) || !same_weighted(transition2, std::pow(hs[k]->GetBinError(b), 2))) return false;
          sum += hs[k]->GetBinContent(b); sum2 += std::pow(hs[k]->GetBinError(b), 2);
        }
        for (int k = 0; k < 4; ++k)
        {
          const auto& h = w ? legacy_weighted[k] : legacy_counts[k]; double row = 0, row2 = 0;
          const auto& m = w ? migration_weighted : migration_count;
          for (std::size_t d = 1; d < keys.size(); ++d) { row += m.GetBinContent(b,k+1,d); row2 += std::pow(m.GetBinError(m.GetBin(b,k+1,d)),2); }
          if (!same_weighted(row,h->GetBinContent(b)) || !same_weighted(row2,std::pow(h->GetBinError(b),2))) return false;
          old += h->GetBinContent(b); old2 += std::pow(h->GetBinError(b), 2);
        }
        if (!same_weighted(sum, hs[0]->GetBinContent(b)) || !same_weighted(sum2, std::pow(hs[0]->GetBinError(b), 2)) ||
            !same_weighted(old, sum) || !same_weighted(old2, sum2)) return false;
      }
    }
    return true;
  }
  static bool same_weighted(double a, double b) { return candidate_composition::same_weighted_integral(a, b); }
  bool write(TDirectory& d) const
  {
    d.cd(); for (std::size_t i = 0; i < keys.size(); ++i) { counts[i]->Write(); weighted[i]->Write(); }
    for (int i = 0; i < 4; ++i) { legacy_counts[i]->Write(); legacy_weighted[i]->Write(); }
    migration_count.Write(); migration_weighted.Write(); migration_reasons_count->Write(); migration_reasons_weighted->Write(); ancestry_count->Write(); ancestry_weighted->Write(); matching_count->Write(); matching_weighted->Write();
    anchor_matching_count->Write(); anchor_matching_weighted->Write();
    unresolved_count->Write(); unresolved_weighted->Write(); return true;
  }
  bool add(TDirectory& d)
  {
    for (std::size_t i = 0; i < keys.size(); ++i) for (auto* h : {counts[i].get(), weighted[i].get()})
    { const auto* v = d.Get<TH1D>(h->GetName()); if (!v || !h->Add(v)) return false; }
    for (int i = 0; i < 4; ++i) for (auto* h : {legacy_counts[i].get(), legacy_weighted[i].get()})
    { const auto* v = d.Get<TH1D>(h->GetName()); if (!v || !h->Add(v)) return false; }
    for (auto* h : {&migration_count, &migration_weighted}) { const auto* v = d.Get<TH3D>(h->GetName()); if (!v || !h->Add(v)) return false; }
    for (auto* h : {ancestry_count.get(), ancestry_weighted.get(), matching_count.get(), matching_weighted.get(), anchor_matching_count.get(), anchor_matching_weighted.get(), migration_reasons_count.get(), migration_reasons_weighted.get()})
    { const auto* v = d.Get<THnSparseD>(h->GetName()); if (!v) return false; h->Add(v); }
    for (auto* h : {unresolved_count.get(), unresolved_weighted.get()}) { const auto* v = d.Get<TH1D>(h->GetName()); if (!v || !h->Add(v)) return false; }
    return true;
  }
};
inline bool valid_populations(const Histograms& common, const candidate_composition::Histograms& composition, const Spectra& pi0, const Spectra& eta)
{
  using namespace candidate_composition;
  const auto matches = [&](const std::vector<int>& indices, const TH1D& target, bool weighted) {
    for (int b = 0; b <= target.GetNbinsX()+1; ++b)
    {
      double sum = 0, sum2 = 0;
      for (int i : indices) { const auto& h = weighted ? common.weighted[i] : common.counts[i]; sum += h->GetBinContent(b); sum2 += std::pow(h->GetBinError(b),2); }
      if (!Histograms::same_weighted(sum, target.GetBinContent(b)) || !Histograms::same_weighted(sum2,std::pow(target.GetBinError(b),2))) return false;
    }
    return true;
  };
  for (bool w : {false,true})
  {
    const auto& c = w ? composition.weighted : composition.counts;
    if (!matches({0}, *c[denominator], w) || !matches({1}, *c[prompt], w) || !matches({7}, *c[pi0_not_applicable], w) ||
        !matches({8}, *c[pi0_unknown], w) || !matches({9,10,11,12,13,14,15,16,26}, *c[candidate_composition::eta], w) ||
        !matches({2,3,4,5,6}, *(w ? pi0.weighted_pb[1] : pi0.counts[1]), w) ||
        !matches({9,10,11,12,13}, *(w ? eta.weighted_pb[1] : eta.counts[1]), w)) return false;
  }
  return true;
}
using Selections = std::array<std::unique_ptr<Histograms>, kSelectionCount>;
}
#endif
