#include "DecayOriginDiagnostics.h"
#include "Pi0ClusterTruthMatcher.h"
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>
#include <g4main/PHG4MCProcessDefs.h>
#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>
#include <HepMC/GenEvent.h>
#include <HepMC/GenParticle.h>
#include <HepMC/GenVertex.h>
#include <HepMC/Units.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include <TTree.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <set>
#include <tuple>
namespace photon_tree
{
namespace
{
using Key = std::tuple<int, int, int>;
struct Node
{
  const PHG4Particle* g4 = nullptr;
  const HepMC::GenParticle* hep = nullptr;
  int embedding = 0;
  Key key() const { return {embedding, hep ? 1 : g4 ? 2 : 0, hep ? hep->barcode() : g4 ? g4->get_track_id() : 0}; }
  int pdg() const { return hep ? hep->pdg_id() : g4 ? g4->get_pid() : 0; }
  explicit operator bool() const { return hep || g4; }
};
bool hadron(int pdg)
{
  const auto* p = TDatabasePDG::Instance()->GetParticle(pdg);
  const std::string type = p ? p->ParticleClass() : "";
  return type.find("Meson") != std::string::npos || type.find("Baryon") != std::string::npos;
}
double proper_length(int pdg)
{
  const auto* p = TDatabasePDG::Instance()->GetParticle(pdg);
  if (!p) return -1;
  if (p->Lifetime() > 0 && std::isfinite(p->Lifetime())) return p->Lifetime() * 2.99792458e10;
  if (p->Width() > 0 && std::isfinite(p->Width())) return 1.973269804e-14 / p->Width();
  return p->Stable() ? 1e30 : -1;
}
bool close(double a, double b) { return std::abs(a - b) <= 1e-5 * std::max({1e-30, std::abs(a), std::abs(b)}); }
}
struct DecayOriginDiagnostics::Impl
{
  PHG4TruthInfoContainer* truth;
  const PHHepMCGenEventMap* events;
  std::map<int, std::vector<const PHG4Particle*>> children;
  std::map<std::pair<int, int>, const PHG4Particle*> tracks;
  std::set<std::pair<int, int>> ambiguous_tracks;
  Impl(PHG4TruthInfoContainer* t, const PHHepMCGenEventMap* e) : truth(t), events(e)
  {
    if (!truth) return;
    const auto range = truth->GetParticleRange();
    for (auto it = range.first; it != range.second; ++it)
    {
      const auto* p = it->second; if (!p) continue;
      children[p->get_parent_id()].push_back(p);
      if (p->get_parent_id() != 0 || p->get_barcode() == 0) continue;
      const auto key = std::make_pair(truth->isEmbeded(p->get_track_id()), p->get_barcode());
      if (!tracks.emplace(key, p).second) ambiguous_tracks.insert(key);
    }
  }
  Node from_g4(const PHG4Particle* p, int embedding) const
  {
    Node n{p, nullptr, embedding};
    if (!p || p->get_parent_id() != 0) return n;
    const auto* sub = events ? events->get(embedding) : nullptr;
    const auto* event = sub ? sub->getEvent() : nullptr;
    const auto key = std::make_pair(embedding, p->get_barcode());
    const auto* hep = event && !ambiguous_tracks.count(key) ? event->barcode_to_particle(p->get_barcode()) : nullptr;
    if (hep && hep->pdg_id() == p->get_pid()) n.hep = hep;
    return n;
  }
  Node from_hep(const HepMC::GenParticle* p, int embedding) const
  {
    Node n{nullptr, p, embedding}; if (!p) return n;
    const auto key = std::make_pair(embedding, p->barcode());
    const auto it = tracks.find(key);
    if (it != tracks.end() && !ambiguous_tracks.count(key) && it->second->get_pid() == p->pdg_id()) n.g4 = it->second;
    return n;
  }
  Node parent(Node n) const
  {
    if (n.g4 && n.g4->get_parent_id() != 0)
    {
      const auto* v = truth->GetVtx(n.g4->get_vtx_id());
      if (!v || v->get_process() != PHG4MCProcess::kPDecay) return {};
      return from_g4(truth->GetParticle(n.g4->get_parent_id()), n.embedding);
    }
    const auto* v = n.hep ? n.hep->production_vertex() : nullptr;
    if (!v || v->particles_in_size() != 1) return {};
    if ((*v->particles_in_const_begin())->status() == 4) return {}; // Beam particles are not decay parents.
    auto p = from_hep(*v->particles_in_const_begin(), n.embedding);
    return hadron(p.pdg()) || (n.pdg() == 22 && p.pdg() == 22) ? p : Node{};
  }
  Node physical(Node n) const
  {
    std::set<Key> seen;
    while (n && seen.insert(n.key()).second)
    {
      const auto p = parent(n);
      if (!p || p.pdg() != n.pdg()) break;
      // Only photon bookkeeping copies are collapsed, never consecutive same-PDG hadrons.
      if (n.pdg() != 22) break;
      n = p;
    }
    return n;
  }
  std::vector<Node> daughters(Node n, int& state) const
  {
    state = 0;
    if (n.g4)
    {
      const auto it = children.find(n.g4->get_track_id());
      if (it != children.end() && !it->second.empty())
      {
        std::vector<Node> result; const PHG4VtxPoint* first = nullptr;
        for (const auto* child : it->second)
        {
          const auto* v = truth->GetVtx(child->get_vtx_id());
          if (!v || v->get_process() != PHG4MCProcess::kPDecay) { state = -2; return {}; }
          if (first && (v->get_x() != first->get_x() || v->get_y() != first->get_y() || v->get_z() != first->get_z() || v->get_t() != first->get_t()))
          { state = -1; return {}; }
          first = v; result.push_back(from_g4(child, n.embedding));
        }
        state = 2; return result; // Choose one representation; never append HepMC and G4 copies of a decay.
      }
    }
    const auto* v = n.hep ? n.hep->end_vertex() : nullptr;
    if (!v) return {};
    if (v->particles_in_size() != 1 || *v->particles_in_const_begin() != n.hep) { state = -1; return {}; }
    std::vector<Node> result;
    for (auto it = v->particles_out_const_begin(); it != v->particles_out_const_end(); ++it) result.push_back(from_hep(*it, n.embedding));
    if (!result.empty()) state = 1;
    return result;
  }
  double flight(Node n) const
  {
    if (n.g4)
    {
      const auto found = children.find(n.g4->get_track_id());
      const auto* begin = truth->GetVtx(n.g4->get_vtx_id());
      const auto* end = found != children.end() && !found->second.empty() ? truth->GetVtx(found->second.front()->get_vtx_id()) : nullptr;
      if (begin && end && end->get_process() == PHG4MCProcess::kPDecay)
      {
        const double distance = std::hypot(std::hypot(end->get_x() - begin->get_x(), end->get_y() - begin->get_y()), end->get_z() - begin->get_z());
        return std::isfinite(distance) ? distance : -1;
      }
    }
    const auto* begin = n.hep ? n.hep->production_vertex() : nullptr;
    const auto* end = n.hep ? n.hep->end_vertex() : nullptr;
    const auto* event = n.hep ? n.hep->parent_event() : nullptr;
    if (!begin || !end || !event) return -1;
    const auto& a = begin->position(); const auto& b = end->position();
    const double distance = std::hypot(std::hypot(b.x() - a.x(), b.y() - a.y()), b.z() - a.z()) * HepMC::Units::conversion_factor(event->length_unit(), HepMC::Units::CM);
    return std::isfinite(distance) ? distance : -1;
  }
  // 1: eligible prompt decay, 0: retain a long-lived/displaced parent, -3: unknown policy inputs.
  int prompt(Node n) const
  {
    const double ctau = proper_length(n.pdg());
    if (!(ctau >= 0)) return -3;
    if (ctau >= kMaxProperLengthCm) return 0;
    const double distance = flight(n);
    if (!(distance >= 0)) return -3;
    return distance <= kMaxFlightCm ? 1 : 0;
  }
  bool diphoton(Node n) const
  {
    if (n.pdg() != 111 && n.pdg() != 221) return false;
    int state = 0; const auto ds = daughters(n, state);
    return state > 0 && ds.size() == 2 && ds[0].pdg() == 22 && ds[1].pdg() == 22;
  }
  Node anchor(Node n) const
  {
    n = physical(n); Node selected = n; bool stopped = false;
    std::set<Key> seen;
    for (auto p = parent(n); p && seen.insert(p.key()).second; p = parent(p))
    {
      if (!hadron(p.pdg())) break;
      const int rule = prompt(p);
      // A non-prompt/unknown ancestor is a barrier in either simulation representation.
      if (rule != 1) { selected = p; stopped = true; }
      else if (!stopped && (p.pdg() == 111 || (p.pdg() == 221 && diphoton(p)))) { selected = p; stopped = true; }
    }
    return selected;
  }
  int ancestry_path(Node input, Node atom) const
  {
    int pathway = 0; std::set<Key> seen;
    for (auto n = input; n && n.key() != atom.key() && seen.insert(n.key()).second; n = parent(n))
      pathway |= n.g4 && n.g4->get_parent_id() != 0 ? 2 : 1;
    return pathway;
  }
  DecayOriginContribution describe(Node n, int root, float fraction, int status, int pathway) const
  {
    DecayOriginContribution r; const auto [embedding, source, id] = n.key();
    r.pdg = n.pdg(); r.embedding = embedding; r.source = source; r.id = id; r.root_track_id = root;
    r.fraction = fraction; r.status = status; r.pathway = pathway; r.diphoton = diphoton(n);
    if (!TDatabasePDG::Instance()->GetParticle(n.pdg()) && r.status >= 0) r.status = -3;
    int child_state = 0; const auto direct = daughters(n, child_state);
    r.daughter_count = direct.size();
    for (const auto& child : direct)
    {
      r.daughter_pi0_count += child.pdg() == 111; r.daughter_gamma_count += child.pdg() == 22;
      const auto* info = TDatabasePDG::Instance()->GetParticle(child.pdg());
      r.daughter_charged_count += info && info->Charge() != 0;
    }
    r.flight_cm = flight(n); r.proper_length_cm = proper_length(n.pdg());
    const auto p = physical(parent(n));
    if (p) { r.parent_pdg = p.pdg(); r.parent_source = std::get<1>(p.key()); r.parent_id = std::get<2>(p.key()); r.parent_flight_cm = flight(p); }
    if (hadron(n.pdg()))
    {
      int state = 0; const auto ds = daughters(n, state);
      if (state < 0 && r.proper_length_cm >= 0 && r.proper_length_cm < kMaxProperLengthCm && r.status >= 0) r.status = state;
      if (!ds.empty())
      {
        const int rule = prompt(n);
        r.long_lived = rule == 0;
        if (rule < 0 && r.status >= 0) r.status = rule;
      }
      else if (r.proper_length_cm >= 0 && r.proper_length_cm < kMaxProperLengthCm && r.status >= 0) r.status = -1;
    }
    return r;
  }
  bool frontier(Node n, std::vector<Node>& leaves, std::set<Key>& active, int& pathway, int& failure) const
  {
    if (!n || !active.insert(n.key()).second) { failure = -1; return false; }
    int state = 0; const auto ds = daughters(n, state);
    if (!hadron(n.pdg()) || n.pdg() == 111 || (n.pdg() == 221 && diphoton(n))) leaves.push_back(n);
    else
    {
      const int rule = prompt(n);
      if (rule == 0) leaves.push_back(n);
      else if (rule < 0 || state <= 0 || ds.empty()) { failure = state < 0 ? state : rule < 0 ? rule : -1; return false; }
      else
      {
        pathway |= state;
        for (const auto& child : ds) if (!frontier(child, leaves, active, pathway, failure)) return false;
      }
    }
    active.erase(n.key()); return true;
  }
  // Measured hit shares always use physical G4 tracks. A HepMC-only pi0 is covered by its G4 photon daughters.
  bool target_tracks(Node n, std::vector<int>& result, std::set<Key>& active) const
  {
    if (!n || !active.insert(n.key()).second) return false;
    if (n.g4) result.push_back(n.g4->get_track_id());
    else
    {
      int state = 0; const auto ds = daughters(n, state);
      if (state <= 0 || ds.empty()) return false;
      for (const auto& child : ds) if (!target_tracks(child, result, active)) return false;
    }
    active.erase(n.key()); return true;
  }
  double measure(Node n, const RawCluster* cluster, TowerInfoContainer* towers, RawTowerContainer* raw, PHG4CellContainer* cells, PHG4HitContainer* hits) const
  {
    std::vector<int> ids; std::set<Key> active;
    if (!target_tracks(n, ids, active)) return -1;
    std::set<int> unique(ids.begin(), ids.end()); if (unique.size() != ids.size()) return -1;
    std::vector<std::array<int, 2>> tracks; for (int id : ids) tracks.push_back({id, 0});
    const auto matches = Pi0ClusterTruthMatcher().match_many(cluster, towers, raw, cells, hits, truth, tracks, true);
    if (matches.size() != ids.size()) return -1;
    double result = 0;
    for (const auto& m : matches) { if (m.status != Pi0ClusterTruthMatchStatus::complete) return -1; result += m.gamma_edep[0]; }
    return result;
  }
  void photon_shares(DecayOriginContribution& r, Node atom, Node input, const RawCluster* cluster, TowerInfoContainer* towers,
                     RawTowerContainer* raw, PHG4CellContainer* cells, PHG4HitContainer* hits) const
  {
    if (!r.diphoton) return;
    int state = 0; const auto ds = daughters(atom, state);
    // Deterministic photon ordering independent of container iteration order.
    auto photons = ds; std::sort(photons.begin(), photons.end(), [](const Node& a, const Node& b) { return a.key() < b.key(); });
    std::set<Key> seen;
    for (auto current = input; current && seen.insert(current.key()).second; current = parent(current))
    {
      for (int i = 0; i < 2; ++i) if (physical(current).key() == physical(photons[i]).key())
      { r.gamma0_fraction = i == 0 ? r.fraction : 0; r.gamma1_fraction = i == 1 ? r.fraction : 0; return; }
      if (current.key() == atom.key()) break;
    }
    const double total = measure(atom, cluster, towers, raw, cells, hits);
    const double a = measure(photons[0], cluster, towers, raw, cells, hits), b = measure(photons[1], cluster, towers, raw, cells, hits);
    if (total > 0 && a >= 0 && b >= 0 && close(a + b, total))
    { r.gamma0_fraction = r.fraction * a / total; r.gamma1_fraction = r.fraction * b / total; }
  }
};
DecayOriginDiagnostics::DecayOriginDiagnostics(PHG4TruthInfoContainer* truth, const PHHepMCGenEventMap* events) : impl_(std::make_unique<Impl>(truth, events)) {}
DecayOriginDiagnostics::~DecayOriginDiagnostics() = default;
DecayOriginMatch DecayOriginDiagnostics::diagnose(const ClusterTruthMatch& original, const RawCluster* cluster, TowerInfoContainer* towers,
                                                 RawTowerContainer* raw, PHG4CellContainer* cells, PHG4HitContainer* hits) const
{
  DecayOriginMatch result; result.valid = original.valid;
  std::map<Key, DecayOriginContribution> grouped;
  const auto append = [&](DecayOriginContribution r) {
    const Key key{r.embedding, r.source, r.id};
    auto [it, inserted] = grouped.emplace(key, r);
    if (!inserted)
    {
      auto& old = it->second; old.fraction += r.fraction; old.pathway |= r.pathway; old.status = old.status < 0 || r.status < 0 ? std::min(old.status, r.status) : std::max(old.status, r.status);
      if (old.root_track_id != r.root_track_id) old.root_track_id = 0;
      if (old.gamma0_fraction >= 0 && r.gamma0_fraction >= 0) { old.gamma0_fraction += r.gamma0_fraction; old.gamma1_fraction += r.gamma1_fraction; }
      else old.gamma0_fraction = old.gamma1_fraction = -1;
    }
  };
  for (const auto& c : original.contributors)
  {
    const auto input = impl_->from_g4(impl_->truth ? impl_->truth->GetParticle(c.g4_track_id) : nullptr, c.embedding_id);
    if (!input)
    {
      DecayOriginContribution r; r.pdg = c.g4_pdg_id; r.embedding = c.embedding_id; r.source = 2; r.id = r.root_track_id = c.g4_track_id;
      r.fraction = c.fraction; r.status = -1; append(r); continue;
    }
    if (input.g4->get_parent_id() == 0 && input.g4->get_barcode() != 0)
    {
      const auto key = std::make_pair(c.embedding_id, input.g4->get_barcode());
      const auto* sub = impl_->events ? impl_->events->get(c.embedding_id) : nullptr;
      const auto* event = sub ? sub->getEvent() : nullptr;
      const auto* hep = event ? event->barcode_to_particle(input.g4->get_barcode()) : nullptr;
      if (impl_->ambiguous_tracks.count(key) || (hep && hep->pdg_id() != input.g4->get_pid()))
      { append(impl_->describe(input, c.g4_track_id, c.fraction, -1, 0)); continue; }
    }
    const auto atom = impl_->anchor(input);
    if (atom.key() != input.key())
    {
      auto r = impl_->describe(atom, c.g4_track_id, c.fraction, 0, impl_->ancestry_path(input, atom));
      impl_->photon_shares(r, atom, input, cluster, towers, raw, cells, hits); append(r); continue;
    }
    std::vector<Node> leaves; std::set<Key> active; int pathway = 0, failure = 0;
    if (!impl_->frontier(atom, leaves, active, pathway, failure))
    { append(impl_->describe(atom, c.g4_track_id, c.fraction, failure, pathway)); continue; }
    if (leaves.size() == 1 && leaves[0].key() == atom.key())
    {
      auto r = impl_->describe(atom, c.g4_track_id, c.fraction, 0, pathway);
      impl_->photon_shares(r, atom, input, cluster, towers, raw, cells, hits); append(r); continue;
    }
    const double root = impl_->measure(input, cluster, towers, raw, cells, hits);
    std::vector<double> measured; double sum = 0; bool complete = root > 0;
    for (const auto& leaf : leaves)
    { const double e = impl_->measure(leaf, cluster, towers, raw, cells, hits); complete = complete && e >= 0; measured.push_back(e); sum += e; }
    if (!complete || !close(sum, root)) { append(impl_->describe(atom, c.g4_track_id, c.fraction, -1, pathway)); continue; }
    for (std::size_t i = 0; i < leaves.size(); ++i)
    {
      if (!(measured[i] > 0)) continue;
      auto r = impl_->describe(leaves[i], c.g4_track_id, c.fraction * measured[i] / sum, 1, pathway);
      impl_->photon_shares(r, leaves[i], leaves[i], cluster, towers, raw, cells, hits); append(r);
    }
  }
  for (const auto& entry : grouped) result.contributors.push_back(entry.second);
  std::sort(result.contributors.begin(), result.contributors.end(), [](const auto& a, const auto& b) {
    return a.fraction == b.fraction ? Key{a.embedding, a.source, a.id} < Key{b.embedding, b.source, b.id} : a.fraction > b.fraction;
  });
  double before = 0, after = 0;
  for (const auto& c : original.contributors) before += c.fraction;
  for (const auto& c : result.contributors) after += c.fraction;
  if (!close(before, after)) result.valid = false;
  return result;
}
void DecayOriginCollection::clear() { *this = DecayOriginCollection{}; offset.push_back(0); }
void DecayOriginCollection::append(const DecayOriginMatch& match)
{
  if (offset.empty()) offset.push_back(0);
  valid.push_back(match.valid);
  for (const auto& c : match.contributors)
  {
#define APPEND(member) member.push_back(c.member)
    APPEND(pdg); APPEND(embedding); APPEND(source); APPEND(id); APPEND(root_track_id); APPEND(parent_pdg); APPEND(parent_source); APPEND(parent_id);
    APPEND(daughter_count); APPEND(daughter_pi0_count); APPEND(daughter_gamma_count); APPEND(daughter_charged_count);
    APPEND(status); APPEND(pathway); APPEND(diphoton); APPEND(long_lived); APPEND(fraction); APPEND(gamma0_fraction); APPEND(gamma1_fraction);
    APPEND(flight_cm); APPEND(parent_flight_cm); APPEND(proper_length_cm);
#undef APPEND
  }
  offset.push_back(pdg.size());
}
void DecayOriginCollection::create_branches(TTree* tree)
{
#define BRANCH(member) tree->Branch((std::string("split_cluster_decay_origin_") + #member).c_str(), &member)
  BRANCH(valid); BRANCH(offset); BRANCH(pdg); BRANCH(embedding); BRANCH(source); BRANCH(id); BRANCH(root_track_id);
  BRANCH(daughter_count); BRANCH(daughter_pi0_count); BRANCH(daughter_gamma_count); BRANCH(daughter_charged_count);
  BRANCH(parent_pdg); BRANCH(parent_source); BRANCH(parent_id); BRANCH(status); BRANCH(pathway); BRANCH(diphoton); BRANCH(long_lived);
  BRANCH(fraction); BRANCH(gamma0_fraction); BRANCH(gamma1_fraction); BRANCH(flight_cm); BRANCH(parent_flight_cm); BRANCH(proper_length_cm);
#undef BRANCH
}
}
