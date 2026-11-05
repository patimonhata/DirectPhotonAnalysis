#include "CommonDecayTruth.h"
#include "HepMCPhotonClassifier.h"
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
#include <functional>
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
struct CommonDecayTruth::Impl
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
    if (n.pdg() == 22 && p.pdg() == 22) return p;
    int state = 0; daughters(p, state);
    return hadron(p.pdg()) && state > 0 ? p : Node{};
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
    // A retained same-species hadron signals bookkeeping/radiation or interaction, not a confirmed decay.
    if (hadron(n.pdg()) && std::any_of(result.begin(), result.end(), [&](Node child) { return child.pdg() == n.pdg(); })) { state = -1; return {}; }
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
  composition::Particle describe(Node n) const
  {
    composition::Particle r; const auto [embedding, source, id] = n.key();
    r.identity = {embedding, source, id}; r.pdg = n.pdg(); r.hadron = hadron(r.pdg);
    const auto* info = TDatabasePDG::Instance()->GetParticle(r.pdg);
    r.species_known = info != nullptr; r.charged = info && info->Charge() != 0;
    if (!info) r.quality |= composition::unknown_species;
    r.flight_cm = flight(n); r.proper_length_cm = proper_length(r.pdg);
    int state = 0; const auto ds = daughters(n, state);
    r.decay = state > 0 ? ds.size() == 2 && ds[0].pdg() == 22 && ds[1].pdg() == 22 ? 1 : 2 : info && info->Stable() ? 3 : 0;
    if (state == -2) r.quality |= composition::interaction;
    if (state == -1) r.quality |= composition::incomplete_split;
    if (r.hadron && r.decay == 0) r.quality |= composition::unknown_decay;
    if (n.hep && r.pdg == 22) { const auto c = HepMCPhotonClassifier().classify(n.hep); if (c.valid) r.prompt = c.category; }
    if (n.g4)
    {
      const auto* v = truth->GetVtx(n.g4->get_vtx_id());
      if (v) r.production_vertex = {v->get_x(), v->get_y(), v->get_z(), v->get_t()};
      const auto found = children.find(n.g4->get_track_id());
      const auto* end = found != children.end() && !found->second.empty() ? truth->GetVtx(found->second.front()->get_vtx_id()) : nullptr;
      if (end && state > 0) r.decay_vertex = {end->get_x(), end->get_y(), end->get_z(), end->get_t()};
    }
    if (n.hep && n.hep->parent_event())
    {
      const double cm = HepMC::Units::conversion_factor(n.hep->parent_event()->length_unit(), HepMC::Units::CM);
      const auto coordinates = [cm](const HepMC::GenVertex* v) {
        return v ? std::vector<double>{v->position().x() * cm, v->position().y() * cm, v->position().z() * cm,
                                      v->position().t() * cm / 29.9792458} : std::vector<double>{};
      };
      if (r.production_vertex.empty()) r.production_vertex = coordinates(n.hep->production_vertex());
      if (r.decay_vertex.empty()) r.decay_vertex = coordinates(n.hep->end_vertex());
    }
    return r;
  }
  bool diphoton(Node n) const { return describe(n).decay == 1; }
  std::vector<composition::Particle> ancestry(Node n, unsigned int& quality) const
  {
    std::vector<composition::Particle> result; std::set<Key> seen{n.key()};
    for (auto a = physical(parent(n)); a; a = physical(parent(a)))
    {
      if (!seen.insert(a.key()).second) { quality |= composition::ancestry_cycle; break; }
      if (!hadron(a.pdg())) break;
      result.push_back(describe(a));
    }
    return result;
  }
  Node atom(Node n) const
  {
    n = physical(n); std::set<Key> seen;
    if (n.pdg() == 111 || (n.pdg() == 221 && diphoton(n))) return n;
    for (auto a = physical(parent(n)); a && seen.insert(a.key()).second; a = physical(parent(a)))
    {
      if (!hadron(a.pdg())) break;
      if (a.pdg() == 111 || (a.pdg() == 221 && diphoton(a))) return a;
    }
    return n;
  }
  void known_photons(Node n, std::vector<Node>& photons, std::set<Key>& seen) const
  {
    n = physical(n); if (!n || !seen.insert(n.key()).second) return;
    if (n.pdg() == 22) { photons.push_back(n); return; }
    if (!hadron(n.pdg())) return; // Never unfold an electromagnetic shower into independent photons.
    int state = 0; const auto ds = daughters(n, state);
    if (state > 0) for (const auto& child : ds) known_photons(child, photons, seen);
  }
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
    if (std::set<int>(ids.begin(), ids.end()).size() != ids.size()) return -1;
    std::vector<std::array<int, 2>> tracks; for (int id : ids) tracks.push_back({id, 0});
    const auto matches = Pi0ClusterTruthMatcher().match_many(cluster, towers, raw, cells, hits, truth, tracks, true);
    if (matches.size() != ids.size()) return -1;
    double sum = 0;
    for (const auto& m : matches) { if (m.status != Pi0ClusterTruthMatchStatus::complete || !m.ancestry_complete) return -1; sum += m.gamma_edep[0]; }
    return sum;
  }
  composition::Piece piece(Node n, double deposit, bool unresolved) const
  {
    composition::Piece p{describe(n), deposit, unresolved, {}};
    p.ancestors = ancestry(n, p.particle.quality); return p;
  }
};
CommonDecayTruth::CommonDecayTruth(PHG4TruthInfoContainer* truth, const PHHepMCGenEventMap* events) : impl_(std::make_unique<Impl>(truth, events)) {}
CommonDecayTruth::~CommonDecayTruth() = default;
composition::Result CommonDecayTruth::resolve(const ClusterTruthMatch& original, const RawCluster* cluster, TowerInfoContainer* towers,
                                             RawTowerContainer* raw, PHG4CellContainer* cells, PHG4HitContainer* hits) const
{
  composition::Result r; r.valid = original.valid; r.total = original.total_edep;
  if (!original.valid) r.quality |= composition::missing_particle;
  if (!(r.total > 0) || !std::isfinite(r.total)) r.quality |= composition::invalid_deposit;
  const auto unknown_photons = [&](Node n) {
    std::vector<Node> photons; std::set<Key> seen; impl_->known_photons(n, photons, seen);
    for (const auto& photon : photons)
    { auto p = impl_->piece(photon, -1, true); p.particle.quality |= composition::incomplete_split; r.photons.push_back(p); }
  };
  const auto add = [&](Node n, Node input, double deposit, bool unresolved) {
    auto p = impl_->piece(n, deposit, unresolved); p.unresolved |= !p.particle.species_known; r.quality |= p.particle.quality;
    r.units.push_back(p); r.attributed += deposit;
    if (p.unresolved) { r.unresolved += deposit; unknown_photons(n); return; }
    if (n.pdg() == 22) r.photons.push_back(p);
    else if (p.particle.charged && n.pdg() != 111 && n.pdg() != 221) r.charged += deposit;
    // Purity components are separate from the composition units, including non-gg pi0.
    if (n.pdg() == 111 || (n.pdg() == 221 && impl_->diphoton(n)))
    {
      if (n.key() != input.key())
      {
        const auto physical = impl_->physical(input);
        if (physical.pdg() == 22) r.photons.push_back(impl_->piece(physical, deposit, false));
        else if (impl_->describe(physical).charged) r.charged += deposit;
        else r.unresolved += deposit;
        return;
      }
      int state = 0; const auto ds = impl_->daughters(n, state);
      const double root = impl_->measure(n, cluster, towers, raw, cells, hits);
      std::vector<double> values; double sum = 0; bool complete = state > 0 && root > 0;
      for (const auto& d : ds) { double v = impl_->measure(d, cluster, towers, raw, cells, hits); complete &= v >= 0; values.push_back(v); sum += v; }
      if (!complete || !close(sum, root)) { r.unresolved += deposit; r.quality |= composition::incomplete_split; unknown_photons(n); return; }
      for (std::size_t i = 0; i < ds.size(); ++i)
      {
        if (!(values[i] > 0)) continue;
        const double share = deposit * values[i] / sum; const auto d = impl_->physical(ds[i]);
        if (d.pdg() == 22) r.photons.push_back(impl_->piece(d, share, false));
        else if (impl_->describe(d).charged) r.charged += share;
        else r.unresolved += share;
      }
    }
  };
  std::function<void(Node, Node, double, std::set<Key>&)> expand;
  expand = [&](Node n, Node input, double deposit, std::set<Key>& active) {
    n = impl_->physical(n);
    if (!n || !active.insert(n.key()).second)
    { r.quality |= composition::ancestry_cycle; add(n, input, deposit, true); r.units.back().particle.quality |= composition::ancestry_cycle; return; }
    const auto info = impl_->describe(n); int state = 0; const auto ds = impl_->daughters(n, state);
    if (!info.hadron || n.pdg() == 111 || (n.pdg() == 221 && info.decay == 1) || info.decay == 3 || state == -2)
      add(n, input, deposit, false);
    else if (state <= 0 || ds.empty())
    { r.quality |= composition::incomplete_split; add(n, input, deposit, true); }
    else
    {
      // Resolve one decay at a time. Failure further downstream must not erase a measured sibling majority.
      const double root = impl_->measure(n, cluster, towers, raw, cells, hits);
      std::vector<double> values; double sum = 0; bool complete = root > 0; std::set<Key> unique;
      for (const auto& child : ds)
      { const double v = impl_->measure(child, cluster, towers, raw, cells, hits); complete &= v >= 0 && unique.insert(impl_->physical(child).key()).second; values.push_back(v); sum += v; }
      if (!complete || !close(sum, root)) { r.quality |= composition::incomplete_split; add(n, input, deposit, true); }
      else for (std::size_t i = 0; i < ds.size(); ++i) if (values[i] > 0) expand(ds[i], ds[i], deposit * values[i] / sum, active);
    }
    active.erase(n.key());
  };
  std::set<int> roots;
  for (const auto& c : original.contributors)
  {
    if (!roots.insert(c.g4_track_id).second) { r.valid = false; r.quality |= composition::ambiguous_identity; continue; }
    const auto input = impl_->from_g4(impl_->truth ? impl_->truth->GetParticle(c.g4_track_id) : nullptr, c.embedding_id);
    if (!input)
    {
      composition::Piece p; p.particle.identity = {c.embedding_id, 2, c.g4_track_id}; p.particle.pdg = c.g4_pdg_id;
      p.deposit = c.edep; p.unresolved = true; r.units.push_back(p); r.attributed += c.edep; r.unresolved += c.edep; r.quality |= composition::missing_particle; continue;
    }
    if (input.g4->get_parent_id() == 0 && input.g4->get_barcode() != 0)
    {
      const auto key = std::make_pair(c.embedding_id, input.g4->get_barcode());
      const auto* sub = impl_->events ? impl_->events->get(c.embedding_id) : nullptr;
      const auto* event = sub ? sub->getEvent() : nullptr;
      const auto* hep = event ? event->barcode_to_particle(input.g4->get_barcode()) : nullptr;
      if (impl_->ambiguous_tracks.count(key) || (hep && hep->pdg_id() != input.g4->get_pid()))
      { r.quality |= composition::ambiguous_identity; add(input, input, c.edep, true); r.units.back().particle.quality |= composition::ambiguous_identity; continue; }
    }
    const auto atom = impl_->atom(input);
    if (atom.key() != input.key()) { add(atom, input, c.edep, false); continue; }
    std::set<Key> active; expand(atom, input, c.edep, active);
  }
  std::map<composition::Identity, double> photons;
  for (const auto& p : r.photons) if (!p.unresolved && p.deposit >= 0) photons[p.particle.identity] += p.deposit;
  double photon_total = 0;
  for (const auto& entry : photons) { r.max_photon = std::max(r.max_photon, entry.second); photon_total += entry.second; }
  r.other_photons = photon_total - r.max_photon;
  if (!close(r.total, r.attributed)) { r.valid = false; r.quality |= composition::invalid_deposit; }
  composition::classify(r);
  if (r.main == composition::Main::pi0 || r.main == composition::Main::eta_gg)
  {
    r.applicability = r.selected.decay == 1 ? composition::Applicability::unknown :
        r.selected.decay == 2 ? composition::Applicability::not_applicable : composition::Applicability::unknown;
    r.topology_reason = r.selected.decay == 2 ? 1 : r.selected.decay != 1 ? 2 : 3;
  }
  else if (r.main == composition::Main::parent || r.main == composition::Main::parent_unresolved)
  {
    const bool unknown_eta_mode = r.selected.pdg == 221 && r.selected.decay == 0;
    r.applicability = unknown_eta_mode ? composition::Applicability::unknown : composition::Applicability::not_applicable;
    r.topology_reason = unknown_eta_mode ? 2 : 1;
  }
  return r;
}
std::vector<CommonDiphoton> CommonDecayTruth::diphotons(int pdg, int embedding) const
{
  std::map<Key, Node> nodes;
  const auto* sub = impl_->events ? impl_->events->get(embedding) : nullptr;
  const auto* event = sub ? sub->getEvent() : nullptr;
  if (event) for (auto it = event->particles_begin(); it != event->particles_end(); ++it)
    if ((*it)->pdg_id() == pdg) { auto n = impl_->from_hep(*it, embedding); nodes.emplace(n.key(), n); }
  if (impl_->truth)
  {
    const auto range = impl_->truth->GetParticleRange();
    for (auto it = range.first; it != range.second; ++it)
    {
      const auto* p = it->second; if (!p || p->get_pid() != pdg) continue;
      const auto* root = p; std::set<int> seen;
      while (root && root->get_parent_id() && seen.insert(root->get_track_id()).second) root = impl_->truth->GetParticle(root->get_parent_id());
      if (!root || root->get_parent_id() || impl_->truth->isEmbeded(root->get_track_id()) != embedding) continue;
      auto n = impl_->from_g4(p, embedding); nodes.emplace(n.key(), n);
    }
  }
  std::vector<CommonDiphoton> result;
  for (const auto& entry : nodes)
  {
    const auto n = entry.second; if (!impl_->diphoton(n)) continue;
    int state = 0; auto ds = impl_->daughters(n, state);
    std::sort(ds.begin(), ds.end(), [](Node a, Node b) { return a.key() < b.key(); });
    CommonDiphoton c; c.parent = impl_->describe(n); c.g4_parent_track = n.g4 ? n.g4->get_track_id() : 0;
    for (int i = 0; i < 2; ++i)
    {
      const auto d = impl_->physical(ds[i]); const auto [e, s, id] = d.key(); c.photon_ids[i] = {e, s, id};
      std::vector<int> tracks; std::set<Key> seen;
      if (impl_->target_tracks(d, tracks, seen) && tracks.size() == 1) c.photons[i] = impl_->truth->GetParticle(tracks[0]);
    }
    if (!c.photons[0] || !c.photons[1] || c.photons[0] == c.photons[1]) continue;
    const double px = n.g4 ? n.g4->get_px() : n.hep->momentum().px();
    const double py = n.g4 ? n.g4->get_py() : n.hep->momentum().py();
    const double pz = n.g4 ? n.g4->get_pz() : n.hep->momentum().pz();
    c.energy = n.g4 ? n.g4->get_e() : n.hep->momentum().e(); c.pt = std::hypot(px, py);
    c.eta = c.pt > 0 ? std::asinh(pz / c.pt) : 0; c.phi = std::atan2(py, px); result.push_back(c);
  }
  return result;
}
void CommonTruthCollection::clear() { *this = CommonTruthCollection{}; offset.push_back(0); ancestor_offset.push_back(0); }
void CommonTruthCollection::append(const composition::Result& r, int legacy)
{
  if (offset.empty()) clear();
  bool long_lived = false;
  for (const auto& p : r.units) for (const auto& a : p.ancestors) long_lived |= a.proper_length_cm >= 1e-4;
  long_lived_ancestor.push_back(long_lived);
  main.push_back(static_cast<int>(r.main)); selected_pdg.push_back(r.selected.pdg);
  selected_embedding.push_back(r.selected.identity.embedding); selected_source.push_back(r.selected.identity.source); selected_id.push_back(r.selected.identity.id);
  applicability.push_back(static_cast<int>(r.applicability)); topology_reason.push_back(r.topology_reason); legacy_origin.push_back(legacy);
#define PUSH(field) field.push_back(r.field)
  quality.push_back(r.quality | (!r.valid ? composition::missing_particle : 0U) | (!(r.total > 0) ? composition::invalid_deposit : 0U));
  PUSH(total); PUSH(attributed); PUSH(unresolved); PUSH(fraction); PUSH(max_photon); PUSH(other_photons); PUSH(charged); PUSH(matching_status); PUSH(minimum_coverage);
  PUSH(anchor_matching_status); PUSH(anchor_matching_coverage); PUSH(anchor_matching_usable);
#undef PUSH
  const std::array<const std::vector<composition::Piece>*, 3> tables = {&r.units, &r.photons, &r.parents};
  for (int t = 0; t < 3; ++t) for (const auto& p : *tables[t])
  {
    const auto& n = p.particle; table.push_back(t); pdg.push_back(n.pdg); embedding.push_back(n.identity.embedding); source.push_back(n.identity.source); id.push_back(n.identity.id);
    decay.push_back(n.decay); prompt.push_back(n.prompt); unresolved_piece.push_back(p.unresolved); piece_quality.push_back(n.quality); deposit.push_back(p.deposit);
    flight_cm.push_back(n.flight_cm); proper_length_cm.push_back(n.proper_length_cm); production_vertex.push_back(n.production_vertex); decay_vertex.push_back(n.decay_vertex);
    for (const auto& a : p.ancestors)
    {
      ancestor_pdg.push_back(a.pdg); ancestor_embedding.push_back(a.identity.embedding); ancestor_source.push_back(a.identity.source); ancestor_id.push_back(a.identity.id);
      ancestor_flight_cm.push_back(a.flight_cm); ancestor_proper_length_cm.push_back(a.proper_length_cm);
    }
    ancestor_offset.push_back(ancestor_id.size());
  }
  offset.push_back(id.size());
}
void CommonTruthCollection::create_branches(TTree* tree)
{
#define B(field) tree->Branch("split_cluster_common_" #field, &field)
  B(main); B(selected_pdg); B(selected_embedding); B(selected_source); B(selected_id); B(applicability); B(topology_reason); B(legacy_origin);
  B(long_lived_ancestor); B(quality); B(total); B(attributed); B(unresolved); B(fraction); B(max_photon); B(other_photons); B(charged); B(matching_status); B(minimum_coverage);
  B(anchor_matching_status); B(anchor_matching_coverage); B(anchor_matching_usable);
  B(offset); B(ancestor_offset); B(table); B(pdg); B(embedding); B(source); B(id); B(decay); B(prompt); B(unresolved_piece); B(piece_quality); B(deposit);
  B(flight_cm); B(proper_length_cm); B(production_vertex); B(decay_vertex);
  B(ancestor_pdg); B(ancestor_embedding); B(ancestor_source); B(ancestor_id); B(ancestor_flight_cm); B(ancestor_proper_length_cm);
#undef B
}
}
