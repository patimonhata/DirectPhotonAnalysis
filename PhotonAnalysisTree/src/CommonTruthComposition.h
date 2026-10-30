#ifndef PHOTON_TREE_COMMON_TRUTH_COMPOSITION_H
#define PHOTON_TREE_COMMON_TRUTH_COMPOSITION_H
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <tuple>
#include <vector>
namespace photon_tree::composition
{
constexpr int kAlgorithmVersion = 1;
// Identities are scoped to the containing event_uid. Source 1 = HepMC, 2 = G4.
struct Identity
{
  int embedding = 0, source = 0, id = 0;
  auto key() const { return std::make_tuple(embedding, source, id); }
  bool operator<(const Identity& b) const { return key() < b.key(); }
  bool operator==(const Identity& b) const { return key() == b.key(); }
  bool operator!=(const Identity& b) const { return !(*this == b); }
  bool valid() const { return (source == 1 || source == 2) && id != 0; }
};
enum class Main : int { unknown = 0, prompt = 1, pi0 = 2, eta_gg = 3, decay_photon = 4, charged = 5,
                        particle = 6, parent = 7, parent_unresolved = 8, mixed = 9, non_signal = 10, photon_unknown = 11 };
enum class Applicability : int { unknown = 0, applicable = 1, not_applicable = 2 };
enum Quality : unsigned int { missing_particle = 1, ambiguous_identity = 2, incomplete_split = 4, interaction = 8,
                              unknown_species = 16, ancestry_cycle = 32, invalid_deposit = 64, unknown_decay = 128 };
struct Particle
{
  Identity identity;
  int pdg = 0, decay = 0, prompt = 0; // decay: 0 unknown, 1 direct gg, 2 confirmed other, 3 stable.
  bool hadron = false, charged = false, species_known = false;
  double proper_length_cm = -1, flight_cm = -1;
  std::vector<double> production_vertex, decay_vertex; // x,y,z,t in cm,ns; empty when unknown.
  unsigned int quality = 0;
};
struct Piece
{
  Particle particle;
  double deposit = 0;
  bool unresolved = false;
  std::vector<Particle> ancestors; // immediate to upstream, confirmed hadron decays only.
};
struct Result
{
  bool valid = false;
  double total = 0, attributed = 0, unresolved = 0, max_photon = 0, other_photons = 0, charged = 0;
  unsigned int quality = 0;
  std::vector<Piece> units, photons;
  std::vector<Piece> parents;
  Main main = Main::unknown;
  Particle selected;
  double fraction = 0;
  Applicability applicability = Applicability::unknown;
  int matching_status = 0;
  double minimum_coverage = 0;
  int topology_reason = 0; // 1 non-gg, 2 unknown mode, 3 missing candidate, 4 matching, 5 projection, 6 first-daughter interaction information.
};
inline int overview(Main main, int pdg)
{
  if (main == Main::prompt) return 0;
  if (main == Main::pi0) return 1;
  if (main == Main::eta_gg || ((main == Main::parent || main == Main::parent_unresolved) && pdg == 221)) return 2;
  return 3;
}
inline void classify(Result& r, int signal_embedding = 1)
{
  r.main = Main::unknown; r.selected = {}; r.fraction = 0; r.parents.clear();
  if (!r.valid || !(r.total > 0) || !std::isfinite(r.total)) return;
  std::map<Identity, Piece> units, parents;
  double non_signal = 0, unknown = 0;
  for (const auto& piece : r.units)
  {
    if (!std::isfinite(piece.deposit) || piece.deposit < 0 || !piece.particle.identity.valid()) { r.quality |= invalid_deposit; return; }
    auto [it, fresh] = units.emplace(piece.particle.identity, piece);
    if (!fresh && it->second.particle.pdg != piece.particle.pdg) { r.quality |= ambiguous_identity; return; }
    if (!fresh) { it->second.deposit += piece.deposit; it->second.unresolved |= piece.unresolved; }
    if (piece.particle.identity.embedding != signal_embedding) non_signal += piece.deposit;
    if (piece.unresolved) unknown += piece.deposit;
    std::set<Identity> seen;
    const bool known_identity = !(piece.particle.quality & (missing_particle | ambiguous_identity | invalid_deposit | ancestry_cycle));
    if (piece.unresolved && piece.particle.hadron && known_identity) parents.emplace(piece.particle.identity, Piece{piece.particle, 0, true, piece.ancestors});
    std::vector<Particle> chain = known_identity ? piece.ancestors : std::vector<Particle>{};
    if (piece.unresolved && piece.particle.hadron && known_identity) chain.insert(chain.begin(), piece.particle);
    for (const auto& parent : chain)
    {
      if (!parent.hadron || !parent.identity.valid() || !seen.insert(parent.identity).second) { r.quality |= ancestry_cycle; return; }
      auto [found, inserted] = parents.emplace(parent.identity, Piece{parent, 0, false, {}});
      found->second.deposit += piece.deposit; found->second.unresolved |= piece.unresolved;
    }
  }
  for (const auto& entry : parents) r.parents.push_back(entry.second);
  const auto select = [&](const Piece& p, Main main) { r.selected = p.particle; r.fraction = p.deposit / r.total; r.main = main; };
  for (const auto& entry : units)
  {
    const auto& p = entry.second; const auto& n = p.particle;
    if (!(p.deposit > r.total * 0.5) || p.unresolved) continue;
    if (n.identity.embedding != signal_embedding) { select(p, Main::non_signal); return; }
    if (!n.species_known) { select(p, Main::unknown); return; }
    if (n.pdg == 111) select(p, Main::pi0);
    else if (n.pdg == 221 && n.decay == 1) select(p, Main::eta_gg);
    else if (n.pdg == 22)
    {
      if (n.prompt == 1 || n.prompt == 2) select(p, Main::prompt);
      else select(p, p.ancestors.empty() ? Main::photon_unknown : Main::decay_photon);
    }
    else select(p, n.charged ? Main::charged : Main::particle);
    return;
  }
  // Majority parents must form one chain. Select the downstream-most member of that chain.
  std::vector<const Piece*> majority;
  for (const auto& p : r.parents) if (p.deposit > r.total * 0.5) majority.push_back(&p);
  for (const auto* p : majority)
  {
    bool upstream = false;
    for (const auto* q : majority)
    {
      if (q == p) continue;
      for (const auto& unit : r.units)
      {
        bool below_q = unit.particle.identity == q->particle.identity;
        for (const auto& a : unit.ancestors)
        {
          if (a.identity == p->particle.identity && below_q) upstream = true;
          if (a.identity == q->particle.identity) below_q = true;
        }
      }
    }
    if (!upstream)
    {
      select(*p, p->particle.identity.embedding != signal_embedding ? Main::non_signal : p->unresolved ? Main::parent_unresolved : Main::parent);
      return;
    }
  }
  if (non_signal > r.total * 0.5) r.main = Main::non_signal;
  else r.main = unknown > 0 || (r.quality & (missing_particle | ambiguous_identity | ancestry_cycle | invalid_deposit | unknown_species)) ? Main::unknown : Main::mixed;
}
}
#endif
