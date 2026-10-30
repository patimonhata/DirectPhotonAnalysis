#include "CommonTruthComposition.h"
#include <cassert>
#include <iostream>
using namespace photon_tree::composition;
Particle node(int pdg, int id, int embedding = 1)
{
  Particle p; p.identity = {embedding, 1, id}; p.pdg = pdg; p.species_known = true;
  p.hadron = pdg == 111 || pdg == 221 || pdg == 223 || pdg == 310 || pdg == 331;
  p.charged = pdg == 211 || pdg == -211 || pdg == 11; p.decay = pdg == 111 || pdg == 221 ? 1 : p.hadron ? 2 : 3;
  return p;
}
Piece piece(Particle p, double deposit, std::vector<Particle> ancestors = {}, bool unresolved = false) { return {p, deposit, unresolved, ancestors}; }
Result run(std::vector<Piece> units)
{
  Result r; r.valid = true; r.total = 100; r.units = std::move(units); classify(r); return r;
}
int main()
{
  auto pi = node(111, 1), pi2 = node(111, 2), eta = node(221, 3), eta2 = node(221, 4), omega = node(223, 5);
  auto gamma = node(22, 6), charged = node(211, 7), ks = node(310, 8), prime = node(331, 9);
  auto prompt = gamma; prompt.prompt = 1;
  assert(run({piece(prompt, 80)}).main == Main::prompt);
  // Multiple deposits of the same atom combine; distinct particles and embeddings never do.
  assert(run({piece(pi, 60), piece(pi, 20)}).main == Main::pi0);
  assert(run({piece(eta, 70), piece(eta, 10)}).main == Main::eta_gg);
  assert(run({piece(pi, 70, {omega}), piece(charged, 20, {omega})}).main == Main::pi0);
  assert(run({piece(gamma, 70, {omega}), piece(pi, 20, {omega})}).main == Main::decay_photon);
  auto r = run({piece(pi, 40, {omega}), piece(gamma, 40, {omega})}); assert(r.main == Main::parent && r.selected.pdg == 223 && r.fraction == .8);
  eta.decay = 2;
  assert(run({piece(pi, 30, {eta}), piece(pi2, 30, {eta}), piece(node(111, 10), 40, {eta})}).main == Main::parent);
  assert(run({piece(pi, 20, {eta}), piece(pi2, 20, {eta}), piece(node(111, 10), 60, {eta})}).main == Main::pi0);
  assert(run({piece(charged, 70, {eta}), piece(node(-211, 11), 10, {eta}), piece(pi, 20, {eta})}).main == Main::charged);
  ks.proper_length_cm = 2.68; ks.flight_cm = 100;
  assert(run({piece(pi, 80, {ks})}).main == Main::pi0);
  r = run({piece(pi, 40, {ks}), piece(pi2, 40, {ks})}); assert(r.main == Main::parent && r.selected.pdg == 310);
  assert(run({piece(eta, 30), piece(eta2, 30)}).main == Main::mixed);
  r = run({piece(pi, 30, {eta, prime}), piece(pi2, 30, {eta, prime})}); assert(r.main == Main::parent && r.selected.identity == eta.identity);
  r = run({piece(pi, 40, {eta, prime}), piece(pi2, 40, {eta2, prime})}); assert(r.main == Main::parent && r.selected.pdg == 331);
  assert(run({piece(pi, 50)}).main == Main::mixed);
  r = run({piece(pi, 50, {omega}), piece(gamma, 10, {omega})}); assert(r.main == Main::parent);
  r = run({piece(eta, 80, {}, true)}); assert(r.main == Main::parent_unresolved && r.selected.pdg == 221);
  assert(run({piece(pi, 60), piece(eta, 40, {}, true)}).main == Main::pi0);
  assert(run({piece(pi, 20), piece(eta, 40, {}, true)}).main == Main::unknown);
  auto negative_id = node(111, -999); negative_id.identity.source = 2;
  assert(run({piece(negative_id, 60)}).main == Main::pi0);
  auto other_pi = pi; other_pi.identity.embedding = 0;
  assert(run({piece(pi, 40), piece(other_pi, 40)}).main == Main::mixed);
  assert(run({piece(other_pi, 60)}).main == Main::non_signal);
  auto gamma2 = node(22, 20); gamma2.prompt = 1;
  assert(run({piece(prompt, 30), piece(gamma2, 30)}).main == Main::mixed);
  assert(run({piece(gamma, 70)}).main == Main::photon_unknown);
  pi.decay = 2; assert(run({piece(pi, 80)}).main == Main::pi0);
  r = run({piece(pi, 30, {eta, eta})}); assert(r.main == Main::unknown && (r.quality & ancestry_cycle));
  std::cout << "Common composition: specified examples, unit priority, deepest parent, identities, strict boundary and incomplete information passed\n";
}
