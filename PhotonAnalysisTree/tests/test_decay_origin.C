#include "../src/DecayOriginDiagnostics.h"
#include <calobase/RawClusterv1.h>
#include <calobase/RawTowerv1.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/TowerInfoContainerv1.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoDefs.h>
#include <g4detectors/PHG4Cellv1.h>
#include <g4detectors/PHG4CellContainer.h>
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particlev2.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPointv2.h>
#include <g4main/PHG4MCProcessDefs.h>
#include <phhepmc/PHHepMCGenEventMap.h>
#include <phhepmc/PHHepMCGenEvent.h>
#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>
#include <HepMC/GenParticle.h>
#include <TTree.h>
#include <TMemFile.h>
#include <algorithm>
#include <cmath>
#include <map>
#undef NDEBUG
#include <cassert>
R__LOAD_LIBRARY(libPhotonAnalysisTree.so)
R__LOAD_LIBRARY(libg4detectors.so)
namespace decay_origin_test
{
using namespace photon_tree;
PHG4Particle* particle(PHG4TruthInfoContainer& truth, int id, int pdg, int parent, int vertex, int barcode = 0)
{
  auto* p = new PHG4Particlev2("test", pdg, 10, 0, 0); p->set_e(10.1);
  p->set_track_id(id); p->set_parent_id(parent); p->set_vtx_id(vertex); p->set_barcode(barcode);
  truth.AddParticle(id, p); if (parent == 0) truth.AddEmbededTrkId(id, 1); return p;
}
TruthContributor contribution(int id, int pdg, float fraction)
{
  TruthContributor c; c.g4_track_id = id; c.g4_pdg_id = pdg; c.embedding_id = 1; c.fraction = fraction; c.edep = fraction * 10; return c;
}
struct Deposits
{
  RawClusterv1 cluster;
  TowerInfoContainerv1 towers{TowerInfoContainer::EMCAL};
  RawTowerContainer raw{RawTowerDefs::CEMC};
  PHG4CellContainer cells;
  PHG4HitContainer hits;
  PHG4Cellv1* cell;
  Deposits()
  {
    const auto key = RawTowerDefs::encode_towerid(RawTowerDefs::CEMC, 0, 0); cluster.addTower(key, 5);
    towers.get_tower_at_key(TowerInfoDefs::encode_emcal(0, 0))->set_energy(10); // 50% split-tower allocation.
    auto* tower = new RawTowerv1(key); tower->add_ecell(1, 10); raw.AddTower(key, tower);
    cell = new PHG4Cellv1(1); cells.AddCell(cell);
  }
  void add(int track, double energy)
  {
    auto* hit = new PHG4Hitv1; hit->set_trkid(track); const auto it = hits.AddHit(1, hit); cell->add_edep(it->first, energy);
  }
};
std::map<int, double> fractions(const DecayOriginMatch& match)
{
  assert(match.valid); std::map<int, double> result; double sum = 0;
  for (const auto& c : match.contributors) { result[c.pdg] += c.fraction; sum += c.fraction; }
  assert(std::abs(sum - 1) < 1e-6); return result;
}
DecayOriginMatch g4_case(bool radiative, bool has_hepmc = false, double omega_flight = 1e-4, bool incomplete = false, bool interaction = false)
{
  PHG4TruthInfoContainer truth; PHHepMCGenEventMap events;
  truth.AddVertex(1, new PHG4VtxPointv2(0, 0, 0, 0, 1, PHG4MCProcess::kPPrimary));
  truth.AddVertex(-1, new PHG4VtxPointv2(omega_flight, 0, 0, 0, -1, interaction ? PHG4MCProcess::kPPrimary : PHG4MCProcess::kPDecay));
  truth.AddVertex(-2, new PHG4VtxPointv2(omega_flight + 1e-6, 0, 0, 0, -2, PHG4MCProcess::kPDecay));
  int barcode = 0;
  if (has_hepmc)
  {
    auto* event = new HepMC::GenEvent; auto* birth = new HepMC::GenVertex; event->add_vertex(birth);
    auto* omega = new HepMC::GenParticle(HepMC::FourVector(10, 0, 0, 10.1), 223, 1); birth->add_particle_out(omega);
    barcode = omega->barcode(); events.insert_event(1)->addEvent(event);
  }
  particle(truth, 1, 223, 0, 1, barcode); particle(truth, 2, 2112, 0, 1);
  particle(truth, -1, 111, 1, -1); particle(truth, -2, radiative ? 22 : 211, 1, -1);
  if (!radiative) particle(truth, -3, -211, 1, -1);
  particle(truth, -4, 22, -1, -2); particle(truth, -5, 22, -1, -2);
  particle(truth, -6, 11, -4, -2); // A shower descendant, not an additional physical photon.
  Deposits d; d.add(-6, radiative ? 1 : 5); d.add(-5, radiative ? 1 : 3); d.add(-2, radiative ? 8 : 1); if (!radiative) d.add(-3, 1);
  ClusterTruthMatch input; input.valid = true; input.total_edep = 10;
  input.contributors = {contribution(1, 223, .8), contribution(2, 2112, .2)};
  const auto original = input.contributors;
  const auto match = DecayOriginDiagnostics(&truth, &events).diagnose(input, &d.cluster, &d.towers, &d.raw, incomplete ? nullptr : &d.cells, &d.hits);
  assert(input.contributors.size() == 2 && input.contributors[0].g4_pdg_id == 223 && input.contributors[0].fraction == original[0].fraction);
  const auto by_pdg = fractions(match);
  if (interaction) { assert(match.contributors[0].pdg == 223 && match.contributors[0].status == -2); return match; }
  if (incomplete) { assert(match.contributors[0].pdg == 223 && match.contributors[0].status == -1); return match; }
  if (omega_flight > DecayOriginDiagnostics::kMaxFlightCm)
  { assert(match.contributors[0].pdg == 223 && match.contributors[0].long_lived); return match; }
  assert(std::abs(by_pdg.at(111) - (radiative ? .16 : .64)) < 1e-6 && !by_pdg.count(223));
  for (const auto& c : match.contributors)
  {
    if (c.root_track_id != 1) continue;
    assert(c.status == 1 && c.parent_pdg == 223 && c.pathway == 2 && std::abs(c.parent_flight_cm - omega_flight) < 1e-7);
    if (c.pdg == 111) assert(c.diphoton && std::abs(c.gamma0_fraction + c.gamma1_fraction - c.fraction) < 1e-6);
    if (c.pdg == 22) assert(std::abs(c.fraction - .64) < 1e-6);
  }
  return match;
}
DecayOriginMatch hepmc_case(bool radiative, bool mixed_boundary = false, int ancestor_pdg = 223, double ancestor_flight_mm = .001)
{
  PHG4TruthInfoContainer truth; PHHepMCGenEventMap events; Deposits d;
  auto* event = new HepMC::GenEvent(HepMC::Units::GEV, HepMC::Units::MM);
  auto* birth = new HepMC::GenVertex; auto* decay = new HepMC::GenVertex(HepMC::FourVector(ancestor_flight_mm, 0, 0, 0));
  auto* pi_decay = new HepMC::GenVertex(HepMC::FourVector(ancestor_flight_mm + .00001, 0, 0, 0));
  event->add_vertex(birth); event->add_vertex(decay); if (!mixed_boundary) event->add_vertex(pi_decay); else delete pi_decay;
  auto* omega = new HepMC::GenParticle(HepMC::FourVector(10, 0, 0, 10.1), ancestor_pdg, 2); birth->add_particle_out(omega); decay->add_particle_in(omega);
  auto* pi0 = new HepMC::GenParticle(HepMC::FourVector(8, 0, 0, 8.1), 111, mixed_boundary ? 1 : 2); decay->add_particle_out(pi0);
  truth.AddVertex(1, new PHG4VtxPointv2(ancestor_flight_mm / 10, 0, 0, 0, 1, PHG4MCProcess::kPPrimary));
  truth.AddVertex(-1, new PHG4VtxPointv2(ancestor_flight_mm / 10 + 1e-6, 0, 0, 0, -1, PHG4MCProcess::kPDecay));
  ClusterTruthMatch input; input.valid = true; input.total_edep = 10;
  int track = 1;
  if (mixed_boundary)
  {
    particle(truth, track, 111, 0, 1, pi0->barcode()); input.contributors.push_back(contribution(track++, 111, radiative ? .16 : .64));
    particle(truth, -1, 22, 1, -1); particle(truth, -2, 22, 1, -1); d.add(-1, 5); d.add(-2, 3);
  }
  else
  {
    pi_decay->add_particle_in(pi0);
    for (int i = 0; i < 2; ++i)
    {
      auto* gamma = new HepMC::GenParticle(HepMC::FourVector(4, 0, 0, 4), 22, 1); pi_decay->add_particle_out(gamma);
      particle(truth, track, 22, 0, 1, gamma->barcode()); input.contributors.push_back(contribution(track++, 22, radiative ? .08 : i == 0 ? .4 : .24));
    }
  }
  for (int pdg : ancestor_pdg == 310 ? std::vector<int>{111} : radiative ? std::vector<int>{22} : std::vector<int>{211, -211})
  {
    auto* child = new HepMC::GenParticle(HepMC::FourVector(1, 0, 0, 1.1), pdg, 1); decay->add_particle_out(child);
    particle(truth, track, pdg, 0, 1, child->barcode()); input.contributors.push_back(contribution(track++, pdg, ancestor_pdg == 310 ? .16 : radiative ? .64 : .08));
  }
  particle(truth, track, 2112, 0, 1); input.contributors.push_back(contribution(track, 2112, .2));
  events.insert_event(1)->addEvent(event);
  const auto result = DecayOriginDiagnostics(&truth, &events).diagnose(input, &d.cluster, &d.towers, &d.raw, &d.cells, &d.hits);
  const auto by_pdg = fractions(result);
  if (ancestor_pdg == 310)
  { assert(std::abs(by_pdg.at(310) - .8) < 1e-6 && result.contributors[0].long_lived && !by_pdg.count(111)); return result; }
  assert(std::abs(by_pdg.at(111) - (radiative ? .16 : .64)) < 1e-6);
  for (const auto& c : result.contributors) if (c.pdg == 111)
  { assert(c.diphoton && c.parent_pdg == 223 && std::abs(c.gamma0_fraction + c.gamma1_fraction - c.fraction) < 1e-6);
    assert(std::abs(c.parent_flight_cm - ancestor_flight_mm / 10) < 1e-6); }
  return result;
}
void ambiguous_identity()
{
  PHG4TruthInfoContainer truth; PHHepMCGenEventMap events;
  auto* event = new HepMC::GenEvent; auto* birth = new HepMC::GenVertex; event->add_vertex(birth);
  auto* gamma = new HepMC::GenParticle(HepMC::FourVector(10, 0, 0, 10), 22, 1); birth->add_particle_out(gamma);
  particle(truth, 1, 22, 0, 1, gamma->barcode()); particle(truth, 2, 22, 0, 1, gamma->barcode());
  events.insert_event(1)->addEvent(event);
  ClusterTruthMatch input; input.valid = true; input.total_edep = 10;
  input.contributors = {contribution(1, 22, .6), contribution(2, 22, .4)};
  const auto result = DecayOriginDiagnostics(&truth, &events).diagnose(input, nullptr, nullptr, nullptr, nullptr, nullptr);
  assert(result.valid && result.contributors.size() == 2);
  for (const auto& c : result.contributors) assert(c.source == 2 && c.status == -1);
}
DecayOriginMatch weak_g4()
{
  PHG4TruthInfoContainer truth;
  truth.AddVertex(1, new PHG4VtxPointv2(0, 0, 0, 0, 1, PHG4MCProcess::kPPrimary));
  truth.AddVertex(-1, new PHG4VtxPointv2(1, 0, 0, 0, -1, PHG4MCProcess::kPDecay));
  particle(truth, 1, 310, 0, 1); particle(truth, 2, 2112, 0, 1);
  particle(truth, -1, 111, 1, -1); particle(truth, -2, 111, 1, -1);
  ClusterTruthMatch input; input.valid = true; input.total_edep = 10; input.contributors = {contribution(1, 310, .8), contribution(2, 2112, .2)};
  const auto result = DecayOriginDiagnostics(&truth, nullptr).diagnose(input, nullptr, nullptr, nullptr, nullptr, nullptr);
  assert(result.valid && result.contributors[0].pdg == 310 && result.contributors[0].long_lived && result.contributors[0].status == 0);
  assert(result.contributors[0].daughter_count == 2 && result.contributors[0].daughter_pi0_count == 2);
  return result;
}
}
void test_decay_origin()
{
  using namespace decay_origin_test;
  for (bool radiative : {false, true})
  {
    const auto g4 = fractions(g4_case(radiative)), hep = fractions(hepmc_case(radiative)), mixed = fractions(hepmc_case(radiative, true));
    for (const auto& p : g4) { assert(std::abs(p.second - hep.at(p.first)) < 1e-6); assert(std::abs(p.second - mixed.at(p.first)) < 1e-6); }
    g4_case(radiative, true); // Canonical HepMC identity at the simulation boundary.
  }
  g4_case(false, false, 1.); g4_case(false, false, 1e-4, true); g4_case(false, false, 1e-4, false, true);
  hepmc_case(false, false, 310, 10.); // Weak parent remains a barrier even if daughters entered G4 as photons.
  weak_g4(); ambiguous_identity();
  const auto match = g4_case(false);
  photon_tree::DecayOriginCollection collection; collection.clear(); collection.append(match);
  TMemFile file("decay_origin_roundtrip.root", "RECREATE"); TTree tree("event_tree", "Diagnostic roundtrip"); collection.create_branches(&tree); tree.Fill(); tree.Write();
  auto* saved = file.Get<TTree>("event_tree"); std::vector<int>* pdg = nullptr; std::vector<float>* fraction = nullptr;
  saved->SetBranchAddress("split_cluster_decay_origin_pdg", &pdg); saved->SetBranchAddress("split_cluster_decay_origin_fraction", &fraction); saved->GetEntry(0);
  assert(pdg && pdg->size() == 4 && fraction && std::abs(fraction->at(0) - .64) < 1e-6); saved->ResetBranchAddresses();
  std::cout << "Unified decay diagnostics: HepMC/G4/mixed omega decays, measured shares, closure, displacement, interaction rejection and branch roundtrip passed" << std::endl;
}
