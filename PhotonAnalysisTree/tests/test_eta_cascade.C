#include "../src/EtaCascadeTruth.h"
#include "../src/EtaOriginFractions.h"
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
#include <phhepmc/PHHepMCGenEventMap.h>
#include <phhepmc/PHHepMCGenEvent.h>
#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>
#include <HepMC/GenParticle.h>
#include <TTree.h>
#include <TFile.h>
#include <TMemFile.h>
#undef NDEBUG
#include <cassert>
#include <iostream>
R__LOAD_LIBRARY(libPhotonAnalysisTree.so)
R__LOAD_LIBRARY(libg4detectors.so)
namespace eta_cascade_test
{
PHG4Particle* particle(PHG4TruthInfoContainer& truth, int track, int pdg, int parent, int barcode = 0)
{
  auto* p = new PHG4Particlev2("test", pdg, 10, 0, 0);
  p->set_track_id(track); p->set_parent_id(parent); p->set_vtx_id(-1); p->set_barcode(barcode); p->set_e(10.1);
  truth.AddParticle(track, p); if (track > 0) truth.AddEmbededTrkId(track, 1); return p;
}
photon_tree::ClusterTruthMatch eta_match()
{
  photon_tree::ClusterTruthMatch m; m.valid = true; m.total_edep = 10;
  photon_tree::TruthContributor eta; eta.g4_track_id = 1; eta.g4_pdg_id = 221; eta.embedding_id = 1;
  eta.edep = 8; eta.fraction = .8; m.contributors.push_back(eta);
  photon_tree::TruthContributor other; other.g4_track_id = 2; other.g4_pdg_id = 211; other.embedding_id = 1;
  other.edep = 2; other.fraction = .2; m.contributors.push_back(other); return m;
}
void g4_case(bool charged)
{
  PHG4TruthInfoContainer truth;
  truth.AddVertex(-1, new PHG4VtxPointv2(0, 0, 0, 0, -1, PHG4MCProcess::kPDecay));
  particle(truth, 1, 221, 0, 100); particle(truth, 2, 211, 0);
  for (int i = 1; i <= 3; ++i) particle(truth, -i, charged && i < 3 ? (i == 1 ? 211 : -211) : 111, 1);
  // Deposits have descendant track IDs, not the pion's own track ID.
  for (int i = 1; i <= 3; ++i) particle(truth, -i-10, 22, -i);
  photon_tree::EtaCascadeTruth index(&truth, nullptr);
  assert(index.is_eta_daughter_pi0(truth.GetParticle(-3)));
  RawClusterv1 cluster; const auto key = RawTowerDefs::encode_towerid(RawTowerDefs::CEMC, 0, 0); cluster.addTower(key, 5);
  TowerInfoContainerv1 towers(TowerInfoContainer::EMCAL); towers.get_tower_at_key(TowerInfoDefs::encode_emcal(0, 0))->set_energy(10);
  RawTowerContainer raw(RawTowerDefs::CEMC); auto* tower = new RawTowerv1(key); tower->add_ecell(1, 10); raw.AddTower(key, tower);
  PHG4CellContainer cells; auto* cell = new PHG4Cellv1(1); cells.AddCell(cell);
  PHG4HitContainer hits;
  for (int i = 1; i <= 3; ++i)
  {
    auto* hit = new PHG4Hitv1; hit->set_trkid(-i-10); const auto inserted = hits.AddHit(1, hit);
    cell->add_edep(inserted->first, i == 3 ? 4 : 3);
  }
  auto m = eta_match(); index.refine(m, &cluster, &towers, &raw, &cells, &hits);
  assert(m.contributors.size() == 4 && m.total_edep == 10);
  double sum = 0, eta_sum = 0, pi0_sum = 0;
  for (const auto& c : m.contributors)
  {
    sum += c.fraction;
    if (c.eta_parent_track_id == 1)
    {
      assert(c.eta_split_status == 1 && c.eta_decay_mode == (charged ? 3 : 2));
      eta_sum += c.fraction; if (c.g4_pdg_id == 111) pi0_sum += c.fraction;
      assert(c.fraction < .5); // .8 * (30%,30%,40%) never creates a false majority.
    }
  }
  assert(std::abs(sum - 1) < 1e-6 && std::abs(eta_sum - .8) < 1e-6);
  assert(std::abs(pi0_sum - (charged ? .32 : .8)) < 1e-6);
  photon_tree::ClusterTruthCollection collection; collection.clear(); collection.append(m);
  assert(std::abs(collection.eta_main_fraction[0] - .8) < 1e-6);
  assert(collection.eta_main_decay_mode[0] == (charged ? 3 : 2) && !collection.eta_main_unresolved[0]);
  assert(std::abs(collection.pi0_total_fraction[0] - pi0_sum) < 1e-6);
  TMemFile file("eta_cascade_roundtrip.root", "RECREATE"); TTree tree("event_tree", "Cascade branch roundtrip");
  collection.create_branches(&tree, "split"); tree.Fill(); tree.Write();
  auto* saved = file.Get<TTree>("event_tree"); std::vector<int>* modes = nullptr; std::vector<float>* eta_fractions = nullptr;
  saved->SetBranchAddress("split_cluster_truth_contributor_eta_decay_mode", &modes);
  saved->SetBranchAddress("split_cluster_truth_eta_main_fraction", &eta_fractions); saved->GetEntry(0);
  assert(modes && modes->size() == 4 && eta_fractions && std::abs(eta_fractions->at(0) - .8) < 1e-6);
  saved->ResetBranchAddresses();

  // Missing hit provenance must retain eta and be marked unresolved, not fabricate daughter shares.
  m = eta_match(); index.refine(m, &cluster, &towers, &raw, nullptr, &hits);
  assert(m.contributors.size() == 2 && m.contributors[0].eta_split_status == -1 && m.contributors[0].g4_pdg_id == 221);
  truth.GetVtx(-1)->set_process(PHG4MCProcess::kPPrimary);
  photon_tree::EtaCascadeTruth interaction(&truth, nullptr);
  assert(!interaction.is_eta_daughter_pi0(truth.GetParticle(-3)));
}
void hepmc_case(bool charged)
{
  auto* event = new HepMC::GenEvent;
  auto* decay = new HepMC::GenVertex; event->add_vertex(decay);
  auto* eta = new HepMC::GenParticle(HepMC::FourVector(30, 0, 0, 31), 221, 2); decay->add_particle_in(eta);
  PHG4TruthInfoContainer truth; photon_tree::ClusterTruthMatch m; m.valid = true; m.total_edep = 10;
  for (int i = 1; i <= 3; ++i)
  {
    const int pdg = charged && i < 3 ? (i == 1 ? 211 : -211) : 111;
    auto* pion = new HepMC::GenParticle(HepMC::FourVector(10, 0, 0, 10.1), pdg, 2); decay->add_particle_out(pion);
    HepMC::GenParticle* primary = pion;
    if (pdg == 111)
    {
      auto* vertex = new HepMC::GenVertex; event->add_vertex(vertex); vertex->add_particle_in(pion);
      auto* gamma = new HepMC::GenParticle(HepMC::FourVector(5, 0, 0, 5), 22, 1); vertex->add_particle_out(gamma);
      vertex->add_particle_out(new HepMC::GenParticle(HepMC::FourVector(5, 0, 0, 5), 22, 1)); primary = gamma;
    }
    particle(truth, i, primary->pdg_id(), 0, primary->barcode());
    photon_tree::TruthContributor c; c.g4_track_id = i; c.g4_pdg_id = primary->pdg_id(); c.embedding_id = 1;
    c.edep = i == 3 ? 4 : 3; c.fraction = c.edep / 10; m.contributors.push_back(c);
  }
  PHHepMCGenEventMap events; events.insert_event(1)->addEvent(event);
  photon_tree::EtaCascadeTruth(&truth, &events).refine(m, nullptr, nullptr, nullptr, nullptr, nullptr);
  for (const auto& c : m.contributors)
  {
    assert(c.eta_parent_barcode == eta->barcode() && c.eta_decay_mode == (charged ? 3 : 2));
    if (c.g4_pdg_id == 22) assert(c.pi0_parent_barcode != 0);
  }
}
void fractions()
{
  std::vector<int> embed = {1, 1, 1}, barcode = {100, 100, 100}, track = {0, 0, 0}, mode = {2, 2, 2}, status = {1, 1, 1};
  std::vector<int> pion = {11, 12, 13}, zero = {0, 0, 0}; std::vector<float> f = {.3, .3, .4};
  auto r = eta_origin::summarize(0, 3, embed, f, barcode, track, mode, status, pion, zero);
  assert(std::abs(r.fraction - 1) < 1e-6 && std::abs(r.pi0_total_fraction - 1) < 1e-6 && r.mode == 2 && !r.unresolved);
  barcode = {100, 101, 102}; r = eta_origin::summarize(0, 3, embed, f, barcode, track, mode, status, pion, zero);
  assert(r.fraction < .5); // No sum across distinct eta parents.
  barcode = {100, 100, 100}; embed[2] = 0;
  r = eta_origin::summarize(0, 3, embed, f, barcode, track, mode, status, pion, zero); assert(std::abs(r.fraction - .6) < 1e-6);
}
}
void test_eta_cascade()
{
  eta_cascade_test::g4_case(false); eta_cascade_test::g4_case(true);
  eta_cascade_test::hepmc_case(false); eta_cascade_test::hepmc_case(true); eta_cascade_test::fractions();
  std::cout << "Eta cascade G4/HepMC ancestry, measured daughter shares, closure, unresolved provenance and parent grouping: passed" << std::endl;
}
