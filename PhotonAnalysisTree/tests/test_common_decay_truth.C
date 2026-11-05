#include "test_decay_origin.C"
#include "../src/CommonDecayTruth.h"
#include "../src/Pi0ClusterTruthMatcher.h"
#include "../analysis/photon_candidate_selection/CommonCompositionHistograms.h"
namespace common_decay_test
{
using namespace photon_tree;
using decay_origin_test::particle;
using decay_origin_test::contribution;
using decay_origin_test::Deposits;
composition::Result g4_case(bool radiative, bool has_hepmc = false, double omega_flight = 1e-4, bool incomplete = false, bool interaction = false)
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
  const auto r = CommonDecayTruth(&truth, &events).resolve(input, &d.cluster, &d.towers, &d.raw, incomplete ? nullptr : &d.cells, &d.hits);
  assert(r.valid && std::abs(r.total - r.attributed) < 1e-8);
  if (interaction) { assert(r.quality & composition::interaction); assert(r.main == composition::Main::particle); return r; }
  if (incomplete) { assert(r.main == composition::Main::parent_unresolved && r.selected.pdg == 223); assert(r.photons.size() == (radiative ? 3U : 2U)); for (const auto& p : r.photons) assert(p.unresolved && p.deposit == -1); return r; }
  assert(r.main == (radiative ? composition::Main::decay_photon : composition::Main::pi0));
  assert(std::abs(r.fraction - .64) < 1e-6);
  assert(CommonDecayTruth(&truth, &events).diphotons(111, 1).size() == 1);
  assert(r.photons.size() == (radiative ? 3U : 2U)); // Includes the direct omega photon, never its shower electrons.
  return r;
}
composition::Result hepmc_case(bool radiative, bool mixed_boundary = false, int ancestor_pdg = 223, double ancestor_flight_mm = .001)
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
  const auto r = CommonDecayTruth(&truth, &events).resolve(input, &d.cluster, &d.towers, &d.raw, &d.cells, &d.hits);
  assert(r.valid && std::abs(r.total - r.attributed) < 1e-5);
  assert(r.main == (radiative ? composition::Main::decay_photon : composition::Main::pi0));
  assert(std::abs(r.fraction - .64) < 1e-6);
  return r;
}

void charged_priority()
{
  PHG4TruthInfoContainer truth;
  truth.AddVertex(1, new PHG4VtxPointv2(0, 0, 0, 0, 1, PHG4MCProcess::kPPrimary));
  truth.AddVertex(-1, new PHG4VtxPointv2(10, 0, 0, 0, -1, PHG4MCProcess::kPDecay));
  particle(truth, 1, 221, 0, 1); particle(truth, -1, 211, 1, -1); particle(truth, -2, -211, 1, -1); particle(truth, -3, 111, 1, -1);
  Deposits d; d.add(-1, 7); d.add(-2, 1); d.add(-3, 2);
  ClusterTruthMatch m; m.valid = true; m.total_edep = 10; m.contributors = {contribution(1, 221, 1)};
  auto r = CommonDecayTruth(&truth, nullptr).resolve(m, &d.cluster, &d.towers, &d.raw, &d.cells, &d.hits);
  assert(r.main == composition::Main::charged && std::abs(r.fraction - .7) < 1e-6);
  assert(r.selected.identity.source == 2 && r.selected.identity.id == -1 && r.parents[0].particle.pdg == 221);
  // Unknown flight never invalidates known ancestry or measured shares.
  truth.GetVtx(1)->set_x(NAN);
  auto no_flight = CommonDecayTruth(&truth, nullptr).resolve(m, &d.cluster, &d.towers, &d.raw, &d.cells, &d.hits);
  assert(no_flight.main == r.main && no_flight.fraction == r.fraction);
}
void weak_parent()
{
  PHG4TruthInfoContainer truth;
  truth.AddVertex(1, new PHG4VtxPointv2(0, 0, 0, 0, 1, PHG4MCProcess::kPPrimary));
  truth.AddVertex(-1, new PHG4VtxPointv2(100, 0, 0, 0, -1, PHG4MCProcess::kPDecay));
  particle(truth, 1, 310, 0, 1); particle(truth, -1, 111, 1, -1); particle(truth, -2, 111, 1, -1); particle(truth, 2, 2112, 0, 1);
  Deposits d; d.add(-1, 5); d.add(-2, 5);
  ClusterTruthMatch m; m.valid = true; m.total_edep = 10; m.contributors = {contribution(1, 310, .8), contribution(2, 2112, .2)};
  const auto r = CommonDecayTruth(&truth, nullptr).resolve(m, &d.cluster, &d.towers, &d.raw, &d.cells, &d.hits);
  assert(r.main == composition::Main::parent && r.selected.pdg == 310 && r.selected.flight_cm == 100);
}
void non_gg_and_duplicate()
{
  PHG4TruthInfoContainer truth; truth.AddVertex(-1, new PHG4VtxPointv2(0, 0, 0, 0, -1, PHG4MCProcess::kPDecay));
  particle(truth, 1, 111, 0, -1); particle(truth, -1, 22, 1, -1); particle(truth, -2, 11, 1, -1); particle(truth, -3, -11, 1, -1);
  ClusterTruthMatch m; m.valid = true; m.total_edep = 10; m.contributors = {contribution(1, 111, 1)};
  auto r = CommonDecayTruth(&truth, nullptr).resolve(m, nullptr, nullptr, nullptr, nullptr, nullptr);
  assert(r.main == composition::Main::pi0 && r.applicability == composition::Applicability::not_applicable);
  m.contributors.push_back(m.contributors[0]); r = CommonDecayTruth(&truth, nullptr).resolve(m, nullptr, nullptr, nullptr, nullptr, nullptr);
  assert(!r.valid && r.main == composition::Main::unknown && (r.quality & composition::ambiguous_identity));
}
void known_majority_with_unresolved_sibling()
{
  PHG4TruthInfoContainer truth; Deposits d;
  truth.AddVertex(1,new PHG4VtxPointv2(0,0,0,0,1,PHG4MCProcess::kPPrimary));
  truth.AddVertex(-1,new PHG4VtxPointv2(0,0,0,0,-1,PHG4MCProcess::kPDecay));
  particle(truth,1,331,0,1); particle(truth,-1,111,1,-1); particle(truth,-2,221,1,-1); particle(truth,2,2112,0,1);
  d.add(-1,7.5); d.add(-2,2.5); // eta has a known identity/deposit but no recorded daughter mode.
  ClusterTruthMatch m; m.valid = true; m.total_edep = 10; m.contributors = {contribution(1,331,.8),contribution(2,2112,.2)};
  const auto r = CommonDecayTruth(&truth,nullptr).resolve(m,&d.cluster,&d.towers,&d.raw,&d.cells,&d.hits);
  assert(r.main == composition::Main::pi0 && std::abs(r.fraction-.6) < 1e-6);
  assert(std::abs(r.attributed-10) < 1e-6 && r.quality & composition::incomplete_split);
  bool unresolved_eta = false; for (const auto& p : r.units) unresolved_eta |= p.particle.pdg == 221 && p.unresolved && std::abs(p.deposit-2) < 1e-6;
  assert(unresolved_eta);
}

void unknown_eta_mode()
{
  PHG4TruthInfoContainer truth; particle(truth,1,221,0,0);
  ClusterTruthMatch m; m.valid = true; m.total_edep = 10; m.contributors = {contribution(1,221,1)};
  const auto r = CommonDecayTruth(&truth,nullptr).resolve(m,nullptr,nullptr,nullptr,nullptr,nullptr);
  assert(r.main == composition::Main::parent_unresolved && r.selected.pdg == 221);
  assert(r.applicability == composition::Applicability::unknown && r.topology_reason == 2);
}

void photon_copies_and_partial_matching()
{
  PHG4TruthInfoContainer truth; PHHepMCGenEventMap events; Deposits d;
  truth.AddVertex(1, new PHG4VtxPointv2(0,0,0,0,1,PHG4MCProcess::kPPrimary));
  auto* event = new HepMC::GenEvent; auto* decay = new HepMC::GenVertex; event->add_vertex(decay);
  auto* pi0 = new HepMC::GenParticle(HepMC::FourVector(10,0,0,10.1),111,2); decay->add_particle_in(pi0);
  for (int i = 0; i < 2; ++i)
  {
    auto* gamma = new HepMC::GenParticle(HepMC::FourVector(5,0,0,5),22,2); decay->add_particle_out(gamma);
    auto* copy_vertex = new HepMC::GenVertex; event->add_vertex(copy_vertex); copy_vertex->add_particle_in(gamma);
    auto* copy = new HepMC::GenParticle(HepMC::FourVector(5,0,0,5),22,1); copy_vertex->add_particle_out(copy);
    particle(truth,i+1,22,0,1,copy->barcode()); d.add(i+1,i == 0 ? 8 : 2);
  }
  events.insert_event(1)->addEvent(event);
  ClusterTruthMatch m; m.valid = true; m.total_edep = 10; m.contributors = {contribution(1,22,.8),contribution(2,22,.2)};
  auto r = CommonDecayTruth(&truth,&events).resolve(m,&d.cluster,&d.towers,&d.raw,&d.cells,&d.hits);
  assert(r.main == composition::Main::pi0 && r.photons.size() == 2 && std::abs(r.max_photon - 8) < 1e-5);
  const auto candidates = CommonDecayTruth(&truth,&events).diphotons(111,1); assert(candidates.size() == 1);
  const std::array<int,2> tracks = {1,2};
  auto complete = Pi0ClusterTruthMatcher().match(&d.cluster,&d.towers,&d.raw,&d.cells,&d.hits,&truth,tracks,true);
  assert(complete.status == Pi0ClusterTruthMatchStatus::complete && complete.ancestry_complete);
  // Half the member energy is artificial/unmatched noise: coverage is .5, not photon purity.
  auto key = RawTowerDefs::encode_towerid(RawTowerDefs::CEMC,0,1); d.cluster.addTower(key,5);
  d.towers.get_tower_at_key(TowerInfoDefs::encode_emcal(0,1))->set_energy(10);
  auto partial = Pi0ClusterTruthMatcher().match(&d.cluster,&d.towers,&d.raw,&d.cells,&d.hits,&truth,tracks,true);
  assert(partial.status == Pi0ClusterTruthMatchStatus::partial && partial.cluster_member_energy_coverage == .5);
  assert(partial.gamma_edep == complete.gamma_edep && partial.total_edep == complete.total_edep);
  r = CommonDecayTruth(&truth,&events).resolve(m,&d.cluster,&d.towers,&d.raw,&d.cells,&d.hits);
  assert(r.main == composition::Main::pi0 && r.total == 10); // Noise adds no truth particle or truth denominator.
  d.add(-12345,1); auto missing = Pi0ClusterTruthMatcher().match(&d.cluster,&d.towers,&d.raw,&d.cells,&d.hits,&truth,tracks,true);
  assert(!missing.ancestry_complete);
}

}
void test_common_decay_truth()
{
  using namespace common_decay_test;
  for (bool radiative : {false, true})
  {
    const auto g4 = g4_case(radiative), hep = hepmc_case(radiative), mixed = hepmc_case(radiative, true);
    assert(g4.main == hep.main && hep.main == mixed.main);
    assert(std::abs(g4.fraction - hep.fraction) < 1e-6 && std::abs(g4.fraction - mixed.fraction) < 1e-6);
    std::map<int, double> a, b, c;
    for (const auto& p : g4.units) a[p.particle.pdg] += p.deposit;
    for (const auto& p : hep.units) b[p.particle.pdg] += p.deposit;
    for (const auto& p : mixed.units) c[p.particle.pdg] += p.deposit;
    for (const auto& p : a) assert(std::abs(p.second - b[p.first]) < 1e-5 && std::abs(p.second - c[p.first]) < 1e-5);
    g4_case(radiative, true); g4_case(radiative, false, 100); hepmc_case(radiative, false, 223, 1000);
  }
  g4_case(false, false, 1e-4, true); g4_case(false, false, 1e-4, false, true);
  charged_priority(); weak_parent(); non_gg_and_duplicate(); unknown_eta_mode(); known_majority_with_unresolved_sibling(); photon_copies_and_partial_matching();
  auto r = g4_case(false);
  r.matching_status = 1; r.minimum_coverage = .1;
  r.anchor_matching_status = 2; r.anchor_matching_coverage = 1; r.anchor_matching_usable = true;
  CommonTruthCollection collection; collection.clear(); collection.append(r, 3);
  TMemFile file("common_roundtrip.root", "RECREATE"); TTree tree("event_tree", "Common truth"); collection.create_branches(&tree); tree.Fill();
  collection.clear(); collection.append({}); tree.Fill(); tree.Write();
  auto* saved = file.Get<TTree>("event_tree"); common_composition::Reader reader; assert(reader.bind(saved)); assert(saved->GetEntry(0) > 0); assert(reader.valid(1));
  assert((*reader.matching_status)[0] == 1 && (*reader.minimum_coverage)[0] == .1);
  assert((*reader.anchor_matching_status)[0] == 2 && (*reader.anchor_matching_coverage)[0] == 1 && (*reader.anchor_matching_usable)[0]);
  common_composition::Histograms h(4, 20);
  for (double et : {-1., 7., 21.}) for (double weight : {-2., 3.}) h.fill(reader, 0, 8, et, weight);
  assert(h.valid()); auto* d = file.mkdir("histograms"); h.write(*d); common_composition::Histograms restored(4,20); assert(restored.add(*d) && restored.valid());
  assert(restored.matching_count->GetEntries() == 6 && restored.anchor_matching_count->GetEntries() == 6);
  assert(restored.matching_count->GetAxis(1)->GetNbins() == 3);
  assert(saved->GetEntry(1) > 0 && reader.valid(1));
  assert((*reader.anchor_matching_status)[0] == -1 && (*reader.anchor_matching_coverage)[0] == -1 && !(*reader.anchor_matching_usable)[0]);
  const auto entries = h.anchor_matching_count->GetEntries(); h.fill(reader, 0, 23, 7, 1);
  assert(h.anchor_matching_count->GetEntries() == entries);
  saved->ResetBranchAddresses();
  std::cout << "Common truth: HepMC/G4/mixed conservation, long-lived and displaced parents, charged priority, secondary mesons, non-gg, missing hits, duplicates, flows/Sumw2 and branch roundtrip passed" << std::endl;
}
