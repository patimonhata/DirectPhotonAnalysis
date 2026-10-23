#include "../workflows/diagnostics/eta_decay_flight/CheckEtaDecayFlight.C"
#include <g4main/PHG4Particlev2.h>
#include <g4main/PHG4VtxPointv1.h>
#include <g4main/PHG4VtxPointv2.h>
#undef NDEBUG
#include <cassert>

namespace eta_flight_test
{
PHG4Particle* particle(PHG4TruthInfoContainer& truth, int track, int pdg, int parent, int vertex)
{
  auto* p = new PHG4Particlev2("test", pdg, 10, 0, 0);
  p->set_e(10.1); p->set_track_id(track); p->set_parent_id(parent); p->set_vtx_id(vertex); p->set_barcode(track);
  truth.AddParticle(track, p); if (track > 0) truth.AddEmbededTrkId(track, 1); return p;
}
void vertex(PHG4TruthInfoContainer& truth, int id, double x, double y, double z, double t, int process)
{
  truth.AddVertex(id, new PHG4VtxPointv2(x, y, z, t, id, static_cast<PHG4MCProcess>(process)));
}
void base(PHG4TruthInfoContainer& truth)
{
  vertex(truth, 1, 12, -4, 37, 2, PHG4MCProcess::kPPrimary);
  particle(truth, 1, 221, 0, 1);
}
void daughters(PHG4TruthInfoContainer& truth, int vtx)
{
  for (int i = 1; i <= 3; ++i) particle(truth, -i, 111, 1, vtx);
}
eta_decay_flight::Record analyze(PHG4TruthInfoContainer& truth)
{
  return eta_decay_flight::analyze(*truth.GetParticle(1), truth, eta_decay_flight::children_by_parent(truth));
}
void cases()
{
  using namespace eta_decay_flight;
  {
    PHG4TruthInfoContainer truth; base(truth); vertex(truth, -1, 12.003, -3.996, 37, 2.01, PHG4MCProcess::kPDecay); daughters(truth, -1);
    const auto r = analyze(truth);
    assert(r.status == 0 && r.usable && r.mode == 2 && r.n_direct_daughters == 3);
    assert(std::abs(r.length_cm - .005) < 1e-12 && std::abs(r.transverse_cm - .005) < 1e-12);
    assert(std::abs(r.parallel_cm - .003) < 1e-12 && std::abs(r.perpendicular_cm - .004) < 1e-12);
    assert(std::abs(r.delta_t - .01) < 1e-12); // Translation of the event vertex does not change the displacement.
  }
  {
    PHG4TruthInfoContainer truth; base(truth); vertex(truth, -1, 12, -4, 37, 2, PHG4MCProcess::kPDecay); daughters(truth, -1);
    auto r = analyze(truth); assert(r.usable && r.zero_displacement && !r.same_vertex_id);
    truth.GetVtx(1)->set_process(PHG4MCProcess::kPDecay);
    for (int i = 1; i <= 3; ++i) truth.GetParticle(-i)->set_vtx_id(1);
    r = analyze(truth); assert(r.usable && r.same_vertex_id && r.zero_displacement);
  }
  {
    PHG4TruthInfoContainer truth; base(truth); daughters(truth, -1);
    assert(analyze(truth).status == 4 && !analyze(truth).distance_valid);
  }
  {
    PHG4TruthInfoContainer truth; particle(truth, 1, 221, 0, 999);
    assert(analyze(truth).status == 2);
  }
  {
    PHG4TruthInfoContainer truth; base(truth); assert(analyze(truth).status == 3);
  }
  {
    PHG4TruthInfoContainer truth; base(truth);
    vertex(truth, -1, 12.01, -4, 37, 2.1, PHG4MCProcess::kPDecay); vertex(truth, -2, 12.02, -4, 37, 2.2, PHG4MCProcess::kPDecay);
    particle(truth, -1, 22, 1, -1); particle(truth, -2, 22, 1, -2);
    assert(analyze(truth).status == 5 && !analyze(truth).distance_valid);
  }
  {
    PHG4TruthInfoContainer truth; base(truth); vertex(truth, -1, 12.1, -4, 37, 2.1, PHG4MCProcess::kPHadronic); daughters(truth, -1);
    assert(analyze(truth).status == 6 && !analyze(truth).usable);
  }
  {
    PHG4TruthInfoContainer truth; base(truth); truth.AddVertex(-1, new PHG4VtxPointv1(12.1, -4, 37, 2.1, -1)); daughters(truth, -1);
    auto r = analyze(truth); assert(r.status == 1 && r.usable && r.mode == 2);
    truth.GetParticle(-1)->set_pid(321);
    assert(analyze(truth).status == 7); // Unknown process and an unrecognized final state are not called a decay.
  }
  {
    PHG4TruthInfoContainer truth; base(truth); vertex(truth, -1, eta_decay_flight::nan, -4, 37, 2.1, PHG4MCProcess::kPDecay); daughters(truth, -1);
    assert(analyze(truth).status == 8);
  }
  {
    PHG4TruthInfoContainer truth; base(truth); vertex(truth, -1, 12.1, -4, 37, 1, PHG4MCProcess::kPDecay); daughters(truth, -1);
    auto r = analyze(truth); assert(r.status == 9 && r.distance_valid && !r.usable);
  }
  {
    PHG4TruthInfoContainer truth; base(truth); vertex(truth, -1, 12.1, -4, 37, 2.1, PHG4MCProcess::kPDecay); daughters(truth, -1);
    vertex(truth, -2, 12.05, -4, 37, 2.05, PHG4MCProcess::kPHadronic); particle(truth, -4, 2112, 1, -2);
    auto r = analyze(truth); assert(r.status == 10 && r.distance_valid && !r.usable);
  }
  {
    PHG4TruthInfoContainer truth; base(truth);
    vertex(truth, -1, 12.1, -4, 37, 2.1, PHG4MCProcess::kPDecay); vertex(truth, -2, 12.1, -4, 37, 2.1, PHG4MCProcess::kPDecay);
    particle(truth, -1, 22, 1, -1); particle(truth, -2, 22, 1, -2);
    auto r = analyze(truth); assert(r.usable && r.mode == 1 && r.n_vertex_groups == 1);
  }
  {
    PHG4TruthInfoContainer truth; base(truth); vertex(truth, -1, 12.1, -4, 37, eta_decay_flight::nan, PHG4MCProcess::kPDecay); daughters(truth, -1);
    auto r = analyze(truth); assert(r.usable && !r.time_valid);
  }
  {
    PHG4TruthInfoContainer truth; base(truth);
    vertex(truth, -1, 12.1, -4, 37, 2.1, PHG4MCProcess::kPDecay); vertex(truth, -2, 12.2, -4, 37, 2.2, PHG4MCProcess::kPDecay);
    truth.GetVtx(-1)->set_id(0); truth.GetVtx(-2)->set_id(0); // Container keys, not stale object IDs, identify distinct vertices.
    particle(truth, -1, 22, 1, -1); particle(truth, -2, 22, 1, -2);
    assert(analyze(truth).status == 5);
  }
  Summary summary;
  Record zero; zero.usable = true; zero.status = 0; zero.length_cm = 0; zero.transverse_cm = 0; zero.zero_displacement = true; zero.weight = 1;
  Record boundary = zero; boundary.length_cm = 1e-4; boundary.zero_displacement = false; boundary.weight = 3;
  Record missing; missing.status = 3;
  summary.add(zero); summary.add(boundary); summary.add(missing);
  assert(summary.total == 3 && summary.usable == 2 && summary.zero == 1);
  assert(fraction_below(summary, 1e-4) == .5 && fraction_below(summary, 1e-4, true) == .25);
  assert(std::abs(fraction_statistics(summary, 1e-4).error - std::sqrt(.125)) < 1e-12);
  assert(std::abs(fraction_statistics(summary, 1e-4, true).error - std::sqrt(.0703125)) < 1e-12);
  assert(fraction_below(summary, std::nextafter(1e-4, 1.)) == 1);
  assert(std::isnan(fraction_below(Summary{}, 1e-4)));
  Record invalid_weight = zero; invalid_weight.weight_valid = false; invalid_weight.weight = -1;
  auto invalid_summary = summary; invalid_summary.add(invalid_weight);
  assert(invalid_summary.usable == 3 && invalid_summary.invalid_weight == 1);
  assert(fraction_below(invalid_summary, 1e-4) == 2. / 3 && fraction_below(invalid_summary, 1e-4, true) == .25);
  assert(decay_mode({111, -211, 211}) == 3 && decay_mode({22, 211, -211}) == 4);
  std::cout << "Eta flight geometry, decay/process states, zeros, missing data, strict CDF and weights: passed" << std::endl;
}

void fixture(const std::string& path)
{
  TFile output(path.c_str(), "RECREATE"); TTree tree("T", "Synthetic G4 truth");
  auto* truth = new PHG4TruthInfoContainer(); tree.Branch("DST#G4TruthInfo", &truth);
  base(*truth); vertex(*truth, -1, 12.001, -4, 37, 2.01, PHG4MCProcess::kPDecay); daughters(*truth, -1); tree.Fill();
  truth->Reset(); base(*truth); tree.Fill(); // A missing decay must remain in the census.
  truth->Reset(); base(*truth); vertex(*truth, -1, 12, -4, 37, 2, PHG4MCProcess::kPDecay); daughters(*truth, -1); tree.Fill();
  tree.Write(); tree.ResetBranchAddresses(); delete truth;
}
}

void test_eta_decay_flight(const std::string output_directory = "")
{
  eta_flight_test::cases();
  const std::string directory = output_directory.empty() ? std::string(gSystem->TempDirectory()) + "/eta_flight_test_" + std::to_string(gSystem->GetPid()) : output_directory;
  gSystem->mkdir(directory.c_str(), true);
  const std::string dst = directory + "/fixture-000123.root", list = directory + "/inputs.list", base = directory + "/inclusive";
  eta_flight_test::fixture(dst);
  { std::ofstream out(list); out << dst << '\n'; }
  assert(CheckEtaDecayFlight(list, base, 0, 1, 0) == 0);
  {
    TFile result((base + ".root").c_str()); auto* tree = result.Get<TTree>("eta_tree");
    assert(tree && tree->GetEntries() == 3);
    assert(tree->GetEntries("usable") == 2 && tree->GetEntries("status==3") == 1 && tree->GetEntries("zero_displacement") == 1);
    auto* metadata = result.Get<TTree>("metadata"); assert(metadata && metadata->GetEntries("completed") == 1);
    assert(result.Get<TH1D>("all/h_log10_length_cm")->GetEntries() == 1);
    assert(result.Get<TTree>("event_tree")->GetEntries("n_unresolved_eta>0") == 1);
  }
  assert(CheckEtaDecayFlight(list, base, 0, 1, 0) != 0); // No overwrite.
  const std::string requests = directory + "/requests.csv";
  { std::ofstream out(requests); out << "source_file_id,event_in_file,eta_track_id\n123,0,1\n123,0,1\n123,2,1\n"; }
  assert(CheckEtaDecayFlight(list, directory + "/requested", 0, 0, 0, requests) == 0);
  { TFile result((directory + "/requested.root").c_str()); assert(result.Get<TTree>("eta_tree")->GetEntries() == 2); }
  { std::ofstream out(requests); out << "source_file_id,event_in_file,eta_track_id\n123,99,1\n"; }
  assert(CheckEtaDecayFlight(list, directory + "/invalid_request", 0, 0, 0, requests) != 0);
  { TFile result((directory + "/invalid_request.root").c_str()); assert(!result.Get<TTree>("metadata")); }
  std::cout << "Eta flight inclusive/requested ROOT pipeline, deduplication, event counts and overwrite rejection: passed" << std::endl;
}
