#include "../analysis/photon_candidate_selection/CandidateMerge.h"
#include "test_decay_origin.C"
namespace
{
decay_origin::Result classify_diagnostic(const photon_tree::DecayOriginMatch& match)
{
  photon_tree::DecayOriginCollection collection; collection.clear(); collection.append(match);
  TMemFile file("diagnostic_classification.root", "RECREATE"); TTree tree("event_tree", "Classification and closure");
  collection.create_branches(&tree); tree.Fill(); tree.Write();
  auto* saved = file.Get<TTree>("event_tree"); decay_origin::Reader reader; assert(reader.bind(saved)); assert(saved->GetEntry(0) > 0);
  assert(reader.consistent(1, {0, 1}, {1})); const auto result = decay_origin::classify(reader, 0);
  (*reader.fraction)[0] += .1; assert(!reader.consistent(1, {0, 1}, {1}));
  saved->ResetBranchAddresses(); return result;
}
}
void test_decay_origin_reduce()
{
  using namespace decay_origin;
  static_assert(kMaxFlightCm == photon_tree::DecayOriginDiagnostics::kMaxFlightCm);
  static_assert(kMaxProperLengthCm == photon_tree::DecayOriginDiagnostics::kMaxProperLengthCm);
  assert(classify_diagnostic(decay_origin_test::g4_case(false)).category == pi0_gg);
  assert(classify_diagnostic(decay_origin_test::hepmc_case(false)).category == pi0_gg);
  assert(classify_diagnostic(decay_origin_test::hepmc_case(false, true)).category == pi0_gg);
  assert(classify_diagnostic(decay_origin_test::g4_case(true)).category == decay_gamma);
  assert(classify_diagnostic(decay_origin_test::hepmc_case(true)).category == decay_gamma);
  assert(classify_diagnostic(decay_origin_test::weak_g4()).category == nonprompt);
  assert(classify_diagnostic(decay_origin_test::hepmc_case(false, false, 310, 10.)).category == nonprompt);
  assert(classify_diagnostic(decay_origin_test::g4_case(false, false, 1e-4, true)).category == unresolved);
  auto half = decay_origin_test::g4_case(false); half.contributors[0].fraction = .5; half.contributors[1].fraction += .14;
  half.contributors[0].gamma0_fraction = half.contributors[0].gamma1_fraction = .25;
  assert(classify_diagnostic(half).category == mixed);
  Histograms histograms(4, 20); candidate_composition::Histograms original(4, 20);
  for (std::size_t i = 0; i < category_count; ++i)
  {
    Result r; r.category = static_cast<Category>(i); r.fraction = .8;
    const double et = i == 0 ? -1 : i == 1 ? 21 : 7, weight = i == 3 ? -2 : 3;
    histograms.fill(r, et, weight, .9, true); original.fill(candidate_composition::other, et, weight);
  }
  assert(histograms.valid(*original.counts[candidate_composition::other], *original.weighted[candidate_composition::other]));
  TMemFile file("diagnostic_histogram.root", "RECREATE"); auto* d = file.mkdir("decay_origin"); assert(histograms.write(*d));
  Histograms restored(4, 20); assert(restored.add(*d));
  assert(restored.valid(*original.counts[candidate_composition::other], *original.weighted[candidate_composition::other]));
  std::cout << "Decay-origin classification, strict majority, branch closure, weighted/flow closure and histogram roundtrip: passed" << std::endl;
}
