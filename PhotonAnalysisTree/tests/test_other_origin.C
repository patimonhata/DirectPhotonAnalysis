#include "../workflows/photon_candidate_selection/MergePythiaPhotonCandidateSelection.C"
#include <TMemFile.h>
#undef NDEBUG
#include <cassert>

void test_other_origin()
{
  using namespace other_origin;
  const auto one = [](int pdg, int parent = 0, float fraction = .8, int embedding = 1) {
    return classify(true, 0, 1, {pdg}, {embedding}, {fraction}, {parent}, {pdg == 22 ? 3 : -1});
  };
  for (int pdg : {211, -211, 321, -321, 2212, -2212}) assert(one(pdg).category == charged_hadron);
  assert(one(11).category == electron && one(-11).category == electron && one(-13).category == muon);
  for (int parent : {223, 331, 113, 553, 441, 3212, -3212}) assert(one(22, parent).category == decay_photon);
  for (int parent : {111, 221, 21, 90, 2101, -999}) assert(one(22, parent).category == other_particle);
  assert(classify(true, 0, 1, {22}, {1}, {.8}, {223}, {0}).category == other_particle);
  assert(one(223).category == other_particle); // No invented G4 omega-daughter split from map-only information.
  assert(one(211, 0, .5).category == mixed && one(22, 223, .5).category == mixed);
  assert(one(211, 0, .9, 0).category == unresolved);
  assert(classify(false, 0, 1, {211}, {1}, {.8}, {0}, {-1}).category == unresolved);
  assert(classify(true, 0, 1, {211}, {1}, {.8}, {0}, {-1}, true).category == unresolved);
  auto contaminated = classify(true, 0, 3, {22, 111, 22}, {1, 1, 1}, {.6, .3, .1}, {223, 0, 111}, {3, -1, 3});
  assert(contaminated.category == decay_photon && std::abs(contaminated.main_fraction - .6) < 1e-6);
  assert(std::abs(contaminated.other_gamma_fraction - .1) < 1e-6);
  assert(classify(true, 0, 3, {22, 111, 211}, {1, 1, 1}, {.45, .35, .2}, {223, 0, 0}, {3, -1, -1}).category == mixed);

  Histograms histograms(4, 20);
  candidate_composition::Histograms composition(4, 20);
  std::vector<Result> cases = {one(211), one(11), one(13), contaminated, one(311), one(211, 0, .5), one(211, 0, .8, 0)};
  for (std::size_t i = 0; i < cases.size(); ++i)
  {
    const double weight = i == 4 ? -1 : 2;
    histograms.fill(cases[i], 7, weight, .9, true);
    composition.fill(candidate_composition::other, 7, weight); composition.fill(candidate_composition::denominator, 7, weight);
  }
  assert(histograms.valid(*composition.counts[candidate_composition::other], *composition.weighted[candidate_composition::other]));
  for (bool detailed : {false, true}) for (bool all : {false, true})
  {
    auto parts = components(composition, histograms, detailed, all);
    double n = 0, w = 0, w2 = 0;
    for (std::size_t i = 0; i < parts.counts.size(); ++i)
    { n += parts.counts[i]->Integral(); w += parts.weighted[i]->Integral(); w2 += std::pow(parts.weighted[i]->GetBinError(2), 2); }
    assert(n == 7 && w == 11 && std::abs(w2 - 25) < 1e-8);
    assert(parts.counts.size() == (detailed ? 7 : 3) + (all ? 3 : 0));
  }
  TMemFile file("other_origin_roundtrip.root", "RECREATE"); auto* directory = file.mkdir("other_origin");
  assert(histograms.write(*directory)); Histograms restored(4, 20); assert(restored.add(*directory));
  assert(restored.valid(*composition.counts[candidate_composition::other], *composition.weighted[candidate_composition::other]));
  assert(restored.purity[decay_photon]->Integral() == 2 && restored.decay_parent->Integral() == 2);
  restored.counts[decay_photon]->Fill(7);
  assert(!restored.valid(*composition.counts[candidate_composition::other], *composition.weighted[candidate_composition::other]));
  std::cout << "Other-origin categories, strict majority, hadron parents, contamination, weighted closure and ROOT roundtrip: passed" << std::endl;
}
