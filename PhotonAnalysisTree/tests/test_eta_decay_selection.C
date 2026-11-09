#include "../analysis/photon_candidate_selection/EtaDecayHistograms.h"
#undef NDEBUG
#include <cassert>

void test_eta_decay_selection()
{
  using photon_tree::composition::Main;
  using photon_tree::composition::Identity;
  using eta_decay::Category;
  // The legacy leading parent is unrelated. Only the common selected parent's mode matters.
  std::vector<int> embedding = {1, 1, 1}, barcode = {999, 100, 100}, track = {-9, -7, -7}, mode = {2, 3, 3};
  const auto classify = [&](Main main, Identity selected = {1, 1, 100}, int pdg = 221) {
    return eta_decay::classify(main, pdg, selected, 0, mode.size(), embedding, barcode, track, mode);
  };
  assert(classify(Main::parent) == Category::charged_pions_pi0);
  assert(classify(Main::parent, {1, 2, -7}) == Category::charged_pions_pi0);
  assert(classify(Main::parent, {1, 1, 999}) == Category::three_pi0);
  assert(classify(Main::parent, {1, 2, 100}) == Category::none); // Never compare a G4 ID with a barcode.
  assert(classify(Main::parent, {2, 1, 100}) == Category::none);
  assert(classify(Main::parent, {1, 1, 777}) == Category::none);
  assert(classify(Main::parent, {1, 0, 100}) == Category::none);
  assert(classify(Main::parent, {1, 1, 100}, 310) == Category::none);
  for (Main main : {Main::eta_gg, Main::pi0, Main::charged, Main::mixed, Main::unknown}) assert(classify(main) == Category::none);
  assert(classify(Main::parent_unresolved) == Category::unresolved);
  mode[2] = 2; assert(classify(Main::parent) == Category::none); // Conflicting annotations.
  mode[2] = 0; assert(classify(Main::parent) == Category::none); // Partially unknown mode.
  mode[1] = -1; mode[2] = 2; assert(classify(Main::parent) == Category::none);
  mode[1] = mode[2] = 4; assert(classify(Main::parent) == Category::none); // Other decay.
  mode[1] = mode[2] = 2; assert(classify(Main::parent) == Category::three_pi0);
  // The observed direct-gg plus unrelated cascade case must never add a cascade entry.
  eta_decay::Histograms modes(10, 10);
  TH1D all("eta_all_test", "", 10, 0, 10), all_w("eta_all_weight_test", "", 10, 0, 10);
  TH1D direct("eta_direct_test", "", 10, 0, 10), direct_w("eta_direct_weight_test", "", 10, 0, 10);
  for (auto* h : {&all, &all_w, &direct, &direct_w}) { h->SetDirectory(nullptr); h->Sumw2(); }
  for (double et : {-1., 5., 11.}) // Include underflow/overflow and signed weights.
  {
    for (const auto& entry : std::vector<std::pair<Main, double>>{{Main::eta_gg, 2.}, {Main::parent, -3.}, {Main::parent_unresolved, 4.}})
    {
      all.Fill(et); all_w.Fill(et, entry.second);
      if (entry.first == Main::eta_gg) { direct.Fill(et); direct_w.Fill(et, entry.second); }
      modes.fill(classify(entry.first), et, entry.second);
    }
  }
  assert(modes.valid(all, all_w, direct, direct_w));
  assert(modes.counts[0]->Integral(0, 11) == 3 && modes.counts[1]->Integral(0, 11) == 0 && modes.counts[2]->Integral(0, 11) == 3);
  // Keep the production integrity check: deliberate double counting must still fail.
  modes.fill(Category::three_pi0, 5., 2.);
  assert(!modes.valid(all, all_w, direct, direct_w));
  std::cout << "Eta selected-parent identity, exclusive subsets, unknown modes and count/Sumw2 guard: passed" << std::endl;
}
