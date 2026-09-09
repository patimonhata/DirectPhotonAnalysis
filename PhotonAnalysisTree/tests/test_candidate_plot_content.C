#include "../workflows/photon_candidate_selection/MergePythiaPhotonCandidateSelection.C"
#include <cassert>

void test_candidate_plot_content(const std::string preview = "")
{
  using namespace candidate_composition;
  plot_caption = {"jet", "Region A + Tagging veto", 0.2, 0.5, 0.7, 0.1, 0.2, 0.45, 0.65};
  const auto lines = caption_lines(plot_caption);
  assert(lines[1] == "Pythia 8 p+p Jet samples");
  assert(lines[4] == "Region A + Tagging veto");
  assert(lines[5] == "Stored/anchor clusters: E > 0.20 GeV");
  assert(lines[6].find("0.50") != std::string::npos && lines[7].find("0.70") != std::string::npos);
  Histograms composition(2, 10.0);
  Spectra topology(2, 10.0);
  const std::array<int, category_count> counts = {0, 40, 15, 10, 5, 10, 5, 10, 5};
  for (std::size_t i = 1; i < category_count; ++i)
    for (int j = 0; j < counts[i]; ++j) { composition.fill(i, 7.0, 2.0); composition.fill(denominator, 7.0, 2.0); }
  const std::array<int, 6> missing_counts = {1, 2, 3, 1, 1, 2};
  for (std::size_t i = 0; i < missing_counts.size(); ++i)
    for (int j = 0; j < missing_counts[i]; ++j) { topology.fill(5, 7.0, 2.0); topology.fill(kMissingSpectrumIndices[i], 7.0, 2.0); }
  auto missing = composition_missing_fractions(composition, topology);
  assert(missing.size() == 6);
  double sum = 0.0;
  for (std::size_t i = 0; i < missing.size(); ++i)
  {
    assert(same_double(missing[i]->GetBinContent(2), missing_counts[i] / 100.0));
    sum += missing[i]->GetBinContent(2);
  }
  assert(same_double(sum, 0.1)); // Denominator is all selected candidates, not only pi0 anchors.
  topology.counts[5]->Fill(7.0);
  assert(composition_missing_fractions(composition, topology).empty());
  if (!preview.empty())
  {
    std::array<std::unique_ptr<TH1D>, category_count> fractions;
    for (std::size_t i = 1; i < category_count; ++i)
      fractions[i] = fraction_histogram(*composition.weighted[i], *composition.weighted[denominator], std::string("qa_fraction_") + candidate_composition::kKeys[i]);
    set_missing_labels(0.2, 0.5);
    draw_superdetailed_stack(fractions, missing, preview);
  }
  std::cout << "plot captions, missing normalization and population guard: passed" << std::endl;
}
