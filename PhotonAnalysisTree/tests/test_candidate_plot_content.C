#include "../workflows/photon_candidate_selection/MergePythiaPhotonCandidateSelection.C"
#include <cassert>

void test_candidate_plot_content(const std::string preview = "")
{
  using namespace candidate_composition;
  plot_caption = {"jet", "Region A + Tagging veto", 0.2, 0.5, 0.7, 0.1, 0.2, 0.45, 0.65};
  const auto lines = caption_lines(plot_caption);
  assert(lines.size() == 5);
  assert(lines[1] == "Pythia 8 p+p Jet samples, |z_{vertex}^{truth}| < 60 cm");
  assert(lines[2] == "Candidate Cluster: after Region A + Tagging veto");
  assert(lines[3].find("0.50") != std::string::npos && lines[4].find("0.70") != std::string::npos);
  plot_caption.selection = kSelectionLabels.front();
  const auto kinematic_lines = caption_lines(plot_caption);
  assert(kinematic_lines.size() == 5);
  assert(kinematic_lines[2] == "Candidate Cluster: 5 < E_{T} < 35 GeV, |#eta| < 0.7");
  assert(std::none_of(kinematic_lines.begin(), kinematic_lines.end(), [](const std::string& line) { return line.find("after") != std::string::npos; }));
  assert(std::none_of(lines.begin(), lines.end(), [](const std::string& line) { return line.find("Stored/anchor") != std::string::npos; }));
  Histograms composition(2, 10.0);
  Spectra topology(2, 10.0);
  const std::array<int, category_count> counts = {0, 40, 15, 10, 5, 10, 5, 10, 5};
  for (std::size_t i = 1; i < category_count; ++i)
    for (int j = 0; j < counts[i]; ++j) { composition.fill(i, 7.0, 2.0); composition.fill(denominator, 7.0, 2.0); }
  const auto origin_density = candidate_composition_merge::make_candidate_origin_density(composition, "qa_candidate_origin_");
  assert(same_double(origin_density[0]->GetBinContent(2), 16.0));
  assert(same_double(origin_density[1]->GetBinContent(2), 18.0));
  assert(same_double(origin_density[2]->GetBinContent(2), 4.0));
  assert(same_double(origin_density[3]->GetBinContent(2), 2.0));
  assert(std::string(origin_density[0]->GetXaxis()->GetTitle()) == "Candidate Cluster E_{T} [GeV]");
  assert(std::string(origin_density[0]->GetYaxis()->GetTitle()) == "Weighted Counts [a.u.]");
  const auto plot_frame = make_plot_frame(*origin_density[0], "qa_plot_frame", "Weighted Counts [a.u.]", 0.0, 1.0);
  assert(same_double(plot_frame->GetXaxis()->GetXmin(), 0.0));
  assert(same_double(plot_frame->GetXaxis()->GetXmax(), 40.0));
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
