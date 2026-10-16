#include "../workflows/photon_candidate_selection/MergePythiaPhotonCandidateSelection.C"
#undef NDEBUG
#include <cassert>

void test_eta_topology(const std::string output_base = "")
{
  using namespace candidate_composition;
  Histograms composition(40, 40);
  Spectra topology(40, 40);
  pi0_bdt::Histograms bdt;
  for (int value : {0, 1, 2, 3, 4})
  {
    const int missing = value == 3 ? 1 : 0;
    composition.fill(eta, 7, 2); composition.fill(denominator, 7, 2);
    assert(eta_topology::fill(topology, value, missing, 7, 2));
    assert(bdt.fill(eta_topology::bdt_category(value), .9, true, true, false, 7, 2));
  }
  for (int missing = 2; missing <= 6; ++missing)
  {
    composition.fill(eta, 7, 3); composition.fill(denominator, 7, 3);
    assert(eta_topology::fill(topology, 3, missing, 7, 3));
    assert(bdt.fill(3, NAN, false, false, true, 7, 3));
  }
  // Multiple eta parents / non-two-photon decays retain their legacy composition category.
  composition.fill(eta, 7, 4); composition.fill(denominator, 7, 4);
  assert(eta_topology::valid(composition, topology, bdt));
  const double before = topology.counts[1]->Integral();
  assert(!eta_topology::fill(topology, 99, 0, 7, 1));
  assert(!eta_topology::fill(topology, 3, 0, 7, 1));
  assert(topology.counts[1]->Integral() == before);
  auto components = eta_topology::composition_components(composition, topology, false);
  auto detailed = eta_topology::composition_components(composition, topology, true);
  assert(components.weighted.size() == 9 && detailed.weighted.size() == 14);
  for (const auto* set : {&components, &detailed})
  {
    double count = 0, weight = 0, variance = 0;
    for (std::size_t i = 0; i < set->weighted.size(); ++i)
    {
      count += set->counts[i]->Integral(); weight += set->weighted[i]->Integral();
      variance += std::pow(set->weighted[i]->GetBinError(8), 2);
    }
    assert(count == composition.counts[denominator]->Integral());
    assert(weight == composition.weighted[denominator]->Integral());
    assert(same_double(variance, std::pow(composition.weighted[denominator]->GetBinError(8), 2)));
  }
  assert(components.weighted[7]->Integral() == 4); // Legacy eta residual.
  assert(components.counts[7]->Integral() == 1);
  assert(same_double(std::pow(components.weighted[7]->GetBinError(8), 2), 16));
  topology.counts[1]->Fill(7);
  assert(!eta_topology::valid(composition, topology, bdt));
  if (!output_base.empty())
  {
    std::array<std::unique_ptr<Histograms>, kSelectionCount> compositions;
    std::array<std::unique_ptr<Spectra>, kSelectionCount> spectra;
    std::array<pi0_bdt::Histograms, kSelectionCount> scores;
    for (std::size_t i = 0; i < kSelectionCount; ++i)
    {
      compositions[i] = std::make_unique<Histograms>(40, 40); spectra[i] = std::make_unique<Spectra>(40, 40);
      compositions[i]->fill(eta, 7, 2); compositions[i]->fill(denominator, 7, 2);
      compositions[i]->fill(eta, 7, 3); compositions[i]->fill(denominator, 7, 3);
      assert(eta_topology::fill(*spectra[i], 1, 0, 7, 2));
      assert(scores[i].fill(0, .9, true, true, false, 7, 2));
    }
    candidate_composition_merge::Metadata metadata;
    metadata.family = "jet"; metadata.schema_version = 7; metadata.source_map_schema_version = 6;
    metadata.eta_topology_algorithm_version = 1;
    plot_caption = {"jet", "Kinematic", .2, .5, .5, .1, .2, .45, .65};
    const auto previous_labels = ::kLabels;
    assert(make_output_directory(output_base + "/placeholder"));
    assert(eta_topology::export_outputs(output_base, metadata, {}, compositions, spectra, scores));
    assert(previous_labels == ::kLabels);
    assert(!gSystem->AccessPathName((output_base + "/kinematic/region_a_eta_anchor_topology.pdf").c_str()));
    assert(gSystem->AccessPathName((output_base + "/kinematic/region_a_pi0_anchor_topology.pdf").c_str()));
    TFile output((output_base + "/candidate_composition_eta_topology.root").c_str());
    assert(output.Get<TH1D>("kinematic/detailed/h_candidate_eta_residual_et_pb")->Integral() == 3);
    TFile partial((output_base + "/partial.root").c_str(), "RECREATE");
    assert(eta_topology::write_partial(partial, spectra, scores)); partial.Close();
    TFile readback((output_base + "/partial.root").c_str());
    assert(readback.Get<TH1D>("eta_anchor_topology/kinematic/h_region_a_eta_anchor_et_count")->Integral() == 1);
    pi0_bdt::Histograms restored;
    assert(restored.add(*readback.GetDirectory("eta_bdt/kinematic")) && restored.valid());
  }
  std::cout << "Eta topology partition, missing categories, BDT population and residual covariance: passed" << std::endl;
}
