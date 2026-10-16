#ifndef ETA_TOPOLOGY_HISTOGRAMS_H
#define ETA_TOPOLOGY_HISTOGRAMS_H

// Included after the existing spectra and composition helpers. Existing pi0 and composition populations remain unchanged.
namespace eta_topology
{
inline const char* spectrum_key(std::size_t index) { return index == 1 ? "eta_anchor" : ::kKeys[index]; }
inline int bdt_category(int topology)
{
  if (topology == 1) return 0;
  if (topology == 2) return 1;
  if (topology == 4) return 2;
  if (topology == 3) return 3;
  return topology == 0 ? 4 : -1;
}
inline bool fill(Spectra& spectra, int topology, int missing, double et, double weight)
{
  constexpr std::array<std::size_t, 5> indices = {2, 3, 4, 5, 12};
  const int category = bdt_category(topology);
  if (category < 0) return false;
  if (topology == 3 && (missing < 1 || missing > 6)) return false;
  spectra.fill(1, et, weight);
  spectra.fill(indices[category], et, weight);
  if (topology == 3)
  {
    constexpr std::array<std::size_t, 7> missing_indices = {0, 6, 8, 11, 7, 9, 10};
    spectra.fill(missing_indices[missing], et, weight);
  }
  return true;
}
inline bool valid(const candidate_composition::Histograms& composition, const Spectra& spectra, const pi0_bdt::Histograms& bdt)
{
  if (!::valid_partition(spectra.counts) || !::valid_partition(spectra.weighted_pb) || !bdt.valid()) return false;
  // Reuse the existing BDT/count/weighted-integral validation with the eta topology population.
  candidate_composition::Histograms selected(spectra.counts[1]->GetNbinsX(), spectra.counts[1]->GetXaxis()->GetXmax());
  constexpr std::array<std::size_t, 5> indices = {2, 3, 4, 5, 12};
  for (std::size_t category = 0; category < indices.size(); ++category)
  {
    selected.counts[candidate_composition::pi0_separated + category]->Add(spectra.counts[indices[category]].get());
    selected.weighted[candidate_composition::pi0_separated + category]->Add(spectra.weighted_pb[indices[category]].get());
  }
  if (!candidate_composition::consistent_pi0_bdt(selected, spectra, bdt)) return false;
  const auto& eta = *composition.counts[candidate_composition::eta];
  const auto& eta_weighted = *composition.weighted[candidate_composition::eta];
  for (int bin = 0; bin <= eta.GetNbinsX() + 1; ++bin)
  {
    if (spectra.counts[1]->GetBinContent(bin) > eta.GetBinContent(bin)) return false;
    const double total_w2 = std::pow(eta_weighted.GetBinError(bin), 2);
    const double anchor_w2 = std::pow(spectra.weighted_pb[1]->GetBinError(bin), 2);
    if (anchor_w2 > total_w2 + 1e-8 * std::max(1.0, total_w2)) return false;
  }
  return true;
}
inline bool write_partial(TFile& output, const std::array<std::unique_ptr<Spectra>, kSelectionCount>& spectra,
                          const std::array<pi0_bdt::Histograms, kSelectionCount>& bdt)
{
  auto* topology_root = output.mkdir("eta_anchor_topology");
  auto* bdt_root = output.mkdir("eta_bdt");
  if (!topology_root || !bdt_root) return false;
  for (std::size_t selection = 0; selection < kSelectionCount; ++selection)
  {
    auto* directory = topology_root->mkdir(kSelectionKeys[selection]);
    auto* scores = bdt_root->mkdir(kSelectionKeys[selection]);
    if (!directory || !scores || !bdt[selection].write(*scores)) return false;
    directory->cd();
    for (std::size_t index = 0; index < kSpectrumCount; ++index)
    {
      spectra[selection]->counts[index]->Write((std::string("h_region_a_") + spectrum_key(index) + "_et_count").c_str());
      spectra[selection]->weighted_pb[index]->Write((std::string("h_region_a_") + spectrum_key(index) + "_et_pb").c_str());
    }
  }
  return true;
}
}
#endif
