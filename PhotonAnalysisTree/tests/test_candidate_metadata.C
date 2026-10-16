#include "../workflows/photon_candidate_selection/MergePythiaPhotonCandidateSelection.C"
#undef NDEBUG
#include <cassert>

void test_candidate_metadata()
{
  MapMetadata map;
  map.schema_version = 6;
  map.eta_topology_algorithm_version = 1;
  map.pi0_topology_algorithm_version = 11;
  map.min_photon_energy_recovery = 0.0;
  assert(compatible(map, map));
  const std::array<double MapMetadata::*, 10> map_fields = {
      &MapMetadata::min_cluster_energy, &MapMetadata::pi0_partner_min_energy, &MapMetadata::eta_partner_min_energy,
      &MapMetadata::pi0_mass_min, &MapMetadata::pi0_mass_max, &MapMetadata::eta_mass_min, &MapMetadata::eta_mass_max,
      &MapMetadata::missing_energy_min, &MapMetadata::missing_energy_max, &MapMetadata::min_photon_energy_recovery};
  for (auto field : map_fields)
  {
    auto changed = map;
    changed.*field += 0.01;
    assert(!compatible(map, changed));
  }
  auto changed = map;
  --changed.pi0_topology_algorithm_version;
  assert(!compatible(map, changed));
  changed = map;
  --changed.eta_topology_algorithm_version;
  assert(!compatible(map, changed));
  candidate_composition_merge::Metadata partial;
  assert(candidate_composition_merge::compatible(partial, partial));
  const std::array<double candidate_composition_merge::Metadata::*, 10> partial_fields = {
      &candidate_composition_merge::Metadata::min_cluster_energy, &candidate_composition_merge::Metadata::pi0_partner_min_energy,
      &candidate_composition_merge::Metadata::eta_partner_min_energy, &candidate_composition_merge::Metadata::pi0_mass_min,
      &candidate_composition_merge::Metadata::pi0_mass_max, &candidate_composition_merge::Metadata::eta_mass_min,
      &candidate_composition_merge::Metadata::eta_mass_max, &candidate_composition_merge::Metadata::missing_energy_min,
      &candidate_composition_merge::Metadata::missing_energy_max, &candidate_composition_merge::Metadata::min_photon_energy_recovery};
  for (auto field : partial_fields)
  {
    auto different = partial;
    different.*field += 0.01;
    assert(!candidate_composition_merge::compatible(partial, different));
  }
  auto different_eta = partial;
  --different_eta.eta_topology_algorithm_version;
  assert(!candidate_composition_merge::compatible(partial, different_eta));
  std::cout << "map and partial metadata compatibility: passed" << std::endl;
}
