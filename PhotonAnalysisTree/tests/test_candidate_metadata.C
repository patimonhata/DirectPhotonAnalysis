#include "../analysis/photon_candidate_selection/CandidateMerge.h"
#undef NDEBUG
#include <cassert>

void test_candidate_metadata()
{
  MapMetadata map;
  map.schema_version = kMbdMapSchema; map.common_composition_algorithm_version = 1;
  map.mbd = {60., 3., .7, .7, 5., 35., "MbdVertexMap", "signal HepMC collision"};
  map.eta_topology_algorithm_version = 4;
  map.pi0_topology_algorithm_version = 14;
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
  for (int field = 0; field < 4; ++field)
  {
    auto policy = map;
    if (field == 0) ++policy.decay_origin_algorithm_version;
    if (field == 1) policy.decay_origin_max_proper_length_cm += .001;
    if (field == 2) policy.decay_origin_max_flight_cm += .001;
    if (field == 3) ++policy.common_composition_algorithm_version;
    assert(!compatible(map, policy));
  }
  candidate_composition_merge::Metadata partial;
  partial.source_map_schema_version = kMbdMapSchema;
  partial.mbd = map.mbd;
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
  for (int field = 0; field < 4; ++field)
  {
    auto policy = partial;
    if (field == 0) ++policy.decay_origin_algorithm_version;
    if (field == 1) policy.decay_origin_max_proper_length_cm += .001;
    if (field == 2) policy.decay_origin_max_flight_cm += .001;
    if (field == 3) ++policy.common_composition_algorithm_version;
    assert(!candidate_composition_merge::compatible(partial, policy));
  }
  for (const auto field : {&MbdMapSettings::max_abs_mbd_vertex_z, &MbdMapSettings::shower_shape_min_cluster_et,
       &MbdMapSettings::candidate_abs_eta_max, &MbdMapSettings::topology_anchor_abs_eta_max, &MbdMapSettings::candidate_et_min, &MbdMapSettings::candidate_et_max})
  {
    auto different_map = map;
    auto different_partial = partial;
    different_map.mbd.*field += .1;
    different_partial.mbd.*field += .1;
    assert(!compatible(map, different_map));
    assert(!candidate_composition_merge::compatible(partial, different_partial));
  }

  auto different_eta = partial;
  --different_eta.eta_topology_algorithm_version;
  assert(!candidate_composition_merge::compatible(partial, different_eta));
  auto different_other = partial;
  --different_other.other_origin_algorithm_version;
  assert(!candidate_composition_merge::compatible(partial, different_other));
  std::cout << "map and partial metadata compatibility: passed" << std::endl;
}
