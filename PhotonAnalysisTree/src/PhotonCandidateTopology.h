#ifndef PHOTON_TREE_PHOTON_CANDIDATE_TOPOLOGY_H
#define PHOTON_TREE_PHOTON_CANDIDATE_TOPOLOGY_H

#include "Pi0AnchorTopologyEvaluator.h"

#include <cmath>
#include <limits>

namespace photon_tree
{
constexpr int kEtaTopologyAlgorithmVersion = 4;

inline Pi0AnchorTopologyConfig photon_candidate_topology_config(Pi0AnchorTopologyConfig config)
{
  // The evaluator includes its energy bound, while the candidate map uses E > min_cluster_energy.
  config.min_cluster_energy = std::nextafter(config.min_cluster_energy, std::numeric_limits<double>::infinity());
  config.enable_missing_diagnostics = true;
  config.eta_pseudo_primary = true;
  config.common_composition = true;
  return config;
}

inline Pi0AnchorTopologyConfig photon_candidate_eta_topology_config(
    Pi0AnchorTopologyConfig config, double partner_min_energy, double mass_min, double mass_max)
{
  config.parent_pdg_id = 221;
  config.require_hepmc_two_photon_decay = true;
  config.tagging_partner_min_cluster_energy = partner_min_energy;
  config.tagging_pi0_mass_min = mass_min;
  config.tagging_pi0_mass_max = mass_max;
  return config;
}

inline void copy_common_topology(composition::Result& target, const composition::Result& evaluated)
{
  target.applicability = evaluated.applicability;
  target.topology_reason = evaluated.topology_reason;
  target.matching_status = evaluated.matching_status;
  target.minimum_coverage = evaluated.minimum_coverage;
  target.anchor_matching_status = evaluated.anchor_matching_status;
  target.anchor_matching_coverage = evaluated.anchor_matching_coverage;
  target.anchor_matching_usable = evaluated.anchor_matching_usable;
}
}

#endif
