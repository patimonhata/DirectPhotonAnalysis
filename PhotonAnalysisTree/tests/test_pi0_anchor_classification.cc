#include "Pi0AnchorTopologyEvaluator.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>

using namespace photon_tree;

int main()
{
  Pi0AnchorTopologyConfig config;
  assert(config.min_photon_energy_recovery == 0.0);
  config.min_photon_energy_recovery = 0.5;
  config.enable_missing_diagnostics = true;
  config.tagging_partner_min_cluster_energy = 0.5;
  Pi0TopologyClusterRecord cluster;
  cluster.cluster_id = 10;
  cluster.energy = 10.0;
  cluster.eta = cluster.phi = 0.0;
  Pi0TopologyCandidateRecord candidate;
  candidate.best_cluster[0] = 0;
  candidate.recovered[0] = true;
  candidate.photon_energy = {10.0, 1.0};
  candidate.photon_projection_valid = {true, true};
  candidate.photon_in_cemc_acceptance = {true, true};
  candidate.photon_cemc_edep[1] = 0.3;
  auto& truth = candidate.truth_partner_clusters[1];
  truth.found = true;
  truth.cluster_id = 20;
  truth.cluster_energy = 0.3;
  truth.direct_edep = 0.1;
  truth.cluster_eta = 0.0;
  const auto set_mass = [&](double mass) { truth.cluster_phi = std::acos(1.0 - mass * mass / (2.0 * cluster.energy * truth.cluster_energy)); };
  const auto classify = [&]() {
    Pi0TopologyAnchorRecord anchor;
    anchor.cluster_index = 0;
    classify_pi0_anchor(config, candidate, cluster, anchor);
    return anchor;
  };
  set_mass(0.135);
  auto anchor = classify();
  assert(anchor.topology == Pi0AnchorTopology::missing);
  assert(anchor.missing_category == Pi0MissingCategory::energy_band_taggable);
  assert(anchor.truth_partner_tag_status == Pi0TruthPartnerTagStatus::below_energy_threshold);
  const double exact_mass = anchor.truth_partner_invariant_mass;
  config.tagging_pi0_mass_min = exact_mass;
  assert(classify().missing_category == Pi0MissingCategory::energy_band_not_taggable);
  config.tagging_pi0_mass_min = 0.10;
  config.tagging_pi0_mass_max = exact_mass;
  assert(classify().missing_category == Pi0MissingCategory::energy_band_not_taggable);
  config.tagging_pi0_mass_max = 0.20;
  set_mass(0.25);
  assert(classify().missing_category == Pi0MissingCategory::energy_band_not_taggable);
  config.tagging_pi0_mass_max = 0.3;
  assert(classify().missing_category == Pi0MissingCategory::energy_band_taggable);
  config.tagging_pi0_mass_max = 0.2;
  for (double energy : {0.05, 0.1, 0.2})
  {
    truth.cluster_energy = energy;
    set_mass(0.135);
    assert(classify().missing_category == Pi0MissingCategory::low_energy);
  }
  truth.cluster_energy = std::nextafter(0.2, 1.0);
  set_mass(0.135);
  assert(classify().missing_category == Pi0MissingCategory::energy_band_taggable);
  truth.cluster_energy = 0.5;
  set_mass(0.135);
  assert(classify().missing_category == Pi0MissingCategory::energy_band_taggable);
  truth.cluster_energy = std::nextafter(0.5, 1.0);
  set_mass(0.135);
  assert(classify().missing_category == Pi0MissingCategory::other);
  truth.cluster_energy = 0.3;
  config.missing_energy_min = 0.35;
  assert(classify().missing_category == Pi0MissingCategory::low_energy);
  config.missing_energy_min = 0.2;
  truth.cluster_phi = std::numeric_limits<double>::quiet_NaN();
  assert(classify().missing_category == Pi0MissingCategory::other);
  set_mass(0.135);
  truth.delta_r = 0.7;
  assert(classify().missing_category == Pi0MissingCategory::energy_band_taggable);
  candidate.photon_in_cemc_acceptance[1] = false;
  assert(classify().missing_category == Pi0MissingCategory::acceptance);
  candidate.photon_projection_valid[1] = false;
  assert(classify().missing_category == Pi0MissingCategory::other);
  candidate.photon_projection_valid[1] = candidate.photon_in_cemc_acceptance[1] = true;
  truth.found = false;
  assert(classify().missing_category == Pi0MissingCategory::unclustered_or_no_cemc_deposit);
  candidate.photon_cemc_edep[1] = 0.0;
  assert(classify().missing_category == Pi0MissingCategory::unclustered_or_no_cemc_deposit);
  candidate.partner_diagnostics[1].found = true;
  candidate.partner_diagnostics[1].match.usable = false;
  assert(classify().missing_category == Pi0MissingCategory::other);
  truth.found = true;
  candidate.topology_partner_clusters[1] = truth;
  auto& partner = candidate.topology_partner_clusters[1];
  partner.cluster_energy = 0.6;
  partner.recovery = 0.5;
  assert(classify().topology == Pi0AnchorTopology::separated);
  partner.recovery = std::nextafter(0.5, 0.0);
  assert(classify().topology == Pi0AnchorTopology::missing);
  config.min_photon_energy_recovery = 0.0;
  assert(classify().topology == Pi0AnchorTopology::separated);
  partner.recovery = 0.0;
  assert(classify().topology == Pi0AnchorTopology::separated);
  partner.recovery = std::numeric_limits<double>::quiet_NaN();
  assert(classify().topology == Pi0AnchorTopology::separated);
  partner.found = false;
  assert(classify().topology == Pi0AnchorTopology::missing);
  partner.found = true;
  partner.recovery = 0.5;
  partner.cluster_id = cluster.cluster_id;
  assert(classify().topology == Pi0AnchorTopology::merged);
  candidate.photon_pre_cemc_interaction[1] = true;
  assert(classify().topology == Pi0AnchorTopology::single_contaminated);
  std::cout << "pi0 anchor classification: passed\n";
}
