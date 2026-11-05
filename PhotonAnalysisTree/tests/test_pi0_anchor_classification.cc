#include "Pi0AnchorTopologyEvaluator.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>

using namespace photon_tree;

int main()
{
  Pi0TruthPartnerClusterRecord best, high_deposit, high_energy;
  high_deposit.found = high_energy.found = true;
  high_deposit.cluster_id = 20; high_deposit.cluster_energy = .3; high_deposit.direct_edep = .2;
  high_energy.cluster_id = 10; high_energy.cluster_energy = 5.; high_energy.direct_edep = .1;
  assert(prefer_pi0_truth_partner(high_energy, best)); best = high_energy;
  assert(prefer_pi0_truth_partner(high_deposit, best)); best = high_deposit;
  assert(!prefer_pi0_truth_partner(high_energy, best)); // No low-deposit replacement after energy threshold failure.
  high_energy.direct_edep = best.direct_edep;
  assert(prefer_pi0_truth_partner(high_energy, best)); best = high_energy;
  assert(!prefer_pi0_truth_partner(high_deposit, best)); // Equal deposits have a deterministic cluster-ID tie break.
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
  candidate.photon_projection_no_intersection[1] = true;
  assert(classify().missing_category == Pi0MissingCategory::acceptance); // Known non-arrival differs from missing vertex information.
  candidate.photon_projection_no_intersection[1] = false;
  candidate.photon_projection_valid[1] = candidate.photon_in_cemc_acceptance[1] = true;
  truth.found = false;
  assert(classify().missing_category == Pi0MissingCategory::unclustered_or_no_cemc_deposit);
  candidate.photon_cemc_edep[1] = 0.0;
  assert(classify().missing_category == Pi0MissingCategory::unclustered_or_no_cemc_deposit);
  candidate.partner_diagnostics[1].found = true;
  candidate.partner_diagnostics[1].match.usable = false;
  assert(classify().missing_category == Pi0MissingCategory::other);
  truth.found = true;
  auto& partner = truth;
  // The representative stays below threshold even if some other cluster could pass it.
  // There is no independently selected threshold-qualified partner anymore.
  config.min_photon_energy_recovery = 0.0;
  partner.cluster_energy = 0.3;
  assert(classify().topology == Pi0AnchorTopology::missing);
  partner.cluster_energy = 0.5;
  assert(classify().topology == Pi0AnchorTopology::missing);
  partner.cluster_energy = std::nextafter(0.5, 1.0);
  assert(classify().topology == Pi0AnchorTopology::separated);
  config.min_photon_energy_recovery = 0.5;
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
  partner.cluster_energy = 0.5;
  assert(classify().topology == Pi0AnchorTopology::missing);
  // A noisy cluster and a failed sibling anchor must not veto the healthy anchor of the same meson.
  candidate.cluster_matches.resize(3);
  candidate.matching_usable = false; candidate.matching_status = 1; candidate.minimum_coverage = .1;
  auto& healthy = candidate.cluster_matches[0];
  healthy.status = Pi0ClusterTruthMatchStatus::partial; healthy.cluster_member_energy_coverage = .5; healthy.usable = true;
  candidate.cluster_matches[1].status = Pi0ClusterTruthMatchStatus::partial;
  candidate.cluster_matches[1].cluster_member_energy_coverage = .45;
  candidate.cluster_matches[2].failure = Pi0ClusterTruthMatchFailure::missing_truth_tower;
  composition::Result good, bad;
  finalize_common_topology(candidate, 0, good); finalize_common_topology(candidate, 1, bad);
  assert(good.applicability == composition::Applicability::applicable && good.topology_reason == 0);
  assert(good.anchor_matching_status == 1 && good.anchor_matching_coverage == .5 && good.anchor_matching_usable);
  assert(bad.applicability == composition::Applicability::unknown && bad.topology_reason == 4 && !bad.anchor_matching_usable);
  assert(!candidate.matching_usable && good.matching_status == 1 && good.minimum_coverage == .1);
  assert(good.minimum_coverage == bad.minimum_coverage && good.matching_status == bad.matching_status);
  // Missing/incomplete partner information stays in the missing diagnostics, not a new applicability veto.
  candidate.photon_pre_cemc_interaction[1] = false;
  truth.found = false; candidate.partner_diagnostics[1].found = true; candidate.partner_diagnostics[1].match.usable = false;
  anchor = classify();
  assert(anchor.topology == Pi0AnchorTopology::missing && anchor.missing_detail == Pi0MissingDetail::partner_direct_match_incomplete);
  finalize_common_topology(candidate, 0, good);
  assert(good.applicability == composition::Applicability::applicable);
  candidate.photon_projection_valid[1] = false;
  finalize_common_topology(candidate, 0, good);
  assert(good.applicability == composition::Applicability::unknown && good.topology_reason == 5);
  candidate.photon_projection_no_intersection[1] = true;
  finalize_common_topology(candidate, 0, good);
  assert(good.applicability == composition::Applicability::applicable);
  candidate.interaction_info_valid = false;
  finalize_common_topology(candidate, 0, good);
  assert(good.applicability == composition::Applicability::unknown && good.topology_reason == 6);
  candidate.interaction_info_valid = true;
  healthy.usable = false; healthy.ancestry_complete = false;
  finalize_common_topology(candidate, 0, good);
  assert(good.applicability == composition::Applicability::unknown && good.topology_reason == 4);
  assert(Pi0AnchorTopologyEvaluator::kAlgorithmVersion == 14);
  std::cout << "pi0 anchor classification: passed\n";
}
