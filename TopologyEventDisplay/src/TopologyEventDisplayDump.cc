#include "TopologyEventDisplayDump.h"

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace
{
constexpr std::size_t invalid_index = std::numeric_limits<std::size_t>::max();
constexpr double display_cemc_inner_radius = 93.0;

int cluster_id_from_index(const photon_tree::Pi0AnchorTopologyEventResult& result, std::size_t index)
{
  return index < result.clusters.size()
      ? static_cast<int>(result.clusters[index].cluster_id) : -999;
}

bool is_selected_candidate_core_track(int track_id, const photon_tree::Pi0AnchorTopologyEventResult& result)
{
  for (const auto& candidate : result.candidates)
  {
    if (track_id == candidate.g4_parent_track_id ||
        track_id == candidate.photon_track_ids[0] ||
        track_id == candidate.photon_track_ids[1])
    {
      return true;
    }
  }
  return false;
}
}

TopologyEventDisplayDump::TopologyEventDisplayDump(const std::string& name)
  : SubsysReco(name)
{
}

TopologyEventDisplayDump::~TopologyEventDisplayDump()
{
  close_output();
}

int TopologyEventDisplayDump::Init(PHCompositeNode*)
{
  const bool valid = !output_file_name_.empty() &&
      std::isfinite(config_.anchor_cluster_eta_max) &&
      config_.anchor_cluster_eta_max > 0.0 &&
      std::isfinite(config_.partner_cluster_eta_max) &&
      std::isfinite(config_.cemc_acceptance_eta_max) &&
      config_.cemc_acceptance_eta_max > 0.0 &&
      std::isfinite(config_.pre_cemc_interaction_radius) &&
      config_.pre_cemc_interaction_radius > 0.0 &&
      std::isfinite(config_.min_cluster_energy) &&
      config_.min_cluster_energy >= 0.0 &&
      config_.dominant_fraction_min >= 0.0 &&
      config_.dominant_fraction_min <= 1.0 &&
      config_.anchor_pi0_fraction_min >= 0.0 &&
      config_.anchor_pi0_fraction_min <= 1.0 &&
      config_.min_energy_contribution_fraction >= 0.0 &&
      first_event_ >= 0 &&
      config_.min_energy_contribution_fraction < 1.0 &&
      config_.min_photon_energy_recovery >= 0.0 &&
      config_.min_photon_energy_recovery <= 1.0 &&
      config_.min_direct_match_cluster_energy_coverage >= 0.0 &&
      config_.min_direct_match_cluster_energy_coverage <= 1.0 &&
      std::isfinite(config_.missing_diagnostic_max_delta_r) &&
      config_.missing_diagnostic_max_delta_r > 0.0;
  if (!valid)
  {
    std::cerr << "TopologyEventDisplayDump::Init - invalid configuration"
              << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  metadata_sample_mode_ = static_cast<int>(config_.sample_mode);
  metadata_write_detail_ = write_detail_ ? 1 : 0;
  config_.evaluate_all_candidates = write_detail_;
  config_.enable_missing_diagnostics = true;
  evaluator_.configure(config_);
  create_output_directory();
  create_output();
  if (!output_file_ || output_file_->IsZombie())
  {
    std::cerr << "TopologyEventDisplayDump::Init - cannot create "
              << output_file_name_ << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int TopologyEventDisplayDump::process_event(PHCompositeNode* topNode)
{
  const int source_event = static_cast<int>(source_events_seen_++);
  if (source_event < first_event_) return Fun4AllReturnCodes::EVENT_OK;
  b_event_ = source_event;
  ++events_processed_;
  const auto result = evaluator_.evaluate(topNode);
  if (result.status != photon_tree::Pi0TopologyEventStatus::accepted)
  {
    ++events_invalid_;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  b_n_truth_particles_ = 0;
  b_n_selected_family_particles_ = 0;
  b_n_g4_secondary_pi0_ = 0;
  b_n_clusters_ = 0;
  if (write_detail_)
  {
    fill_truth(topNode, result);
    fill_clusters(topNode, result);
    fill_candidate_cluster_truth(result);
  }
  else
  {
    auto* truth = findNode::getClass<PHG4TruthInfoContainer>(topNode, config_.truth_node_name);
    if (truth)
    {
      const auto range = truth->GetParticleRange();
      for (auto iterator = range.first; iterator != range.second; ++iterator)
      {
        const PHG4Particle* particle = iterator->second;
        if (!particle)
        {
          continue;
        }
        ++b_n_truth_particles_;
        if (particle->get_pid() == 111 && !truth->is_primary(particle))
        {
          ++b_n_g4_secondary_pi0_;
        }
        if (family_for_track(particle->get_track_id(), truth, result)[0] >= 0)
        {
          ++b_n_selected_family_particles_;
        }
      }
    }
    auto* clusters = findNode::getClass<RawClusterContainer>(topNode, config_.cluster_node_name);
    b_n_clusters_ = clusters ? static_cast<int>(clusters->size()) : 0;
  }
  fill_candidates(result);
  fill_anchors(result);
  fill_event(topNode, result);
  ++events_written_;
  if (verbosity_ > 0)
  {
    std::cout << "TopologyEventDisplayDump - event/candidates/anchors = "
              << b_event_ << "/" << result.candidates.size() << "/"
              << result.anchors.size() << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int TopologyEventDisplayDump::End(PHCompositeNode*)
{
  if (!output_file_ || !metadata_tree_)
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }
  metadata_tree_->Fill();
  output_file_->cd();
  output_file_->Write();
  const bool write_error = output_file_->TestBit(TFile::kWriteError);
  close_output();
  std::cout << "TopologyEventDisplayDump - processed/written/invalid = "
            << events_processed_ << "/" << events_written_ << "/"
            << events_invalid_ << std::endl;
  return write_error ? Fun4AllReturnCodes::ABORTRUN
                     : Fun4AllReturnCodes::EVENT_OK;
}

void TopologyEventDisplayDump::fill_event(PHCompositeNode*, const photon_tree::Pi0AnchorTopologyEventResult& result)
{
  b_collision_x_ = result.collision_vertex[0];
  b_collision_y_ = result.collision_vertex[1];
  b_collision_z_ = result.collision_vertex[2];
  b_n_candidates_ = static_cast<int>(result.candidates.size());
  b_n_g4_primary_pi0_ = 0;
  b_n_generator_pi0_ = 0;
  for (const auto& candidate : result.candidates)
  {
    if (candidate.pathway == photon_tree::Pi0Pathway::generator_decay)
    {
      ++b_n_generator_pi0_;
    }
    else
    {
      ++b_n_g4_primary_pi0_;
    }
  }
  b_n_topology_clusters_ = static_cast<int>(result.clusters.size());
  b_n_anchors_ = static_cast<int>(result.anchors.size());
  b_n_separated_ = 0;
  b_n_merged_ = 0;
  b_n_single_contaminated_ = 0;
  b_n_missing_ = 0;
  b_n_other_ = 0;
  for (const auto& anchor : result.anchors)
  {
    switch (anchor.topology)
    {
    case photon_tree::Pi0AnchorTopology::separated: ++b_n_separated_; break;
    case photon_tree::Pi0AnchorTopology::merged: ++b_n_merged_; break;
    case photon_tree::Pi0AnchorTopology::single_contaminated: ++b_n_single_contaminated_; break;
    case photon_tree::Pi0AnchorTopology::missing: ++b_n_missing_; break;
    case photon_tree::Pi0AnchorTopology::other: ++b_n_other_; break;
    }
  }
  events_tree_->Fill();
}

void TopologyEventDisplayDump::fill_candidates(const photon_tree::Pi0AnchorTopologyEventResult& result)
{
  for (std::size_t index = 0; index < result.candidates.size(); ++index)
  {
    const auto& candidate = result.candidates[index];
    b_candidate_id_ = static_cast<int>(index);
    b_pathway_ = static_cast<int>(candidate.pathway);
    b_pathway_name_ = photon_tree::pi0_pathway_name(candidate.pathway);
    b_parent_barcode_ = candidate.parent_barcode;
    b_g4_parent_track_id_ = candidate.g4_parent_track_id;
    b_energy_ = candidate.energy;
    b_pt_ = candidate.pt;
    b_eta_ = candidate.eta;
    b_phi_ = candidate.phi;
    b_photon0_track_id_ = candidate.photon_track_ids[0];
    b_photon1_track_id_ = candidate.photon_track_ids[1];
    b_photon0_energy_ = candidate.photon_energy[0];
    b_photon1_energy_ = candidate.photon_energy[1];
    b_photon0_eta_ = candidate.photon_eta[0];
    b_photon1_eta_ = candidate.photon_eta[1];
    b_photon0_phi_ = candidate.photon_phi[0];
    b_photon1_phi_ = candidate.photon_phi[1];
    b_photon0_projection_valid_ = candidate.photon_projection_valid[0] ? 1 : 0;
    b_photon1_projection_valid_ = candidate.photon_projection_valid[1] ? 1 : 0;
    b_photon0_projection_eta_ = candidate.photon_projection_eta[0];
    b_photon1_projection_eta_ = candidate.photon_projection_eta[1];
    b_photon0_projection_phi_ = candidate.photon_projection_phi[0];
    b_photon1_projection_phi_ = candidate.photon_projection_phi[1];
    b_photon0_in_cemc_acceptance_ = candidate.photon_in_cemc_acceptance[0] ? 1 : 0;
    b_photon1_in_cemc_acceptance_ = candidate.photon_in_cemc_acceptance[1] ? 1 : 0;
    b_photon0_first_daughter_vertex_valid_ = candidate.photon_first_daughter_vertex_valid[0] ? 1 : 0;
    b_photon1_first_daughter_vertex_valid_ = candidate.photon_first_daughter_vertex_valid[1] ? 1 : 0;
    b_photon0_first_daughter_radius_ = candidate.photon_first_daughter_radius[0];
    b_photon1_first_daughter_radius_ = candidate.photon_first_daughter_radius[1];
    b_photon0_pre_cemc_interaction_ = candidate.photon_pre_cemc_interaction[0] ? 1 : 0;
    b_photon1_pre_cemc_interaction_ = candidate.photon_pre_cemc_interaction[1] ? 1 : 0;
    b_photon0_cemc_edep_ = candidate.photon_cemc_edep[0];
    b_photon1_cemc_edep_ = candidate.photon_cemc_edep[1];
    b_best_cluster0_id_ = cluster_id_from_index(result, candidate.best_cluster[0]);
    b_best_cluster1_id_ = cluster_id_from_index(result, candidate.best_cluster[1]);
    b_maximum_edep0_ = candidate.maximum_edep[0];
    b_maximum_edep1_ = candidate.maximum_edep[1];
    b_reconstructed_photon0_energy_ = candidate.reconstructed_photon_energy[0];
    b_reconstructed_photon1_energy_ = candidate.reconstructed_photon_energy[1];
    b_recovered0_ = candidate.recovered[0] ? 1 : 0;
    b_recovered1_ = candidate.recovered[1] ? 1 : 0;
    b_topology_evaluated_ = candidate.topology_evaluated ? 1 : 0;
    candidates_tree_->Fill();
  }
}

void TopologyEventDisplayDump::fill_anchors(const photon_tree::Pi0AnchorTopologyEventResult& result)
{
  for (std::size_t index = 0; index < result.anchors.size(); ++index)
  {
    const auto& anchor = result.anchors[index];
    if (anchor.cluster_index >= result.clusters.size() ||
        anchor.candidate_index >= result.candidates.size())
    {
      continue;
    }
    const auto& candidate = result.candidates[anchor.candidate_index];
    b_anchor_id_ = static_cast<int>(index);
    b_candidate_id_ = static_cast<int>(anchor.candidate_index);
    b_cluster_id_ = result.clusters[anchor.cluster_index].cluster_id;
    b_energy_ = result.clusters[anchor.cluster_index].energy;
    b_et_ = result.clusters[anchor.cluster_index].et;
    b_topology_ = static_cast<int>(anchor.topology);
    b_topology_name_ = photon_tree::pi0_anchor_topology_name(anchor.topology);
    b_reason_ = static_cast<int>(anchor.reason);
    b_reason_name_ = photon_tree::pi0_anchor_reason_name(anchor.reason);
    b_missing_category_ = static_cast<int>(anchor.missing_category);
    b_missing_category_name_ = photon_tree::pi0_missing_category_name(anchor.missing_category);
    b_missing_detail_ = static_cast<int>(anchor.missing_detail);
    b_missing_detail_name_ = photon_tree::pi0_missing_detail_name(anchor.missing_detail);
    b_partner_photon_index_ = anchor.partner_photon_index;
    b_pre_cemc_photon_index_ = anchor.pre_cemc_photon_index;
    b_partner_cemc_edep_ = anchor.partner_photon_index >= 0 && anchor.partner_photon_index < 2
        ? candidate.photon_cemc_edep[static_cast<std::size_t>(anchor.partner_photon_index)] : 0.0;
    b_partner_diagnostic_invariant_mass_ = anchor.partner_diagnostic_invariant_mass;
    b_main_fraction_ = anchor.main_fraction;
    b_second_fraction_ = anchor.second_fraction;
    b_unmatched_max_fraction_ = anchor.unmatched_max_fraction;
    b_ambiguous_main_ = anchor.ambiguous_main ? 1 : 0;
    b_best_cluster0_id_ = cluster_id_from_index(result, candidate.best_cluster[0]);
    b_best_cluster1_id_ = cluster_id_from_index(result, candidate.best_cluster[1]);
    b_recovered0_ = candidate.recovered[0] ? 1 : 0;
    b_recovered1_ = candidate.recovered[1] ? 1 : 0;
    b_reconstructed_photon0_energy_ = candidate.reconstructed_photon_energy[0];
    b_reconstructed_photon1_energy_ = candidate.reconstructed_photon_energy[1];
    b_photon0_energy_ = candidate.photon_energy[0];
    b_photon1_energy_ = candidate.photon_energy[1];
    b_match_valid_ = 0;
    b_match_usable_ = 0;
    b_match_status_ = static_cast<int>(photon_tree::Pi0ClusterTruthMatchStatus::invalid);
    b_match_status_name_ = "invalid";
    b_match_failure_ = static_cast<int>(photon_tree::Pi0ClusterTruthMatchFailure::none);
    b_match_failure_name_ = "none";
    b_match_failure_ieta_ = b_match_failure_iphi_ = invalid_int_;
    b_match_tower_count_ = b_match_matched_tower_count_ = 0;
    b_match_cluster_member_energy_coverage_ = 0.0;
    b_total_edep_ = b_gamma0_edep_ = b_gamma1_edep_ = b_other_edep_ = 0.0;
    if (anchor.cluster_index < candidate.cluster_matches.size())
    {
      const auto& match = candidate.cluster_matches[anchor.cluster_index];
      b_match_valid_ = match.valid ? 1 : 0;
      b_match_usable_ = match.usable ? 1 : 0;
      b_match_status_ = static_cast<int>(match.status);
      b_match_status_name_ = photon_tree::pi0_cluster_truth_match_status_name(match.status);
      b_match_failure_ = static_cast<int>(match.failure);
      b_match_failure_name_ = photon_tree::pi0_cluster_truth_match_failure_name(match.failure);
      b_match_failure_ieta_ = match.failure_ieta;
      b_match_failure_iphi_ = match.failure_iphi;
      b_match_tower_count_ = match.tower_count;
      b_match_matched_tower_count_ = match.matched_tower_count;
      b_match_cluster_member_energy_coverage_ = match.cluster_member_energy_coverage;
      b_total_edep_ = match.total_edep;
      b_gamma0_edep_ = match.gamma_edep[0];
      b_gamma1_edep_ = match.gamma_edep[1];
      b_other_edep_ = match.other_edep;
    }
    b_partner_diagnostic_found_ = 0;
    b_partner_diagnostic_below_energy_threshold_ = 0;
    b_partner_diagnostic_has_direct_deposit_ = 0;
    b_partner_diagnostic_cluster_id_ = invalid_int_;
    b_partner_diagnostic_cluster_energy_ = b_partner_diagnostic_cluster_eta_ =
        b_partner_diagnostic_cluster_phi_ = b_partner_diagnostic_delta_r_ = invalid_double_;
    b_partner_diagnostic_reconstructed_energy_ = b_partner_diagnostic_recovery_ = 0.0;
    b_partner_diagnostic_match_usable_ = 0;
    b_partner_diagnostic_match_status_ = static_cast<int>(photon_tree::Pi0ClusterTruthMatchStatus::invalid);
    b_partner_diagnostic_match_status_name_ = "invalid";
    b_partner_diagnostic_match_failure_ = static_cast<int>(photon_tree::Pi0ClusterTruthMatchFailure::none);
    b_partner_diagnostic_match_failure_name_ = "none";
    b_partner_diagnostic_failure_ieta_ = b_partner_diagnostic_failure_iphi_ = invalid_int_;
    b_partner_diagnostic_match_coverage_ = 0.0;
    if (anchor.partner_photon_index >= 0 && anchor.partner_photon_index < 2)
    {
      const auto& diagnostic = candidate.partner_diagnostics[static_cast<std::size_t>(anchor.partner_photon_index)];
      b_partner_diagnostic_found_ = diagnostic.found ? 1 : 0;
      b_partner_diagnostic_below_energy_threshold_ = diagnostic.below_energy_threshold ? 1 : 0;
      b_partner_diagnostic_has_direct_deposit_ = diagnostic.has_direct_deposit ? 1 : 0;
      b_partner_diagnostic_cluster_id_ = diagnostic.found ? static_cast<int>(diagnostic.cluster_id) : invalid_int_;
      b_partner_diagnostic_cluster_energy_ = diagnostic.found ? diagnostic.cluster_energy : invalid_double_;
      b_partner_diagnostic_cluster_eta_ = diagnostic.found ? diagnostic.cluster_eta : invalid_double_;
      b_partner_diagnostic_cluster_phi_ = diagnostic.found ? diagnostic.cluster_phi : invalid_double_;
      b_partner_diagnostic_delta_r_ = diagnostic.found ? diagnostic.delta_r : invalid_double_;
      b_partner_diagnostic_reconstructed_energy_ = diagnostic.reconstructed_photon_energy;
      b_partner_diagnostic_recovery_ = diagnostic.recovery;
      b_partner_diagnostic_match_usable_ = diagnostic.match.usable ? 1 : 0;
      b_partner_diagnostic_match_status_ = static_cast<int>(diagnostic.match.status);
      b_partner_diagnostic_match_status_name_ = photon_tree::pi0_cluster_truth_match_status_name(diagnostic.match.status);
      b_partner_diagnostic_match_failure_ = static_cast<int>(diagnostic.match.failure);
      b_partner_diagnostic_match_failure_name_ = photon_tree::pi0_cluster_truth_match_failure_name(diagnostic.match.failure);
      b_partner_diagnostic_failure_ieta_ = diagnostic.match.failure_ieta;
      b_partner_diagnostic_failure_iphi_ = diagnostic.match.failure_iphi;
      b_partner_diagnostic_match_coverage_ = diagnostic.match.cluster_member_energy_coverage;
    }
    anchors_tree_->Fill();
  }
}

void TopologyEventDisplayDump::fill_candidate_cluster_truth(const photon_tree::Pi0AnchorTopologyEventResult& result)
{
  for (std::size_t candidate_index = 0;
       candidate_index < result.candidates.size(); ++candidate_index)
  {
    const auto& candidate = result.candidates[candidate_index];
    if (!candidate.topology_evaluated)
    {
      continue;
    }
    for (std::size_t cluster_index = 0;
         cluster_index < candidate.cluster_matches.size() &&
         cluster_index < result.clusters.size(); ++cluster_index)
    {
      const auto& match = candidate.cluster_matches[cluster_index];
      b_candidate_id_ = static_cast<int>(candidate_index);
      b_cluster_id_ = result.clusters[cluster_index].cluster_id;
      b_match_valid_ = match.valid ? 1 : 0;
      b_match_usable_ = match.usable ? 1 : 0;
      b_match_status_ = static_cast<int>(match.status);
      b_match_status_name_ = photon_tree::pi0_cluster_truth_match_status_name(match.status);
      b_match_failure_ = static_cast<int>(match.failure);
      b_match_failure_name_ = photon_tree::pi0_cluster_truth_match_failure_name(match.failure);
      b_match_failure_ieta_ = match.failure_ieta;
      b_match_failure_iphi_ = match.failure_iphi;
      b_match_tower_count_ = match.tower_count;
      b_match_matched_tower_count_ = match.matched_tower_count;
      b_match_cluster_member_energy_coverage_ = match.cluster_member_energy_coverage;
      b_total_edep_ = match.total_edep;
      b_gamma0_edep_ = match.gamma_edep[0];
      b_gamma1_edep_ = match.gamma_edep[1];
      b_other_edep_ = match.other_edep;
      const double inverse = match.total_edep > 0.0F
          ? 1.0 / match.total_edep : 0.0;
      b_gamma0_fraction_ = match.gamma_edep[0] * inverse;
      b_gamma1_fraction_ = match.gamma_edep[1] * inverse;
      b_other_fraction_ = match.other_edep * inverse;
      b_gamma0_recovery_estimate_ = candidate.photon_energy[0] > 0.0
          ? result.clusters[cluster_index].energy * b_gamma0_fraction_ /
                candidate.photon_energy[0] : 0.0;
      b_gamma1_recovery_estimate_ = candidate.photon_energy[1] > 0.0
          ? result.clusters[cluster_index].energy * b_gamma1_fraction_ /
                candidate.photon_energy[1] : 0.0;
      candidate_cluster_truth_tree_->Fill();
    }
  }
}

void TopologyEventDisplayDump::fill_clusters(PHCompositeNode* topNode, const photon_tree::Pi0AnchorTopologyEventResult& result)
{
  auto* clusters = findNode::getClass<RawClusterContainer>(topNode, config_.cluster_node_name);
  auto* towers = findNode::getClass<TowerInfoContainer>(topNode, config_.tower_node_name);
  if (!clusters || !towers)
  {
    return;
  }
  std::map<unsigned int, std::size_t> considered_by_id;
  for (std::size_t index = 0; index < result.clusters.size(); ++index)
  {
    considered_by_id[result.clusters[index].cluster_id] = index;
  }
  const auto range = clusters->getClusters();
  for (auto iterator = range.first; iterator != range.second; ++iterator)
  {
    const RawCluster* cluster = iterator->second;
    if (!cluster)
    {
      continue;
    }
    ++b_n_clusters_;
    b_cluster_id_ = cluster->get_id();
    b_energy_ = cluster->get_energy();
    b_x_ = cluster->get_x();
    b_y_ = cluster->get_y();
    b_z_ = cluster->get_z();
    b_r_ = std::hypot(b_x_, b_y_);
    const double dx = b_x_ - result.collision_vertex[0];
    const double dy = b_y_ - result.collision_vertex[1];
    const double dz = b_z_ - result.collision_vertex[2];
    const double transverse = std::hypot(dx, dy);
    const double distance = std::sqrt(transverse * transverse + dz * dz);
    b_eta_ = transverse > 0.0 ? std::asinh(dz / transverse) : invalid_double_;
    b_phi_ = std::atan2(dy, dx);
    b_et_ = distance > 0.0 ? b_energy_ * transverse / distance : invalid_double_;
    b_ntowers_ = static_cast<int>(cluster->getNTowers());
    const auto considered = considered_by_id.find(b_cluster_id_);
    b_topology_considered_ = considered != considered_by_id.end() ? 1 : 0;
    b_topology_cluster_index_ = b_topology_considered_
        ? static_cast<int>(considered->second) : invalid_int_;
    b_anchor_acceptance_ = b_topology_considered_
        ? (result.clusters[considered->second].anchor_acceptance ? 1 : 0) : 0;
    b_truth_valid_ = b_topology_considered_
        ? (result.clusters[considered->second].truth.valid ? 1 : 0) : 0;
    b_truth_total_edep_ = b_topology_considered_
        ? result.clusters[considered->second].truth.total_edep : 0.0;
    b_n_contributors_ = b_topology_considered_
        ? static_cast<int>(result.clusters[considered->second].truth.contributors.size()) : 0;
    clusters_tree_->Fill();

    if (b_topology_considered_)
    {
      const auto& contributors = result.clusters[considered->second].truth.contributors;
      for (std::size_t contribution = 0;
           contribution < contributors.size(); ++contribution)
      {
        const auto& value = contributors[contribution];
        b_contributor_index_ = static_cast<int>(contribution);
        b_g4_track_id_ = value.g4_track_id;
        b_g4_pdg_id_ = value.g4_pdg_id;
        b_embedding_id_ = value.embedding_id;
        b_hepmc_barcode_ = value.hepmc_barcode;
        b_contributor_edep_ = value.edep;
        b_contributor_fraction_ = value.fraction;
        cluster_contributors_tree_->Fill();
      }
    }

    const auto tower_range = cluster->get_towers();
    for (auto tower = tower_range.first; tower != tower_range.second; ++tower)
    {
      b_tower_key_ = tower->first;
      b_ieta_ = static_cast<int>(RawTowerDefs::decode_index1(b_tower_key_));
      b_iphi_ = static_cast<int>(RawTowerDefs::decode_index2(b_tower_key_));
      b_cluster_tower_energy_ = tower->second;
      const unsigned int tower_info_key = TowerInfoDefs::encode_emcal(static_cast<unsigned int>(b_ieta_), static_cast<unsigned int>(b_iphi_));
      TowerInfo* tower_info = towers->get_tower_at_key(static_cast<int>(tower_info_key));
      b_tower_energy_ = tower_info ? tower_info->get_energy() : invalid_double_;
      b_allocation_fraction_ = b_tower_energy_ > 0.0
          ? std::clamp(b_cluster_tower_energy_ / b_tower_energy_, 0.0, 1.0)
          : invalid_double_;
      cluster_towers_tree_->Fill();
    }
  }
}

void TopologyEventDisplayDump::fill_truth(PHCompositeNode* topNode, const photon_tree::Pi0AnchorTopologyEventResult& result)
{
  auto* truth = findNode::getClass<PHG4TruthInfoContainer>(topNode, config_.truth_node_name);
  if (!truth)
  {
    return;
  }
  std::map<int, std::vector<const PHG4VtxPoint*>> child_vertices;
  const auto all = truth->GetParticleRange();
  for (auto iterator = all.first; iterator != all.second; ++iterator)
  {
    const PHG4Particle* particle = iterator->second;
    const PHG4VtxPoint* vertex = particle
        ? truth->GetVtx(particle->get_vtx_id()) : nullptr;
    if (particle && vertex && particle->get_parent_id() != 0)
    {
      auto& vertices = child_vertices[particle->get_parent_id()];
      if (std::find(vertices.begin(), vertices.end(), vertex) == vertices.end()) vertices.push_back(vertex);
    }
  }
  for (auto& [parent, vertices] : child_vertices)
  {
    static_cast<void>(parent);
    std::stable_sort(vertices.begin(), vertices.end(), [](const PHG4VtxPoint* left, const PHG4VtxPoint* right) {
      return left->get_t() < right->get_t();
    });
  }

  const auto fill_segment = [&](double x0, double y0, double z0, double x1, double y1, double z1) {
    if (!std::isfinite(x0) || !std::isfinite(y0) || !std::isfinite(z0) ||
        !std::isfinite(x1) || !std::isfinite(y1) || !std::isfinite(z1) ||
        (x0 == x1 && y0 == y1 && z0 == z1)) return;
    b_x0_ = x0; b_y0_ = y0; b_z0_ = z0;
    b_x1_ = x1; b_y1_ = y1; b_z1_ = z1;
    truth_segments_tree_->Fill();
  };
  const auto clip_to_inner_radius = [](double x0, double y0, double z0, double x1, double y1, double z1,
                                       double& clipped_x, double& clipped_y, double& clipped_z) {
    const double dx = x1 - x0;
    const double dy = y1 - y0;
    const double dz = z1 - z0;
    const double a = dx * dx + dy * dy;
    const double b = 2.0 * (x0 * dx + y0 * dy);
    const double c = x0 * x0 + y0 * y0 - display_cemc_inner_radius * display_cemc_inner_radius;
    const double discriminant = b * b - 4.0 * a * c;
    if (!(a > 0.0) || discriminant < 0.0 || !std::isfinite(discriminant)) return false;
    const double root = std::sqrt(discriminant);
    double scale = std::numeric_limits<double>::infinity();
    for (const double candidate : {(-b - root) / (2.0 * a), (-b + root) / (2.0 * a)})
      if (candidate > 0.0 && candidate <= 1.0) scale = std::min(scale, candidate);
    if (!std::isfinite(scale)) return false;
    clipped_x = x0 + scale * dx;
    clipped_y = y0 + scale * dy;
    clipped_z = z0 + scale * dz;
    return std::isfinite(clipped_x) && std::isfinite(clipped_y) && std::isfinite(clipped_z);
  };

  for (auto iterator = all.first; iterator != all.second; ++iterator)
  {
    const PHG4Particle* particle = iterator->second;
    if (!particle)
    {
      continue;
    }
    ++b_n_truth_particles_;
    b_track_id_ = particle->get_track_id();
    b_pid_ = particle->get_pid();
    b_parent_id_ = particle->get_parent_id();
    b_primary_id_ = particle->get_primary_id();
    b_is_primary_ = truth->is_primary(particle) ? 1 : 0;
    b_is_g4_secondary_pi0_ = b_pid_ == 111 && !b_is_primary_ ? 1 : 0;
    if (b_is_g4_secondary_pi0_)
    {
      ++b_n_g4_secondary_pi0_;
    }
    const auto family = family_for_track(b_track_id_, truth, result);
    b_family_candidate_id_ = family[0];
    b_family_gamma_index_ = family[1];
    if (b_family_candidate_id_ >= 0)
    {
      ++b_n_selected_family_particles_;
    }
    b_px_ = particle->get_px();
    b_py_ = particle->get_py();
    b_pz_ = particle->get_pz();
    b_energy_ = particle->get_e();
    b_pt_ = std::hypot(b_px_, b_py_);
    b_eta_ = b_pt_ > 0.0 ? std::asinh(b_pz_ / b_pt_) : invalid_double_;
    b_phi_ = std::atan2(b_py_, b_px_);
    const PHG4VtxPoint* vertex = truth->GetVtx(particle->get_vtx_id());
    b_vx_ = vertex ? vertex->get_x() : invalid_double_;
    b_vy_ = vertex ? vertex->get_y() : invalid_double_;
    b_vz_ = vertex ? vertex->get_z() : invalid_double_;
    truth_particles_tree_->Fill();

    if (!vertex)
    {
      continue;
    }
    const auto children = child_vertices.find(b_track_id_);
    const bool core_track = is_selected_candidate_core_track(b_track_id_, result);
    if (children != child_vertices.end() && !children->second.empty())
    {
      if (b_family_candidate_id_ < 0 || core_track)
      {
        const PHG4VtxPoint* child = children->second.front();
        fill_segment(b_vx_, b_vy_, b_vz_, child->get_x(), child->get_y(), child->get_z());
      }
      else
      {
        double x0 = b_vx_;
        double y0 = b_vy_;
        double z0 = b_vz_;
        for (const PHG4VtxPoint* child : children->second)
        {
          if (std::hypot(x0, y0) >= display_cemc_inner_radius) break;
          const double x1 = child->get_x();
          const double y1 = child->get_y();
          const double z1 = child->get_z();
          if (std::hypot(x1, y1) <= display_cemc_inner_radius)
          {
            fill_segment(x0, y0, z0, x1, y1, z1);
            x0 = x1; y0 = y1; z0 = z1;
            continue;
          }
          double clipped_x = 0.0;
          double clipped_y = 0.0;
          double clipped_z = 0.0;
          if (clip_to_inner_radius(x0, y0, z0, x1, y1, z1, clipped_x, clipped_y, clipped_z))
            fill_segment(x0, y0, z0, clipped_x, clipped_y, clipped_z);
          break;
        }
      }
      continue;
    }
    if (b_family_candidate_id_ < 0 || core_track)
    {
      // Terminal selected-family shower descendants have no stored end vertex.
      // Projecting an inward descendant to r=100 cm would create a false detector-spanning chord.
      double x1 = 0.0;
      double y1 = 0.0;
      double z1 = 0.0;
      if (project_to_radius(b_vx_, b_vy_, b_vz_, b_px_, b_py_, b_pz_, 100.0, x1, y1, z1))
        fill_segment(b_vx_, b_vy_, b_vz_, x1, y1, z1);
    }
  }
}

std::array<int, 2> TopologyEventDisplayDump::family_for_track(
    int track_id, PHG4TruthInfoContainer* truth,
    const photon_tree::Pi0AnchorTopologyEventResult& result)
{
  std::set<int> visited;
  int current = track_id;
  while (current != 0 && visited.insert(current).second)
  {
    for (std::size_t candidate = 0;
         candidate < result.candidates.size(); ++candidate)
    {
      const auto& value = result.candidates[candidate];
      if (current == value.photon_track_ids[0])
      {
        return {static_cast<int>(candidate), 0};
      }
      if (current == value.photon_track_ids[1])
      {
        return {static_cast<int>(candidate), 1};
      }
      if (value.g4_parent_track_id != invalid_int_ &&
          current == value.g4_parent_track_id)
      {
        return {static_cast<int>(candidate), -1};
      }
    }
    const PHG4Particle* particle = truth ? truth->GetParticle(current) : nullptr;
    if (!particle)
    {
      break;
    }
    current = particle->get_parent_id();
  }
  return {-1, -1};
}

bool TopologyEventDisplayDump::project_to_radius(
    double x0, double y0, double z0,
    double px, double py, double pz, double radius,
    double& x1, double& y1, double& z1)
{
  const double a = px * px + py * py;
  if (!(a > 0.0) || !std::isfinite(a))
  {
    return false;
  }
  const double b = 2.0 * (x0 * px + y0 * py);
  const double c = x0 * x0 + y0 * y0 - radius * radius;
  const double discriminant = b * b - 4.0 * a * c;
  if (discriminant < 0.0 || !std::isfinite(discriminant))
  {
    return false;
  }
  const double root = std::sqrt(discriminant);
  const double first = (-b + root) / (2.0 * a);
  const double second = (-b - root) / (2.0 * a);
  double scale = std::numeric_limits<double>::infinity();
  if (first > 0.0) scale = first;
  if (second > 0.0) scale = std::min(scale, second);
  if (!std::isfinite(scale))
  {
    return false;
  }
  x1 = x0 + scale * px;
  y1 = y0 + scale * py;
  z1 = z0 + scale * pz;
  return std::isfinite(x1) && std::isfinite(y1) && std::isfinite(z1);
}

void TopologyEventDisplayDump::create_output_directory() const
{
  const std::size_t slash = output_file_name_.find_last_of('/');
  if (slash != std::string::npos)
  {
    const std::string directory = output_file_name_.substr(0, slash);
    if (!directory.empty())
    {
      gSystem->mkdir(directory.c_str(), true);
    }
  }
}

void TopologyEventDisplayDump::create_output()
{
  output_file_ = TFile::Open(output_file_name_.c_str(), "RECREATE");
  if (!output_file_ || output_file_->IsZombie())
  {
    return;
  }
  static int schema_version = schema_version_;
  static int evaluator_algorithm_version = photon_tree::Pi0AnchorTopologyEvaluator::kAlgorithmVersion;
  static int energy_match_algorithm_version = photon_tree::Pi0ClusterTruthMatcher::kAlgorithmVersion;
  metadata_tree_ = new TTree("metadata", "Topology event-display metadata");
  metadata_tree_->Branch("schema_version", &schema_version);
  metadata_tree_->Branch("evaluator_algorithm_version", &evaluator_algorithm_version);
  metadata_tree_->Branch("energy_match_algorithm_version", &energy_match_algorithm_version);
  metadata_tree_->Branch("sample_mode", &metadata_sample_mode_);
  metadata_tree_->Branch("write_detail", &metadata_write_detail_);
  metadata_tree_->Branch("source_label", &source_label_);
  metadata_tree_->Branch("manifest_path", &manifest_path_);
  metadata_tree_->Branch("manifest_begin", &manifest_begin_);
  metadata_tree_->Branch("manifest_end", &manifest_end_);
  metadata_tree_->Branch("cluster_node", &config_.cluster_node_name);
  metadata_tree_->Branch("tower_geom_node", &config_.tower_geom_node_name);
  metadata_tree_->Branch("first_event", &first_event_);
  metadata_tree_->Branch("anchor_cluster_eta_max", &config_.anchor_cluster_eta_max);
  metadata_tree_->Branch("partner_cluster_eta_max", &config_.partner_cluster_eta_max);
  metadata_tree_->Branch("cemc_acceptance_eta_max", &config_.cemc_acceptance_eta_max);
  metadata_tree_->Branch("pre_cemc_interaction_radius", &config_.pre_cemc_interaction_radius);
  metadata_tree_->Branch("min_cluster_energy", &config_.min_cluster_energy);
  metadata_tree_->Branch("dominant_fraction_min", &config_.dominant_fraction_min);
  metadata_tree_->Branch("anchor_pi0_fraction_min", &config_.anchor_pi0_fraction_min);
  metadata_tree_->Branch("min_energy_contribution_fraction", &config_.min_energy_contribution_fraction);
  metadata_tree_->Branch("min_photon_energy_recovery", &config_.min_photon_energy_recovery);
  metadata_tree_->Branch("min_direct_match_cluster_energy_coverage", &config_.min_direct_match_cluster_energy_coverage);
  metadata_tree_->Branch("missing_diagnostic_max_delta_r", &config_.missing_diagnostic_max_delta_r);
  metadata_tree_->Branch("evaluate_all_candidates", &config_.evaluate_all_candidates);
  metadata_tree_->Branch("enable_missing_diagnostics", &config_.enable_missing_diagnostics);
  metadata_tree_->Branch("events_processed", &events_processed_);
  metadata_tree_->Branch("events_written", &events_written_);
  metadata_tree_->Branch("events_invalid", &events_invalid_);

  events_tree_ = new TTree("events", "Event summary");
#define EVENT_BRANCH(name) events_tree_->Branch(#name, &b_##name##_)
  EVENT_BRANCH(event); EVENT_BRANCH(collision_x); EVENT_BRANCH(collision_y); EVENT_BRANCH(collision_z);
  EVENT_BRANCH(n_candidates); EVENT_BRANCH(n_g4_primary_pi0); EVENT_BRANCH(n_generator_pi0);
  EVENT_BRANCH(n_g4_secondary_pi0); EVENT_BRANCH(n_selected_family_particles); EVENT_BRANCH(n_truth_particles);
  EVENT_BRANCH(n_clusters); EVENT_BRANCH(n_topology_clusters); EVENT_BRANCH(n_anchors);
  EVENT_BRANCH(n_separated); EVENT_BRANCH(n_merged); EVENT_BRANCH(n_single_contaminated); EVENT_BRANCH(n_missing); EVENT_BRANCH(n_other);
#undef EVENT_BRANCH

  candidates_tree_ = new TTree("pi0_candidates", "Selected pi0 candidates");
#define CANDIDATE_BRANCH(name) candidates_tree_->Branch(#name, &b_##name##_)
  CANDIDATE_BRANCH(event); CANDIDATE_BRANCH(candidate_id); CANDIDATE_BRANCH(pathway); CANDIDATE_BRANCH(pathway_name);
  CANDIDATE_BRANCH(parent_barcode); CANDIDATE_BRANCH(g4_parent_track_id);
  CANDIDATE_BRANCH(energy); CANDIDATE_BRANCH(pt); CANDIDATE_BRANCH(eta); CANDIDATE_BRANCH(phi);
  CANDIDATE_BRANCH(photon0_track_id); CANDIDATE_BRANCH(photon1_track_id);
  CANDIDATE_BRANCH(photon0_energy); CANDIDATE_BRANCH(photon1_energy);
  CANDIDATE_BRANCH(photon0_eta); CANDIDATE_BRANCH(photon1_eta);
  CANDIDATE_BRANCH(photon0_phi); CANDIDATE_BRANCH(photon1_phi);
  CANDIDATE_BRANCH(photon0_projection_valid); CANDIDATE_BRANCH(photon1_projection_valid);
  CANDIDATE_BRANCH(photon0_projection_eta); CANDIDATE_BRANCH(photon1_projection_eta);
  CANDIDATE_BRANCH(photon0_projection_phi); CANDIDATE_BRANCH(photon1_projection_phi);
  CANDIDATE_BRANCH(photon0_in_cemc_acceptance); CANDIDATE_BRANCH(photon1_in_cemc_acceptance);
  CANDIDATE_BRANCH(photon0_first_daughter_vertex_valid); CANDIDATE_BRANCH(photon1_first_daughter_vertex_valid);
  CANDIDATE_BRANCH(photon0_first_daughter_radius); CANDIDATE_BRANCH(photon1_first_daughter_radius);
  CANDIDATE_BRANCH(photon0_pre_cemc_interaction); CANDIDATE_BRANCH(photon1_pre_cemc_interaction);
  CANDIDATE_BRANCH(photon0_cemc_edep); CANDIDATE_BRANCH(photon1_cemc_edep);
  CANDIDATE_BRANCH(best_cluster0_id); CANDIDATE_BRANCH(best_cluster1_id);
  CANDIDATE_BRANCH(maximum_edep0); CANDIDATE_BRANCH(maximum_edep1);
  CANDIDATE_BRANCH(reconstructed_photon0_energy); CANDIDATE_BRANCH(reconstructed_photon1_energy);
  CANDIDATE_BRANCH(recovered0); CANDIDATE_BRANCH(recovered1); CANDIDATE_BRANCH(topology_evaluated);
#undef CANDIDATE_BRANCH

  anchors_tree_ = new TTree("anchor_decisions", "Anchor topology decisions");
#define ANCHOR_BRANCH(name) anchors_tree_->Branch(#name, &b_##name##_)
  ANCHOR_BRANCH(event); ANCHOR_BRANCH(anchor_id); ANCHOR_BRANCH(candidate_id); ANCHOR_BRANCH(cluster_id);
  ANCHOR_BRANCH(energy); ANCHOR_BRANCH(et); ANCHOR_BRANCH(topology); ANCHOR_BRANCH(topology_name);
  ANCHOR_BRANCH(reason); ANCHOR_BRANCH(reason_name); ANCHOR_BRANCH(missing_category); ANCHOR_BRANCH(missing_category_name); ANCHOR_BRANCH(missing_detail); ANCHOR_BRANCH(missing_detail_name);
  ANCHOR_BRANCH(partner_photon_index); ANCHOR_BRANCH(pre_cemc_photon_index); ANCHOR_BRANCH(partner_cemc_edep);
  ANCHOR_BRANCH(partner_diagnostic_invariant_mass); ANCHOR_BRANCH(main_fraction);
  ANCHOR_BRANCH(second_fraction); ANCHOR_BRANCH(unmatched_max_fraction); ANCHOR_BRANCH(ambiguous_main);
  ANCHOR_BRANCH(best_cluster0_id); ANCHOR_BRANCH(best_cluster1_id);
  ANCHOR_BRANCH(photon0_energy); ANCHOR_BRANCH(photon1_energy);
  ANCHOR_BRANCH(reconstructed_photon0_energy); ANCHOR_BRANCH(reconstructed_photon1_energy);
  ANCHOR_BRANCH(recovered0); ANCHOR_BRANCH(recovered1);
  ANCHOR_BRANCH(match_valid); ANCHOR_BRANCH(match_usable); ANCHOR_BRANCH(match_status); ANCHOR_BRANCH(match_status_name);
  ANCHOR_BRANCH(match_failure); ANCHOR_BRANCH(match_failure_name);
  ANCHOR_BRANCH(match_failure_ieta); ANCHOR_BRANCH(match_failure_iphi);
  ANCHOR_BRANCH(match_tower_count); ANCHOR_BRANCH(match_matched_tower_count);
  ANCHOR_BRANCH(match_cluster_member_energy_coverage);
  ANCHOR_BRANCH(total_edep); ANCHOR_BRANCH(gamma0_edep); ANCHOR_BRANCH(gamma1_edep); ANCHOR_BRANCH(other_edep);
  ANCHOR_BRANCH(partner_diagnostic_found); ANCHOR_BRANCH(partner_diagnostic_below_energy_threshold);
  ANCHOR_BRANCH(partner_diagnostic_has_direct_deposit); ANCHOR_BRANCH(partner_diagnostic_cluster_id);
  ANCHOR_BRANCH(partner_diagnostic_cluster_energy); ANCHOR_BRANCH(partner_diagnostic_cluster_eta);
  ANCHOR_BRANCH(partner_diagnostic_cluster_phi); ANCHOR_BRANCH(partner_diagnostic_delta_r);
  ANCHOR_BRANCH(partner_diagnostic_reconstructed_energy); ANCHOR_BRANCH(partner_diagnostic_recovery);
  ANCHOR_BRANCH(partner_diagnostic_match_usable); ANCHOR_BRANCH(partner_diagnostic_match_status);
  ANCHOR_BRANCH(partner_diagnostic_match_status_name); ANCHOR_BRANCH(partner_diagnostic_match_failure);
  ANCHOR_BRANCH(partner_diagnostic_match_failure_name);
  ANCHOR_BRANCH(partner_diagnostic_failure_ieta); ANCHOR_BRANCH(partner_diagnostic_failure_iphi);
  ANCHOR_BRANCH(partner_diagnostic_match_coverage);

#undef ANCHOR_BRANCH
  candidate_cluster_truth_tree_ = new TTree("candidate_cluster_truth", "Per-candidate per-cluster direct daughter deposits");
#define MATCH_BRANCH(name) candidate_cluster_truth_tree_->Branch(#name, &b_##name##_)
  MATCH_BRANCH(event); MATCH_BRANCH(candidate_id); MATCH_BRANCH(cluster_id);
  MATCH_BRANCH(match_valid); MATCH_BRANCH(match_usable); MATCH_BRANCH(match_status); MATCH_BRANCH(match_status_name);
  MATCH_BRANCH(match_failure); MATCH_BRANCH(match_failure_name);
  MATCH_BRANCH(match_failure_ieta); MATCH_BRANCH(match_failure_iphi);
  MATCH_BRANCH(match_tower_count); MATCH_BRANCH(match_matched_tower_count);
  MATCH_BRANCH(match_cluster_member_energy_coverage);
  MATCH_BRANCH(total_edep); MATCH_BRANCH(gamma0_edep); MATCH_BRANCH(gamma1_edep); MATCH_BRANCH(other_edep);
  MATCH_BRANCH(gamma0_fraction); MATCH_BRANCH(gamma1_fraction); MATCH_BRANCH(other_fraction);
  MATCH_BRANCH(gamma0_recovery_estimate); MATCH_BRANCH(gamma1_recovery_estimate);
#undef MATCH_BRANCH

  clusters_tree_ = new TTree("clusters", "All reconstructed clusters");
#define CLUSTER_BRANCH(name) clusters_tree_->Branch(#name, &b_##name##_)
  CLUSTER_BRANCH(event); CLUSTER_BRANCH(cluster_id); CLUSTER_BRANCH(energy); CLUSTER_BRANCH(et);
  CLUSTER_BRANCH(x); CLUSTER_BRANCH(y); CLUSTER_BRANCH(z); CLUSTER_BRANCH(r); CLUSTER_BRANCH(eta); CLUSTER_BRANCH(phi);
  CLUSTER_BRANCH(ntowers); CLUSTER_BRANCH(topology_considered); CLUSTER_BRANCH(topology_cluster_index);
  CLUSTER_BRANCH(anchor_acceptance); CLUSTER_BRANCH(truth_valid); CLUSTER_BRANCH(truth_total_edep); CLUSTER_BRANCH(n_contributors);
#undef CLUSTER_BRANCH

  cluster_contributors_tree_ = new TTree("cluster_contributors", "Flat cluster truth contributor table");
#define CONTRIBUTOR_BRANCH(name) cluster_contributors_tree_->Branch(#name, &b_##name##_)
  CONTRIBUTOR_BRANCH(event); CONTRIBUTOR_BRANCH(cluster_id); CONTRIBUTOR_BRANCH(contributor_index);
  CONTRIBUTOR_BRANCH(g4_track_id); CONTRIBUTOR_BRANCH(g4_pdg_id); CONTRIBUTOR_BRANCH(embedding_id);
  CONTRIBUTOR_BRANCH(hepmc_barcode); CONTRIBUTOR_BRANCH(contributor_edep); CONTRIBUTOR_BRANCH(contributor_fraction);
#undef CONTRIBUTOR_BRANCH

  cluster_towers_tree_ = new TTree("cluster_towers", "Cluster member tower allocation");
#define TOWER_BRANCH(name) cluster_towers_tree_->Branch(#name, &b_##name##_)
  TOWER_BRANCH(event); TOWER_BRANCH(cluster_id); TOWER_BRANCH(tower_key); TOWER_BRANCH(ieta); TOWER_BRANCH(iphi);
  TOWER_BRANCH(cluster_tower_energy); TOWER_BRANCH(tower_energy); TOWER_BRANCH(allocation_fraction);
#undef TOWER_BRANCH

  truth_particles_tree_ = new TTree("truth_particles", "G4 truth particles with selected-pi0 family labels");
#define TRUTH_BRANCH(name) truth_particles_tree_->Branch(#name, &b_##name##_)
  TRUTH_BRANCH(event); TRUTH_BRANCH(track_id); TRUTH_BRANCH(pid); TRUTH_BRANCH(parent_id); TRUTH_BRANCH(primary_id);
  TRUTH_BRANCH(is_primary); TRUTH_BRANCH(is_g4_secondary_pi0); TRUTH_BRANCH(family_candidate_id); TRUTH_BRANCH(family_gamma_index);
  TRUTH_BRANCH(px); TRUTH_BRANCH(py); TRUTH_BRANCH(pz); TRUTH_BRANCH(energy); TRUTH_BRANCH(pt); TRUTH_BRANCH(eta); TRUTH_BRANCH(phi);
  TRUTH_BRANCH(vx); TRUTH_BRANCH(vy); TRUTH_BRANCH(vz);
#undef TRUTH_BRANCH

  truth_segments_tree_ = new TTree("truth_segments", "Straight truth display segments; not propagated charged tracks");
#define SEGMENT_BRANCH(name) truth_segments_tree_->Branch(#name, &b_##name##_)
  SEGMENT_BRANCH(event); SEGMENT_BRANCH(track_id); SEGMENT_BRANCH(pid); SEGMENT_BRANCH(parent_id);
  SEGMENT_BRANCH(family_candidate_id); SEGMENT_BRANCH(family_gamma_index);
  SEGMENT_BRANCH(x0); SEGMENT_BRANCH(y0); SEGMENT_BRANCH(z0); SEGMENT_BRANCH(x1); SEGMENT_BRANCH(y1); SEGMENT_BRANCH(z1);
#undef SEGMENT_BRANCH
}

void TopologyEventDisplayDump::close_output()
{
  if (output_file_)
  {
    if (output_file_->IsOpen()) output_file_->Close();
    delete output_file_;
    output_file_ = nullptr;
  }
}
