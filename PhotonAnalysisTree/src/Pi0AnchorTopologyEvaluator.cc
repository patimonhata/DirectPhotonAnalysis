#include "Pi0AnchorTopologyEvaluator.h"

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/TowerInfoContainer.h>
#include <g4detectors/PHG4CellContainer.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>
#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <HepMC/GenEvent.h>
#include <HepMC/GenParticle.h>
#include <HepMC/GenVertex.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <utility>
#include <vector>

namespace
{
constexpr std::size_t invalid_index = std::numeric_limits<std::size_t>::max();

struct Pi0Origin
{
  bool valid = false;
  const HepMC::GenParticle* parent = nullptr;
};

struct CandidateWork
{
  photon_tree::Pi0TopologyCandidateRecord record;
  const HepMC::GenParticle* hepmc_parent = nullptr;
};

std::vector<const HepMC::GenParticle*> incoming(const HepMC::GenVertex* vertex)
{
  std::vector<const HepMC::GenParticle*> result;
  if (!vertex)
  {
    return result;
  }
  for (auto iterator = vertex->particles_in_const_begin();
       iterator != vertex->particles_in_const_end(); ++iterator)
  {
    if (*iterator)
    {
      result.push_back(*iterator);
    }
  }
  return result;
}

Pi0Origin trace_pi0_origin(const HepMC::GenParticle* photon)
{
  Pi0Origin result;
  if (!photon || photon->pdg_id() != 22)
  {
    return result;
  }
  const HepMC::GenParticle* current = photon;
  std::set<int> visited;
  if (current->barcode() != 0)
  {
    visited.insert(current->barcode());
  }
  while (current)
  {
    const auto parents = incoming(current->production_vertex());
    if (parents.size() == 1U && parents.front()->pdg_id() == 22)
    {
      current = parents.front();
      if (current->barcode() != 0 &&
          !visited.insert(current->barcode()).second)
      {
        return result;
      }
      continue;
    }
    if (parents.size() == 1U && parents.front()->pdg_id() == 111)
    {
      result.valid = true;
      result.parent = parents.front();
    }
    return result;
  }
  return result;
}

bool finite_hepmc_kinematics(const HepMC::GenParticle* particle,
                            double& energy, double& pt,
                            double& eta, double& phi)
{
  if (!particle)
  {
    return false;
  }
  const auto& momentum = particle->momentum();
  energy = momentum.e();
  pt = std::hypot(momentum.px(), momentum.py());
  eta = pt > 0.0 ? std::asinh(momentum.pz() / pt) : 0.0;
  phi = std::atan2(momentum.py(), momentum.px());
  return std::isfinite(energy) && pt > 0.0 && std::isfinite(pt) &&
      std::isfinite(eta) && std::isfinite(phi);
}

bool finite_g4_kinematics(const PHG4Particle* particle,
                         double& energy, double& pt,
                         double& eta, double& phi)
{
  if (!particle)
  {
    return false;
  }
  energy = particle->get_e();
  pt = std::hypot(particle->get_px(), particle->get_py());
  eta = pt > 0.0 ? std::asinh(particle->get_pz() / pt) : 0.0;
  phi = std::atan2(particle->get_py(), particle->get_px());
  return std::isfinite(energy) && pt > 0.0 && std::isfinite(pt) &&
      std::isfinite(eta) && std::isfinite(phi);
}

const HepMC::GenParticle* contributor_hepmc_particle(
    const photon_tree::TruthContributor& contributor,
    const PHHepMCGenEventMap* event_map)
{
  const PHHepMCGenEvent* subevent = event_map
      ? event_map->get(contributor.embedding_id) : nullptr;
  const HepMC::GenEvent* event = subevent ? subevent->getEvent() : nullptr;
  return event ? event->barcode_to_particle(contributor.hepmc_barcode) : nullptr;
}

std::size_t contributor_candidate_index(
    const photon_tree::TruthContributor& contributor,
    const photon_tree::Pi0AnchorTopologyConfig& config,
    const PHHepMCGenEventMap* event_map,
    const std::map<int, std::size_t>& g4_candidate_by_barcode,
    const std::map<int, std::size_t>& generator_candidate_by_barcode,
    const std::map<int, std::size_t>& single_candidate_by_track)
{
  if (config.sample_mode == photon_tree::Pi0SampleMode::single_particle)
  {
    const auto found = single_candidate_by_track.find(contributor.g4_track_id);
    return found == single_candidate_by_track.end()
        ? invalid_index : found->second;
  }
  if (contributor.embedding_id != config.signal_embedding_id ||
      !contributor.hepmc_valid)
  {
    return invalid_index;
  }
  const HepMC::GenParticle* particle =
      contributor_hepmc_particle(contributor, event_map);
  if (contributor.g4_pdg_id == 111 && particle && particle->pdg_id() == 111)
  {
    const auto found = g4_candidate_by_barcode.find(particle->barcode());
    return found == g4_candidate_by_barcode.end()
        ? invalid_index : found->second;
  }
  if (contributor.g4_pdg_id != 22)
  {
    return invalid_index;
  }
  const Pi0Origin origin = trace_pi0_origin(particle);
  if (!origin.valid || !origin.parent)
  {
    return invalid_index;
  }
  const auto found =
      generator_candidate_by_barcode.find(origin.parent->barcode());
  return found == generator_candidate_by_barcode.end()
      ? invalid_index : found->second;
}

void fill_photon_kinematics(photon_tree::Pi0TopologyCandidateRecord& record,
                            std::size_t index,
                            const PHG4Particle* photon)
{
  record.photon_track_ids[index] = photon ? photon->get_track_id() : -999;
  double pt = 0.0;
  if (!finite_g4_kinematics(photon, record.photon_energy[index], pt,
                            record.photon_eta[index], record.photon_phi[index]))
  {
    record.photon_energy[index] = photon ? photon->get_e() : 0.0;
    record.photon_eta[index] = 0.0;
    record.photon_phi[index] = 0.0;
  }
}
}

namespace photon_tree
{
const char* pi0_pathway_name(Pi0Pathway value)
{
  switch (value)
  {
  case Pi0Pathway::g4_primary_decay: return "g4_primary_decay";
  case Pi0Pathway::generator_decay: return "generator_decay";
  case Pi0Pathway::single_particle_g4_decay: return "single_particle_g4_decay";
  }
  return "unknown";
}

const char* pi0_anchor_topology_name(Pi0AnchorTopology value)
{
  switch (value)
  {
  case Pi0AnchorTopology::separated: return "separated";
  case Pi0AnchorTopology::merged: return "merged";
  case Pi0AnchorTopology::missing: return "missing";
  case Pi0AnchorTopology::other: return "other";
  }
  return "unknown";
}

const char* pi0_anchor_reason_name(Pi0AnchorReason value)
{
  switch (value)
  {
  case Pi0AnchorReason::other_not_daughter_maximum:
    return "other_not_daughter_maximum";
  case Pi0AnchorReason::separated_distinct_recovered_clusters:
    return "separated_distinct_recovered_clusters";
  case Pi0AnchorReason::merged_shared_recovered_cluster:
    return "merged_shared_recovered_cluster";
  case Pi0AnchorReason::missing_unrecovered_partner:
    return "missing_unrecovered_partner";
  case Pi0AnchorReason::ambiguous_main_contributor:
    return "ambiguous_main_contributor";
  case Pi0AnchorReason::other_best_cluster_below_recovery:
    return "other_best_cluster_below_recovery";
  }
  return "unknown";
}

void Pi0AnchorTopologyEvaluator::configure(
    const Pi0AnchorTopologyConfig& config)
{
  config_ = config;
  truth_matcher_.set_verbosity(config_.verbosity);
}

Pi0AnchorTopologyEventResult Pi0AnchorTopologyEvaluator::evaluate(
    PHCompositeNode* topNode)
{
  Pi0AnchorTopologyEventResult result;
  auto* truth = findNode::getClass<PHG4TruthInfoContainer>(
      topNode, config_.truth_node_name);
  auto* event_map = findNode::getClass<PHHepMCGenEventMap>(
      topNode, config_.hepmc_event_map_node_name);
  auto* towers = findNode::getClass<TowerInfoContainer>(
      topNode, config_.tower_node_name);
  auto* raw_truth_towers = findNode::getClass<RawTowerContainer>(
      topNode, config_.raw_truth_tower_node_name);
  auto* truth_cells = findNode::getClass<PHG4CellContainer>(
      topNode, config_.truth_cell_node_name);
  auto* truth_hits = findNode::getClass<PHG4HitContainer>(
      topNode, config_.truth_hit_node_name);
  auto* clusters = findNode::getClass<RawClusterContainer>(
      topNode, config_.cluster_node_name);

  const bool pythia = config_.sample_mode == Pi0SampleMode::pythia;
  const PHHepMCGenEvent* signal_event = pythia && event_map
      ? event_map->get(config_.signal_embedding_id) : nullptr;
  const HepMC::GenEvent* event = signal_event ? signal_event->getEvent() : nullptr;
  if (!truth || !towers || !raw_truth_towers || !truth_cells || !truth_hits ||
      !clusters || (pythia && (!event_map || !signal_event ||
       !signal_event->is_simulated() || !event)))
  {
    return result;
  }

  if (pythia)
  {
    const auto& collision = signal_event->get_collision_vertex();
    result.collision_vertex = {collision.x(), collision.y(), collision.z()};
  }
  else
  {
    bool found_vertex = false;
    const auto primary_range = truth->GetPrimaryParticleRange();
    for (auto iterator = primary_range.first;
         iterator != primary_range.second && !found_vertex; ++iterator)
    {
      const PHG4Particle* particle = iterator->second;
      const PHG4VtxPoint* vertex = particle
          ? truth->GetVtx(particle->get_vtx_id()) : nullptr;
      if (vertex)
      {
        result.collision_vertex = {vertex->get_x(), vertex->get_y(), vertex->get_z()};
        found_vertex = true;
      }
    }
    if (!found_vertex)
    {
      return result;
    }
  }
  if (!std::isfinite(result.collision_vertex[0]) ||
      !std::isfinite(result.collision_vertex[1]) ||
      !std::isfinite(result.collision_vertex[2]))
  {
    return result;
  }
  if (config_.max_abs_vertex_z > 0.0 &&
      std::abs(result.collision_vertex[2]) >= config_.max_abs_vertex_z)
  {
    result.status = Pi0TopologyEventStatus::vertex_rejected;
    return result;
  }
  if (!truth_matcher_.begin_event(topNode))
  {
    return result;
  }

  std::map<int, std::vector<const PHG4Particle*>> children_by_parent;
  const auto secondary_range = truth->GetSecondaryParticleRange();
  for (auto iterator = secondary_range.first;
       iterator != secondary_range.second; ++iterator)
  {
    const PHG4Particle* particle = iterator->second;
    if (particle)
    {
      children_by_parent[particle->get_parent_id()].push_back(particle);
    }
  }

  std::vector<CandidateWork> candidates;
  if (pythia)
  {
    struct PendingGenerator
    {
      const HepMC::GenParticle* parent = nullptr;
      std::vector<const PHG4Particle*> photons;
    };
    std::map<int, PendingGenerator> pending_generator;
    const auto primary_range = truth->GetPrimaryParticleRange();
    for (auto iterator = primary_range.first;
         iterator != primary_range.second; ++iterator)
    {
      const PHG4Particle* primary = iterator->second;
      if (!primary ||
          truth->isEmbeded(primary->get_track_id()) != config_.signal_embedding_id)
      {
        continue;
      }
      const HepMC::GenParticle* hepmc_particle =
          event->barcode_to_particle(primary->get_barcode());
      if (primary->get_pid() == 111 && hepmc_particle &&
          hepmc_particle->pdg_id() == 111)
      {
        const auto found = children_by_parent.find(primary->get_track_id());
        if (found == children_by_parent.end() || found->second.size() != 2U ||
            found->second[0]->get_pid() != 22 ||
            found->second[1]->get_pid() != 22)
        {
          ++result.malformed_candidate_count;
          continue;
        }
        CandidateWork work;
        work.record.pathway = Pi0Pathway::g4_primary_decay;
        work.record.parent_barcode = hepmc_particle->barcode();
        work.record.g4_parent_track_id = primary->get_track_id();
        if (!finite_g4_kinematics(primary, work.record.energy, work.record.pt,
                                  work.record.eta, work.record.phi) ||
            std::abs(work.record.eta) >= config_.truth_eta_max)
        {
          continue;
        }
        fill_photon_kinematics(work.record, 0, found->second[0]);
        fill_photon_kinematics(work.record, 1, found->second[1]);
        work.hepmc_parent = hepmc_particle;
        candidates.push_back(std::move(work));
        ++result.g4_candidate_count;
        continue;
      }
      if (primary->get_pid() != 22 || !hepmc_particle ||
          hepmc_particle->pdg_id() != 22)
      {
        continue;
      }
      const Pi0Origin origin = trace_pi0_origin(hepmc_particle);
      if (!origin.valid || !origin.parent)
      {
        continue;
      }
      auto& pending = pending_generator[origin.parent->barcode()];
      pending.parent = origin.parent;
      pending.photons.push_back(primary);
    }
    for (const auto& [barcode, pending] : pending_generator)
    {
      if (!pending.parent || pending.photons.size() != 2U)
      {
        ++result.malformed_candidate_count;
        continue;
      }
      CandidateWork work;
      work.record.pathway = Pi0Pathway::generator_decay;
      work.record.parent_barcode = barcode;
      if (!finite_hepmc_kinematics(
              pending.parent, work.record.energy, work.record.pt,
              work.record.eta, work.record.phi) ||
          std::abs(work.record.eta) >= config_.truth_eta_max)
      {
        continue;
      }
      fill_photon_kinematics(work.record, 0, pending.photons[0]);
      fill_photon_kinematics(work.record, 1, pending.photons[1]);
      work.hepmc_parent = pending.parent;
      candidates.push_back(std::move(work));
      ++result.generator_candidate_count;
    }
  }
  else
  {
    const auto primary_range = truth->GetPrimaryParticleRange();
    for (auto iterator = primary_range.first;
         iterator != primary_range.second; ++iterator)
    {
      const PHG4Particle* primary = iterator->second;
      if (!primary || primary->get_pid() != 111)
      {
        continue;
      }
      const auto found = children_by_parent.find(primary->get_track_id());
      if (found == children_by_parent.end() || found->second.size() != 2U ||
          found->second[0]->get_pid() != 22 ||
          found->second[1]->get_pid() != 22)
      {
        ++result.malformed_candidate_count;
        continue;
      }
      CandidateWork work;
      work.record.pathway = Pi0Pathway::single_particle_g4_decay;
      work.record.parent_barcode = primary->get_barcode();
      work.record.g4_parent_track_id = primary->get_track_id();
      if (!finite_g4_kinematics(primary, work.record.energy, work.record.pt,
                                work.record.eta, work.record.phi) ||
          std::abs(work.record.eta) >= config_.truth_eta_max)
      {
        continue;
      }
      fill_photon_kinematics(work.record, 0, found->second[0]);
      fill_photon_kinematics(work.record, 1, found->second[1]);
      candidates.push_back(std::move(work));
      ++result.g4_candidate_count;
    }
  }

  const auto cluster_range = clusters->getClusters();
  for (auto iterator = cluster_range.first; iterator != cluster_range.second;
       ++iterator)
  {
    const RawCluster* cluster = iterator->second;
    if (!cluster || !std::isfinite(cluster->get_energy()) ||
        !std::isfinite(cluster->get_x()) || !std::isfinite(cluster->get_y()) ||
        !std::isfinite(cluster->get_z()) ||
        cluster->get_energy() < config_.min_cluster_energy)
    {
      continue;
    }
    const double dx = cluster->get_x() - result.collision_vertex[0];
    const double dy = cluster->get_y() - result.collision_vertex[1];
    const double dz = cluster->get_z() - result.collision_vertex[2];
    const double distance = std::sqrt(dx * dx + dy * dy + dz * dz);
    const double transverse = std::hypot(dx, dy);
    if (!(distance > 0.0) || !(transverse > 0.0))
    {
      continue;
    }
    Pi0TopologyClusterRecord record;
    record.cluster = cluster;
    record.cluster_id = cluster->get_id();
    record.energy = cluster->get_energy();
    record.et = record.energy * transverse / distance;
    record.eta = std::asinh(dz / transverse);
    if (!std::isfinite(record.et) || !(record.et > 0.0) ||
        !std::isfinite(record.eta))
    {
      continue;
    }
    if (config_.partner_cluster_eta_max > 0.0 &&
        std::abs(record.eta) >= config_.partner_cluster_eta_max)
    {
      continue;
    }
    record.anchor_acceptance =
        std::abs(record.eta) < config_.anchor_cluster_eta_max;
    record.truth = truth_matcher_.match(
        cluster, towers, raw_truth_towers, truth, event_map, true);
    if (!record.truth.valid)
    {
      ++result.cluster_invalid_truth_count;
    }
    result.clusters.push_back(std::move(record));
  }

  result.prompt_cluster.assign(result.clusters.size(), 0U);
  if (pythia)
  {
    for (std::size_t index = 0; index < result.clusters.size(); ++index)
    {
      const auto& cluster = result.clusters[index];
      if (!cluster.anchor_acceptance || !cluster.truth.valid ||
          cluster.truth.contributors.empty())
      {
        continue;
      }
      const auto& dominant = cluster.truth.contributors.front();
      if (dominant.embedding_id != config_.signal_embedding_id ||
          !dominant.hepmc_valid ||
          dominant.fraction < config_.dominant_fraction_min ||
          dominant.g4_pdg_id != 22 ||
          (dominant.photon_category != 1 && dominant.photon_category != 2))
      {
        continue;
      }
      double energy = 0.0;
      double pt = 0.0;
      double eta = 0.0;
      double phi = 0.0;
      if (finite_hepmc_kinematics(
              contributor_hepmc_particle(dominant, event_map),
              energy, pt, eta, phi) && std::abs(eta) < config_.truth_eta_max)
      {
        result.prompt_cluster[index] = 1U;
      }
    }
  }

  result.candidates.reserve(candidates.size());
  for (CandidateWork& candidate : candidates)
  {
    result.candidates.push_back(std::move(candidate.record));
  }

  std::map<int, std::size_t> g4_candidate_by_barcode;
  std::map<int, std::size_t> generator_candidate_by_barcode;
  std::map<int, std::size_t> single_candidate_by_track;
  for (std::size_t index = 0; index < result.candidates.size(); ++index)
  {
    const auto& candidate = result.candidates[index];
    if (candidate.pathway == Pi0Pathway::g4_primary_decay)
    {
      g4_candidate_by_barcode[candidate.parent_barcode] = index;
    }
    else if (candidate.pathway == Pi0Pathway::generator_decay)
    {
      generator_candidate_by_barcode[candidate.parent_barcode] = index;
    }
    else
    {
      single_candidate_by_track[candidate.g4_parent_track_id] = index;
    }
  }

  std::vector<std::vector<std::size_t>> anchors_by_candidate(
      result.candidates.size());
  constexpr double tie_tolerance = 1e-6;
  for (std::size_t cluster_index = 0;
       cluster_index < result.clusters.size(); ++cluster_index)
  {
    const auto& cluster = result.clusters[cluster_index];
    if (!cluster.anchor_acceptance || !cluster.truth.valid ||
        cluster.truth.contributors.empty())
    {
      continue;
    }
    std::map<std::size_t, double> fraction_by_candidate;
    double unmatched_max_fraction = 0.0;
    for (const auto& contributor : cluster.truth.contributors)
    {
      const std::size_t candidate_index = contributor_candidate_index(
          contributor, config_, event_map, g4_candidate_by_barcode,
          generator_candidate_by_barcode, single_candidate_by_track);
      if (candidate_index == invalid_index)
      {
        unmatched_max_fraction = std::max(
            unmatched_max_fraction, static_cast<double>(contributor.fraction));
      }
      else
      {
        fraction_by_candidate[candidate_index] += contributor.fraction;
      }
    }
    std::size_t best_candidate = invalid_index;
    double best_fraction = -1.0;
    double second_fraction = -1.0;
    for (const auto& [candidate_index, fraction] : fraction_by_candidate)
    {
      if (fraction > best_fraction)
      {
        second_fraction = best_fraction;
        best_fraction = fraction;
        best_candidate = candidate_index;
      }
      else
      {
        second_fraction = std::max(second_fraction, fraction);
      }
    }
    if (best_candidate == invalid_index ||
        best_fraction < config_.anchor_pi0_fraction_min)
    {
      continue;
    }
    Pi0TopologyAnchorRecord anchor;
    anchor.cluster_index = cluster_index;
    anchor.candidate_index = best_candidate;
    anchor.main_fraction = best_fraction;
    anchor.second_fraction = second_fraction;
    anchor.unmatched_max_fraction = unmatched_max_fraction;
    anchor.ambiguous_main =
        second_fraction >= best_fraction - tie_tolerance ||
        unmatched_max_fraction >= best_fraction - tie_tolerance;
    if (anchor.ambiguous_main)
    {
      anchor.reason = Pi0AnchorReason::ambiguous_main_contributor;
    }
    result.anchors.push_back(anchor);
    anchors_by_candidate[best_candidate].push_back(result.anchors.size() - 1U);
  }

  for (std::size_t candidate_index = 0;
       candidate_index < result.candidates.size(); ++candidate_index)
  {
    if (anchors_by_candidate[candidate_index].empty() &&
        !config_.evaluate_all_candidates)
    {
      continue;
    }
    auto& candidate = result.candidates[candidate_index];
    candidate.topology_evaluated = true;
    candidate.cluster_matches.assign(
        result.clusters.size(), Pi0ClusterTruthMatch{});
    for (std::size_t cluster_index = 0;
         cluster_index < result.clusters.size(); ++cluster_index)
    {
      const auto match = pi0_truth_matcher_.match(
          result.clusters[cluster_index].cluster, towers, raw_truth_towers,
          truth_cells, truth_hits, truth, candidate.photon_track_ids, true);
      candidate.cluster_matches[cluster_index] = match;
      if (!match.valid)
      {
        ++result.energy_match_invalid_count;
        continue;
      }
      for (std::size_t photon = 0; photon < 2U; ++photon)
      {
        const double deposit = match.gamma_edep[photon];
        const double fraction = match.total_edep > 0.0F
            ? deposit / match.total_edep : 0.0;
        if (deposit > 0.0 &&
            fraction > config_.min_energy_contribution_fraction &&
            deposit > candidate.maximum_edep[photon])
        {
          candidate.best_cluster[photon] = cluster_index;
          candidate.maximum_edep[photon] = deposit;
          candidate.reconstructed_photon_energy[photon] =
              result.clusters[cluster_index].energy * fraction;
        }
      }
    }
    for (std::size_t photon = 0; photon < 2U; ++photon)
    {
      candidate.recovered[photon] =
          candidate.best_cluster[photon] != invalid_index &&
          std::isfinite(candidate.photon_energy[photon]) &&
          std::isfinite(candidate.reconstructed_photon_energy[photon]) &&
          candidate.photon_energy[photon] > 0.0 &&
          candidate.reconstructed_photon_energy[photon] /
              candidate.photon_energy[photon] >=
                  config_.min_photon_energy_recovery;
    }

    for (const std::size_t anchor_position :
         anchors_by_candidate[candidate_index])
    {
      auto& anchor = result.anchors[anchor_position];
      if (anchor.ambiguous_main)
      {
        anchor.topology = Pi0AnchorTopology::other;
        anchor.reason = Pi0AnchorReason::ambiguous_main_contributor;
        continue;
      }
      const bool is_best0 = candidate.recovered[0] &&
          candidate.best_cluster[0] == anchor.cluster_index;
      const bool is_best1 = candidate.recovered[1] &&
          candidate.best_cluster[1] == anchor.cluster_index;
      if (is_best0 && is_best1)
      {
        anchor.topology = Pi0AnchorTopology::merged;
        anchor.reason = Pi0AnchorReason::merged_shared_recovered_cluster;
      }
      else if ((is_best0 && candidate.recovered[1] &&
                candidate.best_cluster[1] != anchor.cluster_index) ||
               (is_best1 && candidate.recovered[0] &&
                candidate.best_cluster[0] != anchor.cluster_index))
      {
        anchor.topology = Pi0AnchorTopology::separated;
        anchor.reason = Pi0AnchorReason::separated_distinct_recovered_clusters;
      }
      else if ((is_best0 && !candidate.recovered[1]) ||
               (is_best1 && !candidate.recovered[0]))
      {
        anchor.topology = Pi0AnchorTopology::missing;
        anchor.reason = Pi0AnchorReason::missing_unrecovered_partner;
      }
      else
      {
        anchor.topology = Pi0AnchorTopology::other;
        const bool anchor_is_unrecovered_best =
            (!candidate.recovered[0] &&
             candidate.best_cluster[0] == anchor.cluster_index) ||
            (!candidate.recovered[1] &&
             candidate.best_cluster[1] == anchor.cluster_index);
        anchor.reason = anchor_is_unrecovered_best
            ? Pi0AnchorReason::other_best_cluster_below_recovery
            : Pi0AnchorReason::other_not_daughter_maximum;
      }
    }
  }

  result.status = Pi0TopologyEventStatus::accepted;
  return result;
}
}
