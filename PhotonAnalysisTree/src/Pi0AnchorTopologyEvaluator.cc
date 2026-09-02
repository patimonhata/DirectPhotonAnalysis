#include "Pi0AnchorTopologyEvaluator.h"

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfoContainer.h>
#include <g4detectors/PHG4CellContainer.h>
#include <g4main/PHG4Hit.h>
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

struct Projection
{
  bool valid = false;
  double eta = -999.0;
  double phi = -999.0;
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

const HepMC::GenParticle* contributor_hepmc_particle(const photon_tree::TruthContributor& contributor, const PHHepMCGenEventMap* event_map)
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
  const HepMC::GenParticle* particle = contributor_hepmc_particle(contributor, event_map);
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
  const auto found = generator_candidate_by_barcode.find(origin.parent->barcode());
  return found == generator_candidate_by_barcode.end()
      ? invalid_index : found->second;
}

double cemc_radius(RawTowerGeomContainer* geometry)
{
  constexpr double fallback = 95.0;
  if (!geometry) return fallback;
  for (int ieta = 0; ieta < 96; ++ieta)
  {
    for (int iphi = 0; iphi < 256; ++iphi)
    {
      const unsigned int key = RawTowerDefs::encode_towerid(RawTowerDefs::CEMC, static_cast<unsigned int>(ieta), static_cast<unsigned int>(iphi));
      const RawTowerGeom* tower = geometry->get_tower_geometry(key);
      if (tower)
      {
        const double radius = std::hypot(tower->get_center_x(), tower->get_center_y());
        return radius > 0.0 && std::isfinite(radius) ? radius : fallback;
      }
    }
  }
  return fallback;
}

Projection project_photon(const PHG4Particle* photon, PHG4TruthInfoContainer* truth, double target_radius)
{
  Projection result;
  if (!photon || !truth || !(target_radius > 0.0)) return result;
  const PHG4VtxPoint* vertex = truth->GetVtx(photon->get_vtx_id());
  if (!vertex) return result;
  const double px = photon->get_px();
  const double py = photon->get_py();
  const double pz = photon->get_pz();
  const double momentum = std::sqrt(px * px + py * py + pz * pz);
  if (!(momentum > 0.0) || !std::isfinite(vertex->get_x()) || !std::isfinite(vertex->get_y()) || !std::isfinite(vertex->get_z())) return result;
  const double ux = px / momentum;
  const double uy = py / momentum;
  const double uz = pz / momentum;
  const double a = ux * ux + uy * uy;
  const double b = 2.0 * (vertex->get_x() * ux + vertex->get_y() * uy);
  const double c = vertex->get_x() * vertex->get_x() + vertex->get_y() * vertex->get_y() - target_radius * target_radius;
  const double discriminant = b * b - 4.0 * a * c;
  if (a <= std::numeric_limits<double>::epsilon() || discriminant < 0.0) return result;
  const double root = std::sqrt(discriminant);
  const double first = (-b + root) / (2.0 * a);
  const double second = (-b - root) / (2.0 * a);
  double distance = std::numeric_limits<double>::infinity();
  if (first > 0.0) distance = std::min(distance, first);
  if (second > 0.0) distance = std::min(distance, second);
  if (!std::isfinite(distance)) return result;
  const double x = vertex->get_x() + distance * ux;
  const double y = vertex->get_y() + distance * uy;
  const double z = vertex->get_z() + distance * uz;
  const double radius = std::hypot(x, y);
  if (!(radius > 0.0) || !std::isfinite(z)) return result;
  result.eta = std::asinh(z / radius);
  result.phi = std::atan2(y, x);
  result.valid = std::isfinite(result.eta) && std::isfinite(result.phi);
  return result;
}

void fill_photon_kinematics(photon_tree::Pi0TopologyCandidateRecord& record, std::size_t index, const PHG4Particle* photon,
                            PHG4TruthInfoContainer* truth, double target_radius, double acceptance_eta_max)
{
  record.photon_track_ids[index] = photon ? photon->get_track_id() : -999;
  double pt = 0.0;
  if (!finite_g4_kinematics(photon, record.photon_energy[index], pt, record.photon_eta[index], record.photon_phi[index]))
  {
    record.photon_energy[index] = photon ? photon->get_e() : 0.0;
    record.photon_eta[index] = 0.0;
    record.photon_phi[index] = 0.0;
  }
  const Projection projection = project_photon(photon, truth, target_radius);
  record.photon_projection_valid[index] = projection.valid;
  if (projection.valid)
  {
    record.photon_projection_eta[index] = projection.eta;
    record.photon_projection_phi[index] = projection.phi;
    record.photon_in_cemc_acceptance[index] = std::abs(projection.eta) < acceptance_eta_max;
  }
}

void fill_first_daughter_diagnostic(photon_tree::Pi0TopologyCandidateRecord& record, std::size_t index,
                                    PHG4TruthInfoContainer* truth,
                                    const std::map<int, std::vector<const PHG4Particle*>>& children_by_parent,
                                    double pre_cemc_interaction_radius)
{
  const auto found = children_by_parent.find(record.photon_track_ids[index]);
  if (found == children_by_parent.end()) return;
  const PHG4VtxPoint* first_vertex = nullptr;
  double first_time = std::numeric_limits<double>::infinity();
  for (const PHG4Particle* child : found->second)
  {
    const PHG4VtxPoint* vertex = child && truth ? truth->GetVtx(child->get_vtx_id()) : nullptr;
    if (!vertex || !std::isfinite(vertex->get_t()) || vertex->get_t() >= first_time) continue;
    first_vertex = vertex;
    first_time = vertex->get_t();
  }
  if (!first_vertex) return;
  const double radius = std::hypot(first_vertex->get_x(), first_vertex->get_y());
  if (!std::isfinite(radius)) return;
  record.photon_first_daughter_vertex_valid[index] = true;
  record.photon_first_daughter_radius[index] = radius;
  record.photon_pre_cemc_interaction[index] = radius < pre_cemc_interaction_radius;
}
void fill_cemc_edep(std::vector<photon_tree::Pi0TopologyCandidateRecord>& candidates,
                    PHG4HitContainer* hits, PHG4TruthInfoContainer* truth)
{
  if (!hits || !truth || candidates.empty()) return;
  std::map<int, std::vector<std::pair<std::size_t, std::size_t>>> owners;
  for (std::size_t candidate = 0; candidate < candidates.size(); ++candidate)
    for (std::size_t photon = 0; photon < 2U; ++photon)
    {
      const int track_id = candidates[candidate].photon_track_ids[photon];
      if (track_id != 0 && track_id != -999) owners[track_id].push_back({candidate, photon});
    }

  const PHG4HitContainer::ConstRange range = hits->getHits();
  for (auto iterator = range.first; iterator != range.second; ++iterator)
  {
    const PHG4Hit* hit = iterator->second;
    const double edep = hit ? hit->get_edep() : 0.0;
    if (!(edep > 0.0) || !std::isfinite(edep)) continue;
    int track_id = hit->get_trkid();
    std::set<int> visited;
    while (track_id != 0 && visited.insert(track_id).second)
    {
      const auto owner = owners.find(track_id);
      if (owner != owners.end())
      {
        for (const auto& [candidate, photon] : owner->second) candidates[candidate].photon_cemc_edep[photon] += edep;
        break;
      }
      const PHG4Particle* particle = truth->GetParticle(track_id);
      if (!particle) break;
      track_id = particle->get_parent_id();
    }
  }
}

double diphoton_invariant_mass(const photon_tree::Pi0TopologyClusterRecord& anchor,
                               const photon_tree::Pi0PartnerDiagnosticRecord& partner)
{
  if (!(anchor.energy > 0.0) || !(partner.cluster_energy > 0.0)) return -1.0;
  const double delta_phi = std::atan2(std::sin(anchor.phi - partner.cluster_phi), std::cos(anchor.phi - partner.cluster_phi));
  const double anchor_pt = anchor.energy / std::cosh(anchor.eta);
  const double partner_pt = partner.cluster_energy / std::cosh(partner.cluster_eta);
  const double mass_squared = 2.0 * anchor_pt * partner_pt *
      (std::cosh(anchor.eta - partner.cluster_eta) - std::cos(delta_phi));
  return mass_squared >= 0.0 && std::isfinite(mass_squared) ? std::sqrt(mass_squared) : -1.0;
}

double diphoton_invariant_mass(const photon_tree::Pi0TopologyClusterRecord& anchor,
                               const photon_tree::Pi0TruthPartnerClusterRecord& partner)
{
  if (!(anchor.energy > 0.0) || !(partner.cluster_energy > 0.0)) return -1.0;
  const double delta_phi = std::atan2(std::sin(anchor.phi - partner.cluster_phi), std::cos(anchor.phi - partner.cluster_phi));
  const double anchor_pt = anchor.energy / std::cosh(anchor.eta);
  const double partner_pt = partner.cluster_energy / std::cosh(partner.cluster_eta);
  const double mass_squared = 2.0 * anchor_pt * partner_pt *
      (std::cosh(anchor.eta - partner.cluster_eta) - std::cos(delta_phi));
  return mass_squared >= 0.0 && std::isfinite(mass_squared) ? std::sqrt(mass_squared) : -1.0;
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
  case Pi0AnchorTopology::single_contaminated: return "single_contaminated";
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
  case Pi0AnchorReason::single_contaminated_pre_cemc_partner:
    return "single_contaminated_pre_cemc_partner";
  }
  return "unknown";
}

const char* pi0_missing_category_name(Pi0MissingCategory value)
{
  switch (value)
  {
  case Pi0MissingCategory::not_missing: return "not_missing";
  case Pi0MissingCategory::energy_threshold: return "energy_threshold";
  case Pi0MissingCategory::acceptance: return "acceptance";
  case Pi0MissingCategory::other: return "other";
  case Pi0MissingCategory::displaced_partner_cluster: return "displaced_partner_cluster";
  case Pi0MissingCategory::no_cemc_deposit: return "no_cemc_deposit";
  case Pi0MissingCategory::unclustered_deposit: return "unclustered_deposit";
  case Pi0MissingCategory::match_incomplete: return "match_incomplete";
  }
  return "unknown";
}

const char* pi0_partner_alignment_name(Pi0PartnerAlignment value)
{
  switch (value)
  {
  case Pi0PartnerAlignment::not_applicable: return "not_applicable";
  case Pi0PartnerAlignment::near: return "near";
  case Pi0PartnerAlignment::displaced: return "displaced";
  case Pi0PartnerAlignment::projection_invalid: return "projection_invalid";
  case Pi0PartnerAlignment::cluster_unavailable: return "cluster_unavailable";
  }
  return "unknown";
}

const char* pi0_truth_partner_tag_status_name(Pi0TruthPartnerTagStatus value)
{
  switch (value)
  {
  case Pi0TruthPartnerTagStatus::not_applicable: return "not_applicable";
  case Pi0TruthPartnerTagStatus::taggable: return "taggable";
  case Pi0TruthPartnerTagStatus::cluster_unavailable: return "cluster_unavailable";
  case Pi0TruthPartnerTagStatus::same_as_anchor: return "same_as_anchor";
  case Pi0TruthPartnerTagStatus::below_energy_threshold: return "below_energy_threshold";
  case Pi0TruthPartnerTagStatus::mass_outside_window: return "mass_outside_window";
  case Pi0TruthPartnerTagStatus::invalid_mass: return "invalid_mass";
  }
  return "unknown";
}

const char* pi0_anchor_tag_result_name(Pi0AnchorTagResult value)
{
  switch (value)
  {
  case Pi0AnchorTagResult::not_applicable: return "not_applicable";
  case Pi0AnchorTagResult::survived: return "survived";
  case Pi0AnchorTagResult::truth_pair_taggable_veto: return "truth_pair_taggable_veto";
  case Pi0AnchorTagResult::combinatorial_only_veto: return "combinatorial_only_veto";
  }
  return "unknown";
}

const char* pi0_missing_detail_name(Pi0MissingDetail value)
{
  switch (value)
  {
  case Pi0MissingDetail::not_missing: return "not_missing";
  case Pi0MissingDetail::partner_best_below_recovery: return "partner_best_below_recovery";
  case Pi0MissingDetail::partner_cluster_below_energy_threshold_recovered: return "partner_cluster_below_energy_threshold_recovered";
  case Pi0MissingDetail::partner_cluster_below_energy_threshold_below_recovery: return "partner_cluster_below_energy_threshold_below_recovery";
  case Pi0MissingDetail::partner_direct_match_incomplete: return "partner_direct_match_incomplete";
  case Pi0MissingDetail::partner_no_cemc_deposit: return "partner_no_cemc_deposit";
  case Pi0MissingDetail::partner_outside_cemc_acceptance: return "partner_outside_cemc_acceptance";
  case Pi0MissingDetail::partner_projection_invalid: return "partner_projection_invalid";
  case Pi0MissingDetail::partner_displaced_cluster_below_energy_threshold_recovered: return "partner_displaced_cluster_below_energy_threshold_recovered";
  case Pi0MissingDetail::partner_displaced_cluster_below_energy_threshold_below_recovery: return "partner_displaced_cluster_below_energy_threshold_below_recovery";
  case Pi0MissingDetail::partner_unclustered_cemc_deposit: return "partner_unclustered_cemc_deposit";
  case Pi0MissingDetail::partner_diagnostics_disabled: return "partner_diagnostics_disabled";
  }
  return "unknown";
}

void Pi0AnchorTopologyEvaluator::configure(const Pi0AnchorTopologyConfig& config)
{
  config_ = config;
  truth_matcher_.set_verbosity(config_.verbosity);
}

Pi0AnchorTopologyEventResult Pi0AnchorTopologyEvaluator::evaluate(PHCompositeNode* topNode)
{
  Pi0AnchorTopologyEventResult result;
  auto* truth = findNode::getClass<PHG4TruthInfoContainer>(topNode, config_.truth_node_name);
  auto* event_map = findNode::getClass<PHHepMCGenEventMap>(topNode, config_.hepmc_event_map_node_name);
  auto* towers = findNode::getClass<TowerInfoContainer>(topNode, config_.tower_node_name);
  auto* raw_truth_towers = findNode::getClass<RawTowerContainer>(topNode, config_.raw_truth_tower_node_name);
  auto* truth_cells = findNode::getClass<PHG4CellContainer>(topNode, config_.truth_cell_node_name);
  auto* truth_hits = findNode::getClass<PHG4HitContainer>(topNode, config_.truth_hit_node_name);
  auto* clusters = findNode::getClass<RawClusterContainer>(topNode, config_.cluster_node_name);
  auto* geometry = findNode::getClass<RawTowerGeomContainer>(topNode, config_.tower_geom_node_name);

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
  const double projection_radius = cemc_radius(geometry);

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
      const HepMC::GenParticle* hepmc_particle = event->barcode_to_particle(primary->get_barcode());
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
                                  work.record.eta, work.record.phi))
        {
          continue;
        }
        fill_photon_kinematics(work.record, 0, found->second[0], truth, projection_radius, config_.cemc_acceptance_eta_max);
        fill_photon_kinematics(work.record, 1, found->second[1], truth, projection_radius, config_.cemc_acceptance_eta_max);
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
              work.record.eta, work.record.phi))
      {
        continue;
      }
      fill_photon_kinematics(work.record, 0, pending.photons[0], truth, projection_radius, config_.cemc_acceptance_eta_max);
      fill_photon_kinematics(work.record, 1, pending.photons[1], truth, projection_radius, config_.cemc_acceptance_eta_max);
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
                                work.record.eta, work.record.phi))
      {
        continue;
      }
      fill_photon_kinematics(work.record, 0, found->second[0], truth, projection_radius, config_.cemc_acceptance_eta_max);
      fill_photon_kinematics(work.record, 1, found->second[1], truth, projection_radius, config_.cemc_acceptance_eta_max);
      candidates.push_back(std::move(work));
      ++result.g4_candidate_count;
    }
  }

  struct DiagnosticCluster
  {
    const RawCluster* cluster = nullptr;
    unsigned int id = 0;
    double energy = 0.0;
    double eta = 0.0;
    double phi = 0.0;
  };
  std::vector<DiagnosticCluster> below_threshold_clusters;
  const auto cluster_range = clusters->getClusters();
  for (auto iterator = cluster_range.first; iterator != cluster_range.second;
       ++iterator)
  {
    const RawCluster* cluster = iterator->second;
    if (!cluster || !std::isfinite(cluster->get_energy()) ||
        !std::isfinite(cluster->get_x()) || !std::isfinite(cluster->get_y()) ||
        !std::isfinite(cluster->get_z()))
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
    record.phi = std::atan2(dy, dx);
    if (!std::isfinite(record.et) || !(record.et > 0.0) ||
        !std::isfinite(record.eta) || !std::isfinite(record.phi))
    {
      continue;
    }
    if (config_.partner_cluster_eta_max > 0.0 &&
        std::abs(record.eta) >= config_.partner_cluster_eta_max)
    {
      continue;
    }
    if (record.energy < config_.min_cluster_energy)
    {
      if (config_.enable_missing_diagnostics && record.energy > config_.partner_diagnostic_min_cluster_energy)
      {
        below_threshold_clusters.push_back({cluster, record.cluster_id, record.energy, record.eta, record.phi});
      }
      continue;
    }
    record.anchor_acceptance = std::abs(record.eta) < config_.anchor_cluster_eta_max;
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
              energy, pt, eta, phi))
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
  for (auto& candidate : result.candidates)
  {
    for (std::size_t photon = 0; photon < 2U; ++photon)
      fill_first_daughter_diagnostic(candidate, photon, truth, children_by_parent, config_.pre_cemc_interaction_radius);
  }
  if (config_.enable_missing_diagnostics)
  {
    fill_cemc_edep(result.candidates, truth_hits, truth);
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
        unmatched_max_fraction = std::max(unmatched_max_fraction, static_cast<double>(contributor.fraction));
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

  const auto update_truth_partner = [&](Pi0TopologyCandidateRecord& candidate, std::size_t photon,
                                        unsigned int cluster_id, double cluster_energy, double cluster_eta,
                                        double cluster_phi, bool below_topology_threshold, const Pi0ClusterTruthMatch& match) {
    if (!(cluster_energy > config_.partner_diagnostic_min_cluster_energy) || !match.usable || !(match.total_edep > 0.0F)) return;
    const double direct_edep = match.gamma_edep[photon];
    const double fraction = direct_edep / match.total_edep;
    auto& partner = candidate.truth_partner_clusters[photon];
    if (!(direct_edep > 0.0) || !(fraction > config_.min_energy_contribution_fraction) ||
        (partner.found && !(direct_edep > partner.direct_edep))) return;
    partner.found = true;
    partner.below_topology_threshold = below_topology_threshold;
    partner.cluster_id = cluster_id;
    partner.cluster_energy = cluster_energy;
    partner.cluster_eta = cluster_eta;
    partner.cluster_phi = cluster_phi;
    partner.direct_edep = direct_edep;
    partner.reconstructed_photon_energy = cluster_energy * fraction;
    partner.recovery = candidate.photon_energy[photon] > 0.0
        ? partner.reconstructed_photon_energy / candidate.photon_energy[photon] : 0.0;
    if (candidate.photon_projection_valid[photon])
    {
      const double delta_phi = std::atan2(std::sin(cluster_phi - candidate.photon_projection_phi[photon]),
                                          std::cos(cluster_phi - candidate.photon_projection_phi[photon]));
      partner.delta_r = std::hypot(cluster_eta - candidate.photon_projection_eta[photon], delta_phi);
    }
    else
    {
      partner.delta_r = -1.0;
    }
  };

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
    candidate.cluster_matches.assign(result.clusters.size(), Pi0ClusterTruthMatch{});
    for (std::size_t cluster_index = 0;
         cluster_index < result.clusters.size(); ++cluster_index)
    {
      auto match = pi0_truth_matcher_.match(
          result.clusters[cluster_index].cluster, towers, raw_truth_towers,
          truth_cells, truth_hits, truth, candidate.photon_track_ids, true);
      if (match.status == Pi0ClusterTruthMatchStatus::partial &&
          match.cluster_member_energy_coverage >= config_.min_direct_match_cluster_energy_coverage)
      {
        match.usable = true;
      }
      candidate.cluster_matches[cluster_index] = match;
      if (!match.usable)
      {
        ++result.energy_match_invalid_count;
        continue;
      }
      for (std::size_t photon = 0; photon < 2U; ++photon)
      {
        const double deposit = match.gamma_edep[photon];
        const auto& cluster = result.clusters[cluster_index];
        update_truth_partner(candidate, photon, cluster.cluster_id, cluster.energy, cluster.eta, cluster.phi, false, match);
        const double fraction = match.total_edep > 0.0F
            ? deposit / match.total_edep : 0.0;
        if (deposit > 0.0 &&
            fraction > config_.min_energy_contribution_fraction &&
            deposit > candidate.maximum_edep[photon])
        {
          candidate.best_cluster[photon] = cluster_index;
          candidate.maximum_edep[photon] = deposit;
          candidate.reconstructed_photon_energy[photon] = result.clusters[cluster_index].energy * fraction;
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
              candidate.photon_energy[photon] >= config_.min_photon_energy_recovery;
    }

  }

  const auto update_diagnostic = [&](Pi0TopologyCandidateRecord& candidate, std::size_t photon,
                                     unsigned int cluster_id, double cluster_energy,
                                     double cluster_eta, double cluster_phi, bool below_threshold,
                                     const Pi0ClusterTruthMatch& match) {
    if (candidate.recovered[photon]) return;
    const double reference_eta = candidate.photon_projection_valid[photon]
        ? candidate.photon_projection_eta[photon] : candidate.photon_eta[photon];
    const double reference_phi = candidate.photon_projection_valid[photon]
        ? candidate.photon_projection_phi[photon] : candidate.photon_phi[photon];
    const double delta_phi = std::atan2(std::sin(cluster_phi - reference_phi), std::cos(cluster_phi - reference_phi));
    const double delta_r = std::hypot(cluster_eta - reference_eta, delta_phi);
    const double fraction = match.total_edep > 0.0F ? match.gamma_edep[photon] / match.total_edep : 0.0;
    const bool has_direct_deposit = match.usable && match.gamma_edep[photon] > 0.0F &&
        fraction > config_.min_energy_contribution_fraction;
    if (!has_direct_deposit && delta_r > config_.missing_diagnostic_max_delta_r) return;
    auto& diagnostic = candidate.partner_diagnostics[photon];
    const int priority = has_direct_deposit ? 2 : (!match.usable ? 1 : 0);
    const int current_priority = !diagnostic.found ? -1
        : (diagnostic.has_direct_deposit ? 2 : (!diagnostic.match.usable ? 1 : 0));
    const bool replace = !diagnostic.found || priority > current_priority ||
        (priority == 2 && current_priority == 2 && match.gamma_edep[photon] > diagnostic.match.gamma_edep[photon]) ||
        (priority < 2 && priority == current_priority && delta_r < diagnostic.delta_r);
    if (!replace) return;
    diagnostic.found = true;
    diagnostic.below_energy_threshold = below_threshold;
    diagnostic.has_direct_deposit = has_direct_deposit;
    diagnostic.cluster_id = cluster_id;
    diagnostic.cluster_energy = cluster_energy;
    diagnostic.cluster_eta = cluster_eta;
    diagnostic.cluster_phi = cluster_phi;
    diagnostic.delta_r = delta_r;
    diagnostic.match = match;
    diagnostic.reconstructed_photon_energy = has_direct_deposit ? cluster_energy * fraction : 0.0;
    diagnostic.recovery = candidate.photon_energy[photon] > 0.0
        ? diagnostic.reconstructed_photon_energy / candidate.photon_energy[photon] : 0.0;
  };

  if (config_.enable_missing_diagnostics)
  {
    std::vector<std::size_t> diagnostic_candidate_indices;
    std::vector<std::array<int, 2>> diagnostic_track_ids;
    for (std::size_t candidate_index = 0; candidate_index < result.candidates.size(); ++candidate_index)
    {
      const auto& candidate = result.candidates[candidate_index];
      if (!candidate.topology_evaluated) continue;
      diagnostic_candidate_indices.push_back(candidate_index);
      diagnostic_track_ids.push_back(candidate.photon_track_ids);
    }

    for (const std::size_t candidate_index : diagnostic_candidate_indices)
    {
      auto& candidate = result.candidates[candidate_index];
      for (std::size_t cluster_index = 0; cluster_index < result.clusters.size(); ++cluster_index)
      {
        const auto& match = candidate.cluster_matches[cluster_index];
        if (match.usable) continue;
        const auto& cluster = result.clusters[cluster_index];
        for (std::size_t photon = 0; photon < 2U; ++photon)
          update_diagnostic(candidate, photon, cluster.cluster_id, cluster.energy, cluster.eta, cluster.phi, false, match);
      }
    }

    for (const auto& cluster : below_threshold_clusters)
    {
      auto matches = pi0_truth_matcher_.match_many(
          cluster.cluster, towers, raw_truth_towers, truth_cells, truth_hits, truth, diagnostic_track_ids, true);
      for (std::size_t position = 0; position < matches.size(); ++position)
      {
        auto& match = matches[position];
        if (match.status == Pi0ClusterTruthMatchStatus::partial &&
            match.cluster_member_energy_coverage >= config_.min_direct_match_cluster_energy_coverage)
        {
          match.usable = true;
        }
        auto& candidate = result.candidates[diagnostic_candidate_indices[position]];
        for (std::size_t photon = 0; photon < 2U; ++photon)
        {
          update_truth_partner(candidate, photon, cluster.id, cluster.energy, cluster.eta, cluster.phi, true, match);
          update_diagnostic(candidate, photon, cluster.id, cluster.energy, cluster.eta, cluster.phi, true, match);
        }
      }
    }
  }

  for (std::size_t candidate_index = 0; candidate_index < result.candidates.size(); ++candidate_index)
  {
    if (anchors_by_candidate[candidate_index].empty() && !config_.evaluate_all_candidates) continue;
    auto& candidate = result.candidates[candidate_index];
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
        const bool pre_cemc0 = candidate.photon_pre_cemc_interaction[0];
        const bool pre_cemc1 = candidate.photon_pre_cemc_interaction[1];
        if (pre_cemc0 != pre_cemc1)
        {
          anchor.topology = Pi0AnchorTopology::single_contaminated;
          anchor.reason = Pi0AnchorReason::single_contaminated_pre_cemc_partner;
          anchor.pre_cemc_photon_index = pre_cemc0 ? 0 : 1;
        }
        else
        {
          anchor.topology = Pi0AnchorTopology::merged;
          anchor.reason = Pi0AnchorReason::merged_shared_recovered_cluster;
        }
      }
      else if ((is_best0 && candidate.recovered[1] &&
                candidate.best_cluster[1] != anchor.cluster_index) ||
               (is_best1 && candidate.recovered[0] &&
                candidate.best_cluster[0] != anchor.cluster_index))
      {
        anchor.topology = Pi0AnchorTopology::separated;
        anchor.reason = Pi0AnchorReason::separated_distinct_recovered_clusters;
        anchor.partner_photon_index = is_best0 ? 1 : 0;
      }
      else if ((is_best0 && !candidate.recovered[1]) ||
               (is_best1 && !candidate.recovered[0]))
      {
        anchor.topology = Pi0AnchorTopology::missing;
        anchor.reason = Pi0AnchorReason::missing_unrecovered_partner;
        anchor.partner_photon_index = is_best0 ? 1 : 0;
        const std::size_t partner = static_cast<std::size_t>(anchor.partner_photon_index);
        const auto& diagnostic = candidate.partner_diagnostics[partner];
        const bool projection_valid = candidate.photon_projection_valid[partner];
        const bool has_cemc_deposit = candidate.photon_cemc_edep[partner] > 0.0 ||
            candidate.best_cluster[partner] != invalid_index || diagnostic.has_direct_deposit;
        if (diagnostic.found && diagnostic.has_direct_deposit)
        {
          anchor.partner_diagnostic_invariant_mass =
              diphoton_invariant_mass(result.clusters[anchor.cluster_index], diagnostic);
        }
        if (!projection_valid)
        {
          anchor.missing_category = Pi0MissingCategory::other;
          anchor.missing_detail = Pi0MissingDetail::partner_projection_invalid;
        }
        else if (!candidate.photon_in_cemc_acceptance[partner])
        {
          anchor.missing_category = Pi0MissingCategory::acceptance;
          anchor.missing_detail = Pi0MissingDetail::partner_outside_cemc_acceptance;
        }
        else if (!config_.enable_missing_diagnostics)
        {
          anchor.missing_category = Pi0MissingCategory::other;
          anchor.missing_detail = Pi0MissingDetail::partner_diagnostics_disabled;
        }
        else if (!has_cemc_deposit)
        {
          anchor.missing_category = Pi0MissingCategory::no_cemc_deposit;
          anchor.missing_detail = Pi0MissingDetail::partner_no_cemc_deposit;
        }
        else if (diagnostic.found && diagnostic.below_energy_threshold && diagnostic.has_direct_deposit)
        {
          const bool displaced = diagnostic.delta_r > config_.missing_diagnostic_max_delta_r;
          anchor.missing_category = displaced
              ? Pi0MissingCategory::displaced_partner_cluster
              : Pi0MissingCategory::energy_threshold;
          if (displaced)
          {
            anchor.missing_detail = diagnostic.recovery >= config_.min_photon_energy_recovery
                ? Pi0MissingDetail::partner_displaced_cluster_below_energy_threshold_recovered
                : Pi0MissingDetail::partner_displaced_cluster_below_energy_threshold_below_recovery;
          }
          else
          {
            anchor.missing_detail = diagnostic.recovery >= config_.min_photon_energy_recovery
                ? Pi0MissingDetail::partner_cluster_below_energy_threshold_recovered
                : Pi0MissingDetail::partner_cluster_below_energy_threshold_below_recovery;
          }
        }
        else if (candidate.best_cluster[partner] != invalid_index)
        {
          anchor.missing_category = Pi0MissingCategory::other;
          anchor.missing_detail = Pi0MissingDetail::partner_best_below_recovery;
        }
        else if (diagnostic.found && !diagnostic.match.usable)
        {
          anchor.missing_category = Pi0MissingCategory::match_incomplete;
          anchor.missing_detail = Pi0MissingDetail::partner_direct_match_incomplete;
        }
        else
        {
          anchor.missing_category = Pi0MissingCategory::unclustered_deposit;
          anchor.missing_detail = Pi0MissingDetail::partner_unclustered_cemc_deposit;
        }
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

      if (anchor.partner_photon_index >= 0 && anchor.partner_photon_index < 2)
      {
        const std::size_t partner_index = static_cast<std::size_t>(anchor.partner_photon_index);
        const auto& partner = candidate.truth_partner_clusters[partner_index];
        if (!candidate.photon_projection_valid[partner_index])
        {
          anchor.partner_alignment = Pi0PartnerAlignment::projection_invalid;
        }
        else if (!partner.found)
        {
          anchor.partner_alignment = Pi0PartnerAlignment::cluster_unavailable;
        }
        else
        {
          anchor.partner_alignment = partner.delta_r > config_.missing_diagnostic_max_delta_r
              ? Pi0PartnerAlignment::displaced : Pi0PartnerAlignment::near;
        }

        if (!partner.found)
        {
          anchor.truth_partner_tag_status = Pi0TruthPartnerTagStatus::cluster_unavailable;
        }
        else
        {
          anchor.truth_partner_cluster_id = static_cast<int>(partner.cluster_id);
          anchor.truth_partner_cluster_energy = partner.cluster_energy;
          anchor.truth_partner_cluster_eta = partner.cluster_eta;
          anchor.truth_partner_cluster_phi = partner.cluster_phi;
          anchor.truth_partner_delta_r = partner.delta_r;
          anchor.truth_partner_direct_edep = partner.direct_edep;
          anchor.truth_partner_reconstructed_photon_energy = partner.reconstructed_photon_energy;
          anchor.truth_partner_recovery = partner.recovery;
          anchor.truth_partner_invariant_mass = diphoton_invariant_mass(result.clusters[anchor.cluster_index], partner);
          if (partner.cluster_id == result.clusters[anchor.cluster_index].cluster_id)
            anchor.truth_partner_tag_status = Pi0TruthPartnerTagStatus::same_as_anchor;
          else if (!(partner.cluster_energy > config_.tagging_partner_min_cluster_energy))
            anchor.truth_partner_tag_status = Pi0TruthPartnerTagStatus::below_energy_threshold;
          else if (!(anchor.truth_partner_invariant_mass >= 0.0) || !std::isfinite(anchor.truth_partner_invariant_mass))
            anchor.truth_partner_tag_status = Pi0TruthPartnerTagStatus::invalid_mass;
          else if (!(anchor.truth_partner_invariant_mass > config_.tagging_pi0_mass_min &&
                     anchor.truth_partner_invariant_mass < config_.tagging_pi0_mass_max))
            anchor.truth_partner_tag_status = Pi0TruthPartnerTagStatus::mass_outside_window;
          else
            anchor.truth_partner_tag_status = Pi0TruthPartnerTagStatus::taggable;
        }
      }
    }
  }

  result.status = Pi0TopologyEventStatus::accepted;
  return result;
}
}
