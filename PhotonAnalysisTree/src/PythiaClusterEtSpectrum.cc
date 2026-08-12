#include "PythiaClusterEtSpectrum.h"

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/TowerInfoContainer.h>
#include <fun4all/Fun4AllReturnCodes.h>
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
#include <HepMC/SimpleVector.h>

#include <TFile.h>
#include <TH1D.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <utility>
#include <vector>

namespace
{
enum class Pi0Pathway : int
{
  g4_primary_decay = 1,
  g4_secondary_decay = 2,
  generator_decay = 3
};

struct Pi0Origin
{
  bool valid = false;
  const HepMC::GenParticle* parent = nullptr;
};

struct Pi0Candidate
{
  Pi0Pathway pathway = Pi0Pathway::g4_primary_decay;
  int parent_barcode = 0;
  const HepMC::GenParticle* parent = nullptr;
  const PHG4Particle* g4_parent = nullptr;
  std::vector<const PHG4Particle*> photons;
};

struct ClusterRecord
{
  const RawCluster* cluster = nullptr;
  double et = 0.0;
  double eta = 0.0;
  double surface_eta = 0.0;
  double surface_phi = 0.0;
  photon_tree::ClusterTruthMatch truth;
};

struct Projection
{
  bool valid = false;
  double eta = 0.0;
  double phi = 0.0;
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
    const std::vector<const HepMC::GenParticle*> parents =
        incoming(current->production_vertex());
    if (parents.size() == 1U && parents.front()->pdg_id() == 22)
    {
      current = parents.front();
      if (current->barcode() != 0 && !visited.insert(current->barcode()).second)
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

bool finite_eta(const HepMC::GenParticle* particle, double& eta)
{
  if (!particle)
  {
    return false;
  }
  const HepMC::FourVector& momentum = particle->momentum();
  const double pt = std::hypot(momentum.px(), momentum.py());
  if (!(pt > 0.0) || !std::isfinite(momentum.pz()))
  {
    return false;
  }
  eta = std::asinh(momentum.pz() / pt);
  return std::isfinite(eta);
}

double particle_pt(const PHG4Particle* particle)
{
  return particle ? std::hypot(particle->get_px(), particle->get_py()) : 0.0;
}

double hepmc_pt(const HepMC::GenParticle* particle)
{
  return particle ? std::hypot(particle->momentum().px(), particle->momentum().py()) : 0.0;
}

double wrap_delta_phi(double value)
{
  constexpr double pi = 3.14159265358979323846;
  while (value > pi)
  {
    value -= 2.0 * pi;
  }
  while (value <= -pi)
  {
    value += 2.0 * pi;
  }
  return value;
}

double delta_r(double eta0, double phi0, double eta1, double phi1)
{
  return std::hypot(eta0 - eta1, wrap_delta_phi(phi0 - phi1));
}

double cemc_radius(RawTowerGeomContainer* geometry)
{
  constexpr double fallback = 95.0;
  if (!geometry)
  {
    return fallback;
  }
  for (int ieta = 0; ieta < 96; ++ieta)
  {
    for (int iphi = 0; iphi < 256; ++iphi)
    {
      const unsigned int key = RawTowerDefs::encode_towerid(
          RawTowerDefs::CEMC, static_cast<unsigned int>(ieta),
          static_cast<unsigned int>(iphi));
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

Projection project_photon(const PHG4Particle* photon,
                          PHG4TruthInfoContainer* truth,
                          double target_radius)
{
  Projection result;
  if (!photon || !truth || !(target_radius > 0.0))
  {
    return result;
  }
  const PHG4VtxPoint* vertex = truth->GetVtx(photon->get_vtx_id());
  if (!vertex)
  {
    return result;
  }
  const double px = photon->get_px();
  const double py = photon->get_py();
  const double pz = photon->get_pz();
  const double momentum = std::sqrt(px * px + py * py + pz * pz);
  if (!(momentum > 0.0) || !std::isfinite(vertex->get_x()) ||
      !std::isfinite(vertex->get_y()) || !std::isfinite(vertex->get_z()))
  {
    return result;
  }
  const double ux = px / momentum;
  const double uy = py / momentum;
  const double uz = pz / momentum;
  const double a = ux * ux + uy * uy;
  const double b = 2.0 * (vertex->get_x() * ux + vertex->get_y() * uy);
  const double c = vertex->get_x() * vertex->get_x() +
      vertex->get_y() * vertex->get_y() - target_radius * target_radius;
  const double discriminant = b * b - 4.0 * a * c;
  if (a <= std::numeric_limits<double>::epsilon() || discriminant < 0.0)
  {
    return result;
  }
  const double root = std::sqrt(discriminant);
  const double first = (-b + root) / (2.0 * a);
  const double second = (-b - root) / (2.0 * a);
  double distance = std::numeric_limits<double>::infinity();
  if (first > 0.0)
  {
    distance = std::min(distance, first);
  }
  if (second > 0.0)
  {
    distance = std::min(distance, second);
  }
  if (!std::isfinite(distance))
  {
    return result;
  }
  const double x = vertex->get_x() + distance * ux;
  const double y = vertex->get_y() + distance * uy;
  const double z = vertex->get_z() + distance * uz;
  const double radius = std::hypot(x, y);
  if (!(radius > 0.0) || !std::isfinite(z))
  {
    return result;
  }
  result.eta = std::asinh(z / radius);
  result.phi = std::atan2(y, x);
  result.valid = std::isfinite(result.eta) && std::isfinite(result.phi);
  return result;
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

bool contributor_matches_pi0(const photon_tree::TruthContributor& contributor,
                             Pi0Pathway pathway,
                             int parent_barcode,
                             const PHHepMCGenEventMap* event_map,
                             int signal_embedding_id)
{
  if (contributor.embedding_id != signal_embedding_id || !contributor.hepmc_valid)
  {
    return false;
  }
  const HepMC::GenParticle* particle =
      contributor_hepmc_particle(contributor, event_map);
  if (pathway == Pi0Pathway::g4_secondary_decay)
  {
    return contributor.g4_track_id == parent_barcode;
  }
  if (pathway == Pi0Pathway::g4_primary_decay)
  {
    return contributor.g4_pdg_id == 111 && particle &&
        particle->pdg_id() == 111 && particle->barcode() == parent_barcode;
  }
  const Pi0Origin origin = trace_pi0_origin(particle);
  return contributor.g4_pdg_id == 22 && origin.valid && origin.parent &&
      origin.parent->barcode() == parent_barcode;
}
}

PythiaClusterEtSpectrum::PythiaClusterEtSpectrum(const std::string& name)
  : SubsysReco(name)
{
}

PythiaClusterEtSpectrum::~PythiaClusterEtSpectrum()
{
  close_output();
}

int PythiaClusterEtSpectrum::Init(PHCompositeNode* /*topNode*/)
{
  const bool valid = !output_file_name_.empty() && !manifest_path_.empty() &&
      manifest_begin_ >= 0 && manifest_end_ > manifest_begin_ &&
      !first_suffix_.empty() && !last_suffix_.empty() && signal_embedding_id_ > 0 &&
      n_bins_ > 0 && std::isfinite(et_max_) && et_max_ > 0.0 &&
      std::isfinite(truth_eta_max_) && truth_eta_max_ > 0.0 &&
      std::isfinite(cluster_eta_max_) && cluster_eta_max_ > 0.0 &&
      std::isfinite(min_cluster_energy_) && min_cluster_energy_ >= 0.0 &&
      dominant_fraction_min_ >= 0.0 && dominant_fraction_min_ <= 1.0 &&
      pi0_contributor_fraction_min_ >= 0.0 && pi0_contributor_fraction_min_ <= 1.0 &&
      separated_delta_r_cut_ > 0.0 && merged_delta_r_cut_ >= separated_delta_r_cut_ &&
      response_min_ >= 0.0 && response_max_ > response_min_;
  if (!valid)
  {
    std::cerr << "PythiaClusterEtSpectrum::Init - invalid configuration" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  truth_matcher_.set_verbosity(verbosity_);
  create_output_directory();
  create_output();
  if (!output_file_ || output_file_->IsZombie())
  {
    std::cerr << "PythiaClusterEtSpectrum::Init - failed to create "
              << output_file_name_ << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int PythiaClusterEtSpectrum::process_event(PHCompositeNode* topNode)
{
  ++n_events_processed_;
  auto* truth = findNode::getClass<PHG4TruthInfoContainer>(topNode, truth_node_name_);
  auto* event_map = findNode::getClass<PHHepMCGenEventMap>(topNode, hepmc_event_map_node_name_);
  auto* towers = findNode::getClass<TowerInfoContainer>(topNode, tower_node_name_);
  auto* raw_truth_towers = findNode::getClass<RawTowerContainer>(topNode, raw_truth_tower_node_name_);
  auto* geometry = findNode::getClass<RawTowerGeomContainer>(topNode, tower_geom_node_name_);
  auto* clusters = findNode::getClass<RawClusterContainer>(topNode, split_cluster_node_name_);
  const PHHepMCGenEvent* signal_event = event_map ? event_map->get(signal_embedding_id_) : nullptr;
  const HepMC::GenEvent* event = signal_event ? signal_event->getEvent() : nullptr;
  if (!truth || !event_map || !towers || !geometry || !clusters || !signal_event || !signal_event->is_simulated() || !event || !truth_matcher_.begin_event(topNode)) {
    ++n_events_invalid_;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  const HepMC::FourVector& collision = signal_event->get_collision_vertex();
  if (!std::isfinite(collision.x()) || !std::isfinite(collision.y()) || !std::isfinite(collision.z())) {
    ++n_events_invalid_;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  std::vector<ClusterRecord> cluster_records;
  const auto cluster_range = clusters->getClusters();
  for (auto iterator = cluster_range.first; iterator != cluster_range.second; ++iterator) {
    const RawCluster* cluster = iterator->second;
    if (!cluster || !std::isfinite(cluster->get_energy()) || !std::isfinite(cluster->get_x()) || !std::isfinite(cluster->get_y()) || !std::isfinite(cluster->get_z()) ||
        cluster->get_energy() < min_cluster_energy_) {
      continue;
    }
    const double dx = cluster->get_x() - collision.x();
    const double dy = cluster->get_y() - collision.y();
    const double dz = cluster->get_z() - collision.z();
    const double distance = std::sqrt(dx * dx + dy * dy + dz * dz);
    const double transverse = std::hypot(dx, dy);
    const double surface_radius = std::hypot(cluster->get_x(), cluster->get_y());
    if (!(distance > 0.0) || !(transverse > 0.0) || !(surface_radius > 0.0)) {
      continue;
    }
    ClusterRecord record;
    record.cluster = cluster;
    record.et = cluster->get_energy() * transverse / distance;
    record.eta = std::asinh(dz / transverse);
    record.surface_eta = std::asinh(cluster->get_z() / surface_radius);
    record.surface_phi = std::atan2(cluster->get_y(), cluster->get_x());
    if (!std::isfinite(record.et) || !(record.et > 0.0) || !std::isfinite(record.eta) || std::abs(record.eta) >= cluster_eta_max_) {
      continue;
    }
    record.truth = truth_matcher_.match(cluster, towers, raw_truth_towers, truth, event_map, true);
    ++n_cluster_considered_;
    if (!record.truth.valid) {
      ++n_cluster_invalid_truth_;
    }
    cluster_records.push_back(std::move(record));
  }

  std::set<const RawCluster*> pi0_clusters_filled;
  for (const ClusterRecord& cluster : cluster_records) {
    if (!cluster.truth.valid || cluster.truth.contributors.empty()) {
      continue;
    }
    const photon_tree::TruthContributor& dominant = cluster.truth.contributors.front();
    if (dominant.embedding_id != signal_embedding_id_ || !dominant.hepmc_valid || dominant.fraction < dominant_fraction_min_) {
      continue;
    }
    const HepMC::GenParticle* particle = contributor_hepmc_particle(dominant, event_map);
    double truth_eta = 0.0;
    if (dominant.g4_pdg_id == 22 &&
        (dominant.photon_category == 1 || dominant.photon_category == 2) &&
        finite_eta(particle, truth_eta) && std::abs(truth_eta) < truth_eta_max_)
    {
      h_prompt_->Fill(cluster.et);
      ++n_prompt_cluster_;
    }
    bool from_pi0 = false;
    if (dominant.g4_pdg_id == 111 && particle && particle->pdg_id() == 111) {
      double parent_eta = 0.0;
      from_pi0 = finite_eta(particle, parent_eta) && std::abs(parent_eta) < truth_eta_max_;
      if (from_pi0) {
        ++n_pi0_cluster_g4_decay_;
      }
    }
    else if (dominant.g4_pdg_id == 22)
    {
      const Pi0Origin origin = trace_pi0_origin(particle);
      if (origin.valid && origin.parent) {
        double parent_eta = 0.0;
        from_pi0 = finite_eta(origin.parent, parent_eta) && std::abs(parent_eta) < truth_eta_max_;
        if (from_pi0) {
          ++n_pi0_cluster_generator_decay_;
        }
      }
    }
    if (from_pi0) {
      h_pi0_->Fill(cluster.et);
      pi0_clusters_filled.insert(cluster.cluster);
      ++n_pi0_cluster_;
    }
  }

  std::map<int, std::vector<const PHG4Particle*>> children_by_parent;
  const auto secondary_range = truth->GetSecondaryParticleRange();
  for (auto iterator = secondary_range.first; iterator != secondary_range.second; ++iterator) {
    const PHG4Particle* particle = iterator->second;
    if (particle) {
      children_by_parent[particle->get_parent_id()].push_back(particle);
    }
  }

  std::vector<Pi0Candidate> candidates;
  std::map<int, std::size_t> generator_candidate_by_barcode;
  const auto primary_range = truth->GetPrimaryParticleRange();
  for (auto iterator = primary_range.first; iterator != primary_range.second; ++iterator) {
    const PHG4Particle* primary = iterator->second;
    if (!primary || truth->isEmbeded(primary->get_track_id()) != signal_embedding_id_) {
      continue;
    }
    const HepMC::GenParticle* hepmc_particle = event->barcode_to_particle(primary->get_barcode());
    if (primary->get_pid() == 111 && hepmc_particle && hepmc_particle->pdg_id() == 111) {
      const auto found = children_by_parent.find(primary->get_track_id());
      if (found == children_by_parent.end() || found->second.size() != 2U || found->second[0]->get_pid() != 22 || found->second[1]->get_pid() != 22) {
        ++n_pi0_malformed_daughters_;
        continue;
      }
      candidates.push_back({Pi0Pathway::g4_primary_decay, hepmc_particle->barcode(), hepmc_particle, primary, found->second});
      continue;
    }
    if (primary->get_pid() != 22 || !hepmc_particle || hepmc_particle->pdg_id() != 22) {
      continue;
    }
    const Pi0Origin origin = trace_pi0_origin(hepmc_particle);
    if (!origin.valid || !origin.parent) {
      continue;
    }
    const int barcode = origin.parent->barcode();
    auto [position, inserted] = generator_candidate_by_barcode.emplace(barcode, candidates.size());
    if (inserted) {
      candidates.push_back({Pi0Pathway::generator_decay, barcode, origin.parent, nullptr, {}});
    }
    candidates[position->second].photons.push_back(primary);
  }

  for (auto iterator = secondary_range.first; iterator != secondary_range.second; ++iterator) {
    const PHG4Particle* pi0 = iterator->second;
    if (!pi0 || pi0->get_pid() != 111) {
      continue;
    }
    const auto found = children_by_parent.find(pi0->get_track_id());
    if (found == children_by_parent.end() || found->second.size() != 2U || found->second[0]->get_pid() != 22 || found->second[1]->get_pid() != 22) {
      ++n_pi0_malformed_daughters_;
      continue;
    }
    const PHG4Particle* ancestor = pi0;
    std::set<int> visited;
    while (ancestor && !truth->is_primary(ancestor) && visited.insert(ancestor->get_track_id()).second) {
      ancestor = truth->GetParticle(ancestor->get_parent_id());
    }
    if (!ancestor || !truth->is_primary(ancestor) || truth->isEmbeded(ancestor->get_track_id()) != signal_embedding_id_) {
      continue;
    }
    candidates.push_back({Pi0Pathway::g4_secondary_decay, ancestor->get_track_id(), nullptr, pi0, found->second});
  }

  const double target_radius = cemc_radius(geometry);
  for (const Pi0Candidate& candidate : candidates) {
    double parent_eta = 0.0;
    bool valid_parent_eta = finite_eta(candidate.parent, parent_eta);
    if (candidate.g4_parent) {
      const double parent_pt = particle_pt(candidate.g4_parent);
      valid_parent_eta = parent_pt > 0.0;
      parent_eta = valid_parent_eta
          ? std::asinh(candidate.g4_parent->get_pz() / parent_pt) : 0.0;
    }
    if (!valid_parent_eta || !std::isfinite(parent_eta) || std::abs(parent_eta) >= truth_eta_max_) {
      continue;
    }
    if (candidate.photons.size() != 2U) {
      ++n_pi0_malformed_daughters_;
      continue;
    }
    if (candidate.pathway != Pi0Pathway::generator_decay) {
      ++n_pi0_candidate_g4_decay_;
    } else {
      ++n_pi0_candidate_generator_decay_;
    }
    const std::array<Projection, 2> projections = {
        project_photon(candidate.photons[0], truth, target_radius),
        project_photon(candidate.photons[1], truth, target_radius)};
    if (!projections[0].valid || !projections[1].valid) {
      ++n_pi0_projection_failure_;
      continue;
    }

    std::vector<std::size_t> eligible;
    for (std::size_t index = 0; index < cluster_records.size(); ++index) {
      const ClusterRecord& cluster = cluster_records[index];
      if (!cluster.truth.valid) {
        continue;
      }
      const bool matching_contributor = std::any_of(
          cluster.truth.contributors.begin(), cluster.truth.contributors.end(),
          [&](const photon_tree::TruthContributor& contributor) {
            return contributor.fraction >= pi0_contributor_fraction_min_ &&
                contributor_matches_pi0(contributor, candidate.pathway,
                    candidate.parent_barcode, event_map, signal_embedding_id_);
          });
      if (matching_contributor) {
        eligible.push_back(index);
      }
    }

    if (candidate.pathway == Pi0Pathway::g4_secondary_decay) {
      for (const std::size_t index : eligible) {
        const ClusterRecord& cluster = cluster_records[index];
        if (pi0_clusters_filled.insert(cluster.cluster).second) {
          h_pi0_->Fill(cluster.et);
          ++n_pi0_cluster_;
          ++n_pi0_cluster_g4_decay_;
        }
      }
    }

    constexpr std::size_t invalid_index = std::numeric_limits<std::size_t>::max();
    std::array<std::size_t, 2> individual = {invalid_index, invalid_index};
    std::array<double, 2> individual_distance = {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
    for (std::size_t photon = 0; photon < 2U; ++photon) {
      for (const std::size_t index : eligible) {
        const ClusterRecord& cluster = cluster_records[index];
        const double distance = delta_r(projections[photon].eta, projections[photon].phi, cluster.surface_eta, cluster.surface_phi);
        if (distance < separated_delta_r_cut_ && distance < individual_distance[photon]) {
          individual[photon] = index;
          individual_distance[photon] = distance;
        }
      }
    }

    if (individual[0] != invalid_index && individual[1] != invalid_index && individual[0] != individual[1]) {
      h_pi0_separated_->Fill(cluster_records[individual[0]].et);
      h_pi0_separated_->Fill(cluster_records[individual[1]].et);
      ++n_pi0_separated_;
      n_pi0_separated_cluster_fill_ += 2ULL;
      continue;
    }

    std::size_t merged = invalid_index;
    double merged_distance = std::numeric_limits<double>::infinity();
    const double parent_pt = candidate.g4_parent
        ? particle_pt(candidate.g4_parent) : hepmc_pt(candidate.parent);
    for (const std::size_t index : eligible) {
      const ClusterRecord& cluster = cluster_records[index];
      const double distance0 = delta_r(projections[0].eta, projections[0].phi, cluster.surface_eta, cluster.surface_phi);
      const double distance1 = delta_r(projections[1].eta, projections[1].phi, cluster.surface_eta, cluster.surface_phi);
      const double maximum_distance = std::max(distance0, distance1);
      const double response = parent_pt > 0.0 ? cluster.et / parent_pt : -1.0;
      if (maximum_distance < merged_delta_r_cut_ && response >= response_min_ && response <= response_max_ && maximum_distance < merged_distance) {
        merged = index;
        merged_distance = maximum_distance;
      }
    }
    if (merged != invalid_index) {
      h_pi0_merged_->Fill(cluster_records[merged].et);
      ++n_pi0_merged_;
      ++n_pi0_merged_cluster_fill_;
      continue;
    }

    const bool matched0 = individual[0] != invalid_index;
    const bool matched1 = individual[1] != invalid_index;
    if (matched0 != matched1) {
      const std::size_t photon = matched0 ? 0U : 1U;
      const std::size_t index = individual[photon];
      const double daughter_pt = particle_pt(candidate.photons[photon]);
      const double response = daughter_pt > 0.0
          ? cluster_records[index].et / daughter_pt : -1.0;
      if (response >= response_min_ && response <= response_max_) {
        h_pi0_missing_->Fill(cluster_records[index].et);
        ++n_pi0_missing_;
        ++n_pi0_missing_cluster_fill_;
      } else {
        ++n_pi0_none_;
      }
      continue;
    }
    if (matched0 && matched1 && individual[0] == individual[1]) {
      ++n_pi0_ambiguous_;
    } else {
      ++n_pi0_none_;
    }
  }

  ++n_events_written_;
  if (verbosity_ > 0 && n_events_processed_ % 1000ULL == 0ULL) {
    std::cout << "PythiaClusterEtSpectrum - processed " << n_events_processed_ << " events" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int PythiaClusterEtSpectrum::End(PHCompositeNode* /*topNode*/)
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
  std::cout << "PythiaClusterEtSpectrum - events/prompt/pi0/separated/merged/missing = "
            << n_events_written_ << "/" << n_prompt_cluster_ << "/"
            << n_pi0_cluster_ << "/" << n_pi0_separated_cluster_fill_ << "/"
            << n_pi0_merged_cluster_fill_ << "/" << n_pi0_missing_cluster_fill_
            << std::endl;
  return write_error ? Fun4AllReturnCodes::ABORTRUN : Fun4AllReturnCodes::EVENT_OK;
}

void PythiaClusterEtSpectrum::create_output_directory() const
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

void PythiaClusterEtSpectrum::create_output()
{
  output_file_ = TFile::Open(output_file_name_.c_str(), "RECREATE");
  if (!output_file_ || output_file_->IsZombie())
  {
    return;
  }
  output_file_->cd();
  h_prompt_ = new TH1D("h_prompt_cluster_et_raw", "", n_bins_, 0.0, et_max_);
  h_pi0_ = new TH1D("h_pi0_cluster_et_raw", "", n_bins_, 0.0, et_max_);
  h_pi0_separated_ = new TH1D("h_pi0_separated_cluster_et_raw", "", n_bins_, 0.0, et_max_);
  h_pi0_merged_ = new TH1D("h_pi0_merged_cluster_et_raw", "", n_bins_, 0.0, et_max_);
  h_pi0_missing_ = new TH1D("h_pi0_missing_cluster_et_raw", "", n_bins_, 0.0, et_max_);
  for (TH1D* histogram : {h_prompt_, h_pi0_, h_pi0_separated_,
                          h_pi0_merged_, h_pi0_missing_})
  {
    histogram->Sumw2();
  }

  metadata_tree_ = new TTree("metadata", "Pythia cluster ET partial metadata");
  static int schema_version = schema_version_;
  static unsigned char bin_width_normalized = 0U;
  metadata_tree_->Branch("schema_version", &schema_version);
  metadata_tree_->Branch("manifest_path", &manifest_path_);
  metadata_tree_->Branch("manifest_begin", &manifest_begin_);
  metadata_tree_->Branch("manifest_end", &manifest_end_);
  metadata_tree_->Branch("first_suffix", &first_suffix_);
  metadata_tree_->Branch("last_suffix", &last_suffix_);
  metadata_tree_->Branch("cluster_collection", &cluster_collection_);
  metadata_tree_->Branch("prompt_selection", &prompt_selection_);
  metadata_tree_->Branch("pi0_selection", &pi0_selection_);
  metadata_tree_->Branch("topology_priority", &topology_priority_);
  metadata_tree_->Branch("projection_scheme", &projection_scheme_);
  metadata_tree_->Branch("signal_embedding_id", &signal_embedding_id_);
  metadata_tree_->Branch("n_bins", &n_bins_);
  metadata_tree_->Branch("et_max", &et_max_);
  metadata_tree_->Branch("truth_eta_max", &truth_eta_max_);
  metadata_tree_->Branch("cluster_eta_max", &cluster_eta_max_);
  metadata_tree_->Branch("min_cluster_energy", &min_cluster_energy_);
  metadata_tree_->Branch("dominant_fraction_min", &dominant_fraction_min_);
  metadata_tree_->Branch("pi0_contributor_fraction_min", &pi0_contributor_fraction_min_);
  metadata_tree_->Branch("separated_delta_r_cut", &separated_delta_r_cut_);
  metadata_tree_->Branch("merged_delta_r_cut", &merged_delta_r_cut_);
  metadata_tree_->Branch("response_min", &response_min_);
  metadata_tree_->Branch("response_max", &response_max_);
  metadata_tree_->Branch("bin_width_normalized", &bin_width_normalized);
  metadata_tree_->Branch("events_processed", &n_events_processed_);
  metadata_tree_->Branch("events_written", &n_events_written_);
  metadata_tree_->Branch("events_invalid", &n_events_invalid_);
  metadata_tree_->Branch("cluster_considered_count", &n_cluster_considered_);
  metadata_tree_->Branch("cluster_invalid_truth_count", &n_cluster_invalid_truth_);
  metadata_tree_->Branch("prompt_cluster_count", &n_prompt_cluster_);
  metadata_tree_->Branch("pi0_cluster_count", &n_pi0_cluster_);
  metadata_tree_->Branch("pi0_cluster_g4_decay_count", &n_pi0_cluster_g4_decay_);
  metadata_tree_->Branch("pi0_cluster_generator_decay_count", &n_pi0_cluster_generator_decay_);
  metadata_tree_->Branch("pi0_candidate_g4_decay_count", &n_pi0_candidate_g4_decay_);
  metadata_tree_->Branch("pi0_candidate_generator_decay_count", &n_pi0_candidate_generator_decay_);
  metadata_tree_->Branch("pi0_malformed_daughters_count", &n_pi0_malformed_daughters_);
  metadata_tree_->Branch("pi0_projection_failure_count", &n_pi0_projection_failure_);
  metadata_tree_->Branch("pi0_separated_count", &n_pi0_separated_);
  metadata_tree_->Branch("pi0_merged_count", &n_pi0_merged_);
  metadata_tree_->Branch("pi0_missing_count", &n_pi0_missing_);
  metadata_tree_->Branch("pi0_none_count", &n_pi0_none_);
  metadata_tree_->Branch("pi0_ambiguous_count", &n_pi0_ambiguous_);
  metadata_tree_->Branch("pi0_separated_cluster_fill_count", &n_pi0_separated_cluster_fill_);
  metadata_tree_->Branch("pi0_merged_cluster_fill_count", &n_pi0_merged_cluster_fill_);
  metadata_tree_->Branch("pi0_missing_cluster_fill_count", &n_pi0_missing_cluster_fill_);
}

void PythiaClusterEtSpectrum::close_output()
{
  if (output_file_)
  {
    if (output_file_->IsOpen())
    {
      output_file_->Close();
    }
    delete output_file_;
    output_file_ = nullptr;
  }
  h_prompt_ = nullptr;
  h_pi0_ = nullptr;
  h_pi0_separated_ = nullptr;
  h_pi0_merged_ = nullptr;
  h_pi0_missing_ = nullptr;
  metadata_tree_ = nullptr;
}
