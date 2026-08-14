#include "PythiaPi0AnchorClusterSpectrum.h"

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/TowerInfoContainer.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <g4detectors/PHG4CellContainer.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
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
#include <string>
#include <utility>
#include <vector>

namespace
{
constexpr std::size_t invalid_index = std::numeric_limits<std::size_t>::max();

enum class Pi0Pathway : int
{
  g4_primary_decay = 1,
  generator_decay = 2
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
  std::array<const PHG4Particle*, 2> photons = {nullptr, nullptr};
};

struct PendingGeneratorCandidate
{
  const HepMC::GenParticle* parent = nullptr;
  std::vector<const PHG4Particle*> photons;
};

struct ClusterRecord
{
  const RawCluster* cluster = nullptr;
  double et = 0.0;
  double eta = 0.0;
  bool anchor_acceptance = false;
  photon_tree::ClusterTruthMatch truth;
};

struct AnchorRecord
{
  std::size_t cluster_index = invalid_index;
  std::size_t candidate_index = invalid_index;
  bool ambiguous_main = false;
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

bool finite_eta(const HepMC::GenParticle* particle, double& eta)
{
  if (!particle)
  {
    return false;
  }
  const auto& momentum = particle->momentum();
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
    const PHHepMCGenEventMap* event_map,
    int signal_embedding_id,
    const std::map<int, std::size_t>& g4_candidate_by_barcode,
    const std::map<int, std::size_t>& generator_candidate_by_barcode)
{
  if (contributor.embedding_id != signal_embedding_id ||
      !contributor.hepmc_valid)
  {
    return invalid_index;
  }
  const HepMC::GenParticle* particle =
      contributor_hepmc_particle(contributor, event_map);
  if (contributor.g4_pdg_id == 111 && particle &&
      particle->pdg_id() == 111)
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
}

PythiaPi0AnchorClusterSpectrum::PythiaPi0AnchorClusterSpectrum(
    const std::string& name)
  : SubsysReco(name)
{
}

PythiaPi0AnchorClusterSpectrum::~PythiaPi0AnchorClusterSpectrum()
{
  close_output();
}

int PythiaPi0AnchorClusterSpectrum::Init(PHCompositeNode* /*topNode*/)
{
  const bool valid = !output_file_name_.empty() && !manifest_path_.empty() &&
      manifest_begin_ >= 0 && manifest_end_ > manifest_begin_ &&
      !first_suffix_.empty() && !last_suffix_.empty() &&
      signal_embedding_id_ > 0 && n_bins_ > 0 &&
      std::isfinite(et_max_) && et_max_ > 0.0 &&
      std::isfinite(truth_eta_max_) && truth_eta_max_ > 0.0 &&
      std::isfinite(anchor_cluster_eta_max_) &&
      anchor_cluster_eta_max_ > 0.0 &&
      std::isfinite(partner_cluster_eta_max_) &&
      std::isfinite(min_cluster_energy_) && min_cluster_energy_ >= 0.0 &&
      dominant_fraction_min_ >= 0.0 && dominant_fraction_min_ <= 1.0 &&
      anchor_pi0_fraction_min_ >= 0.0 &&
      anchor_pi0_fraction_min_ <= 1.0 &&
      min_energy_contribution_fraction_ >= 0.0 &&
      min_energy_contribution_fraction_ < 1.0;
  if (!valid)
  {
    std::cerr << "PythiaPi0AnchorClusterSpectrum::Init - invalid configuration"
              << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  truth_matcher_.set_verbosity(verbosity_);
  create_output_directory();
  create_output();
  if (!output_file_ || output_file_->IsZombie())
  {
    std::cerr << "PythiaPi0AnchorClusterSpectrum::Init - failed to create "
              << output_file_name_ << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int PythiaPi0AnchorClusterSpectrum::process_event(PHCompositeNode* topNode)
{
  ++n_events_processed_;
  auto* truth =
      findNode::getClass<PHG4TruthInfoContainer>(topNode, truth_node_name_);
  auto* event_map = findNode::getClass<PHHepMCGenEventMap>(
      topNode, hepmc_event_map_node_name_);
  auto* towers =
      findNode::getClass<TowerInfoContainer>(topNode, tower_node_name_);
  auto* raw_truth_towers = findNode::getClass<RawTowerContainer>(
      topNode, raw_truth_tower_node_name_);
  auto* truth_cells = findNode::getClass<PHG4CellContainer>(
      topNode, truth_cell_node_name_);
  auto* truth_hits = findNode::getClass<PHG4HitContainer>(
      topNode, truth_hit_node_name_);
  auto* clusters = findNode::getClass<RawClusterContainer>(
      topNode, split_cluster_node_name_);
  const PHHepMCGenEvent* signal_event = event_map
      ? event_map->get(signal_embedding_id_) : nullptr;
  const HepMC::GenEvent* event =
      signal_event ? signal_event->getEvent() : nullptr;
  if (!truth || !event_map || !towers || !raw_truth_towers || !truth_cells ||
      !truth_hits || !clusters || !signal_event ||
      !signal_event->is_simulated() || !event ||
      !truth_matcher_.begin_event(topNode))
  {
    ++n_events_invalid_;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  const HepMC::FourVector& collision = signal_event->get_collision_vertex();
  if (!std::isfinite(collision.x()) || !std::isfinite(collision.y()) ||
      !std::isfinite(collision.z()))
  {
    ++n_events_invalid_;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  std::vector<ClusterRecord> cluster_records;
  const auto cluster_range = clusters->getClusters();
  for (auto iterator = cluster_range.first; iterator != cluster_range.second;
       ++iterator)
  {
    const RawCluster* cluster = iterator->second;
    if (!cluster || !std::isfinite(cluster->get_energy()) ||
        !std::isfinite(cluster->get_x()) || !std::isfinite(cluster->get_y()) ||
        !std::isfinite(cluster->get_z()) ||
        cluster->get_energy() < min_cluster_energy_)
    {
      continue;
    }
    const double dx = cluster->get_x() - collision.x();
    const double dy = cluster->get_y() - collision.y();
    const double dz = cluster->get_z() - collision.z();
    const double distance = std::sqrt(dx * dx + dy * dy + dz * dz);
    const double transverse = std::hypot(dx, dy);
    if (!(distance > 0.0) || !(transverse > 0.0))
    {
      continue;
    }
    ClusterRecord record;
    record.cluster = cluster;
    record.et = cluster->get_energy() * transverse / distance;
    record.eta = std::asinh(dz / transverse);
    if (!std::isfinite(record.et) || !(record.et > 0.0) ||
        !std::isfinite(record.eta))
    {
      continue;
    }
    if (partner_cluster_eta_max_ > 0.0 &&
        std::abs(record.eta) >= partner_cluster_eta_max_)
    {
      continue;
    }
    record.anchor_acceptance =
        std::abs(record.eta) < anchor_cluster_eta_max_;
    record.truth = truth_matcher_.match(
        cluster, towers, raw_truth_towers, truth, event_map, true);
    ++n_cluster_considered_;
    if (!record.truth.valid)
    {
      ++n_cluster_invalid_truth_;
    }
    cluster_records.push_back(std::move(record));
  }

  // Prompt is retained only as an external spectrum reference. It is not part
  // of the pi0-anchor partition.
  for (const ClusterRecord& cluster : cluster_records)
  {
    if (!cluster.anchor_acceptance || !cluster.truth.valid ||
        cluster.truth.contributors.empty())
    {
      continue;
    }
    const auto& dominant = cluster.truth.contributors.front();
    if (dominant.embedding_id != signal_embedding_id_ ||
        !dominant.hepmc_valid ||
        dominant.fraction < dominant_fraction_min_ ||
        dominant.g4_pdg_id != 22 ||
        (dominant.photon_category != 1 && dominant.photon_category != 2))
    {
      continue;
    }
    double truth_eta = 0.0;
    const HepMC::GenParticle* particle =
        contributor_hepmc_particle(dominant, event_map);
    if (finite_eta(particle, truth_eta) &&
        std::abs(truth_eta) < truth_eta_max_)
    {
      h_prompt_->Fill(cluster.et);
      ++n_prompt_cluster_;
    }
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

  std::vector<Pi0Candidate> candidates;
  std::map<int, PendingGeneratorCandidate> pending_generator;
  const auto primary_range = truth->GetPrimaryParticleRange();
  for (auto iterator = primary_range.first; iterator != primary_range.second;
       ++iterator)
  {
    const PHG4Particle* primary = iterator->second;
    if (!primary ||
        truth->isEmbeded(primary->get_track_id()) != signal_embedding_id_)
    {
      continue;
    }
    const HepMC::GenParticle* hepmc_particle =
        event->barcode_to_particle(primary->get_barcode());
    if (primary->get_pid() == 111 && hepmc_particle &&
        hepmc_particle->pdg_id() == 111)
    {
      double eta = 0.0;
      const auto found = children_by_parent.find(primary->get_track_id());
      if (found == children_by_parent.end() || found->second.size() != 2U ||
          found->second[0]->get_pid() != 22 ||
          found->second[1]->get_pid() != 22)
      {
        ++n_pi0_malformed_daughters_;
        continue;
      }
      const double pt = particle_pt(primary);
      eta = pt > 0.0 ? std::asinh(primary->get_pz() / pt) : 0.0;
      if (!(pt > 0.0) || !std::isfinite(eta) ||
          std::abs(eta) >= truth_eta_max_)
      {
        continue;
      }
      candidates.push_back({Pi0Pathway::g4_primary_decay,
          hepmc_particle->barcode(), hepmc_particle, primary,
          {found->second[0], found->second[1]}});
      ++n_pi0_candidate_g4_decay_;
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
    double eta = 0.0;
    if (!pending.parent || pending.photons.size() != 2U)
    {
      ++n_pi0_malformed_daughters_;
      continue;
    }
    if (!finite_eta(pending.parent, eta) || std::abs(eta) >= truth_eta_max_)
    {
      continue;
    }
    candidates.push_back({Pi0Pathway::generator_decay, barcode,
        pending.parent, nullptr,
        {pending.photons[0], pending.photons[1]}});
    ++n_pi0_candidate_generator_decay_;
  }

  std::map<int, std::size_t> g4_candidate_by_barcode;
  std::map<int, std::size_t> generator_candidate_by_barcode;
  for (std::size_t index = 0; index < candidates.size(); ++index)
  {
    const Pi0Candidate& candidate = candidates[index];
    if (candidate.pathway == Pi0Pathway::g4_primary_decay)
    {
      g4_candidate_by_barcode[candidate.parent_barcode] = index;
    }
    else
    {
      generator_candidate_by_barcode[candidate.parent_barcode] = index;
    }
  }

  std::vector<AnchorRecord> anchors;
  std::vector<std::vector<std::size_t>> anchors_by_candidate(
      candidates.size());
  constexpr double tie_tolerance = 1e-6;
  for (std::size_t cluster_index = 0;
       cluster_index < cluster_records.size(); ++cluster_index)
  {
    const ClusterRecord& cluster = cluster_records[cluster_index];
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
          contributor, event_map, signal_embedding_id_,
          g4_candidate_by_barcode, generator_candidate_by_barcode);
      if (candidate_index == invalid_index)
      {
        unmatched_max_fraction =
            std::max(unmatched_max_fraction,
                     static_cast<double>(contributor.fraction));
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
        best_fraction < anchor_pi0_fraction_min_)
    {
      continue;
    }

    const bool ambiguous_main =
        second_fraction >= best_fraction - tie_tolerance ||
        unmatched_max_fraction >= best_fraction - tie_tolerance;
    anchors.push_back(
        {cluster_index, best_candidate, ambiguous_main});
    anchors_by_candidate[best_candidate].push_back(anchors.size() - 1U);
    h_anchor_->Fill(cluster.et);
    ++n_anchor_cluster_;
    if (candidates[best_candidate].pathway ==
        Pi0Pathway::g4_primary_decay)
    {
      ++n_anchor_g4_decay_;
    }
    else
    {
      ++n_anchor_generator_decay_;
    }
    if (ambiguous_main)
    {
      ++n_anchor_ambiguous_main_;
    }
  }

  for (std::size_t candidate_index = 0;
       candidate_index < candidates.size(); ++candidate_index)
  {
    if (anchors_by_candidate[candidate_index].empty())
    {
      continue;
    }
    const Pi0Candidate& candidate = candidates[candidate_index];
    const std::array<int, 2> direct_gamma_track_ids = {
        candidate.photons[0]->get_track_id(),
        candidate.photons[1]->get_track_id()};
    std::array<std::size_t, 2> best_cluster = {
        invalid_index, invalid_index};
    std::array<double, 2> maximum_edep = {-1.0, -1.0};

    for (std::size_t cluster_index = 0;
         cluster_index < cluster_records.size(); ++cluster_index)
    {
      const ClusterRecord& cluster = cluster_records[cluster_index];
      const photon_tree::Pi0ClusterTruthMatch match =
          pi0_truth_matcher_.match(
              cluster.cluster, towers, raw_truth_towers, truth_cells,
              truth_hits, truth, direct_gamma_track_ids, true);
      if (!match.valid)
      {
        ++n_energy_match_invalid_;
        continue;
      }
      for (std::size_t photon = 0; photon < 2U; ++photon)
      {
        const double deposit = match.gamma_edep[photon];
        const double fraction = match.total_edep > 0.0F
            ? deposit / match.total_edep : 0.0;
        if (deposit > 0.0 &&
            fraction > min_energy_contribution_fraction_ &&
            deposit > maximum_edep[photon])
        {
          best_cluster[photon] = cluster_index;
          maximum_edep[photon] = deposit;
        }
      }
    }

    for (const std::size_t anchor_position :
         anchors_by_candidate[candidate_index])
    {
      const AnchorRecord& anchor = anchors[anchor_position];
      const ClusterRecord& cluster =
          cluster_records[anchor.cluster_index];
      if (anchor.ambiguous_main)
      {
        h_other_->Fill(cluster.et);
        ++n_other_;
        continue;
      }

      const bool is_best0 =
          best_cluster[0] == anchor.cluster_index;
      const bool is_best1 =
          best_cluster[1] == anchor.cluster_index;
      const bool found0 = best_cluster[0] != invalid_index;
      const bool found1 = best_cluster[1] != invalid_index;

      if (is_best0 && is_best1)
      {
        h_merged_->Fill(cluster.et);
        ++n_merged_;
      }
      else if ((is_best0 && found1 &&
                best_cluster[1] != anchor.cluster_index) ||
               (is_best1 && found0 &&
                best_cluster[0] != anchor.cluster_index))
      {
        h_separated_->Fill(cluster.et);
        ++n_separated_;
      }
      else if ((is_best0 && !found1) || (is_best1 && !found0))
      {
        h_missing_->Fill(cluster.et);
        ++n_missing_;
      }
      else
      {
        h_other_->Fill(cluster.et);
        ++n_other_;
      }
    }
  }

  ++n_events_written_;
  if (verbosity_ > 0 && n_events_processed_ % 1000ULL == 0ULL)
  {
    std::cout << "PythiaPi0AnchorClusterSpectrum - processed "
              << n_events_processed_ << " events" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int PythiaPi0AnchorClusterSpectrum::End(PHCompositeNode* /*topNode*/)
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
  std::cout
      << "PythiaPi0AnchorClusterSpectrum - events/anchors/separated/merged/missing/other = "
      << n_events_written_ << "/" << n_anchor_cluster_ << "/"
      << n_separated_ << "/" << n_merged_ << "/" << n_missing_ << "/"
      << n_other_ << std::endl;
  return write_error ? Fun4AllReturnCodes::ABORTRUN
                     : Fun4AllReturnCodes::EVENT_OK;
}

void PythiaPi0AnchorClusterSpectrum::create_output_directory() const
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

void PythiaPi0AnchorClusterSpectrum::create_output()
{
  output_file_ = TFile::Open(output_file_name_.c_str(), "RECREATE");
  if (!output_file_ || output_file_->IsZombie())
  {
    return;
  }
  output_file_->cd();
  h_prompt_ = new TH1D(
      "h_prompt_cluster_et_raw", "", n_bins_, 0.0, et_max_);
  h_anchor_ = new TH1D(
      "h_pi0_anchor_cluster_et_raw", "", n_bins_, 0.0, et_max_);
  h_separated_ = new TH1D(
      "h_pi0_anchor_separated_cluster_et_raw", "",
      n_bins_, 0.0, et_max_);
  h_merged_ = new TH1D(
      "h_pi0_anchor_merged_cluster_et_raw", "",
      n_bins_, 0.0, et_max_);
  h_missing_ = new TH1D(
      "h_pi0_anchor_missing_cluster_et_raw", "",
      n_bins_, 0.0, et_max_);
  h_other_ = new TH1D(
      "h_pi0_anchor_other_cluster_et_raw", "",
      n_bins_, 0.0, et_max_);
  for (TH1D* histogram : {
           h_prompt_, h_anchor_, h_separated_, h_merged_, h_missing_,
           h_other_})
  {
    histogram->Sumw2();
  }

  metadata_tree_ = new TTree(
      "metadata", "Pythia pi0 anchor-cluster partial metadata");
  static int schema_version = schema_version_;
  static unsigned char bin_width_normalized = 0U;
  metadata_tree_->Branch("schema_version", &schema_version);
  metadata_tree_->Branch("manifest_path", &manifest_path_);
  metadata_tree_->Branch("manifest_begin", &manifest_begin_);
  metadata_tree_->Branch("manifest_end", &manifest_end_);
  metadata_tree_->Branch("first_suffix", &first_suffix_);
  metadata_tree_->Branch("last_suffix", &last_suffix_);
  metadata_tree_->Branch("cluster_collection", &cluster_collection_);
  metadata_tree_->Branch("classification_unit", &classification_unit_);
  metadata_tree_->Branch("pi0_selection", &pi0_selection_);
  metadata_tree_->Branch("partner_selection", &partner_selection_);
  metadata_tree_->Branch("topology_definition", &topology_definition_);
  metadata_tree_->Branch("topology_priority", &topology_priority_);
  metadata_tree_->Branch("response_policy", &response_policy_);
  metadata_tree_->Branch("signal_embedding_id", &signal_embedding_id_);
  metadata_tree_->Branch("n_bins", &n_bins_);
  metadata_tree_->Branch("et_max", &et_max_);
  metadata_tree_->Branch("truth_eta_max", &truth_eta_max_);
  metadata_tree_->Branch(
      "anchor_cluster_eta_max", &anchor_cluster_eta_max_);
  metadata_tree_->Branch(
      "partner_cluster_eta_max", &partner_cluster_eta_max_);
  metadata_tree_->Branch("min_cluster_energy", &min_cluster_energy_);
  metadata_tree_->Branch(
      "dominant_fraction_min", &dominant_fraction_min_);
  metadata_tree_->Branch(
      "anchor_pi0_fraction_min", &anchor_pi0_fraction_min_);
  metadata_tree_->Branch(
      "min_energy_contribution_fraction",
      &min_energy_contribution_fraction_);
  metadata_tree_->Branch(
      "pi0_truth_matching_algorithm_version",
      &pi0_truth_matching_algorithm_version_);
  metadata_tree_->Branch(
      "bin_width_normalized", &bin_width_normalized);
  metadata_tree_->Branch("events_processed", &n_events_processed_);
  metadata_tree_->Branch("events_written", &n_events_written_);
  metadata_tree_->Branch("events_invalid", &n_events_invalid_);
  metadata_tree_->Branch(
      "cluster_considered_count", &n_cluster_considered_);
  metadata_tree_->Branch(
      "cluster_invalid_truth_count", &n_cluster_invalid_truth_);
  metadata_tree_->Branch(
      "prompt_cluster_count", &n_prompt_cluster_);
  metadata_tree_->Branch(
      "pi0_candidate_g4_decay_count", &n_pi0_candidate_g4_decay_);
  metadata_tree_->Branch(
      "pi0_candidate_generator_decay_count",
      &n_pi0_candidate_generator_decay_);
  metadata_tree_->Branch(
      "pi0_malformed_daughters_count", &n_pi0_malformed_daughters_);
  metadata_tree_->Branch(
      "anchor_cluster_count", &n_anchor_cluster_);
  metadata_tree_->Branch(
      "anchor_g4_decay_count", &n_anchor_g4_decay_);
  metadata_tree_->Branch(
      "anchor_generator_decay_count", &n_anchor_generator_decay_);
  metadata_tree_->Branch(
      "anchor_ambiguous_main_count", &n_anchor_ambiguous_main_);
  metadata_tree_->Branch(
      "energy_match_invalid_count", &n_energy_match_invalid_);
  metadata_tree_->Branch("separated_count", &n_separated_);
  metadata_tree_->Branch("merged_count", &n_merged_);
  metadata_tree_->Branch("missing_count", &n_missing_);
  metadata_tree_->Branch("other_count", &n_other_);
}

void PythiaPi0AnchorClusterSpectrum::close_output()
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
  h_anchor_ = nullptr;
  h_separated_ = nullptr;
  h_merged_ = nullptr;
  h_missing_ = nullptr;
  h_other_ = nullptr;
  metadata_tree_ = nullptr;
}
