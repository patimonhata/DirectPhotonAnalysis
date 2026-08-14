#include "PythiaPi0ClusterMultiplicity.h"

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/TowerInfoContainer.h>
#include <fun4all/Fun4AllReturnCodes.h>
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
#include <TH2D.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace
{
enum class Pi0Pathway : int
{
  g4_primary_decay = 0,
  generator_decay = 1
};

struct Pi0Origin
{
  bool valid = false;
  const HepMC::GenParticle* parent = nullptr;
};

struct Pi0Candidate
{
  Pi0Pathway pathway = Pi0Pathway::g4_primary_decay;
  int parent_key = 0;
  const HepMC::GenParticle* parent = nullptr;
  const PHG4Particle* g4_parent = nullptr;
  std::vector<const PHG4Particle*> photons;
};

struct ClusterRecord
{
  double energy = 0.0;
  photon_tree::ClusterTruthMatch truth;
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

double hepmc_pt(const HepMC::GenParticle* particle)
{
  return particle
      ? std::hypot(particle->momentum().px(), particle->momentum().py())
      : 0.0;
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
                             const Pi0Candidate& candidate,
                             const PHHepMCGenEventMap* event_map,
                             int signal_embedding_id)
{
  if (contributor.embedding_id != signal_embedding_id ||
      !contributor.hepmc_valid)
  {
    return false;
  }
  const HepMC::GenParticle* particle =
      contributor_hepmc_particle(contributor, event_map);
  if (candidate.pathway == Pi0Pathway::g4_primary_decay)
  {
    return contributor.g4_pdg_id == 111 && particle &&
        particle->pdg_id() == 111 && particle->barcode() == candidate.parent_key;
  }
  const Pi0Origin origin = trace_pi0_origin(particle);
  return contributor.g4_pdg_id == 22 && origin.valid && origin.parent &&
      origin.parent->barcode() == candidate.parent_key;
}

const char* threshold_tag(std::size_t index)
{
  static constexpr std::array<const char*, 4> tags = {
      "0p0", "0p1", "0p3", "0p5"};
  return tags.at(index);
}

const char* pathway_tag(std::size_t index)
{
  static constexpr std::array<const char*, 2> tags = {
      "g4_primary", "generator"};
  return tags.at(index);
}
}

PythiaPi0ClusterMultiplicity::PythiaPi0ClusterMultiplicity(
    const std::string& name)
  : SubsysReco(name)
{
}

PythiaPi0ClusterMultiplicity::~PythiaPi0ClusterMultiplicity()
{
  close_output();
}

int PythiaPi0ClusterMultiplicity::Init(PHCompositeNode* /*topNode*/)
{
  const bool valid = !output_file_name_.empty() && !manifest_path_.empty() &&
      manifest_begin_ >= 0 && manifest_end_ > manifest_begin_ &&
      !first_suffix_.empty() && !last_suffix_.empty() &&
      signal_embedding_id_ > 0 && pt_bins_ > 0 && pt_max_ > 0.0 &&
      multiplicity_max_ > 0 && cluster_energy_bins_ > 0 &&
      cluster_energy_max_ > 0.0 && truth_eta_max_ > 0.0 &&
      cluster_eta_max_ > 0.0;
  if (!valid)
  {
    std::cerr << "PythiaPi0ClusterMultiplicity::Init - invalid configuration"
              << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  truth_matcher_.set_verbosity(verbosity_);
  create_output_directory();
  create_output();
  if (!output_file_ || output_file_->IsZombie())
  {
    std::cerr << "PythiaPi0ClusterMultiplicity::Init - failed to create "
              << output_file_name_ << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int PythiaPi0ClusterMultiplicity::process_event(PHCompositeNode* topNode)
{
  ++n_events_processed_;
  auto* truth = findNode::getClass<PHG4TruthInfoContainer>(topNode, truth_node_name_);
  auto* event_map = findNode::getClass<PHHepMCGenEventMap>(
      topNode, hepmc_event_map_node_name_);
  auto* towers = findNode::getClass<TowerInfoContainer>(topNode, tower_node_name_);
  auto* raw_truth_towers = findNode::getClass<RawTowerContainer>(
      topNode, raw_truth_tower_node_name_);
  auto* clusters = findNode::getClass<RawClusterContainer>(
      topNode, split_cluster_node_name_);
  const PHHepMCGenEvent* signal_event = event_map
      ? event_map->get(signal_embedding_id_) : nullptr;
  const HepMC::GenEvent* event = signal_event ? signal_event->getEvent() : nullptr;
  if (!truth || !event_map || !towers || !clusters || !signal_event ||
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
        !std::isfinite(cluster->get_z()))
    {
      continue;
    }
    const double dx = cluster->get_x() - collision.x();
    const double dy = cluster->get_y() - collision.y();
    const double dz = cluster->get_z() - collision.z();
    const double transverse = std::hypot(dx, dy);
    if (!(transverse > 0.0))
    {
      continue;
    }
    const double eta = std::asinh(dz / transverse);
    if (!std::isfinite(eta) || std::abs(eta) >= cluster_eta_max_)
    {
      continue;
    }
    ClusterRecord record;
    record.energy = cluster->get_energy();
    record.truth = truth_matcher_.match(
        cluster, towers, raw_truth_towers, truth, event_map, true);
    ++n_cluster_considered_;
    if (!record.truth.valid)
    {
      ++n_cluster_invalid_truth_;
    }
    cluster_records.push_back(std::move(record));
  }

  // Secondary particles are inspected only as daughters of a G4-primary pi0;
  // G4-secondary pi0s themselves are not candidates.
  std::map<int, std::vector<const PHG4Particle*>> children_by_parent;
  const auto secondary_range = truth->GetSecondaryParticleRange();
  for (auto iterator = secondary_range.first; iterator != secondary_range.second;
       ++iterator)
  {
    const PHG4Particle* particle = iterator->second;
    if (particle)
    {
      children_by_parent[particle->get_parent_id()].push_back(particle);
    }
  }

  std::vector<Pi0Candidate> candidates;
  std::map<int, std::size_t> generator_candidate_by_barcode;
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
      const auto found = children_by_parent.find(primary->get_track_id());
      if (found == children_by_parent.end() || found->second.size() != 2U ||
          found->second[0]->get_pid() != 22 ||
          found->second[1]->get_pid() != 22)
      {
        ++n_pi0_malformed_daughters_;
        continue;
      }
      candidates.push_back({Pi0Pathway::g4_primary_decay,
          hepmc_particle->barcode(), hepmc_particle, primary, found->second});
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
    const int barcode = origin.parent->barcode();
    const auto insertion = generator_candidate_by_barcode.emplace(
        barcode, candidates.size());
    if (insertion.second)
    {
      candidates.push_back({Pi0Pathway::generator_decay, barcode,
          origin.parent, nullptr, {}});
    }
    candidates[insertion.first->second].photons.push_back(primary);
  }

  for (const Pi0Candidate& candidate : candidates)
  {
    double parent_eta = 0.0;
    bool valid_parent_eta = finite_eta(candidate.parent, parent_eta);
    if (candidate.g4_parent)
    {
      const double parent_pt = particle_pt(candidate.g4_parent);
      valid_parent_eta = parent_pt > 0.0;
      parent_eta = valid_parent_eta
          ? std::asinh(candidate.g4_parent->get_pz() / parent_pt) : 0.0;
    }
    if (!valid_parent_eta || !std::isfinite(parent_eta) ||
        std::abs(parent_eta) >= truth_eta_max_)
    {
      continue;
    }
    if (candidate.photons.size() != 2U)
    {
      ++n_pi0_malformed_daughters_;
      continue;
    }

    const double parent_pt = candidate.g4_parent
        ? particle_pt(candidate.g4_parent) : hepmc_pt(candidate.parent);
    if (!(parent_pt > 0.0) || !std::isfinite(parent_pt))
    {
      continue;
    }

    ++n_pi0_candidate_;
    const std::size_t pathway = static_cast<std::size_t>(candidate.pathway);
    if (candidate.pathway == Pi0Pathway::g4_primary_decay)
    {
      ++n_pi0_candidate_g4_primary_;
    }
    else
    {
      ++n_pi0_candidate_generator_;
    }

    std::array<int, threshold_count_> multiplicity{};
    std::vector<double> positive_fractions;
    for (const ClusterRecord& cluster : cluster_records)
    {
      if (!cluster.truth.valid)
      {
        continue;
      }
      ++n_pi0_cluster_pair_evaluated_;
      double matching_fraction = 0.0;
      for (const auto& contributor : cluster.truth.contributors)
      {
        if (contributor_matches_pi0(
                contributor, candidate, event_map, signal_embedding_id_))
        {
          matching_fraction += contributor.fraction;
        }
      }
      if (!(matching_fraction > 0.0) || !std::isfinite(matching_fraction))
      {
        continue;
      }
      ++n_pi0_cluster_pair_positive_;
      positive_fractions.push_back(matching_fraction);
      h_fraction_vs_cluster_energy_->Fill(cluster.energy, matching_fraction);
      ++multiplicity[0];
      for (std::size_t threshold = 1; threshold < threshold_count_; ++threshold)
      {
        if (matching_fraction >= thresholds_[threshold])
        {
          ++multiplicity[threshold];
        }
      }
    }

    std::sort(positive_fractions.begin(), positive_fractions.end(),
        std::greater<double>());
    h_maximum_fraction_->Fill(
        positive_fractions.empty() ? 0.0 : positive_fractions[0]);
    h_second_fraction_->Fill(
        positive_fractions.size() < 2U ? 0.0 : positive_fractions[1]);
    for (std::size_t threshold = 0; threshold < threshold_count_; ++threshold)
    {
      h_multiplicity_[threshold]->Fill(multiplicity[threshold]);
      h_multiplicity_vs_pt_[threshold]->Fill(
          parent_pt, multiplicity[threshold]);
      h_pathway_multiplicity_[pathway][threshold]->Fill(
          multiplicity[threshold]);
    }
  }

  ++n_events_written_;
  if (verbosity_ > 0 && n_events_processed_ % 1000ULL == 0ULL)
  {
    std::cout << "PythiaPi0ClusterMultiplicity - processed "
              << n_events_processed_ << " events" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int PythiaPi0ClusterMultiplicity::End(PHCompositeNode* /*topNode*/)
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
  std::cout << "PythiaPi0ClusterMultiplicity - events/pi0/positive pairs = "
            << n_events_written_ << "/" << n_pi0_candidate_ << "/"
            << n_pi0_cluster_pair_positive_ << std::endl;
  return write_error ? Fun4AllReturnCodes::ABORTRUN
                     : Fun4AllReturnCodes::EVENT_OK;
}

void PythiaPi0ClusterMultiplicity::create_output_directory() const
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

void PythiaPi0ClusterMultiplicity::create_output()
{
  output_file_ = TFile::Open(output_file_name_.c_str(), "RECREATE");
  if (!output_file_ || output_file_->IsZombie())
  {
    return;
  }
  output_file_->cd();
  for (std::size_t threshold = 0; threshold < threshold_count_; ++threshold)
  {
    const std::string tag = threshold_tag(threshold);
    h_multiplicity_[threshold] = new TH1D(
        ("h_pi0_cluster_multiplicity_fmin_" + tag + "_raw").c_str(), "",
        multiplicity_max_ + 1, -0.5, multiplicity_max_ + 0.5);
    h_multiplicity_vs_pt_[threshold] = new TH2D(
        ("h_pi0_cluster_multiplicity_vs_truth_pt_fmin_" + tag + "_raw").c_str(),
        "", pt_bins_, 0.0, pt_max_, multiplicity_max_ + 1, -0.5,
        multiplicity_max_ + 0.5);
    h_multiplicity_[threshold]->Sumw2();
    h_multiplicity_vs_pt_[threshold]->Sumw2();
    for (std::size_t pathway = 0; pathway < pathway_count_; ++pathway)
    {
      h_pathway_multiplicity_[pathway][threshold] = new TH1D(
          ("h_pi0_cluster_multiplicity_fmin_" + tag + "_" +
           pathway_tag(pathway) + "_raw").c_str(), "", multiplicity_max_ + 1,
          -0.5, multiplicity_max_ + 0.5);
      h_pathway_multiplicity_[pathway][threshold]->Sumw2();
    }
  }
  constexpr double fraction_axis_max = 1.000001;
  h_maximum_fraction_ = new TH1D(
      "h_pi0_maximum_compatible_fraction_raw", "", 100, 0.0,
      fraction_axis_max);
  h_second_fraction_ = new TH1D(
      "h_pi0_second_compatible_fraction_raw", "", 100, 0.0,
      fraction_axis_max);
  h_fraction_vs_cluster_energy_ = new TH2D(
      "h_pi0_compatible_fraction_vs_cluster_energy_raw", "",
      cluster_energy_bins_, 0.0, cluster_energy_max_, 100, 0.0,
      fraction_axis_max);
  h_maximum_fraction_->Sumw2();
  h_second_fraction_->Sumw2();
  h_fraction_vs_cluster_energy_->Sumw2();

  metadata_tree_ = new TTree(
      "metadata", "Pythia pi0 cluster multiplicity partial metadata");
  static int schema_version = schema_version_;
  static unsigned char cluster_energy_cut_applied = 0U;
  static int truth_matcher_algorithm_version =
      photon_tree::PythiaClusterTruthMatcher::kAlgorithmVersion;
  metadata_tree_->Branch("schema_version", &schema_version);
  metadata_tree_->Branch("manifest_path", &manifest_path_);
  metadata_tree_->Branch("manifest_begin", &manifest_begin_);
  metadata_tree_->Branch("manifest_end", &manifest_end_);
  metadata_tree_->Branch("first_suffix", &first_suffix_);
  metadata_tree_->Branch("last_suffix", &last_suffix_);
  metadata_tree_->Branch("cluster_collection", &cluster_collection_);
  metadata_tree_->Branch("pi0_selection", &pi0_selection_);
  metadata_tree_->Branch("cluster_selection", &cluster_selection_);
  metadata_tree_->Branch("fraction_definition", &fraction_definition_);
  metadata_tree_->Branch("zero_threshold_definition", &zero_threshold_definition_);
  metadata_tree_->Branch("signal_embedding_id", &signal_embedding_id_);
  metadata_tree_->Branch("truth_matcher_algorithm_version",
                         &truth_matcher_algorithm_version);
  metadata_tree_->Branch("pt_bins", &pt_bins_);
  metadata_tree_->Branch("pt_max", &pt_max_);
  metadata_tree_->Branch("multiplicity_max", &multiplicity_max_);
  metadata_tree_->Branch("cluster_energy_bins", &cluster_energy_bins_);
  metadata_tree_->Branch("cluster_energy_max", &cluster_energy_max_);
  metadata_tree_->Branch("truth_eta_max", &truth_eta_max_);
  metadata_tree_->Branch("cluster_eta_max", &cluster_eta_max_);
  metadata_tree_->Branch("cluster_energy_cut_applied",
                         &cluster_energy_cut_applied);
  metadata_tree_->Branch("fraction_thresholds", thresholds_.data(),
                         "fraction_thresholds[4]/D");
  metadata_tree_->Branch("events_processed", &n_events_processed_);
  metadata_tree_->Branch("events_written", &n_events_written_);
  metadata_tree_->Branch("events_invalid", &n_events_invalid_);
  metadata_tree_->Branch("cluster_considered_count", &n_cluster_considered_);
  metadata_tree_->Branch("cluster_invalid_truth_count",
                         &n_cluster_invalid_truth_);
  metadata_tree_->Branch("pi0_candidate_count", &n_pi0_candidate_);
  metadata_tree_->Branch("pi0_candidate_g4_primary_count",
                         &n_pi0_candidate_g4_primary_);
  metadata_tree_->Branch("pi0_candidate_generator_count",
                         &n_pi0_candidate_generator_);
  metadata_tree_->Branch("pi0_malformed_daughters_count",
                         &n_pi0_malformed_daughters_);
  metadata_tree_->Branch("pi0_cluster_pair_evaluated_count",
                         &n_pi0_cluster_pair_evaluated_);
  metadata_tree_->Branch("pi0_cluster_pair_positive_count",
                         &n_pi0_cluster_pair_positive_);
}

void PythiaPi0ClusterMultiplicity::close_output()
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
  h_multiplicity_.fill(nullptr);
  h_multiplicity_vs_pt_.fill(nullptr);
  for (auto& pathway : h_pathway_multiplicity_)
  {
    pathway.fill(nullptr);
  }
  h_maximum_fraction_ = nullptr;
  h_second_fraction_ = nullptr;
  h_fraction_vs_cluster_energy_ = nullptr;
  metadata_tree_ = nullptr;
}
