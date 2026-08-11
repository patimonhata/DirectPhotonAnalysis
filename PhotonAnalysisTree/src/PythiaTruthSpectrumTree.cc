#include "PythiaTruthSpectrumTree.h"

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
#include <HepMC/WeightContainer.h>

#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <set>
#include <vector>

namespace
{
constexpr float invalid_float = -999.0F;

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

bool has_same_pdg_daughter(const HepMC::GenParticle* particle)
{
  if (!particle || !particle->end_vertex())
  {
    return false;
  }
  for (auto iterator = particle->end_vertex()->particles_out_const_begin();
       iterator != particle->end_vertex()->particles_out_const_end(); ++iterator)
  {
    if (*iterator && (*iterator)->pdg_id() == particle->pdg_id())
    {
      return true;
    }
  }
  return false;
}

unsigned int count_direct_photon_daughters(const HepMC::GenParticle* particle)
{
  if (!particle || !particle->end_vertex())
  {
    return 0U;
  }
  unsigned int count = 0U;
  for (auto iterator = particle->end_vertex()->particles_out_const_begin();
       iterator != particle->end_vertex()->particles_out_const_end(); ++iterator)
  {
    count += (*iterator && (*iterator)->pdg_id() == 22) ? 1U : 0U;
  }
  return count;
}

struct Kinematics
{
  bool valid = false;
  float e = invalid_float;
  float pt = invalid_float;
  float eta = invalid_float;
  float phi = invalid_float;
};

Kinematics kinematics(const HepMC::GenParticle* particle)
{
  Kinematics result;
  if (!particle)
  {
    return result;
  }
  const HepMC::FourVector& momentum = particle->momentum();
  const double px = momentum.px();
  const double py = momentum.py();
  const double pz = momentum.pz();
  const double energy = momentum.e();
  const double pt = std::hypot(px, py);
  if (!std::isfinite(px) || !std::isfinite(py) || !std::isfinite(pz) ||
      !std::isfinite(energy) || !std::isfinite(pt) || pt <= 0.0)
  {
    return result;
  }
  const double eta = std::asinh(pz / pt);
  const double phi = std::atan2(py, px);
  if (!std::isfinite(eta) || !std::isfinite(phi))
  {
    return result;
  }
  result.valid = true;
  result.e = static_cast<float>(energy);
  result.pt = static_cast<float>(pt);
  result.eta = static_cast<float>(eta);
  result.phi = static_cast<float>(phi);
  return result;
}

Kinematics kinematics(const PHG4Particle* particle)
{
  Kinematics result;
  if (!particle)
  {
    return result;
  }
  const double px = particle->get_px();
  const double py = particle->get_py();
  const double pz = particle->get_pz();
  const double energy = particle->get_e();
  const double pt = std::hypot(px, py);
  if (!std::isfinite(px) || !std::isfinite(py) || !std::isfinite(pz) ||
      !std::isfinite(energy) || !std::isfinite(pt) || pt <= 0.0)
  {
    return result;
  }
  const double eta = std::asinh(pz / pt);
  const double phi = std::atan2(py, px);
  if (!std::isfinite(eta) || !std::isfinite(phi))
  {
    return result;
  }
  result.valid = true;
  result.e = static_cast<float>(energy);
  result.pt = static_cast<float>(pt);
  result.eta = static_cast<float>(eta);
  result.phi = static_cast<float>(phi);
  return result;
}

struct PhotonOrigin
{
  bool valid = true;
  unsigned int copy_depth = 0U;
  int parent_count = 0;
  int parent_pdg = 0;
  int parent_barcode = 0;
};

PhotonOrigin trace_photon_origin(const HepMC::GenParticle* photon)
{
  PhotonOrigin result;
  const HepMC::GenParticle* current = photon;
  std::set<int> visited_barcodes;
  if (current && current->barcode() != 0)
  {
    visited_barcodes.insert(current->barcode());
  }
  while (current)
  {
    const auto parents = incoming(current->production_vertex());
    if (parents.size() != 1U || parents.front()->pdg_id() != 22)
    {
      result.parent_count = static_cast<int>(parents.size());
      if (parents.size() == 1U)
      {
        result.parent_pdg = parents.front()->pdg_id();
        result.parent_barcode = parents.front()->barcode();
      }
      return result;
    }
    const HepMC::GenParticle* parent = parents.front();
    if (parent->barcode() != 0 && !visited_barcodes.insert(parent->barcode()).second)
    {
      result.valid = false;
      result.parent_count = 0;
      result.parent_pdg = 0;
      result.parent_barcode = 0;
      return result;
    }
    current = parent;
    ++result.copy_depth;
  }
  result.valid = false;
  return result;
}
}

PythiaTruthSpectrumTree::PythiaTruthSpectrumTree(const std::string& name)
  : SubsysReco(name)
{
}

PythiaTruthSpectrumTree::~PythiaTruthSpectrumTree()
{
  close_output_file();
}

int PythiaTruthSpectrumTree::Init(PHCompositeNode* /*topNode*/)
{
  if (signal_embedding_id_ <= 0 || output_file_name_.empty())
  {
    std::cerr << "PythiaTruthSpectrumTree::Init - invalid configuration" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  create_output_directory();
  output_file_ = TFile::Open(output_file_name_.c_str(), "RECREATE");
  if (!output_file_ || output_file_->IsZombie())
  {
    std::cerr << "PythiaTruthSpectrumTree::Init - failed to create "
              << output_file_name_ << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  create_trees();
  std::cout << "PythiaTruthSpectrumTree::Init - input: " << input_file_name_ << '\n'
            << "PythiaTruthSpectrumTree::Init - output: " << output_file_name_ << '\n'
            << "PythiaTruthSpectrumTree::Init - signal embedding ID: "
            << signal_embedding_id_ << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

int PythiaTruthSpectrumTree::process_event(PHCompositeNode* topNode)
{
  reset_event();
  b_source_file_id_ = source_file_id_;
  b_event_in_file_ = static_cast<unsigned int>(n_events_processed_);
  b_event_uid_ = (static_cast<unsigned long long>(source_file_id_) << 32U) |
                 static_cast<unsigned long long>(b_event_in_file_);
  ++n_events_processed_;

  auto* event_map = findNode::getClass<PHHepMCGenEventMap>(
      topNode, hepmc_event_map_node_name_);
  auto* truth = findNode::getClass<PHG4TruthInfoContainer>(topNode, truth_node_name_);
  if (!fill_truth(event_map, truth))
  {
    ++n_events_invalid_truth_;
    if (verbosity_ > 0)
    {
      std::cerr << "PythiaTruthSpectrumTree::process_event - missing/invalid signal HepMC event "
                << signal_embedding_id_ << " in event " << b_event_in_file_ << std::endl;
    }
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  event_tree_->Fill();
  ++n_events_written_;
  if (verbosity_ > 0 && b_event_in_file_ < 5U)
  {
    std::cout << "PythiaTruthSpectrumTree::process_event - event " << b_event_in_file_
              << " particles/final photons/terminal pi0 = "
              << b_hepmc_n_particle_record_ << "/" << b_truth_photon_n_
              << "/" << b_truth_pi0_n_ << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int PythiaTruthSpectrumTree::End(PHCompositeNode* /*topNode*/)
{
  close_output_file();
  std::cout << "PythiaTruthSpectrumTree::End - processed/written/invalid_truth = "
            << n_events_processed_ << "/" << n_events_written_ << "/"
            << n_events_invalid_truth_ << '\n'
            << "PythiaTruthSpectrumTree::End - records/final photons/terminal pi0 = "
            << n_hepmc_particle_record_ << "/" << n_final_photon_ << "/"
            << n_terminal_pi0_ << '\n'
            << "PythiaTruthSpectrumTree::End - G4 pi0 decay photons/photon copy edges/"
               "nonterminal pi0 copies/invalid ancestry/invalid kinematics/unmatched G4 photons = "
            << n_g4_pi0_decay_photon_ << "/"
            << n_photon_copy_edge_ << "/" << n_nonterminal_pi0_copy_ << "/"
            << n_invalid_photon_ancestry_ << "/" << n_invalid_kinematics_ << "/"
            << n_unmatched_g4_pi0_decay_photon_ << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

bool PythiaTruthSpectrumTree::fill_truth(const PHHepMCGenEventMap* event_map,
                                         const PHG4TruthInfoContainer* truth)
{
  const PHHepMCGenEvent* signal_event = event_map ? event_map->get(signal_embedding_id_) : nullptr;
  const HepMC::GenEvent* event = signal_event ? signal_event->getEvent() : nullptr;
  if (!signal_event || !event || !signal_event->is_simulated() || !truth)
  {
    return false;
  }

  b_hepmc_event_number_ = event->event_number();
  b_event_weights_ = event->weights().weights();
  const bool finite_weights = std::all_of(
      b_event_weights_.begin(), b_event_weights_.end(),
      [](const double weight) { return std::isfinite(weight); });
  b_event_weight_valid_ = finite_weights ? 1U : 0U;
  b_event_weight_ = (!b_event_weights_.empty() && std::isfinite(b_event_weights_.front()))
      ? b_event_weights_.front() : 1.0;

  for (auto iterator = event->particles_begin(); iterator != event->particles_end(); ++iterator)
  {
    const HepMC::GenParticle* particle = *iterator;
    if (!particle)
    {
      continue;
    }
    ++b_hepmc_n_particle_record_;
    ++n_hepmc_particle_record_;

    if (particle->pdg_id() == 22 && particle->status() == 1 && !particle->end_vertex())
    {
      const Kinematics momentum = kinematics(particle);
      const auto classification = photon_classifier_.classify(particle);
      const auto immediate_parents = incoming(particle->production_vertex());
      const PhotonOrigin origin = trace_photon_origin(particle);

      b_truth_photon_barcode_.push_back(particle->barcode());
      b_truth_photon_status_.push_back(particle->status());
      b_truth_photon_kinematics_valid_.push_back(momentum.valid ? 1U : 0U);
      b_truth_photon_e_.push_back(momentum.e);
      b_truth_photon_pt_.push_back(momentum.pt);
      b_truth_photon_eta_.push_back(momentum.eta);
      b_truth_photon_phi_.push_back(momentum.phi);
      b_truth_photon_classification_valid_.push_back(classification.valid ? 1U : 0U);
      b_truth_photon_category_.push_back(classification.category);
      b_truth_photon_immediate_parent_count_.push_back(
          static_cast<int>(immediate_parents.size()));
      b_truth_photon_immediate_parent_pdg_.push_back(
          immediate_parents.size() == 1U ? immediate_parents.front()->pdg_id() : 0);
      b_truth_photon_copy_chain_valid_.push_back(origin.valid ? 1U : 0U);
      b_truth_photon_copy_depth_.push_back(origin.copy_depth);
      b_truth_photon_origin_parent_count_.push_back(origin.parent_count);
      b_truth_photon_origin_parent_pdg_.push_back(origin.parent_pdg);
      b_truth_photon_origin_parent_barcode_.push_back(origin.parent_barcode);
      ++b_truth_photon_n_;
      ++n_final_photon_;
      n_photon_copy_edge_ += origin.copy_depth;
      n_invalid_photon_ancestry_ += origin.valid ? 0U : 1U;
      n_invalid_kinematics_ += momentum.valid ? 0U : 1U;
    }

    if (particle->pdg_id() == 111)
    {
      if (has_same_pdg_daughter(particle))
      {
        ++n_nonterminal_pi0_copy_;
        continue;
      }
      const Kinematics momentum = kinematics(particle);
      b_truth_pi0_barcode_.push_back(particle->barcode());
      b_truth_pi0_status_.push_back(particle->status());
      b_truth_pi0_kinematics_valid_.push_back(momentum.valid ? 1U : 0U);
      b_truth_pi0_e_.push_back(momentum.e);
      b_truth_pi0_pt_.push_back(momentum.pt);
      b_truth_pi0_eta_.push_back(momentum.eta);
      b_truth_pi0_phi_.push_back(momentum.phi);
      b_truth_pi0_hepmc_direct_photon_count_.push_back(count_direct_photon_daughters(particle));
      ++b_truth_pi0_n_;
      ++n_terminal_pi0_;
      n_invalid_kinematics_ += momentum.valid ? 0U : 1U;
    }
  }

  const auto secondary_range = truth->GetSecondaryParticleRange();
  for (auto iterator = secondary_range.first; iterator != secondary_range.second; ++iterator)
  {
    const PHG4Particle* photon = iterator->second;
    if (!photon || photon->get_pid() != 22)
    {
      continue;
    }
    const PHG4Particle* parent = truth->GetParticle(photon->get_parent_id());
    if (!parent || parent->get_pid() != 111 || !truth->is_primary(parent))
    {
      continue;
    }
    const int parent_embedding_id = truth->isEmbeded(parent->get_track_id());
    const HepMC::GenParticle* hepmc_parent = event->barcode_to_particle(parent->get_barcode());
    if (parent_embedding_id != signal_embedding_id_ || !hepmc_parent ||
        hepmc_parent->pdg_id() != 111)
    {
      ++n_unmatched_g4_pi0_decay_photon_;
      continue;
    }
    const Kinematics momentum = kinematics(photon);
    b_truth_pi0_decay_photon_g4_track_id_.push_back(photon->get_track_id());
    b_truth_pi0_decay_photon_parent_g4_track_id_.push_back(parent->get_track_id());
    b_truth_pi0_decay_photon_parent_hepmc_barcode_.push_back(parent->get_barcode());
    b_truth_pi0_decay_photon_kinematics_valid_.push_back(momentum.valid ? 1U : 0U);
    b_truth_pi0_decay_photon_e_.push_back(momentum.e);
    b_truth_pi0_decay_photon_pt_.push_back(momentum.pt);
    b_truth_pi0_decay_photon_eta_.push_back(momentum.eta);
    b_truth_pi0_decay_photon_phi_.push_back(momentum.phi);
    ++b_truth_pi0_decay_photon_n_;
    ++n_g4_pi0_decay_photon_;
    n_invalid_kinematics_ += momentum.valid ? 0U : 1U;
  }
  return true;
}

void PythiaTruthSpectrumTree::reset_event()
{
  b_source_file_id_ = 0U;
  b_event_in_file_ = 0U;
  b_event_uid_ = 0ULL;
  b_hepmc_event_number_ = -999;
  b_event_weight_valid_ = 0U;
  b_event_weight_ = 1.0;
  b_event_weights_.clear();
  b_hepmc_n_particle_record_ = 0U;

  b_truth_photon_n_ = 0U;
  b_truth_photon_barcode_.clear();
  b_truth_photon_status_.clear();
  b_truth_photon_kinematics_valid_.clear();
  b_truth_photon_e_.clear();
  b_truth_photon_pt_.clear();
  b_truth_photon_eta_.clear();
  b_truth_photon_phi_.clear();
  b_truth_photon_classification_valid_.clear();
  b_truth_photon_category_.clear();
  b_truth_photon_immediate_parent_count_.clear();
  b_truth_photon_immediate_parent_pdg_.clear();
  b_truth_photon_copy_chain_valid_.clear();
  b_truth_photon_copy_depth_.clear();
  b_truth_photon_origin_parent_count_.clear();
  b_truth_photon_origin_parent_pdg_.clear();
  b_truth_photon_origin_parent_barcode_.clear();

  b_truth_pi0_n_ = 0U;
  b_truth_pi0_barcode_.clear();
  b_truth_pi0_status_.clear();
  b_truth_pi0_kinematics_valid_.clear();
  b_truth_pi0_e_.clear();
  b_truth_pi0_pt_.clear();
  b_truth_pi0_eta_.clear();
  b_truth_pi0_phi_.clear();
  b_truth_pi0_hepmc_direct_photon_count_.clear();

  b_truth_pi0_decay_photon_n_ = 0U;
  b_truth_pi0_decay_photon_g4_track_id_.clear();
  b_truth_pi0_decay_photon_parent_g4_track_id_.clear();
  b_truth_pi0_decay_photon_parent_hepmc_barcode_.clear();
  b_truth_pi0_decay_photon_kinematics_valid_.clear();
  b_truth_pi0_decay_photon_e_.clear();
  b_truth_pi0_decay_photon_pt_.clear();
  b_truth_pi0_decay_photon_eta_.clear();
  b_truth_pi0_decay_photon_phi_.clear();
}

void PythiaTruthSpectrumTree::create_trees()
{
  output_file_->cd();
  event_tree_ = new TTree("event_tree", "Compact Pythia HepMC truth spectrum tree");
  event_tree_->Branch("source_file_id", &b_source_file_id_);
  event_tree_->Branch("event_in_file", &b_event_in_file_);
  event_tree_->Branch("event_uid", &b_event_uid_);
  event_tree_->Branch("hepmc_event_number", &b_hepmc_event_number_);
  event_tree_->Branch("event_weight_valid", &b_event_weight_valid_);
  event_tree_->Branch("event_weight", &b_event_weight_);
  event_tree_->Branch("event_weights", &b_event_weights_);
  event_tree_->Branch("hepmc_n_particle_record", &b_hepmc_n_particle_record_);

  event_tree_->Branch("truth_photon_n", &b_truth_photon_n_);
  event_tree_->Branch("truth_photon_barcode", &b_truth_photon_barcode_);
  event_tree_->Branch("truth_photon_status", &b_truth_photon_status_);
  event_tree_->Branch("truth_photon_kinematics_valid", &b_truth_photon_kinematics_valid_);
  event_tree_->Branch("truth_photon_e", &b_truth_photon_e_);
  event_tree_->Branch("truth_photon_pt", &b_truth_photon_pt_);
  event_tree_->Branch("truth_photon_eta", &b_truth_photon_eta_);
  event_tree_->Branch("truth_photon_phi", &b_truth_photon_phi_);
  event_tree_->Branch("truth_photon_classification_valid", &b_truth_photon_classification_valid_);
  event_tree_->Branch("truth_photon_category", &b_truth_photon_category_);
  event_tree_->Branch("truth_photon_immediate_parent_count", &b_truth_photon_immediate_parent_count_);
  event_tree_->Branch("truth_photon_immediate_parent_pdg", &b_truth_photon_immediate_parent_pdg_);
  event_tree_->Branch("truth_photon_copy_chain_valid", &b_truth_photon_copy_chain_valid_);
  event_tree_->Branch("truth_photon_copy_depth", &b_truth_photon_copy_depth_);
  event_tree_->Branch("truth_photon_origin_parent_count", &b_truth_photon_origin_parent_count_);
  event_tree_->Branch("truth_photon_origin_parent_pdg", &b_truth_photon_origin_parent_pdg_);
  event_tree_->Branch("truth_photon_origin_parent_barcode", &b_truth_photon_origin_parent_barcode_);

  event_tree_->Branch("truth_pi0_n", &b_truth_pi0_n_);
  event_tree_->Branch("truth_pi0_barcode", &b_truth_pi0_barcode_);
  event_tree_->Branch("truth_pi0_status", &b_truth_pi0_status_);
  event_tree_->Branch("truth_pi0_kinematics_valid", &b_truth_pi0_kinematics_valid_);
  event_tree_->Branch("truth_pi0_e", &b_truth_pi0_e_);
  event_tree_->Branch("truth_pi0_pt", &b_truth_pi0_pt_);
  event_tree_->Branch("truth_pi0_eta", &b_truth_pi0_eta_);
  event_tree_->Branch("truth_pi0_phi", &b_truth_pi0_phi_);
  event_tree_->Branch("truth_pi0_hepmc_direct_photon_count", &b_truth_pi0_hepmc_direct_photon_count_);
  event_tree_->Branch("truth_pi0_decay_photon_n", &b_truth_pi0_decay_photon_n_);
  event_tree_->Branch("truth_pi0_decay_photon_g4_track_id", &b_truth_pi0_decay_photon_g4_track_id_);
  event_tree_->Branch("truth_pi0_decay_photon_parent_g4_track_id", &b_truth_pi0_decay_photon_parent_g4_track_id_);
  event_tree_->Branch("truth_pi0_decay_photon_parent_hepmc_barcode", &b_truth_pi0_decay_photon_parent_hepmc_barcode_);
  event_tree_->Branch("truth_pi0_decay_photon_kinematics_valid", &b_truth_pi0_decay_photon_kinematics_valid_);
  event_tree_->Branch("truth_pi0_decay_photon_e", &b_truth_pi0_decay_photon_e_);
  event_tree_->Branch("truth_pi0_decay_photon_pt", &b_truth_pi0_decay_photon_pt_);
  event_tree_->Branch("truth_pi0_decay_photon_eta", &b_truth_pi0_decay_photon_eta_);
  event_tree_->Branch("truth_pi0_decay_photon_phi", &b_truth_pi0_decay_photon_phi_);

  metadata_tree_ = new TTree("metadata", "One entry per Pythia truth source DST");
  metadata_tree_->Branch("schema_version", &metadata_schema_version_);
  metadata_tree_->Branch("sample_type", &sample_type_);
  metadata_tree_->Branch("input_file", &input_file_name_);
  metadata_tree_->Branch("output_file", &output_file_name_);
  metadata_tree_->Branch("source_file_id", &source_file_id_);
  metadata_tree_->Branch("signal_embedding_id", &signal_embedding_id_);
  metadata_tree_->Branch("hepmc_event_map_node", &hepmc_event_map_node_name_);
  metadata_tree_->Branch("truth_node", &truth_node_name_);
  metadata_tree_->Branch("photon_classifier_version", &metadata_classifier_version_);
  metadata_tree_->Branch("photon_selection", &photon_selection_);
  metadata_tree_->Branch("pi0_selection", &pi0_selection_);
  metadata_tree_->Branch("photon_origin_scheme", &photon_origin_scheme_);
  metadata_tree_->Branch("pi0_decay_photon_scheme", &pi0_decay_photon_scheme_);
  metadata_tree_->Branch("event_weight_scheme", &event_weight_scheme_);
  metadata_tree_->Branch("n_events_processed", &n_events_processed_);
  metadata_tree_->Branch("n_events_written", &n_events_written_);
  metadata_tree_->Branch("n_events_invalid_truth", &n_events_invalid_truth_);
  metadata_tree_->Branch("n_hepmc_particle_record", &n_hepmc_particle_record_);
  metadata_tree_->Branch("n_final_photon", &n_final_photon_);
  metadata_tree_->Branch("n_terminal_pi0", &n_terminal_pi0_);
  metadata_tree_->Branch("n_g4_pi0_decay_photon", &n_g4_pi0_decay_photon_);
  metadata_tree_->Branch("n_photon_copy_edge", &n_photon_copy_edge_);
  metadata_tree_->Branch("n_nonterminal_pi0_copy", &n_nonterminal_pi0_copy_);
  metadata_tree_->Branch("n_invalid_photon_ancestry", &n_invalid_photon_ancestry_);
  metadata_tree_->Branch("n_invalid_kinematics", &n_invalid_kinematics_);
  metadata_tree_->Branch("n_unmatched_g4_pi0_decay_photon", &n_unmatched_g4_pi0_decay_photon_);
}

void PythiaTruthSpectrumTree::create_output_directory() const
{
  const auto slash = output_file_name_.find_last_of('/');
  if (slash == std::string::npos)
  {
    return;
  }
  const std::string directory = output_file_name_.substr(0, slash);
  if (!directory.empty())
  {
    gSystem->mkdir(directory.c_str(), true);
  }
}

void PythiaTruthSpectrumTree::close_output_file()
{
  if (!output_file_)
  {
    return;
  }
  output_file_->cd();
  if (metadata_tree_ && !metadata_filled_)
  {
    metadata_tree_->Fill();
    metadata_filled_ = true;
  }
  if (event_tree_)
  {
    event_tree_->Write();
  }
  if (metadata_tree_)
  {
    metadata_tree_->Write();
  }
  output_file_->Close();
  delete output_file_;
  output_file_ = nullptr;
  event_tree_ = nullptr;
  metadata_tree_ = nullptr;
}
