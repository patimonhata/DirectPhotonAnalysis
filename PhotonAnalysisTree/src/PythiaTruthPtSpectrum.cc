#include "PythiaTruthPtSpectrum.h"

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
#include <TH1D.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <set>
#include <vector>

namespace
{
struct Kinematics
{
  bool valid = false;
  float pt = 0.0F;
  float eta = 0.0F;
};

Kinematics make_kinematics(const double px, const double py, const double pz, const double energy)
{
  Kinematics result;
  const double pt = std::hypot(px, py);
  if (!std::isfinite(px) || !std::isfinite(py) || !std::isfinite(pz) || !std::isfinite(energy) ||
      !std::isfinite(pt) || pt <= 0.0)
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
  result.pt = static_cast<float>(pt);
  result.eta = static_cast<float>(eta);
  return result;
}

Kinematics kinematics(const HepMC::GenParticle* particle)
{
  if (!particle)
  {
    return {};
  }
  const HepMC::FourVector& momentum = particle->momentum();
  return make_kinematics(momentum.px(), momentum.py(), momentum.pz(), momentum.e());
}

Kinematics kinematics(const PHG4Particle* particle)
{
  return particle ? make_kinematics(particle->get_px(), particle->get_py(), particle->get_pz(), particle->get_e())
                  : Kinematics{};
}

bool in_acceptance(const Kinematics& value, const double max_abs_eta)
{
  return value.valid && (max_abs_eta < 0.0 || std::abs(value.eta) < max_abs_eta);
}

std::vector<const HepMC::GenParticle*> incoming(const HepMC::GenVertex* vertex)
{
  std::vector<const HepMC::GenParticle*> result;
  if (!vertex)
  {
    return result;
  }
  for (auto iterator = vertex->particles_in_const_begin(); iterator != vertex->particles_in_const_end(); ++iterator)
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
  for (auto iterator = particle->end_vertex()->particles_out_const_begin(); iterator != particle->end_vertex()->particles_out_const_end(); ++iterator)
  {
    if (*iterator && (*iterator)->pdg_id() == particle->pdg_id())
    {
      return true;
    }
  }
  return false;
}

struct PhotonOrigin
{
  bool valid = true;
  int parent_count = 0;
  int parent_pdg = 0;
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
      result.parent_pdg = parents.size() == 1U ? parents.front()->pdg_id() : 0;
      return result;
    }
    const HepMC::GenParticle* parent = parents.front();
    if (parent->barcode() != 0 && !visited_barcodes.insert(parent->barcode()).second)
    {
      result.valid = false;
      result.parent_count = 0;
      result.parent_pdg = 0;
      return result;
    }
    current = parent;
  }
  result.valid = false;
  return result;
}
}

PythiaTruthPtSpectrum::PythiaTruthPtSpectrum(const std::string& name)
  : SubsysReco(name)
{
}

PythiaTruthPtSpectrum::~PythiaTruthPtSpectrum()
{
  close_output();
}

int PythiaTruthPtSpectrum::Init(PHCompositeNode* /*topNode*/)
{
  const bool valid = !output_file_name_.empty() && !manifest_path_.empty() && !input_file_prefix_.empty() &&
      manifest_begin_ >= 0 && manifest_end_ > manifest_begin_ && !first_suffix_.empty() && !last_suffix_.empty() &&
      signal_embedding_id_ > 0 && n_bins_ > 0 && std::isfinite(pt_max_) && pt_max_ > 0.0 && std::isfinite(max_abs_eta_);
  if (!valid)
  {
    std::cerr << "PythiaTruthPtSpectrum::Init - invalid configuration" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  files_added_ = manifest_end_ - manifest_begin_;
  metadata_use_event_weight_ = use_event_weight_ ? 1U : 0U;
  create_output_directory();
  create_output();
  if (!output_file_ || output_file_->IsZombie())
  {
    std::cerr << "PythiaTruthPtSpectrum::Init - failed to create " << output_file_name_ << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int PythiaTruthPtSpectrum::process_event(PHCompositeNode* topNode)
{
  ++n_events_processed_;
  const auto* event_map = findNode::getClass<PHHepMCGenEventMap>(topNode, hepmc_event_map_node_name_);
  const auto* truth = findNode::getClass<PHG4TruthInfoContainer>(topNode, truth_node_name_);
  const PHHepMCGenEvent* signal_event = event_map ? event_map->get(signal_embedding_id_) : nullptr;
  const HepMC::GenEvent* event = signal_event ? signal_event->getEvent() : nullptr;
  if (!signal_event || !event || !signal_event->is_simulated() || !truth)
  {
    ++n_malformed_events_;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  const auto& event_weights = event->weights().weights();
  const bool valid_weight = std::all_of(event_weights.begin(), event_weights.end(), [](const double weight) { return std::isfinite(weight); });
  if (use_event_weight_ && !valid_weight)
  {
    ++n_invalid_weight_events_;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  const double weight = use_event_weight_ && !event_weights.empty() ? event_weights.front() : 1.0;

  for (auto iterator = event->particles_begin(); iterator != event->particles_end(); ++iterator)
  {
    const HepMC::GenParticle* particle = *iterator;
    if (!particle)
    {
      continue;
    }
    if (particle->pdg_id() == 22 && particle->status() == 1 && !particle->end_vertex())
    {
      const Kinematics momentum = kinematics(particle);
      if (!in_acceptance(momentum, max_abs_eta_))
      {
        continue;
      }
      const auto classification = photon_classifier_.classify(particle);
      if (classification.valid && (classification.category == 1 || classification.category == 2))
      {
        h_prompt_photon_->Fill(momentum.pt, weight);
        ++n_prompt_photon_;
      }
      const PhotonOrigin origin = trace_photon_origin(particle);
      if (origin.valid && origin.parent_count == 1 && origin.parent_pdg == 111)
      {
        h_hepmc_pi0_decay_photon_->Fill(momentum.pt, weight);
        h_pi0_decay_photon_->Fill(momentum.pt, weight);
        ++n_hepmc_pi0_decay_photon_;
        ++n_pi0_decay_photon_;
      }
    }
    if (particle->pdg_id() == 111 && !has_same_pdg_daughter(particle))
    {
      const Kinematics momentum = kinematics(particle);
      if (in_acceptance(momentum, max_abs_eta_))
      {
        h_pi0_->Fill(momentum.pt, weight);
        ++n_pi0_;
      }
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
    if (!parent || parent->get_pid() != 111 || !truth->is_primary(parent) || truth->isEmbeded(parent->get_track_id()) != signal_embedding_id_)
    {
      continue;
    }
    const HepMC::GenParticle* hepmc_parent = event->barcode_to_particle(parent->get_barcode());
    if (!hepmc_parent || hepmc_parent->pdg_id() != 111)
    {
      continue;
    }
    const Kinematics momentum = kinematics(photon);
    if (in_acceptance(momentum, max_abs_eta_))
    {
      h_g4_pi0_decay_photon_->Fill(momentum.pt, weight);
      h_pi0_decay_photon_->Fill(momentum.pt, weight);
      ++n_g4_pi0_decay_photon_;
      ++n_pi0_decay_photon_;
    }
  }

  ++n_events_written_;
  if (verbosity_ > 0 && n_events_processed_ % 100000 == 0)
  {
    std::cout << "PythiaTruthPtSpectrum - processed " << n_events_processed_ << " events" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int PythiaTruthPtSpectrum::End(PHCompositeNode* /*topNode*/)
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
  std::cout << "PythiaTruthPtSpectrum - files/events/written/malformed/invalid-weight/prompt/pi0/pi0 decay (total/HepMC/G4) = "
            << files_added_ << "/" << n_events_processed_ << "/" << n_events_written_ << "/" << n_malformed_events_ << "/"
            << n_invalid_weight_events_ << "/" << n_prompt_photon_ << "/" << n_pi0_ << "/" << n_pi0_decay_photon_ << "/"
            << n_hepmc_pi0_decay_photon_ << "/" << n_g4_pi0_decay_photon_ << std::endl;
  return write_error ? Fun4AllReturnCodes::ABORTRUN : Fun4AllReturnCodes::EVENT_OK;
}

void PythiaTruthPtSpectrum::create_output_directory() const
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

void PythiaTruthPtSpectrum::create_output()
{
  output_file_ = TFile::Open(output_file_name_.c_str(), "RECREATE");
  if (!output_file_ || output_file_->IsZombie())
  {
    return;
  }
  output_file_->cd();
  h_prompt_photon_ = new TH1D("h_prompt_photon_truth_pt_raw", "", n_bins_, 0.0, pt_max_);
  h_pi0_ = new TH1D("h_pi0_truth_pt_raw", "", n_bins_, 0.0, pt_max_);
  h_pi0_decay_photon_ = new TH1D("h_pi0_decay_photon_truth_pt_raw", "", n_bins_, 0.0, pt_max_);
  h_hepmc_pi0_decay_photon_ = new TH1D("h_hepmc_pi0_decay_photon_truth_pt_raw", "", n_bins_, 0.0, pt_max_);
  h_g4_pi0_decay_photon_ = new TH1D("h_g4_pi0_decay_photon_truth_pt_raw", "", n_bins_, 0.0, pt_max_);
  for (TH1D* histogram : {h_prompt_photon_, h_pi0_, h_pi0_decay_photon_, h_hepmc_pi0_decay_photon_, h_g4_pi0_decay_photon_})
  {
    histogram->Sumw2();
  }

  metadata_tree_ = new TTree("metadata", "Pythia truth pT direct-DST partial metadata");
  metadata_tree_->Branch("schema_version", &metadata_schema_version_);
  metadata_tree_->Branch("photon_selection", &photon_selection_);
  metadata_tree_->Branch("pi0_decay_photon_selection", &pi0_decay_photon_selection_);
  metadata_tree_->Branch("manifest_path", &manifest_path_);
  metadata_tree_->Branch("input_file_prefix", &input_file_prefix_);
  metadata_tree_->Branch("manifest_begin", &manifest_begin_);
  metadata_tree_->Branch("manifest_end", &manifest_end_);
  metadata_tree_->Branch("files_added", &files_added_);
  metadata_tree_->Branch("first_suffix", &first_suffix_);
  metadata_tree_->Branch("last_suffix", &last_suffix_);
  metadata_tree_->Branch("signal_embedding_id", &signal_embedding_id_);
  metadata_tree_->Branch("hepmc_event_map_node", &hepmc_event_map_node_name_);
  metadata_tree_->Branch("truth_node", &truth_node_name_);
  metadata_tree_->Branch("events_processed", &n_events_processed_);
  metadata_tree_->Branch("events_written", &n_events_written_);
  metadata_tree_->Branch("n_bins", &n_bins_);
  metadata_tree_->Branch("pt_max", &pt_max_);
  metadata_tree_->Branch("max_abs_eta", &max_abs_eta_);
  metadata_tree_->Branch("use_event_weight", &metadata_use_event_weight_);
  metadata_tree_->Branch("bin_width_normalized", &metadata_bin_width_normalized_);
  metadata_tree_->Branch("prompt_photon_count", &n_prompt_photon_);
  metadata_tree_->Branch("pi0_count", &n_pi0_);
  metadata_tree_->Branch("pi0_decay_photon_count", &n_pi0_decay_photon_);
  metadata_tree_->Branch("hepmc_pi0_decay_photon_count", &n_hepmc_pi0_decay_photon_);
  metadata_tree_->Branch("g4_pi0_decay_photon_count", &n_g4_pi0_decay_photon_);
  metadata_tree_->Branch("malformed_event_count", &n_malformed_events_);
  metadata_tree_->Branch("invalid_weight_event_count", &n_invalid_weight_events_);
}

void PythiaTruthPtSpectrum::close_output()
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
  h_prompt_photon_ = nullptr;
  h_pi0_ = nullptr;
  h_pi0_decay_photon_ = nullptr;
  h_hepmc_pi0_decay_photon_ = nullptr;
  h_g4_pi0_decay_photon_ = nullptr;
  metadata_tree_ = nullptr;
}
