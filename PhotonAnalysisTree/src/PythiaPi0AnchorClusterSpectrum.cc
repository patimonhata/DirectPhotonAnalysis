#include "PythiaPi0AnchorClusterSpectrum.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <TFile.h>
#include <TH1D.h>
#include <TSystem.h>
#include <TTree.h>

#include <cmath>
#include <cstddef>
#include <iostream>
#include <string>

PythiaPi0AnchorClusterSpectrum::PythiaPi0AnchorClusterSpectrum(const std::string& name)
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
      std::isfinite(anchor_cluster_eta_max_) &&
      anchor_cluster_eta_max_ > 0.0 &&
      std::isfinite(partner_cluster_eta_max_) &&
      std::isfinite(cemc_acceptance_eta_max_) && cemc_acceptance_eta_max_ > 0.0 &&
      std::isfinite(min_cluster_energy_) && min_cluster_energy_ >= 0.0 &&
      dominant_fraction_min_ >= 0.0 && dominant_fraction_min_ <= 1.0 &&
      anchor_pi0_fraction_min_ >= 0.0 &&
      anchor_pi0_fraction_min_ <= 1.0 &&
      min_energy_contribution_fraction_ >= 0.0 &&
      min_energy_contribution_fraction_ < 1.0 &&
      std::isfinite(min_photon_energy_recovery_) &&
      min_photon_energy_recovery_ >= 0.0 &&
      min_photon_energy_recovery_ <= 1.0 &&
      min_direct_match_cluster_energy_coverage_ >= 0.0 && min_direct_match_cluster_energy_coverage_ <= 1.0 &&
      std::isfinite(missing_diagnostic_max_delta_r_) && missing_diagnostic_max_delta_r_ > 0.0 &&
      std::isfinite(max_abs_vertex_z_) && max_abs_vertex_z_ > 0.0;
  if (!valid)
  {
    std::cerr << "PythiaPi0AnchorClusterSpectrum::Init - invalid configuration"
              << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  photon_tree::Pi0AnchorTopologyConfig topology_config;
  topology_config.sample_mode = photon_tree::Pi0SampleMode::pythia;
  topology_config.truth_node_name = truth_node_name_;
  topology_config.hepmc_event_map_node_name = hepmc_event_map_node_name_;
  topology_config.tower_node_name = tower_node_name_;
  topology_config.raw_truth_tower_node_name = raw_truth_tower_node_name_;
  topology_config.truth_cell_node_name = truth_cell_node_name_;
  topology_config.truth_hit_node_name = truth_hit_node_name_;
  topology_config.cluster_node_name = split_cluster_node_name_;
  topology_config.tower_geom_node_name = tower_geom_node_name_;
  topology_config.signal_embedding_id = signal_embedding_id_;
  topology_config.anchor_cluster_eta_max = anchor_cluster_eta_max_;
  topology_config.partner_cluster_eta_max = partner_cluster_eta_max_;
  topology_config.cemc_acceptance_eta_max = cemc_acceptance_eta_max_;
  topology_config.min_cluster_energy = min_cluster_energy_;
  topology_config.dominant_fraction_min = dominant_fraction_min_;
  topology_config.anchor_pi0_fraction_min = anchor_pi0_fraction_min_;
  topology_config.min_energy_contribution_fraction = min_energy_contribution_fraction_;
  topology_config.min_photon_energy_recovery = min_photon_energy_recovery_;
  topology_config.min_direct_match_cluster_energy_coverage = min_direct_match_cluster_energy_coverage_;
  topology_config.missing_diagnostic_max_delta_r = missing_diagnostic_max_delta_r_;
  topology_config.enable_missing_diagnostics = enable_missing_diagnostics_;
  topology_config.max_abs_vertex_z = max_abs_vertex_z_;
  topology_config.verbosity = verbosity_;
  topology_evaluator_.configure(topology_config);
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
  const photon_tree::Pi0AnchorTopologyEventResult result = topology_evaluator_.evaluate(topNode);
  if (result.status == photon_tree::Pi0TopologyEventStatus::vertex_rejected)
  {
    ++n_events_vertex_rejected_;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  if (result.status != photon_tree::Pi0TopologyEventStatus::accepted)
  {
    ++n_events_invalid_;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  n_cluster_considered_ += result.clusters.size();
  n_cluster_invalid_truth_ += result.cluster_invalid_truth_count;
  n_pi0_candidate_g4_decay_ += result.g4_candidate_count;
  n_pi0_candidate_generator_decay_ += result.generator_candidate_count;
  n_pi0_malformed_daughters_ += result.malformed_candidate_count;
  n_energy_match_invalid_ += result.energy_match_invalid_count;

  for (std::size_t index = 0; index < result.clusters.size(); ++index)
  {
    if (index < result.prompt_cluster.size() && result.prompt_cluster[index])
    {
      h_prompt_->Fill(result.clusters[index].et);
      ++n_prompt_cluster_;
    }
  }

  for (const photon_tree::Pi0TopologyAnchorRecord& anchor : result.anchors)
  {
    if (anchor.cluster_index >= result.clusters.size() ||
        anchor.candidate_index >= result.candidates.size())
    {
      ++n_events_invalid_;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    const double et = result.clusters[anchor.cluster_index].et;
    h_anchor_->Fill(et);
    ++n_anchor_cluster_;
    if (result.candidates[anchor.candidate_index].pathway ==
        photon_tree::Pi0Pathway::g4_primary_decay)
    {
      ++n_anchor_g4_decay_;
    }
    else
    {
      ++n_anchor_generator_decay_;
    }
    if (anchor.ambiguous_main)
    {
      ++n_anchor_ambiguous_main_;
    }
    switch (anchor.topology)
    {
    case photon_tree::Pi0AnchorTopology::separated:
      h_separated_->Fill(et);
      ++n_separated_;
      break;
    case photon_tree::Pi0AnchorTopology::merged:
      h_merged_->Fill(et);
      ++n_merged_;
      break;
    case photon_tree::Pi0AnchorTopology::missing:
      h_missing_->Fill(et);
      ++n_missing_;
      switch (anchor.missing_category)
      {
      case photon_tree::Pi0MissingCategory::energy_threshold:
        h_missing_energy_threshold_->Fill(et);
        ++n_missing_energy_threshold_;
        break;
      case photon_tree::Pi0MissingCategory::acceptance:
        h_missing_acceptance_->Fill(et);
        ++n_missing_acceptance_;
        break;
      case photon_tree::Pi0MissingCategory::other:
        h_missing_other_->Fill(et);
        ++n_missing_other_;
        break;
      case photon_tree::Pi0MissingCategory::not_missing:
        ++n_events_invalid_;
        return Fun4AllReturnCodes::ABORTEVENT;
      }
      break;
    case photon_tree::Pi0AnchorTopology::other:
      h_other_->Fill(et);
      ++n_other_;
      break;
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
      << "PythiaPi0AnchorClusterSpectrum - processed/written/vertex-rejected/invalid/anchors/separated/merged/missing(energy/acceptance/other)/other = "
      << n_events_processed_ << "/" << n_events_written_ << "/"
      << n_events_vertex_rejected_ << "/" << n_events_invalid_ << "/" << n_anchor_cluster_ << "/"
      << n_separated_ << "/" << n_merged_ << "/" << n_missing_ << "("
      << n_missing_energy_threshold_ << "/" << n_missing_acceptance_ << "/" << n_missing_other_ << ")/"
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
  h_missing_energy_threshold_ = new TH1D(
      "h_pi0_anchor_missing_energy_threshold_cluster_et_raw", "", n_bins_, 0.0, et_max_);
  h_missing_acceptance_ = new TH1D(
      "h_pi0_anchor_missing_acceptance_cluster_et_raw", "", n_bins_, 0.0, et_max_);
  h_missing_other_ = new TH1D(
      "h_pi0_anchor_missing_other_cluster_et_raw", "", n_bins_, 0.0, et_max_);
  h_other_ = new TH1D(
      "h_pi0_anchor_other_cluster_et_raw", "",
      n_bins_, 0.0, et_max_);
  for (TH1D* histogram : {
           h_prompt_, h_anchor_, h_separated_, h_merged_, h_missing_,
           h_missing_energy_threshold_, h_missing_acceptance_, h_missing_other_, h_other_})
  {
    histogram->Sumw2();
  }

  metadata_tree_ = new TTree("metadata", "Pythia pi0 anchor-cluster partial metadata");
  static int schema_version = schema_version_;
  static unsigned char bin_width_normalized = 0U;
  metadata_tree_->Branch("schema_version", &schema_version);
  metadata_tree_->Branch("manifest_path", &manifest_path_);
  metadata_tree_->Branch("manifest_begin", &manifest_begin_);
  metadata_tree_->Branch("manifest_end", &manifest_end_);
  metadata_tree_->Branch("first_suffix", &first_suffix_);
  metadata_tree_->Branch("last_suffix", &last_suffix_);
  metadata_tree_->Branch("cluster_collection", &cluster_collection_);
  metadata_tree_->Branch("tower_geom_node", &tower_geom_node_name_);
  metadata_tree_->Branch("classification_unit", &classification_unit_);
  metadata_tree_->Branch("pi0_selection", &pi0_selection_);
  metadata_tree_->Branch("partner_selection", &partner_selection_);
  metadata_tree_->Branch("topology_definition", &topology_definition_);
  metadata_tree_->Branch("topology_priority", &topology_priority_);
  metadata_tree_->Branch("missing_category_priority", &missing_category_priority_);
  metadata_tree_->Branch("response_policy", &response_policy_);
  metadata_tree_->Branch("photon_recovery_policy", &photon_recovery_policy_);
  metadata_tree_->Branch("vertex_selection", &vertex_selection_);
  metadata_tree_->Branch("signal_embedding_id", &signal_embedding_id_);
  metadata_tree_->Branch("n_bins", &n_bins_);
  metadata_tree_->Branch("et_max", &et_max_);
  metadata_tree_->Branch("anchor_cluster_eta_max", &anchor_cluster_eta_max_);
  metadata_tree_->Branch("partner_cluster_eta_max", &partner_cluster_eta_max_);
  metadata_tree_->Branch("cemc_acceptance_eta_max", &cemc_acceptance_eta_max_);
  metadata_tree_->Branch("min_cluster_energy", &min_cluster_energy_);
  metadata_tree_->Branch("dominant_fraction_min", &dominant_fraction_min_);
  metadata_tree_->Branch("anchor_pi0_fraction_min", &anchor_pi0_fraction_min_);
  metadata_tree_->Branch("min_energy_contribution_fraction", &min_energy_contribution_fraction_);
  metadata_tree_->Branch("min_photon_energy_recovery", &min_photon_energy_recovery_);
  metadata_tree_->Branch("min_direct_match_cluster_energy_coverage", &min_direct_match_cluster_energy_coverage_);
  metadata_tree_->Branch("missing_diagnostic_max_delta_r", &missing_diagnostic_max_delta_r_);
  metadata_tree_->Branch("enable_missing_diagnostics", &enable_missing_diagnostics_);
  metadata_tree_->Branch("max_abs_vertex_z", &max_abs_vertex_z_);
  metadata_tree_->Branch("pi0_truth_matching_algorithm_version", &pi0_truth_matching_algorithm_version_);
  metadata_tree_->Branch("pi0_topology_algorithm_version", &pi0_topology_algorithm_version_);
  metadata_tree_->Branch("bin_width_normalized", &bin_width_normalized);
  metadata_tree_->Branch("events_processed", &n_events_processed_);
  metadata_tree_->Branch("events_written", &n_events_written_);
  metadata_tree_->Branch("events_invalid", &n_events_invalid_);
  metadata_tree_->Branch("events_vertex_rejected", &n_events_vertex_rejected_);
  metadata_tree_->Branch("cluster_considered_count", &n_cluster_considered_);
  metadata_tree_->Branch("cluster_invalid_truth_count", &n_cluster_invalid_truth_);
  metadata_tree_->Branch("prompt_cluster_count", &n_prompt_cluster_);
  metadata_tree_->Branch("pi0_candidate_g4_decay_count", &n_pi0_candidate_g4_decay_);
  metadata_tree_->Branch("pi0_candidate_generator_decay_count", &n_pi0_candidate_generator_decay_);
  metadata_tree_->Branch("pi0_malformed_daughters_count", &n_pi0_malformed_daughters_);
  metadata_tree_->Branch("anchor_cluster_count", &n_anchor_cluster_);
  metadata_tree_->Branch("anchor_g4_decay_count", &n_anchor_g4_decay_);
  metadata_tree_->Branch("anchor_generator_decay_count", &n_anchor_generator_decay_);
  metadata_tree_->Branch("anchor_ambiguous_main_count", &n_anchor_ambiguous_main_);
  metadata_tree_->Branch("energy_match_invalid_count", &n_energy_match_invalid_);
  metadata_tree_->Branch("separated_count", &n_separated_);
  metadata_tree_->Branch("merged_count", &n_merged_);
  metadata_tree_->Branch("missing_count", &n_missing_);
  metadata_tree_->Branch("missing_energy_threshold_count", &n_missing_energy_threshold_);
  metadata_tree_->Branch("missing_acceptance_count", &n_missing_acceptance_);
  metadata_tree_->Branch("missing_other_count", &n_missing_other_);
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
  h_missing_energy_threshold_ = nullptr;
  h_missing_acceptance_ = nullptr;
  h_missing_other_ = nullptr;
  h_other_ = nullptr;
  metadata_tree_ = nullptr;
}
