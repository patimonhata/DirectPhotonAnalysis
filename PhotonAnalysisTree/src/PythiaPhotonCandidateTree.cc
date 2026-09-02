#include "PythiaPhotonCandidateTree.h"

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfoContainer.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllInputManager.h>
#include <jetbase/Jet.h>
#include <jetbase/JetContainer.h>
#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <HepMC/GenEvent.h>
#include <HepMC/GenParticle.h>
#include <HepMC/WeightContainer.h>
#include <TMVA/RBDT.hxx>
#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
#include <limits>
#include <map>
#include <utility>

namespace
{
constexpr float invalid_float = -999.0F;

double delta_phi(double first, double second)
{
  return std::atan2(std::sin(first - second), std::cos(first - second));
}

bool raw_cluster_kinematics(const RawCluster* cluster, double vx, double vy, double vz,
                            double& energy, double& et, double& eta, double& phi)
{
  if (!cluster || !std::isfinite(cluster->get_energy()) || !std::isfinite(cluster->get_x()) ||
      !std::isfinite(cluster->get_y()) || !std::isfinite(cluster->get_z()))
  {
    return false;
  }
  const double dx = cluster->get_x() - vx;
  const double dy = cluster->get_y() - vy;
  const double dz = cluster->get_z() - vz;
  const double distance = std::sqrt(dx * dx + dy * dy + dz * dz);
  const double transverse = std::hypot(dx, dy);
  if (!(distance > 0.0) || !(transverse > 0.0))
  {
    return false;
  }
  energy = cluster->get_energy();
  et = energy * transverse / distance;
  eta = std::asinh(dz / transverse);
  phi = std::atan2(dy, dx);
  return std::isfinite(et) && std::isfinite(eta) && std::isfinite(phi);
}

double diphoton_mass(double e1, double eta1, double phi1, double e2, double eta2, double phi2)
{
  const double pt1 = e1 / std::cosh(eta1);
  const double pt2 = e2 / std::cosh(eta2);
  const double mass2 = 2.0 * pt1 * pt2 * (std::cosh(eta1 - eta2) - std::cos(delta_phi(phi1, phi2)));
  return mass2 >= 0.0 && std::isfinite(mass2) ? std::sqrt(mass2) : -1.0;
}

int cluster_id_from_index(const photon_tree::Pi0AnchorTopologyEventResult& result, std::size_t index)
{
  return index < result.clusters.size() ? static_cast<int>(result.clusters[index].cluster_id) : -999;
}
}

namespace photon_tree
{
void PhotonCandidateSelectionBranches::clear()
{
  *this = PhotonCandidateSelectionBranches{};
}

void PhotonCandidateSelectionBranches::create_branches(TTree* tree)
{
#define BRANCH(member) tree->Branch("split_cluster_" #member, &member)
  BRANCH(bdt_score);
  BRANCH(bdt_valid);
  BRANCH(bdt_tight_boundary);
  BRANCH(bdt_nontight_lower);
  BRANCH(bdt_nontight_upper);
  BRANCH(pass_kinematics);
  BRANCH(pass_preselection);
  BRANCH(pass_tight);
  BRANCH(pass_nontight);
  BRANCH(iso_raw_et);
  BRANCH(iso_corrected_et);
  BRANCH(iso_boundary);
  BRANCH(noniso_boundary);
  BRANCH(iso_topocluster_count);
  BRANCH(pass_isolated);
  BRANCH(pass_nonisolated);
  BRANCH(pass_region_a);
  BRANCH(pass_region_b);
  BRANCH(pass_region_c);
  BRANCH(pass_region_d);
  BRANCH(pass_final_photon);
  BRANCH(pi0_tag);
  BRANCH(eta_tag);
  BRANCH(pi0_partner_cluster_id);
  BRANCH(eta_partner_cluster_id);
  BRANCH(pi0_partner_mass);
  BRANCH(eta_partner_mass);
  BRANCH(truth_prompt_cluster);
  BRANCH(pi0_anchor_valid);
  BRANCH(pi0_anchor_candidate_index);
  BRANCH(pi0_anchor_main_fraction);
  BRANCH(pi0_anchor_second_fraction);
  BRANCH(pi0_anchor_unmatched_max_fraction);
  BRANCH(pi0_anchor_ambiguous_main);
  BRANCH(pi0_anchor_topology);
  BRANCH(pi0_anchor_reason);
  BRANCH(pi0_anchor_missing_category);
  BRANCH(pi0_anchor_missing_detail);
  BRANCH(pi0_anchor_partner_photon_index);
  BRANCH(pi0_anchor_partner_diagnostic_mass);
  BRANCH(pi0_anchor_partner_alignment);
  BRANCH(pi0_anchor_truth_partner_tag_status);
  BRANCH(pi0_anchor_tag_result);
  BRANCH(pi0_anchor_truth_partner_cluster_id);
  BRANCH(pi0_anchor_truth_partner_cluster_e);
  BRANCH(pi0_anchor_truth_partner_cluster_eta);
  BRANCH(pi0_anchor_truth_partner_cluster_phi);
  BRANCH(pi0_anchor_truth_partner_delta_r);
  BRANCH(pi0_anchor_truth_partner_direct_edep);
  BRANCH(pi0_anchor_truth_partner_reconstructed_e);
  BRANCH(pi0_anchor_truth_partner_recovery);
  BRANCH(pi0_anchor_truth_partner_mass);
  BRANCH(pi0_anchor_selected_tag_partner_matches_truth_partner);
#undef BRANCH
}

void Pi0TopologyTreeBranches::clear()
{
  *this = Pi0TopologyTreeBranches{};
}

void Pi0TopologyTreeBranches::fill(const Pi0AnchorTopologyEventResult& result)
{
  clear();
  candidate_count = static_cast<unsigned int>(result.candidates.size());
  for (const auto& candidate : result.candidates)
  {
    candidate_pathway.push_back(static_cast<int>(candidate.pathway));
    candidate_parent_barcode.push_back(candidate.parent_barcode);
    candidate_g4_parent_track_id.push_back(candidate.g4_parent_track_id);
    candidate_e.push_back(candidate.energy);
    candidate_pt.push_back(candidate.pt);
    candidate_eta.push_back(candidate.eta);
    candidate_phi.push_back(candidate.phi);
    candidate_evaluated.push_back(candidate.topology_evaluated ? 1U : 0U);

    photon0_track_id.push_back(candidate.photon_track_ids[0]);
    photon1_track_id.push_back(candidate.photon_track_ids[1]);
    photon0_e.push_back(candidate.photon_energy[0]);
    photon1_e.push_back(candidate.photon_energy[1]);
    photon0_eta.push_back(candidate.photon_eta[0]);
    photon1_eta.push_back(candidate.photon_eta[1]);
    photon0_phi.push_back(candidate.photon_phi[0]);
    photon1_phi.push_back(candidate.photon_phi[1]);
    photon0_projection_valid.push_back(candidate.photon_projection_valid[0] ? 1U : 0U);
    photon1_projection_valid.push_back(candidate.photon_projection_valid[1] ? 1U : 0U);
    photon0_projection_eta.push_back(candidate.photon_projection_eta[0]);
    photon1_projection_eta.push_back(candidate.photon_projection_eta[1]);
    photon0_projection_phi.push_back(candidate.photon_projection_phi[0]);
    photon1_projection_phi.push_back(candidate.photon_projection_phi[1]);
    photon0_in_cemc_acceptance.push_back(candidate.photon_in_cemc_acceptance[0] ? 1U : 0U);
    photon1_in_cemc_acceptance.push_back(candidate.photon_in_cemc_acceptance[1] ? 1U : 0U);
    photon0_pre_cemc_interaction.push_back(candidate.photon_pre_cemc_interaction[0] ? 1U : 0U);
    photon1_pre_cemc_interaction.push_back(candidate.photon_pre_cemc_interaction[1] ? 1U : 0U);
    photon0_first_daughter_radius.push_back(candidate.photon_first_daughter_radius[0]);
    photon1_first_daughter_radius.push_back(candidate.photon_first_daughter_radius[1]);
    photon0_cemc_edep.push_back(candidate.photon_cemc_edep[0]);
    photon1_cemc_edep.push_back(candidate.photon_cemc_edep[1]);
    photon0_best_cluster_id.push_back(cluster_id_from_index(result, candidate.best_cluster[0]));
    photon1_best_cluster_id.push_back(cluster_id_from_index(result, candidate.best_cluster[1]));
    photon0_maximum_edep.push_back(candidate.maximum_edep[0]);
    photon1_maximum_edep.push_back(candidate.maximum_edep[1]);
    photon0_reconstructed_e.push_back(candidate.reconstructed_photon_energy[0]);
    photon1_reconstructed_e.push_back(candidate.reconstructed_photon_energy[1]);
    photon0_recovered.push_back(candidate.recovered[0] ? 1U : 0U);
    photon1_recovered.push_back(candidate.recovered[1] ? 1U : 0U);

    const auto& diagnostic0 = candidate.partner_diagnostics[0];
    const auto& diagnostic1 = candidate.partner_diagnostics[1];
    photon0_diagnostic_found.push_back(diagnostic0.found ? 1U : 0U);
    photon1_diagnostic_found.push_back(diagnostic1.found ? 1U : 0U);
    photon0_diagnostic_cluster_id.push_back(diagnostic0.found ? static_cast<int>(diagnostic0.cluster_id) : -999);
    photon1_diagnostic_cluster_id.push_back(diagnostic1.found ? static_cast<int>(diagnostic1.cluster_id) : -999);
    photon0_diagnostic_cluster_e.push_back(diagnostic0.found ? diagnostic0.cluster_energy : invalid_float);
    photon1_diagnostic_cluster_e.push_back(diagnostic1.found ? diagnostic1.cluster_energy : invalid_float);
    photon0_diagnostic_delta_r.push_back(diagnostic0.found ? diagnostic0.delta_r : invalid_float);
    photon1_diagnostic_delta_r.push_back(diagnostic1.found ? diagnostic1.delta_r : invalid_float);
    photon0_diagnostic_recovery.push_back(diagnostic0.found ? diagnostic0.recovery : invalid_float);
    photon1_diagnostic_recovery.push_back(diagnostic1.found ? diagnostic1.recovery : invalid_float);
    photon0_diagnostic_below_threshold.push_back(diagnostic0.below_energy_threshold ? 1U : 0U);
    photon1_diagnostic_below_threshold.push_back(diagnostic1.below_energy_threshold ? 1U : 0U);
    photon0_diagnostic_direct_deposit.push_back(diagnostic0.has_direct_deposit ? 1U : 0U);
    photon1_diagnostic_direct_deposit.push_back(diagnostic1.has_direct_deposit ? 1U : 0U);
  }

  anchor_count = static_cast<unsigned int>(result.anchors.size());
  for (const auto& anchor : result.anchors)
  {
    anchor_cluster_id.push_back(anchor.cluster_index < result.clusters.size() ? result.clusters[anchor.cluster_index].cluster_id : 0U);
    anchor_candidate_index.push_back(static_cast<unsigned int>(anchor.candidate_index));
    anchor_main_fraction.push_back(anchor.main_fraction);
    anchor_second_fraction.push_back(anchor.second_fraction);
    anchor_unmatched_max_fraction.push_back(anchor.unmatched_max_fraction);
    anchor_ambiguous_main.push_back(anchor.ambiguous_main ? 1U : 0U);
    anchor_topology.push_back(static_cast<int>(anchor.topology));
    anchor_reason.push_back(static_cast<int>(anchor.reason));
    anchor_missing_category.push_back(static_cast<int>(anchor.missing_category));
    anchor_missing_detail.push_back(static_cast<int>(anchor.missing_detail));
    anchor_partner_photon_index.push_back(anchor.partner_photon_index);
    anchor_pre_cemc_photon_index.push_back(anchor.pre_cemc_photon_index);
    anchor_partner_diagnostic_mass.push_back(anchor.partner_diagnostic_invariant_mass);
    anchor_partner_alignment.push_back(static_cast<int>(anchor.partner_alignment));
    anchor_truth_partner_tag_status.push_back(static_cast<int>(anchor.truth_partner_tag_status));
    anchor_truth_partner_cluster_id.push_back(anchor.truth_partner_cluster_id);
    anchor_truth_partner_cluster_e.push_back(anchor.truth_partner_cluster_energy);
    anchor_truth_partner_cluster_eta.push_back(anchor.truth_partner_cluster_eta);
    anchor_truth_partner_cluster_phi.push_back(anchor.truth_partner_cluster_phi);
    anchor_truth_partner_delta_r.push_back(anchor.truth_partner_delta_r);
    anchor_truth_partner_direct_edep.push_back(anchor.truth_partner_direct_edep);
    anchor_truth_partner_reconstructed_e.push_back(anchor.truth_partner_reconstructed_photon_energy);
    anchor_truth_partner_recovery.push_back(anchor.truth_partner_recovery);
    anchor_truth_partner_mass.push_back(anchor.truth_partner_invariant_mass);
  }
}

void Pi0TopologyTreeBranches::create_branches(TTree* tree)
{
#define BRANCH(member) tree->Branch("pi0_topology_" #member, &member)
  BRANCH(candidate_count);
  BRANCH(candidate_pathway);
  BRANCH(candidate_parent_barcode);
  BRANCH(candidate_g4_parent_track_id);
  BRANCH(candidate_e);
  BRANCH(candidate_pt);
  BRANCH(candidate_eta);
  BRANCH(candidate_phi);
  BRANCH(candidate_evaluated);
  BRANCH(photon0_track_id);
  BRANCH(photon1_track_id);
  BRANCH(photon0_e);
  BRANCH(photon1_e);
  BRANCH(photon0_eta);
  BRANCH(photon1_eta);
  BRANCH(photon0_phi);
  BRANCH(photon1_phi);
  BRANCH(photon0_projection_valid);
  BRANCH(photon1_projection_valid);
  BRANCH(photon0_projection_eta);
  BRANCH(photon1_projection_eta);
  BRANCH(photon0_projection_phi);
  BRANCH(photon1_projection_phi);
  BRANCH(photon0_in_cemc_acceptance);
  BRANCH(photon1_in_cemc_acceptance);
  BRANCH(photon0_pre_cemc_interaction);
  BRANCH(photon1_pre_cemc_interaction);
  BRANCH(photon0_first_daughter_radius);
  BRANCH(photon1_first_daughter_radius);
  BRANCH(photon0_cemc_edep);
  BRANCH(photon1_cemc_edep);
  BRANCH(photon0_best_cluster_id);
  BRANCH(photon1_best_cluster_id);
  BRANCH(photon0_maximum_edep);
  BRANCH(photon1_maximum_edep);
  BRANCH(photon0_reconstructed_e);
  BRANCH(photon1_reconstructed_e);
  BRANCH(photon0_recovered);
  BRANCH(photon1_recovered);
  BRANCH(photon0_diagnostic_found);
  BRANCH(photon1_diagnostic_found);
  BRANCH(photon0_diagnostic_cluster_id);
  BRANCH(photon1_diagnostic_cluster_id);
  BRANCH(photon0_diagnostic_cluster_e);
  BRANCH(photon1_diagnostic_cluster_e);
  BRANCH(photon0_diagnostic_delta_r);
  BRANCH(photon1_diagnostic_delta_r);
  BRANCH(photon0_diagnostic_recovery);
  BRANCH(photon1_diagnostic_recovery);
  BRANCH(photon0_diagnostic_below_threshold);
  BRANCH(photon1_diagnostic_below_threshold);
  BRANCH(photon0_diagnostic_direct_deposit);
  BRANCH(photon1_diagnostic_direct_deposit);
  BRANCH(anchor_count);
  BRANCH(anchor_cluster_id);
  BRANCH(anchor_candidate_index);
  BRANCH(anchor_main_fraction);
  BRANCH(anchor_second_fraction);
  BRANCH(anchor_unmatched_max_fraction);
  BRANCH(anchor_ambiguous_main);
  BRANCH(anchor_topology);
  BRANCH(anchor_reason);
  BRANCH(anchor_missing_category);
  BRANCH(anchor_missing_detail);
  BRANCH(anchor_partner_photon_index);
  BRANCH(anchor_pre_cemc_photon_index);
  BRANCH(anchor_partner_diagnostic_mass);
  BRANCH(anchor_partner_alignment);
  BRANCH(anchor_truth_partner_tag_status);
  BRANCH(anchor_truth_partner_cluster_id);
  BRANCH(anchor_truth_partner_cluster_e);
  BRANCH(anchor_truth_partner_cluster_eta);
  BRANCH(anchor_truth_partner_cluster_phi);
  BRANCH(anchor_truth_partner_delta_r);
  BRANCH(anchor_truth_partner_direct_edep);
  BRANCH(anchor_truth_partner_reconstructed_e);
  BRANCH(anchor_truth_partner_recovery);
  BRANCH(anchor_truth_partner_mass);
#undef BRANCH
}
}

PythiaPhotonCandidateTree::PythiaPhotonCandidateTree(const std::string& name)
  : SubsysReco(name)
{
}

PythiaPhotonCandidateTree::~PythiaPhotonCandidateTree()
{
  close_output();
}

bool PythiaPhotonCandidateTree::configure_sample()
{
  struct Sample
  {
    const char* name;
    double minimum;
    double maximum;
    double cross_section;
    bool photonjet;
    bool unbounded;
  };
  static const Sample samples[] = {
      {"photonjet3", 0.0, 5.0, 1.0348e6, true, false},
      {"photonjet5", 5.0, 14.0, 1.4636e5, true, false},
      {"photonjet10", 14.0, 22.0, 6.9447e3, true, false},
      {"photonjet20", 22.0, 200.0, 1.3045e2, true, false},
      {"jet3", 0.0, 5.0, 1.2147e9, false, false},
      {"jet5", 5.0, 9.0, 1.3878e8, false, false},
      {"jet8", 9.0, 14.0, 1.3013e7, false, false},
      {"jet12", 14.0, 21.0, 1.4903e6, false, false},
      {"jet20", 21.0, 32.0, 6.2623e4, false, false},
      {"jet30", 32.0, 42.0, 2.5298e3, false, false},
      {"jet40", 42.0, 0.0, 1.3553e2, false, true}};
  std::transform(sample_name_.begin(), sample_name_.end(), sample_name_.begin(),
                 [](unsigned char value) { return static_cast<char>(std::tolower(value)); });
  for (const Sample& sample : samples)
  {
    if (sample_name_ != sample.name)
    {
      continue;
    }
    sample_window_min_ = sample.minimum;
    sample_window_max_ = sample.maximum;
    sample_cross_section_pb_ = sample.cross_section;
    sample_is_photonjet_ = sample.photonjet;
    sample_upper_unbounded_ = sample.unbounded;
    return true;
  }
  return false;
}

bool PythiaPhotonCandidateTree::update_input_provenance()
{
  const std::string input_file = primary_input_manager_ ? primary_input_manager_->FileName() : std::string{};
  if (primary_input_manager_ && input_file.empty())
  {
    std::cerr << "PythiaPhotonCandidateTree::process_event - primary input filename is unavailable" << std::endl;
    return false;
  }
  if (!input_file.empty() && input_file != current_input_file_)
  {
    const std::size_t basename_begin = input_file.find_last_of("/\\");
    const std::string basename = input_file.substr(basename_begin == std::string::npos ? 0U : basename_begin + 1U);
    const std::size_t extension_position = basename.rfind(".root");
    const std::size_t segment_separator = extension_position == std::string::npos ? std::string::npos : basename.rfind("-", extension_position);
    if (extension_position == std::string::npos || extension_position + 5U != basename.size() ||
        segment_separator == std::string::npos || segment_separator + 1U >= extension_position)
    {
      std::cerr << "PythiaPhotonCandidateTree::process_event - cannot parse input segment: " << input_file << std::endl;
      return false;
    }
    const std::string segment_text = basename.substr(segment_separator + 1U, extension_position - segment_separator - 1U);
    std::size_t parsed = 0U;
    unsigned long segment_id = 0UL;
    try
    {
      segment_id = std::stoul(segment_text, &parsed);
    }
    catch (const std::exception&)
    {
      parsed = 0U;
    }
    if (parsed != segment_text.size() || segment_id > std::numeric_limits<unsigned int>::max())
    {
      std::cerr << "PythiaPhotonCandidateTree::process_event - invalid input segment: " << input_file << std::endl;
      return false;
    }
    current_input_file_ = input_file;
    source_file_id_ = static_cast<unsigned int>(segment_id);
    next_event_in_file_ = 0U;
  }
  if (next_event_in_file_ == std::numeric_limits<unsigned int>::max())
  {
    std::cerr << "PythiaPhotonCandidateTree::process_event - event index overflow" << std::endl;
    return false;
  }
  b_source_file_id_ = source_file_id_;
  b_event_in_file_ = next_event_in_file_++;
  b_event_uid_ = (static_cast<unsigned long long>(b_source_file_id_) << 32U) | b_event_in_file_;
  return true;
}

int PythiaPhotonCandidateTree::Init(PHCompositeNode* /*topNode*/)
{
  const bool manifest_valid = manifest_path_.empty() || (manifest_begin_ >= 0 && manifest_end_ > manifest_begin_ && primary_input_manager_);
  if (output_file_name_.empty() || model_file_name_.empty() || signal_embedding_id_ <= 0 || !manifest_valid ||
      !std::isfinite(min_cluster_energy_) || min_cluster_energy_ < 0.0 || !std::isfinite(meson_partner_min_energy_) ||
      meson_partner_min_energy_ < 0.0 || !configure_sample())
  {
    std::cerr << "PythiaPhotonCandidateTree::Init - invalid output/model/sample configuration" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  input_file_count_ = manifest_path_.empty() ? 1 : manifest_end_ - manifest_begin_;
  common_.set_min_cluster_energy(min_cluster_energy_);
  common_.set_min_cluster_energy_inclusive(false);
  common_.set_shower_shape_min_tower_energy(shower_shape_min_tower_energy_);
  common_.set_store_shower_shape_tower_patch(false);
  common_.set_store_cluster_pairs(false);
  if (!common_.initialize())
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }

  photon_tree::Pi0AnchorTopologyConfig config;
  config.sample_mode = photon_tree::Pi0SampleMode::pythia;
  config.truth_node_name = truth_node_name_;
  config.hepmc_event_map_node_name = hepmc_event_map_node_name_;
  config.tower_node_name = tower_node_name_;
  config.raw_truth_tower_node_name = raw_truth_tower_node_name_;
  config.truth_cell_node_name = truth_cell_node_name_;
  config.truth_hit_node_name = truth_hit_node_name_;
  config.cluster_node_name = split_cluster_node_name_;
  config.tower_geom_node_name = tower_geom_node_name_;
  config.signal_embedding_id = signal_embedding_id_;
  config.anchor_cluster_eta_max = candidate_abs_eta_max_;
  config.partner_cluster_eta_max = -1.0;
  config.min_cluster_energy = std::nextafter(min_cluster_energy_, std::numeric_limits<double>::infinity());
  config.partner_diagnostic_min_cluster_energy = partner_diagnostic_min_cluster_energy_;
  config.tagging_partner_min_cluster_energy = meson_partner_min_energy_;
  config.tagging_pi0_mass_min = pi0_mass_min_;
  config.tagging_pi0_mass_max = pi0_mass_max_;
  config.max_abs_vertex_z = max_abs_vertex_z_;
  config.enable_missing_diagnostics = true;
  config.verbosity = verbosity_;
  topology_evaluator_.configure(config);

  try
  {
    bdt_ = std::make_unique<TMVA::Experimental::RBDT>(model_key_, model_file_name_);
  }
  catch (const std::exception& error)
  {
    std::cerr << "PythiaPhotonCandidateTree::Init - cannot load BDT: " << error.what() << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  create_output();
  if (!output_file_ || output_file_->IsZombie())
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }
  std::cout << "PythiaPhotonCandidateTree::Init - sample/model/output = "
            << sample_name_ << "/" << model_file_name_ << "/" << output_file_name_ << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

int PythiaPhotonCandidateTree::process_event(PHCompositeNode* topNode)
{
  reset_event();
  ++n_events_processed_;
  if (!update_input_provenance())
  {
    ++n_events_invalid_;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  auto* event_map = findNode::getClass<PHHepMCGenEventMap>(topNode, hepmc_event_map_node_name_);
  if (!event_map || !fill_event_truth(event_map, topNode))
  {
    ++n_events_invalid_;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  sum_generator_weight_processed_ += b_generator_weight_;
  if (b_sample_stitching_pass_)
  {
    ++n_events_stitch_pass_;
    sum_generator_weight_stitch_pass_ += b_generator_weight_;
  }

  const photon_tree::Pi0AnchorTopologyEventResult topology_result = topology_evaluator_.evaluate(topNode);
  if (topology_result.status == photon_tree::Pi0TopologyEventStatus::vertex_rejected)
  {
    ++n_events_vertex_rejected_;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  if (topology_result.status != photon_tree::Pi0TopologyEventStatus::accepted)
  {
    ++n_events_invalid_;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  auto* split_clusters = findNode::getClass<RawClusterContainer>(topNode, split_cluster_node_name_);
  auto* topo_clusters = findNode::getClass<RawClusterContainer>(topNode, topo_cluster_node_name_);
  auto* towers = findNode::getClass<TowerInfoContainer>(topNode, tower_node_name_);
  auto* geometry = findNode::getClass<RawTowerGeomContainer>(topNode, tower_geom_node_name_);
  if (!split_clusters || !topo_clusters || !towers || !geometry)
  {
    ++n_events_invalid_;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  photon_tree::EventVertex vertex;
  vertex.valid = true;
  vertex.x = b_vertex_x_;
  vertex.y = b_vertex_y_;
  vertex.z = b_vertex_z_;
  vertex.source = signal_embedding_id_;
  if (!common_.fill_collection(split_clusters, towers, geometry, vertex, false, false, common_.split()))
  {
    ++n_events_invalid_;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  topology_.fill(topology_result);
  fill_topology_cluster_links(topology_result);
  if (!fill_candidate_selection(split_clusters, topo_clusters))
  {
    ++n_events_invalid_;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  const auto topo_range = topo_clusters->getClusters();
  for (auto iterator = topo_range.first; iterator != topo_range.second; ++iterator)
  {
    b_topocluster_count_ += iterator->second ? 1U : 0U;
  }
  event_tree_->Fill();
  ++n_events_written_;
  return Fun4AllReturnCodes::EVENT_OK;
}

bool PythiaPhotonCandidateTree::fill_event_truth(const PHHepMCGenEventMap* event_map, PHCompositeNode* topNode)
{
  const PHHepMCGenEvent* signal_event = event_map ? event_map->get(signal_embedding_id_) : nullptr;
  const HepMC::GenEvent* event = signal_event ? signal_event->getEvent() : nullptr;
  if (!signal_event || !event || !signal_event->is_simulated())
  {
    return false;
  }

  const auto& vertex = signal_event->get_collision_vertex();
  if (!std::isfinite(vertex.x()) || !std::isfinite(vertex.y()) || !std::isfinite(vertex.z()) ||
      !std::isfinite(vertex.t()))
  {
    return false;
  }
  b_vertex_x_ = vertex.x();
  b_vertex_y_ = vertex.y();
  b_vertex_z_ = vertex.z();
  b_vertex_t_ = vertex.t();
  b_hepmc_event_number_ = event->event_number();

  const auto weights = event->weights().weights();
  b_event_weight_valid_ = std::all_of(weights.begin(), weights.end(), [](double value) { return std::isfinite(value); }) ? 1U : 0U;
  b_generator_weight_ = weights.empty() ? 1.0 : weights.front();
  if (!b_event_weight_valid_ || !std::isfinite(b_generator_weight_))
  {
    return false;
  }
  b_weight_numerator_pb_ = sample_cross_section_pb_ * b_generator_weight_;

  bool found_prompt = false;
  double leading_prompt_pt = 0.0;
  for (auto iterator = event->particles_begin(); iterator != event->particles_end(); ++iterator)
  {
    const HepMC::GenParticle* particle = *iterator;
    if (!particle || particle->pdg_id() != 22 || particle->status() != 1 || particle->end_vertex())
    {
      continue;
    }
    const auto classification = photon_classifier_.classify(particle);
    if (!classification.valid || (classification.category != 1 && classification.category != 2))
    {
      continue;
    }
    const double pt = std::hypot(particle->momentum().px(), particle->momentum().py());
    if (std::isfinite(pt))
    {
      leading_prompt_pt = std::max(leading_prompt_pt, pt);
      found_prompt = true;
    }
  }
  b_leading_truth_photon_pt_ = found_prompt ? leading_prompt_pt : photon_tree::kInvalidDouble;

  bool stitch_valid = false;
  double stitch_value = photon_tree::kInvalidDouble;
  if (sample_is_photonjet_)
  {
    stitch_valid = found_prompt;
    stitch_value = b_leading_truth_photon_pt_;
  }
  else
  {
    auto* jets = findNode::getClass<JetContainer>(topNode, truth_jet_node_name_);
    if (jets)
    {
      double leading_jet_pt = -1.0;
      for (Jet* jet : *jets)
      {
        if (jet && std::isfinite(jet->get_pt()))
        {
          leading_jet_pt = std::max(leading_jet_pt, static_cast<double>(jet->get_pt()));
        }
      }
      if (leading_jet_pt >= 0.0)
      {
        b_leading_truth_jet_pt_ = leading_jet_pt;
        stitch_value = leading_jet_pt;
        stitch_valid = true;
      }
    }
  }
  b_sample_stitching_valid_ = stitch_valid ? 1U : 0U;
  b_sample_stitching_pass_ = stitch_valid && stitch_value >= sample_window_min_ &&
      (sample_upper_unbounded_ || stitch_value < sample_window_max_) ? 1U : 0U;
  return true;
}

void PythiaPhotonCandidateTree::fill_topology_cluster_links(const photon_tree::Pi0AnchorTopologyEventResult& result)
{
  const auto& ids = common_.split().data.cluster_id;
  std::map<unsigned int, std::size_t> result_index_by_id;
  for (std::size_t index = 0; index < result.clusters.size(); ++index)
  {
    result_index_by_id[result.clusters[index].cluster_id] = index;
  }
  std::map<unsigned int, const photon_tree::Pi0TopologyAnchorRecord*> anchor_by_id;
  for (const auto& anchor : result.anchors)
  {
    if (anchor.cluster_index < result.clusters.size())
    {
      anchor_by_id[result.clusters[anchor.cluster_index].cluster_id] = &anchor;
    }
  }

  cluster_truth_.clear();
  for (unsigned int id : ids)
  {
    const auto result_found = result_index_by_id.find(id);
    if (result_found == result_index_by_id.end())
    {
      cluster_truth_.append({});
      selection_.truth_prompt_cluster.push_back(0U);
    }
    else
    {
      const std::size_t index = result_found->second;
      cluster_truth_.append(result.clusters[index].truth);
      selection_.truth_prompt_cluster.push_back(index < result.prompt_cluster.size() ? result.prompt_cluster[index] : 0U);
    }

    const auto anchor_found = anchor_by_id.find(id);
    const auto* anchor = anchor_found == anchor_by_id.end() ? nullptr : anchor_found->second;
    selection_.pi0_anchor_valid.push_back(anchor ? 1U : 0U);
    selection_.pi0_anchor_candidate_index.push_back(anchor ? static_cast<int>(anchor->candidate_index) : -999);
    selection_.pi0_anchor_main_fraction.push_back(anchor ? anchor->main_fraction : invalid_float);
    selection_.pi0_anchor_second_fraction.push_back(anchor ? anchor->second_fraction : invalid_float);
    selection_.pi0_anchor_unmatched_max_fraction.push_back(anchor ? anchor->unmatched_max_fraction : invalid_float);
    selection_.pi0_anchor_ambiguous_main.push_back(anchor && anchor->ambiguous_main ? 1U : 0U);
    selection_.pi0_anchor_topology.push_back(anchor ? static_cast<int>(anchor->topology) : -999);
    selection_.pi0_anchor_reason.push_back(anchor ? static_cast<int>(anchor->reason) : -999);
    selection_.pi0_anchor_missing_category.push_back(anchor ? static_cast<int>(anchor->missing_category) : -999);
    selection_.pi0_anchor_missing_detail.push_back(anchor ? static_cast<int>(anchor->missing_detail) : -999);
    selection_.pi0_anchor_partner_photon_index.push_back(anchor ? anchor->partner_photon_index : -999);
    selection_.pi0_anchor_partner_diagnostic_mass.push_back(anchor ? anchor->partner_diagnostic_invariant_mass : invalid_float);
    selection_.pi0_anchor_partner_alignment.push_back(anchor ? static_cast<int>(anchor->partner_alignment) : -999);
    selection_.pi0_anchor_truth_partner_tag_status.push_back(anchor ? static_cast<int>(anchor->truth_partner_tag_status) : -999);
    selection_.pi0_anchor_tag_result.push_back(static_cast<int>(photon_tree::Pi0AnchorTagResult::not_applicable));
    selection_.pi0_anchor_truth_partner_cluster_id.push_back(anchor ? anchor->truth_partner_cluster_id : -999);
    selection_.pi0_anchor_truth_partner_cluster_e.push_back(anchor ? anchor->truth_partner_cluster_energy : invalid_float);
    selection_.pi0_anchor_truth_partner_cluster_eta.push_back(anchor ? anchor->truth_partner_cluster_eta : invalid_float);
    selection_.pi0_anchor_truth_partner_cluster_phi.push_back(anchor ? anchor->truth_partner_cluster_phi : invalid_float);
    selection_.pi0_anchor_truth_partner_delta_r.push_back(anchor ? anchor->truth_partner_delta_r : invalid_float);
    selection_.pi0_anchor_truth_partner_direct_edep.push_back(anchor ? anchor->truth_partner_direct_edep : invalid_float);
    selection_.pi0_anchor_truth_partner_reconstructed_e.push_back(anchor ? anchor->truth_partner_reconstructed_photon_energy : invalid_float);
    selection_.pi0_anchor_truth_partner_recovery.push_back(anchor ? anchor->truth_partner_recovery : invalid_float);
    selection_.pi0_anchor_truth_partner_mass.push_back(anchor ? anchor->truth_partner_invariant_mass : invalid_float);
    selection_.pi0_anchor_selected_tag_partner_matches_truth_partner.push_back(0U);
  }
}

bool PythiaPhotonCandidateTree::fill_candidate_selection(RawClusterContainer* split_clusters, RawClusterContainer* topo_clusters)
{
  struct KinematicCluster
  {
    const RawCluster* cluster = nullptr;
    double energy = 0.0;
    double et = 0.0;
    double eta = 0.0;
    double phi = 0.0;
  };
  std::vector<KinematicCluster> all_split;
  const auto split_range = split_clusters->getClusters();
  for (auto iterator = split_range.first; iterator != split_range.second; ++iterator)
  {
    KinematicCluster value;
    value.cluster = iterator->second;
    if (raw_cluster_kinematics(value.cluster, b_vertex_x_, b_vertex_y_, b_vertex_z_, value.energy, value.et, value.eta, value.phi))
    {
      all_split.push_back(value);
    }
  }
  std::vector<KinematicCluster> all_topo;
  const auto topo_range = topo_clusters->getClusters();
  for (auto iterator = topo_range.first; iterator != topo_range.second; ++iterator)
  {
    KinematicCluster value;
    value.cluster = iterator->second;
    if (raw_cluster_kinematics(value.cluster, b_vertex_x_, b_vertex_y_, b_vertex_z_, value.energy, value.et, value.eta, value.phi) &&
        value.energy > 0.0 && value.et > 0.0)
    {
      all_topo.push_back(value);
    }
  }

  const auto& data = common_.split().data;
  for (std::size_t index = 0; index < data.ncluster; ++index)
  {
    const double et = data.cluster_et[index];
    const double eta = data.cluster_eta[index];
    const double phi = data.cluster_phi[index];
    std::vector<float> features = {
        static_cast<float>(et), data.shower_w_eta_cogx[index], data.shower_w_phi_cogx[index],
        static_cast<float>(b_vertex_z_), static_cast<float>(eta), data.shower_e11_over_e33[index],
        data.shower_et1[index], data.shower_et2[index], data.shower_et3[index], data.shower_et4[index],
        data.shower_e32_over_e35[index]};
    const bool finite_features = std::all_of(features.begin(), features.end(), [](float value) { return std::isfinite(value); });
    float score = invalid_float;
    bool bdt_valid = false;
    if (finite_features && data.shower_valid[index])
    {
      const auto result = bdt_->Compute(features);
      bdt_valid = !result.empty() && std::isfinite(result.front());
      score = bdt_valid ? result.front() : invalid_float;
    }

    const float tight_boundary = static_cast<float>(0.8156 - 0.00156 * et);
    const float nontight_lower = static_cast<float>(0.7333 - 0.01333 * et);
    const float nontight_upper = static_cast<float>(0.6844 + 0.00156 * et);
    const bool kinematics = et > candidate_et_min_ && et < candidate_et_max_ && std::abs(eta) < candidate_abs_eta_max_;
    const bool preselection = data.shower_valid[index] && data.shower_e11_over_e33[index] < 0.98F &&
        data.shower_et1[index] > 0.6F && data.shower_et1[index] <= 1.0F &&
        data.shower_e32_over_e35[index] > 0.8F && data.shower_e32_over_e35[index] <= 1.0F;
    const bool tight = bdt_valid && score > tight_boundary;
    const bool nontight = bdt_valid && score > nontight_lower && score < nontight_upper;

    double iso_sum = 0.0;
    unsigned int iso_count = 0U;
    for (const auto& topo : all_topo)
    {
      if (std::hypot(topo.eta - eta, delta_phi(topo.phi, phi)) < isolation_radius_)
      {
        iso_sum += topo.et;
        ++iso_count;
      }
    }
    const double iso_raw = iso_sum - et;
    const double iso_corrected = isolation_scale_ * iso_raw + isolation_offset_;
    const double iso_boundary = 0.490 + 0.037 * et;
    const double noniso_boundary = iso_boundary + nonisolation_gap_;
    const bool isolated = iso_corrected < iso_boundary;
    const bool nonisolated = iso_corrected > noniso_boundary;

    bool pi0_tag = false;
    bool eta_tag = false;
    int pi0_partner = -999;
    int eta_partner = -999;
    double pi0_mass = invalid_float;
    double eta_mass = invalid_float;
    double pi0_distance = std::numeric_limits<double>::infinity();
    double eta_distance = std::numeric_limits<double>::infinity();
    for (const auto& partner : all_split)
    {
      if (partner.cluster->get_id() == data.cluster_id[index] || !(partner.energy > meson_partner_min_energy_))
      {
        continue;
      }
      const double mass = diphoton_mass(data.cluster_e[index], eta, phi, partner.energy, partner.eta, partner.phi);
      if (mass > pi0_mass_min_ && mass < pi0_mass_max_ && std::abs(mass - 0.134977) < pi0_distance)
      {
        pi0_tag = true;
        pi0_distance = std::abs(mass - 0.134977);
        pi0_mass = mass;
        pi0_partner = static_cast<int>(partner.cluster->get_id());
      }
      if (mass > eta_mass_min_ && mass < eta_mass_max_ && std::abs(mass - 0.547862) < eta_distance)
      {
        eta_tag = true;
        eta_distance = std::abs(mass - 0.547862);
        eta_mass = mass;
        eta_partner = static_cast<int>(partner.cluster->get_id());
      }
    }

    if (selection_.pi0_anchor_valid[index])
    {
      const bool truth_partner_taggable = selection_.pi0_anchor_truth_partner_tag_status[index] ==
          static_cast<int>(photon_tree::Pi0TruthPartnerTagStatus::taggable);
      const photon_tree::Pi0AnchorTagResult tag_result = !pi0_tag
          ? photon_tree::Pi0AnchorTagResult::survived
          : (truth_partner_taggable ? photon_tree::Pi0AnchorTagResult::truth_pair_taggable_veto
                                    : photon_tree::Pi0AnchorTagResult::combinatorial_only_veto);
      selection_.pi0_anchor_tag_result[index] = static_cast<int>(tag_result);
      selection_.pi0_anchor_selected_tag_partner_matches_truth_partner[index] = pi0_tag &&
          pi0_partner == selection_.pi0_anchor_truth_partner_cluster_id[index] ? 1U : 0U;
    }

    const bool common_selection = kinematics && preselection;
    const bool region_a = common_selection && isolated && tight;
    const bool region_b = common_selection && nonisolated && tight;
    const bool region_c = common_selection && isolated && nontight;
    const bool region_d = common_selection && nonisolated && nontight;
    const bool final_photon = region_a && !pi0_tag && !eta_tag;
    selection_.bdt_score.push_back(score);
    selection_.bdt_valid.push_back(bdt_valid ? 1U : 0U);
    selection_.bdt_tight_boundary.push_back(tight_boundary);
    selection_.bdt_nontight_lower.push_back(nontight_lower);
    selection_.bdt_nontight_upper.push_back(nontight_upper);
    selection_.pass_kinematics.push_back(kinematics ? 1U : 0U);
    selection_.pass_preselection.push_back(preselection ? 1U : 0U);
    selection_.pass_tight.push_back(tight ? 1U : 0U);
    selection_.pass_nontight.push_back(nontight ? 1U : 0U);
    selection_.iso_raw_et.push_back(iso_raw);
    selection_.iso_corrected_et.push_back(iso_corrected);
    selection_.iso_boundary.push_back(iso_boundary);
    selection_.noniso_boundary.push_back(noniso_boundary);
    selection_.iso_topocluster_count.push_back(iso_count);
    selection_.pass_isolated.push_back(isolated ? 1U : 0U);
    selection_.pass_nonisolated.push_back(nonisolated ? 1U : 0U);
    selection_.pass_region_a.push_back(region_a ? 1U : 0U);
    selection_.pass_region_b.push_back(region_b ? 1U : 0U);
    selection_.pass_region_c.push_back(region_c ? 1U : 0U);
    selection_.pass_region_d.push_back(region_d ? 1U : 0U);
    selection_.pass_final_photon.push_back(final_photon ? 1U : 0U);
    b_region_a_count_ += region_a ? 1U : 0U;
    b_region_b_count_ += region_b ? 1U : 0U;
    b_region_c_count_ += region_c ? 1U : 0U;
    b_region_d_count_ += region_d ? 1U : 0U;
    b_final_photon_count_ += final_photon ? 1U : 0U;
    selection_.pi0_tag.push_back(pi0_tag ? 1U : 0U);
    selection_.eta_tag.push_back(eta_tag ? 1U : 0U);
    selection_.pi0_partner_cluster_id.push_back(pi0_partner);
    selection_.eta_partner_cluster_id.push_back(eta_partner);
    selection_.pi0_partner_mass.push_back(pi0_mass);
    selection_.eta_partner_mass.push_back(eta_mass);
  }
  const bool valid = selection_.bdt_score.size() == data.ncluster;
  if (valid)
  {
    n_clusters_region_a_ += b_region_a_count_;
    n_clusters_region_b_ += b_region_b_count_;
    n_clusters_region_c_ += b_region_c_count_;
    n_clusters_region_d_ += b_region_d_count_;
    n_clusters_final_photon_ += b_final_photon_count_;
  }
  return valid;
}

void PythiaPhotonCandidateTree::reset_event()
{
  b_source_file_id_ = 0U;
  b_event_in_file_ = 0U;
  b_event_uid_ = 0ULL;
  b_hepmc_event_number_ = -999;
  b_vertex_x_ = b_vertex_y_ = b_vertex_z_ = b_vertex_t_ = photon_tree::kInvalidDouble;
  b_event_weight_valid_ = 0U;
  b_generator_weight_ = 1.0;
  b_weight_numerator_pb_ = 0.0;
  b_leading_truth_photon_pt_ = photon_tree::kInvalidDouble;
  b_leading_truth_jet_pt_ = photon_tree::kInvalidDouble;
  b_sample_stitching_valid_ = 0U;
  b_sample_stitching_pass_ = 0U;
  b_topocluster_count_ = 0U;
  b_region_a_count_ = b_region_b_count_ = b_region_c_count_ = b_region_d_count_ = b_final_photon_count_ = 0U;
  common_.clear_event();
  cluster_truth_.clear();
  selection_.clear();
  topology_.clear();
}

void PythiaPhotonCandidateTree::create_output()
{
  const auto slash = output_file_name_.find_last_of('/');
  if (slash != std::string::npos && slash > 0U)
  {
    gSystem->mkdir(output_file_name_.substr(0, slash).c_str(), true);
  }
  output_file_ = TFile::Open(output_file_name_.c_str(), "RECREATE");
  if (!output_file_ || output_file_->IsZombie())
  {
    return;
  }

  output_file_->cd();
  event_tree_ = new TTree("event_tree", "Pythia split-cluster photon-candidate selection tree");
#define EVENT_BRANCH(member) event_tree_->Branch(#member, &b_##member##_)
  EVENT_BRANCH(source_file_id);
  EVENT_BRANCH(event_in_file);
  EVENT_BRANCH(event_uid);
  EVENT_BRANCH(hepmc_event_number);
  EVENT_BRANCH(vertex_x);
  EVENT_BRANCH(vertex_y);
  EVENT_BRANCH(vertex_z);
  EVENT_BRANCH(vertex_t);
  EVENT_BRANCH(event_weight_valid);
  EVENT_BRANCH(generator_weight);
  EVENT_BRANCH(weight_numerator_pb);
  EVENT_BRANCH(leading_truth_photon_pt);
  EVENT_BRANCH(leading_truth_jet_pt);
  EVENT_BRANCH(sample_stitching_valid);
  EVENT_BRANCH(sample_stitching_pass);
  EVENT_BRANCH(topocluster_count);
  EVENT_BRANCH(region_a_count);
  EVENT_BRANCH(region_b_count);
  EVENT_BRANCH(region_c_count);
  EVENT_BRANCH(region_d_count);
  EVENT_BRANCH(final_photon_count);
#undef EVENT_BRANCH
  common_.create_collection_branches(event_tree_, "split", common_.split().data, false);
  cluster_truth_.create_branches(event_tree_, "split");
  selection_.create_branches(event_tree_);
  topology_.create_branches(event_tree_);

  metadata_tree_ = new TTree("metadata", "Pythia photon-candidate workflow metadata");
  metadata_tree_->Branch("schema_version", &metadata_schema_version_);
  metadata_tree_->Branch("analysis_release", &analysis_release_);
  metadata_tree_->Branch("input_file", &input_file_name_);
  metadata_tree_->Branch("input_manifest", &manifest_path_);
  metadata_tree_->Branch("manifest_begin", &manifest_begin_);
  metadata_tree_->Branch("manifest_end", &manifest_end_);
  metadata_tree_->Branch("input_file_count", &input_file_count_);
  metadata_tree_->Branch("first_input_suffix", &first_input_suffix_);
  metadata_tree_->Branch("last_input_suffix", &last_input_suffix_);
  metadata_tree_->Branch("map_chunk_id", &map_chunk_id_);
  metadata_tree_->Branch("output_file", &output_file_name_);
  metadata_tree_->Branch("sample_name", &sample_name_);
  metadata_tree_->Branch("sample_cross_section_pb", &sample_cross_section_pb_);
  metadata_tree_->Branch("sample_window_min", &sample_window_min_);
  metadata_tree_->Branch("sample_window_max", &sample_window_max_);
  metadata_tree_->Branch("sample_upper_unbounded", &sample_upper_unbounded_);
  metadata_tree_->Branch("event_weight_definition", &event_weight_definition_);
  metadata_tree_->Branch("photon_stitch_definition", &photon_stitch_definition_);
  metadata_tree_->Branch("jet_stitch_definition", &jet_stitch_definition_);
  metadata_tree_->Branch("model_file", &model_file_name_);
  metadata_tree_->Branch("model_key", &model_key_);
  metadata_tree_->Branch("model_sha256", &model_sha256_);
  metadata_tree_->Branch("split_cluster_node", &split_cluster_node_name_);
  metadata_tree_->Branch("topo_cluster_node", &topo_cluster_node_name_);
  metadata_tree_->Branch("truth_jet_node", &truth_jet_node_name_);
  metadata_tree_->Branch("topocluster_configuration", &topocluster_configuration_);
  metadata_tree_->Branch("isolation_definition", &isolation_definition_);
  metadata_tree_->Branch("min_cluster_energy", &min_cluster_energy_);
  metadata_tree_->Branch("partner_diagnostic_min_cluster_energy", &partner_diagnostic_min_cluster_energy_);
  metadata_tree_->Branch("pi0_topology_algorithm_version", &pi0_topology_algorithm_version_);
  metadata_tree_->Branch("shower_shape_min_tower_energy", &shower_shape_min_tower_energy_);
  metadata_tree_->Branch("candidate_et_min", &candidate_et_min_);
  metadata_tree_->Branch("candidate_et_max", &candidate_et_max_);
  metadata_tree_->Branch("candidate_abs_eta_max", &candidate_abs_eta_max_);
  metadata_tree_->Branch("max_abs_vertex_z", &max_abs_vertex_z_);
  metadata_tree_->Branch("isolation_radius", &isolation_radius_);
  metadata_tree_->Branch("isolation_scale", &isolation_scale_);
  metadata_tree_->Branch("isolation_offset", &isolation_offset_);
  metadata_tree_->Branch("nonisolation_gap", &nonisolation_gap_);
  metadata_tree_->Branch("meson_partner_min_energy", &meson_partner_min_energy_);
  metadata_tree_->Branch("pi0_mass_min", &pi0_mass_min_);
  metadata_tree_->Branch("pi0_mass_max", &pi0_mass_max_);
  metadata_tree_->Branch("eta_mass_min", &eta_mass_min_);
  metadata_tree_->Branch("eta_mass_max", &eta_mass_max_);
  metadata_tree_->Branch("n_events_processed", &n_events_processed_);
  metadata_tree_->Branch("n_events_written", &n_events_written_);
  metadata_tree_->Branch("n_events_vertex_rejected", &n_events_vertex_rejected_);
  metadata_tree_->Branch("n_events_invalid", &n_events_invalid_);
  metadata_tree_->Branch("n_events_stitch_pass", &n_events_stitch_pass_);
  metadata_tree_->Branch("n_clusters_region_a", &n_clusters_region_a_);
  metadata_tree_->Branch("n_clusters_region_b", &n_clusters_region_b_);
  metadata_tree_->Branch("n_clusters_region_c", &n_clusters_region_c_);
  metadata_tree_->Branch("n_clusters_region_d", &n_clusters_region_d_);
  metadata_tree_->Branch("n_clusters_final_photon", &n_clusters_final_photon_);
  metadata_tree_->Branch("sum_generator_weight_processed", &sum_generator_weight_processed_);
  metadata_tree_->Branch("sum_generator_weight_stitch_pass", &sum_generator_weight_stitch_pass_);
}

int PythiaPhotonCandidateTree::End(PHCompositeNode* /*topNode*/)
{
  close_output();
  std::cout << "PythiaPhotonCandidateTree::End - processed/written/vertex-rejected/invalid/stitch-pass = "
            << n_events_processed_ << "/" << n_events_written_ << "/" << n_events_vertex_rejected_ << "/"
            << n_events_invalid_ << "/" << n_events_stitch_pass_ << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

void PythiaPhotonCandidateTree::close_output()
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
  bdt_.reset();
}
