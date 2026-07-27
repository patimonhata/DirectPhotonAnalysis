#include "PhotonAnalysisTree.h"

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
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
#include <set>
#include <utility>

PhotonAnalysisTree::PhotonAnalysisTree(const std::string& name)
  : SubsysReco(name)
{
}

PhotonAnalysisTree::~PhotonAnalysisTree()
{
  close_output_file();
}

void PhotonAnalysisTree::ClusterCollection::clear()
{
  *this = ClusterCollection{};
}

int PhotonAnalysisTree::Init(PHCompositeNode* /*topNode*/)
{
  if (expected_primary_pdg_ != 22 && expected_primary_pdg_ != 111 && expected_primary_pdg_ != 221)
  {
    std::cout << "PhotonAnalysisTree::Init - expected primary PDG must be 22 or 111" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  if (!(acceptance_eta_max_ > 0.0) || min_cluster_energy_ < 0.0 || shower_shape_min_tower_energy_ < 0.0)
  {
    std::cout << "PhotonAnalysisTree::Init - invalid numeric configuration" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  ShowerShapeCalculator::Config shower_config;
  shower_config.min_tower_energy = static_cast<float>(shower_shape_min_tower_energy_);
  shower_shape_calculator_ = ShowerShapeCalculator(shower_config);

  create_output_directory();
  output_file_ = TFile::Open(output_file_name_.c_str(), "RECREATE");
  if (!output_file_ || output_file_->IsZombie())
  {
    std::cout << "PhotonAnalysisTree::Init - failed to create " << output_file_name_ << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  create_trees();

  std::cout << "PhotonAnalysisTree::Init - input: " << input_file_name_ << '\n'
            << "PhotonAnalysisTree::Init - output: " << output_file_name_ << '\n'
            << "PhotonAnalysisTree::Init - split/no-split nodes: "
            << split_cluster_node_name_ << "/" << nosplit_cluster_node_name_ << '\n'
            << "PhotonAnalysisTree::Init - truth/no-split required: "
            << require_truth_node_ << "/" << require_nosplit_cluster_node_ << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

int PhotonAnalysisTree::process_event(PHCompositeNode* topNode)
{
  reset_event();
  b_source_file_id_ = source_file_id_;
  b_event_in_file_ = static_cast<unsigned int>(n_events_processed_);
  b_event_uid_ = (static_cast<unsigned long long>(source_file_id_) << 32U) |
                 static_cast<unsigned long long>(b_event_in_file_);
  ++n_events_processed_;

  auto* truth = findNode::getClass<PHG4TruthInfoContainer>(topNode, truth_node_name_);
  auto* towers = findNode::getClass<TowerInfoContainer>(topNode, tower_node_name_);
  auto* geometry = findNode::getClass<RawTowerGeomContainer>(topNode, tower_geom_node_name_);
  auto* split_clusters = findNode::getClass<RawClusterContainer>(topNode, split_cluster_node_name_);
  auto* nosplit_clusters = findNode::getClass<RawClusterContainer>(topNode, nosplit_cluster_node_name_);

  const bool missing_required_node =
      !towers || !geometry || !split_clusters ||
      (require_truth_node_ && !truth) ||
      (require_nosplit_cluster_node_ && !nosplit_clusters);
  if (missing_required_node)
  {
    ++n_events_invalid_detector_;
    std::cout << "PhotonAnalysisTree::process_event - missing node in event " << b_event_in_file_
              << ": truth=" << static_cast<bool>(truth)
              << " towers=" << static_cast<bool>(towers)
              << " geometry=" << static_cast<bool>(geometry)
              << " split=" << static_cast<bool>(split_clusters)
              << " nosplit=" << static_cast<bool>(nosplit_clusters) << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if (truth && !fill_truth(truth, geometry))
  {
    ++n_events_invalid_truth_;
  }
  const bool split_valid =
      fill_collection(split_clusters, towers, geometry, true, false, split_);
  const bool nosplit_valid =
      !nosplit_clusters ||
      fill_collection(nosplit_clusters, towers, geometry, true, true, nosplit_);
  if (!split_valid || !nosplit_valid)
  {
    ++n_events_invalid_detector_;
    std::cout << "PhotonAnalysisTree::process_event - invalid cluster/tower content in event "
              << b_event_in_file_ << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  event_tree_->Fill();
  ++n_events_written_;
  if (verbosity_ > 0 && b_event_in_file_ < 5)
  {
    std::cout << "PhotonAnalysisTree::process_event - event " << b_event_in_file_
              << " truth_valid=" << static_cast<int>(b_truth_valid_)
              << " split/nosplit=" << split_.ncluster << "/" << nosplit_.ncluster << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int PhotonAnalysisTree::End(PHCompositeNode* /*topNode*/)
{
  close_output_file();
  std::cout << "PhotonAnalysisTree::End - processed/written/invalid_truth/invalid_detector = "
            << n_events_processed_ << "/" << n_events_written_ << "/"
            << n_events_invalid_truth_ << "/" << n_events_invalid_detector_ << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

bool PhotonAnalysisTree::fill_truth(PHG4TruthInfoContainer* truth, RawTowerGeomContainer* geometry)
{
  if (!truth)
  {
    return false;
  }

  PHG4Particle* primary = nullptr;
  unsigned int primary_count = 0;
  const auto primary_range = truth->GetPrimaryParticleRange();
  for (auto iter = primary_range.first; iter != primary_range.second; ++iter)
  {
    if (iter->second)
    {
      primary = iter->second;
      ++primary_count;
    }
  }
  if (primary_count != 1 || !primary || primary->get_pid() != expected_primary_pdg_)
  {
    return false;
  }

  PHG4VtxPoint* primary_vertex = truth->GetVtx(primary->get_vtx_id());
  const double pt = std::hypot(primary->get_px(), primary->get_py());
  if (!primary_vertex || pt <= 0.0 ||
      !std::isfinite(primary->get_e()) || !std::isfinite(primary->get_px()) ||
      !std::isfinite(primary->get_py()) || !std::isfinite(primary->get_pz()) ||
      !std::isfinite(primary_vertex->get_x()) || !std::isfinite(primary_vertex->get_y()) ||
      !std::isfinite(primary_vertex->get_z()))
  {
    return false;
  }

  b_label_ = primary->get_pid() == 22 ? 1 : 0;
  b_truth_primary_pdg_id_ = primary->get_pid();
  b_truth_primary_track_id_ = primary->get_track_id();
  b_truth_e_ = primary->get_e();
  b_truth_px_ = primary->get_px();
  b_truth_py_ = primary->get_py();
  b_truth_pz_ = primary->get_pz();
  b_truth_pt_ = pt;
  b_truth_eta_ = std::asinh(primary->get_pz() / pt);
  b_truth_phi_ = std::atan2(primary->get_py(), primary->get_px());
  b_truth_vx_ = primary_vertex->get_x();
  b_truth_vy_ = primary_vertex->get_y();
  b_truth_vz_ = primary_vertex->get_z();

  std::vector<PHG4Particle*> daughters;
  const auto particle_range = truth->GetParticleRange();
  for (auto iter = particle_range.first; iter != particle_range.second; ++iter)
  {
    PHG4Particle* particle = iter->second;
    if (particle && particle->get_parent_id() == primary->get_track_id())
    {
      daughters.push_back(particle);
    }
  }
  std::sort(daughters.begin(), daughters.end(), [](const PHG4Particle* lhs, const PHG4Particle* rhs) {
    return lhs->get_track_id() < rhs->get_track_id();
  });
  b_truth_n_direct_daughter_ = static_cast<unsigned int>(daughters.size());
  b_truth_is_pi0_to_2gamma_ = daughters.size() == 2 && daughters[0]->get_pid() == 22 &&
                             daughters[1]->get_pid() == 22 ? 1U : 0U;

  const double projection_radius = cemc_radius(geometry);
  for (PHG4Particle* daughter : daughters)
  {
    const double daughter_pt = std::hypot(daughter->get_px(), daughter->get_py());
    if (daughter_pt <= 0.0 || !std::isfinite(daughter->get_e()) ||
        !std::isfinite(daughter->get_px()) || !std::isfinite(daughter->get_py()) ||
        !std::isfinite(daughter->get_pz()))
    {
      return false;
    }

    b_truth_daughter_track_id_.push_back(daughter->get_track_id());
    b_truth_daughter_pdg_id_.push_back(daughter->get_pid());
    b_truth_daughter_e_.push_back(daughter->get_e());
    b_truth_daughter_px_.push_back(daughter->get_px());
    b_truth_daughter_py_.push_back(daughter->get_py());
    b_truth_daughter_pz_.push_back(daughter->get_pz());
    b_truth_daughter_pt_.push_back(daughter_pt);
    b_truth_daughter_eta_.push_back(std::asinh(daughter->get_pz() / daughter_pt));
    b_truth_daughter_phi_.push_back(std::atan2(daughter->get_py(), daughter->get_px()));

    double projection_eta = invalid_double_;
    double projection_phi = invalid_double_;
    bool projection_valid = false;
    PHG4VtxPoint* vertex = truth->GetVtx(daughter->get_vtx_id());
    if (vertex)
    {
      double x = 0.0;
      double y = 0.0;
      double z = 0.0;
      projection_valid = project_to_radius(vertex->get_x(), vertex->get_y(), vertex->get_z(),
                                           daughter->get_px(), daughter->get_py(), daughter->get_pz(),
                                           projection_radius, x, y, z);
      if (projection_valid)
      {
        projection_eta = eta_from_xyz(x, y, z);
        projection_phi = phi_from_xy(x, y);
        projection_valid = std::isfinite(projection_eta) && std::isfinite(projection_phi) &&
                           projection_eta != invalid_double_ && projection_phi != invalid_double_;
      }
    }
    b_truth_daughter_projection_eta_.push_back(projection_valid ? projection_eta : invalid_double_);
    b_truth_daughter_projection_phi_.push_back(projection_valid ? projection_phi : invalid_double_);
    b_truth_daughter_projection_valid_.push_back(projection_valid ? 1U : 0U);
    b_truth_daughter_in_acceptance_.push_back(
        projection_valid && std::abs(projection_eta) < acceptance_eta_max_ ? 1U : 0U);
  }

  if (b_truth_is_pi0_to_2gamma_)
  {
    b_truth_both_gamma_in_acceptance_ =
        b_truth_daughter_in_acceptance_[0] && b_truth_daughter_in_acceptance_[1] ? 1U : 0U;
    b_truth_at_least_one_gamma_out_acceptance_ = b_truth_both_gamma_in_acceptance_ ? 0U : 1U;
    b_truth_missing_gamma_projection_ =
        b_truth_daughter_projection_valid_[0] && b_truth_daughter_projection_valid_[1] ? 0U : 1U;
    b_truth_m_gg_ = invariant_mass(daughters[0], daughters[1]);
    b_truth_pair_e_asym_ = energy_asymmetry(daughters[0], daughters[1]);
  }

  b_truth_valid_ = 1U;
  return true;
}

bool PhotonAnalysisTree::fill_collection(RawClusterContainer* clusters,
                                         TowerInfoContainer* towers,
                                         RawTowerGeomContainer* geometry,
                                         bool include_towers,
                                         bool require_unique_tower_keys,
                                         ClusterCollection& output)
{
  std::vector<const RawCluster*> ordered_clusters;
  const auto range = clusters->getClusters();
  for (auto iter = range.first; iter != range.second; ++iter)
  {
    const RawCluster* cluster = iter->second;
    if (!cluster)
    {
      continue;
    }
    if (!std::isfinite(cluster->get_energy()) || !std::isfinite(cluster->get_x()) ||
        !std::isfinite(cluster->get_y()) || !std::isfinite(cluster->get_z()))
    {
      return false;
    }
    if (cluster->get_energy() >= min_cluster_energy_)
    {
      ordered_clusters.push_back(cluster);
    }
  }
  std::sort(ordered_clusters.begin(), ordered_clusters.end(), [](const RawCluster* lhs, const RawCluster* rhs) {
    if (lhs->get_energy() != rhs->get_energy())
    {
      return lhs->get_energy() > rhs->get_energy();
    }
    return lhs->get_id() < rhs->get_id();
  });

  std::set<unsigned int> seen_tower_keys;
  for (std::size_t cluster_index = 0; cluster_index < ordered_clusters.size(); ++cluster_index)
  {
    const RawCluster* cluster = ordered_clusters[cluster_index];
    const double dx = cluster->get_x() - b_vertex_x_;
    const double dy = cluster->get_y() - b_vertex_y_;
    const double dz = cluster->get_z() - b_vertex_z_;
    const double distance = std::sqrt(dx * dx + dy * dy + dz * dz);
    if (!(distance > std::numeric_limits<double>::epsilon()))
    {
      return false;
    }

    const double energy = cluster->get_energy();
    const double px = energy * dx / distance;
    const double py = energy * dy / distance;
    const double pz = energy * dz / distance;
    const double pt = std::hypot(px, py);
    if (!(pt > 0.0))
    {
      return false;
    }

    output.cluster_id.push_back(cluster->get_id());
    output.cluster_ntower.push_back(static_cast<int>(cluster->getNTowers()));
    output.cluster_e.push_back(energy);
    output.cluster_et.push_back(pt);
    output.cluster_eta.push_back(std::asinh(pz / pt));
    output.cluster_phi.push_back(std::atan2(py, px));
    output.cluster_x.push_back(cluster->get_x());
    output.cluster_y.push_back(cluster->get_y());
    output.cluster_z.push_back(cluster->get_z());
    output.cluster_px.push_back(px);
    output.cluster_py.push_back(py);
    output.cluster_pz.push_back(pz);
    append_shower_shape(shower_shape_calculator_.calculate(*cluster, *towers), output);

    if (!include_towers)
    {
      continue;
    }
    const auto tower_range = cluster->get_towers();
    for (auto tower_iter = tower_range.first; tower_iter != tower_range.second; ++tower_iter)
    {
      const unsigned int raw_key = tower_iter->first;
      if (require_unique_tower_keys && !seen_tower_keys.insert(raw_key).second)
      {
        std::cout << "PhotonAnalysisTree::fill_collection - duplicate tower key "
                  << raw_key << std::endl;
        return false;
      }
      const int ieta = static_cast<int>(RawTowerDefs::decode_index1(raw_key));
      const int iphi = static_cast<int>(RawTowerDefs::decode_index2(raw_key));
      const unsigned int tower_info_key = TowerInfoDefs::encode_emcal(
          static_cast<unsigned int>(ieta), static_cast<unsigned int>(iphi));
      TowerInfo* tower = towers->get_tower_at_key(static_cast<int>(tower_info_key));
      RawTowerGeom* tower_geometry = geometry->get_tower_geometry(raw_key);
      if (!tower || !tower_geometry || !std::isfinite(tower->get_energy()))
      {
        return false;
      }
      const double x = tower_geometry->get_center_x();
      const double y = tower_geometry->get_center_y();
      const double z = tower_geometry->get_center_z();
      const double tower_radius = radius(x, y);
      const double eta = eta_from_xyz(x, y, z);
      const double phi = phi_from_xy(x, y);
      if (tower_radius == invalid_double_ || eta == invalid_double_ || phi == invalid_double_)
      {
        return false;
      }

      output.tower_cluster_index.push_back(static_cast<int>(cluster_index));
      output.tower_key.push_back(raw_key);
      output.tower_ieta.push_back(ieta);
      output.tower_iphi.push_back(iphi);
      output.tower_x.push_back(x);
      output.tower_y.push_back(y);
      output.tower_z.push_back(z);
      output.tower_r.push_back(tower_radius);
      output.tower_eta.push_back(eta);
      output.tower_phi.push_back(phi);
      output.tower_energy.push_back(tower->get_energy());
      output.tower_cluster_value.push_back(
          std::isfinite(tower_iter->second) ? tower_iter->second : invalid_double_);
      output.tower_time.push_back(std::isfinite(tower->get_time()) ? tower->get_time() : invalid_double_);
      output.tower_is_good.push_back(tower->get_isGood() ? 1 : 0);
      output.tower_status.push_back(static_cast<int>(tower->get_status()));
    }
  }

  output.ncluster = static_cast<unsigned int>(output.cluster_id.size());
  output.ntower = static_cast<unsigned int>(output.tower_key.size());
  for (std::size_t i = 0; i < output.ncluster; ++i)
  {
    for (std::size_t j = i + 1; j < output.ncluster; ++j)
    {
      const double total_e = output.cluster_e[i] + output.cluster_e[j];
      const double px = output.cluster_px[i] + output.cluster_px[j];
      const double py = output.cluster_py[i] + output.cluster_py[j];
      const double pz = output.cluster_pz[i] + output.cluster_pz[j];
      const double mass2 = total_e * total_e - px * px - py * py - pz * pz;
      output.pair_cluster_i.push_back(static_cast<unsigned int>(i));
      output.pair_cluster_j.push_back(static_cast<unsigned int>(j));
      output.pair_m_gg.push_back(std::sqrt(std::max(0.0, mass2)));
      output.pair_e_asym.push_back(total_e > 0.0
                                      ? std::abs(output.cluster_e[i] - output.cluster_e[j]) / total_e
                                      : invalid_double_);
    }
  }
  return true;
}

void PhotonAnalysisTree::append_shower_shape(const ShowerShapeCalculator::Result& result,
                                             ClusterCollection& output)
{
  output.shower_valid.push_back(result.valid ? 1U : 0U);
  output.shower_full_containment.push_back(result.full_containment ? 1U : 0U);
  output.shower_edge_padded.push_back(result.edge_padded ? 1U : 0U);
  output.shower_tower_data_complete.push_back(result.tower_data_complete ? 1U : 0U);
  output.shower_cog_ieta.push_back(result.cog_ieta);
  output.shower_cog_iphi.push_back(result.cog_iphi);
  output.shower_cluster_e_thresholded.push_back(result.cluster_energy_above_threshold);
  output.shower_owned_patch_e.push_back(result.owned_patch_energy);
  output.shower_w_eta_cogx.push_back(result.w_eta_cogx);
  output.shower_w_phi_cogx.push_back(result.w_phi_cogx);
  output.shower_e11.push_back(result.e11);
  output.shower_e33.push_back(result.e33);
  output.shower_e32.push_back(result.e32);
  output.shower_e35.push_back(result.e35);
  output.shower_e11_over_e33.push_back(result.e11_over_e33);
  output.shower_e32_over_e35.push_back(result.e32_over_e35);
  output.shower_et1.push_back(result.et1);
  output.shower_et2.push_back(result.et2);
  output.shower_et3.push_back(result.et3);
  output.shower_et4.push_back(result.et4);
  if (store_shower_shape_tower_patch_)
  {
    output.shower_patch_e.insert(output.shower_patch_e.end(), result.patch_energy.begin(), result.patch_energy.end());
    output.shower_patch_good.insert(output.shower_patch_good.end(), result.patch_good.begin(), result.patch_good.end());
    output.shower_patch_owned.insert(output.shower_patch_owned.end(), result.patch_owned.begin(), result.patch_owned.end());
  }
}

void PhotonAnalysisTree::reset_event()
{
  b_source_file_id_ = 0;
  b_event_in_file_ = 0;
  b_event_uid_ = 0;
  b_vertex_x_ = 0.0;
  b_vertex_y_ = 0.0;
  b_vertex_z_ = 0.0;
  b_label_ = invalid_int_;
  b_truth_valid_ = 0;
  b_truth_primary_pdg_id_ = invalid_int_;
  b_truth_primary_track_id_ = invalid_int_;
  b_truth_e_ = b_truth_px_ = b_truth_py_ = b_truth_pz_ = invalid_double_;
  b_truth_pt_ = b_truth_eta_ = b_truth_phi_ = invalid_double_;
  b_truth_vx_ = b_truth_vy_ = b_truth_vz_ = invalid_double_;
  b_truth_n_direct_daughter_ = 0;
  b_truth_is_pi0_to_2gamma_ = 0;
  b_truth_daughter_track_id_.clear();
  b_truth_daughter_pdg_id_.clear();
  b_truth_daughter_e_.clear();
  b_truth_daughter_px_.clear();
  b_truth_daughter_py_.clear();
  b_truth_daughter_pz_.clear();
  b_truth_daughter_pt_.clear();
  b_truth_daughter_eta_.clear();
  b_truth_daughter_phi_.clear();
  b_truth_daughter_projection_eta_.clear();
  b_truth_daughter_projection_phi_.clear();
  b_truth_daughter_projection_valid_.clear();
  b_truth_daughter_in_acceptance_.clear();
  b_truth_both_gamma_in_acceptance_ = 0;
  b_truth_at_least_one_gamma_out_acceptance_ = 0;
  b_truth_missing_gamma_projection_ = 0;
  b_truth_m_gg_ = invalid_double_;
  b_truth_pair_e_asym_ = invalid_double_;
  split_.clear();
  nosplit_.clear();
}

void PhotonAnalysisTree::create_trees()
{
  output_file_->cd();
  event_tree_ = new TTree("event_tree", "Unified truth, split, and no-split photon analysis tree");
  event_tree_->Branch("source_file_id", &b_source_file_id_);
  event_tree_->Branch("event_in_file", &b_event_in_file_);
  event_tree_->Branch("event_uid", &b_event_uid_);
  event_tree_->Branch("vertex_x", &b_vertex_x_);
  event_tree_->Branch("vertex_y", &b_vertex_y_);
  event_tree_->Branch("vertex_z", &b_vertex_z_);
  event_tree_->Branch("label", &b_label_);
  event_tree_->Branch("truth_valid", &b_truth_valid_);
  event_tree_->Branch("truth_primary_pdg_id", &b_truth_primary_pdg_id_);
  event_tree_->Branch("truth_primary_track_id", &b_truth_primary_track_id_);
  event_tree_->Branch("truth_e", &b_truth_e_);
  event_tree_->Branch("truth_px", &b_truth_px_);
  event_tree_->Branch("truth_py", &b_truth_py_);
  event_tree_->Branch("truth_pz", &b_truth_pz_);
  event_tree_->Branch("truth_pt", &b_truth_pt_);
  event_tree_->Branch("truth_eta", &b_truth_eta_);
  event_tree_->Branch("truth_phi", &b_truth_phi_);
  event_tree_->Branch("truth_vx", &b_truth_vx_);
  event_tree_->Branch("truth_vy", &b_truth_vy_);
  event_tree_->Branch("truth_vz", &b_truth_vz_);
  event_tree_->Branch("truth_n_direct_daughter", &b_truth_n_direct_daughter_);
  event_tree_->Branch("truth_is_pi0_to_2gamma", &b_truth_is_pi0_to_2gamma_);
  event_tree_->Branch("truth_daughter_track_id", &b_truth_daughter_track_id_);
  event_tree_->Branch("truth_daughter_pdg_id", &b_truth_daughter_pdg_id_);
  event_tree_->Branch("truth_daughter_e", &b_truth_daughter_e_);
  event_tree_->Branch("truth_daughter_px", &b_truth_daughter_px_);
  event_tree_->Branch("truth_daughter_py", &b_truth_daughter_py_);
  event_tree_->Branch("truth_daughter_pz", &b_truth_daughter_pz_);
  event_tree_->Branch("truth_daughter_pt", &b_truth_daughter_pt_);
  event_tree_->Branch("truth_daughter_eta", &b_truth_daughter_eta_);
  event_tree_->Branch("truth_daughter_phi", &b_truth_daughter_phi_);
  event_tree_->Branch("truth_daughter_projection_eta", &b_truth_daughter_projection_eta_);
  event_tree_->Branch("truth_daughter_projection_phi", &b_truth_daughter_projection_phi_);
  event_tree_->Branch("truth_daughter_projection_valid", &b_truth_daughter_projection_valid_);
  event_tree_->Branch("truth_daughter_in_acceptance", &b_truth_daughter_in_acceptance_);
  event_tree_->Branch("truth_both_gamma_in_acceptance", &b_truth_both_gamma_in_acceptance_);
  event_tree_->Branch("truth_at_least_one_gamma_out_acceptance", &b_truth_at_least_one_gamma_out_acceptance_);
  event_tree_->Branch("truth_missing_gamma_projection", &b_truth_missing_gamma_projection_);
  event_tree_->Branch("truth_m_gg", &b_truth_m_gg_);
  event_tree_->Branch("truth_pair_e_asym", &b_truth_pair_e_asym_);

  event_tree_->Branch("min_cluster_energy", &min_cluster_energy_);
  event_tree_->Branch("shower_shape_min_tower_energy", &shower_shape_min_tower_energy_);
  event_tree_->Branch("shower_shape_algorithm_version", &shower_shape_algorithm_version_);
  event_tree_->Branch("shower_shape_patch_side", &shower_shape_patch_side_);
  event_tree_->Branch("store_shower_shape_tower_patch", &store_shower_shape_tower_patch_);
  create_collection_branches("split", split_, true);
  create_collection_branches("nosplit", nosplit_, true);

  metadata_tree_ = new TTree("metadata", "One entry per source DST");
  metadata_tree_->Branch("schema_version", &metadata_schema_version_);
  metadata_tree_->Branch("input_file", &input_file_name_);
  metadata_tree_->Branch("output_file", &output_file_name_);
  metadata_tree_->Branch("source_file_id", &source_file_id_);
  metadata_tree_->Branch("expected_primary_pdg", &expected_primary_pdg_);
  metadata_tree_->Branch("truth_node", &truth_node_name_);
  metadata_tree_->Branch("tower_node", &tower_node_name_);
  metadata_tree_->Branch("tower_geom_node", &tower_geom_node_name_);
  metadata_tree_->Branch("split_cluster_node", &split_cluster_node_name_);
  metadata_tree_->Branch("nosplit_cluster_node", &nosplit_cluster_node_name_);
  metadata_tree_->Branch("require_truth_node", &require_truth_node_);
  metadata_tree_->Branch("require_nosplit_cluster_node", &require_nosplit_cluster_node_);
  metadata_tree_->Branch("cluster_ordering", &cluster_ordering_);
  metadata_tree_->Branch("acceptance_eta_max", &acceptance_eta_max_);
  metadata_tree_->Branch("n_events_processed", &n_events_processed_);
  metadata_tree_->Branch("n_events_written", &n_events_written_);
  metadata_tree_->Branch("n_events_invalid_truth", &n_events_invalid_truth_);
  metadata_tree_->Branch("n_events_invalid_detector", &n_events_invalid_detector_);
}

void PhotonAnalysisTree::create_collection_branches(const std::string& prefix,
                                                    ClusterCollection& c,
                                                    bool include_towers)
{
  const auto name = [&prefix](const char* suffix) { return prefix + "_" + suffix; };
  event_tree_->Branch(name("ncluster").c_str(), &c.ncluster);
  event_tree_->Branch(name("cluster_id").c_str(), &c.cluster_id);
  event_tree_->Branch(name("cluster_ntower").c_str(), &c.cluster_ntower);
  event_tree_->Branch(name("cluster_e").c_str(), &c.cluster_e);
  event_tree_->Branch(name("cluster_et").c_str(), &c.cluster_et);
  event_tree_->Branch(name("cluster_eta").c_str(), &c.cluster_eta);
  event_tree_->Branch(name("cluster_phi").c_str(), &c.cluster_phi);
  event_tree_->Branch(name("cluster_x").c_str(), &c.cluster_x);
  event_tree_->Branch(name("cluster_y").c_str(), &c.cluster_y);
  event_tree_->Branch(name("cluster_z").c_str(), &c.cluster_z);
  event_tree_->Branch(name("cluster_px").c_str(), &c.cluster_px);
  event_tree_->Branch(name("cluster_py").c_str(), &c.cluster_py);
  event_tree_->Branch(name("cluster_pz").c_str(), &c.cluster_pz);
  event_tree_->Branch(name("cluster_shower_valid").c_str(), &c.shower_valid);
  event_tree_->Branch(name("cluster_shower_full_containment").c_str(), &c.shower_full_containment);
  event_tree_->Branch(name("cluster_shower_edge_padded").c_str(), &c.shower_edge_padded);
  event_tree_->Branch(name("cluster_shower_tower_data_complete").c_str(), &c.shower_tower_data_complete);
  event_tree_->Branch(name("cluster_shower_cog_ieta").c_str(), &c.shower_cog_ieta);
  event_tree_->Branch(name("cluster_shower_cog_iphi").c_str(), &c.shower_cog_iphi);
  event_tree_->Branch(name("cluster_shower_cluster_e_thresholded").c_str(), &c.shower_cluster_e_thresholded);
  event_tree_->Branch(name("cluster_shower_owned_patch_e").c_str(), &c.shower_owned_patch_e);
  event_tree_->Branch(name("cluster_shower_w_eta_cogx").c_str(), &c.shower_w_eta_cogx);
  event_tree_->Branch(name("cluster_shower_w_phi_cogx").c_str(), &c.shower_w_phi_cogx);
  event_tree_->Branch(name("cluster_shower_e11").c_str(), &c.shower_e11);
  event_tree_->Branch(name("cluster_shower_e33").c_str(), &c.shower_e33);
  event_tree_->Branch(name("cluster_shower_e32").c_str(), &c.shower_e32);
  event_tree_->Branch(name("cluster_shower_e35").c_str(), &c.shower_e35);
  event_tree_->Branch(name("cluster_shower_e11_over_e33").c_str(), &c.shower_e11_over_e33);
  event_tree_->Branch(name("cluster_shower_e32_over_e35").c_str(), &c.shower_e32_over_e35);
  event_tree_->Branch(name("cluster_shower_et1").c_str(), &c.shower_et1);
  event_tree_->Branch(name("cluster_shower_et2").c_str(), &c.shower_et2);
  event_tree_->Branch(name("cluster_shower_et3").c_str(), &c.shower_et3);
  event_tree_->Branch(name("cluster_shower_et4").c_str(), &c.shower_et4);
  event_tree_->Branch(name("cluster_shower_patch_e").c_str(), &c.shower_patch_e);
  event_tree_->Branch(name("cluster_shower_patch_good").c_str(), &c.shower_patch_good);
  event_tree_->Branch(name("cluster_shower_patch_owned").c_str(), &c.shower_patch_owned);
  event_tree_->Branch(name("pair_cluster_i").c_str(), &c.pair_cluster_i);
  event_tree_->Branch(name("pair_cluster_j").c_str(), &c.pair_cluster_j);
  event_tree_->Branch(name("pair_m_gg").c_str(), &c.pair_m_gg);
  event_tree_->Branch(name("pair_e_asym").c_str(), &c.pair_e_asym);

  if (include_towers)
  {
    event_tree_->Branch(name("ntower").c_str(), &c.ntower);
    event_tree_->Branch(name("tower_cluster_index").c_str(), &c.tower_cluster_index);
    event_tree_->Branch(name("tower_key").c_str(), &c.tower_key);
    event_tree_->Branch(name("tower_ieta").c_str(), &c.tower_ieta);
    event_tree_->Branch(name("tower_iphi").c_str(), &c.tower_iphi);
    event_tree_->Branch(name("tower_x").c_str(), &c.tower_x);
    event_tree_->Branch(name("tower_y").c_str(), &c.tower_y);
    event_tree_->Branch(name("tower_z").c_str(), &c.tower_z);
    event_tree_->Branch(name("tower_r").c_str(), &c.tower_r);
    event_tree_->Branch(name("tower_eta").c_str(), &c.tower_eta);
    event_tree_->Branch(name("tower_phi").c_str(), &c.tower_phi);
    event_tree_->Branch(name("tower_energy").c_str(), &c.tower_energy);
    event_tree_->Branch(name("tower_cluster_value").c_str(), &c.tower_cluster_value);
    event_tree_->Branch(name("tower_time").c_str(), &c.tower_time);
    event_tree_->Branch(name("tower_is_good").c_str(), &c.tower_is_good);
    event_tree_->Branch(name("tower_status").c_str(), &c.tower_status);
  }
}

void PhotonAnalysisTree::create_output_directory() const
{
  const auto slash = output_file_name_.find_last_of('/');
  if (slash != std::string::npos)
  {
    const std::string directory = output_file_name_.substr(0, slash);
    if (!directory.empty())
    {
      gSystem->mkdir(directory.c_str(), true);
    }
  }
}

void PhotonAnalysisTree::close_output_file()
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

double PhotonAnalysisTree::radius(double x, double y)
{
  return std::isfinite(x) && std::isfinite(y) ? std::hypot(x, y) : invalid_double_;
}

double PhotonAnalysisTree::eta_from_xyz(double x, double y, double z)
{
  const double r = radius(x, y);
  return r > 0.0 && std::isfinite(z) ? std::asinh(z / r) : invalid_double_;
}

double PhotonAnalysisTree::phi_from_xy(double x, double y)
{
  return std::isfinite(x) && std::isfinite(y) ? std::atan2(y, x) : invalid_double_;
}

double PhotonAnalysisTree::wrap_delta_phi(double value)
{
  constexpr double pi = 3.14159265358979323846;
  while (value > pi) value -= 2.0 * pi;
  while (value <= -pi) value += 2.0 * pi;
  return value;
}

bool PhotonAnalysisTree::project_to_radius(double x0, double y0, double z0,
                                           double px, double py, double pz,
                                           double target_radius,
                                           double& x1, double& y1, double& z1)
{
  const double momentum = std::sqrt(px * px + py * py + pz * pz);
  if (!(momentum > 0.0) || !(target_radius > 0.0))
  {
    return false;
  }
  const double ux = px / momentum;
  const double uy = py / momentum;
  const double uz = pz / momentum;
  const double a = ux * ux + uy * uy;
  const double b = 2.0 * (x0 * ux + y0 * uy);
  const double c = x0 * x0 + y0 * y0 - target_radius * target_radius;
  const double discriminant = b * b - 4.0 * a * c;
  if (a <= std::numeric_limits<double>::epsilon() || discriminant < 0.0)
  {
    return false;
  }
  const double root = std::sqrt(discriminant);
  const double s1 = (-b + root) / (2.0 * a);
  const double s2 = (-b - root) / (2.0 * a);
  double distance = std::numeric_limits<double>::max();
  if (s1 > 0.0) distance = std::min(distance, s1);
  if (s2 > 0.0) distance = std::min(distance, s2);
  if (distance == std::numeric_limits<double>::max())
  {
    return false;
  }
  x1 = x0 + distance * ux;
  y1 = y0 + distance * uy;
  z1 = z0 + distance * uz;
  return std::isfinite(x1) && std::isfinite(y1) && std::isfinite(z1);
}

double PhotonAnalysisTree::invariant_mass(const PHG4Particle* first, const PHG4Particle* second)
{
  const double energy = first->get_e() + second->get_e();
  const double px = first->get_px() + second->get_px();
  const double py = first->get_py() + second->get_py();
  const double pz = first->get_pz() + second->get_pz();
  return std::sqrt(std::max(0.0, energy * energy - px * px - py * py - pz * pz));
}

double PhotonAnalysisTree::energy_asymmetry(const PHG4Particle* first, const PHG4Particle* second)
{
  const double total = first->get_e() + second->get_e();
  return total > 0.0 ? std::abs(first->get_e() - second->get_e()) / total : invalid_double_;
}

double PhotonAnalysisTree::cemc_radius(RawTowerGeomContainer* geometry) const
{
  if (!geometry)
  {
    return default_cemc_radius_;
  }
  for (int ieta = 0; ieta < 96; ++ieta)
  {
    for (int iphi = 0; iphi < 256; ++iphi)
    {
      const unsigned int key = RawTowerDefs::encode_towerid(
          RawTowerDefs::CEMC, static_cast<unsigned int>(ieta), static_cast<unsigned int>(iphi));
      RawTowerGeom* tower = geometry->get_tower_geometry(key);
      if (tower)
      {
        const double value = radius(tower->get_center_x(), tower->get_center_y());
        return value > 0.0 ? value : default_cemc_radius_;
      }
    }
  }
  return default_cemc_radius_;
}
