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

  common_.set_min_cluster_energy(min_cluster_energy_);
  common_.set_shower_shape_min_tower_energy(shower_shape_min_tower_energy_);
  common_.set_store_shower_shape_tower_patch(store_shower_shape_tower_patch_);
  if (!common_.initialize())
  {
    std::cout << "PhotonAnalysisTree::Init - failed to initialize common reco filler" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

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
  photon_tree::EventVertex event_vertex;
  event_vertex.valid = true;
  event_vertex.x = b_vertex_x_;
  event_vertex.y = b_vertex_y_;
  event_vertex.z = b_vertex_z_;
  event_vertex.source = b_truth_valid_ ? 1 : 0;
  const bool split_valid =
      common_.fill_collection(split_clusters, towers, geometry, event_vertex, true, false, common_.split());
  const bool nosplit_valid =
      !nosplit_clusters ||
      common_.fill_collection(nosplit_clusters, towers, geometry, event_vertex, true, true, common_.nosplit());
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
              << " split/nosplit=" << common_.split().data.ncluster << "/"
              << common_.nosplit().data.ncluster << std::endl;
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
  b_vertex_x_ = b_truth_vx_;
  b_vertex_y_ = b_truth_vy_;
  b_vertex_z_ = b_truth_vz_;

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
  common_.clear_event();
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
  common_.create_collection_branches(event_tree_, "split", common_.split().data, true);
  common_.create_collection_branches(event_tree_, "nosplit", common_.nosplit().data, true);

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
