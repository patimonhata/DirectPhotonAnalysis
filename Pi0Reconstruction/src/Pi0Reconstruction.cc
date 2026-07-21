#include "Pi0Reconstruction.h"

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/TowerInfoContainer.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <globalvertex/GlobalVertexv3.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/MbdVertexv3.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

Pi0Reconstruction::Pi0Reconstruction(const std::string &name)
  : SubsysReco(name)
{
}

Pi0Reconstruction::~Pi0Reconstruction()
{
  if (output_file_)
  {
    output_file_->Close();
    delete output_file_;
    output_file_ = nullptr;
  }
}

void Pi0Reconstruction::set_output_file_name(const std::string &output_file_name)
{
  output_file_name_ = output_file_name;
}

void Pi0Reconstruction::set_cluster_node_name(const std::string &cluster_node_name)
{
  cluster_node_name_ = cluster_node_name;
}

void Pi0Reconstruction::set_process_id(unsigned int process_id)
{
  process_id_ = process_id;
}

void Pi0Reconstruction::set_vertex_node_name(const std::string &vertex_node_name)
{
  vertex_node_name_ = vertex_node_name;
}

void Pi0Reconstruction::set_vertex_mode(VertexMode vertex_mode)
{
  vertex_mode_ = vertex_mode;
}

void Pi0Reconstruction::set_abort_on_missing_cluster_node(bool abort_on_missing_cluster_node)
{
  abort_on_missing_cluster_node_ = abort_on_missing_cluster_node;
}

void Pi0Reconstruction::set_abort_on_missing_vertex_node(bool abort_on_missing_vertex_node)
{
  abort_on_missing_vertex_node_ = abort_on_missing_vertex_node;
}

void Pi0Reconstruction::set_cemc_tower_node_name(const std::string &cemc_tower_node_name)
{
  cemc_tower_node_name_ = cemc_tower_node_name;
}

void Pi0Reconstruction::set_abort_on_missing_cemc_tower_node(bool abort_on_missing_cemc_tower_node)
{
  abort_on_missing_cemc_tower_node_ = abort_on_missing_cemc_tower_node;
}

void Pi0Reconstruction::set_shower_shape_min_tower_energy(double min_tower_energy)
{
  if (!std::isfinite(min_tower_energy) || min_tower_energy < 0.0)
  {
    std::cout << "Pi0Reconstruction::set_shower_shape_min_tower_energy - invalid threshold requested; keeping current setting" << std::endl;
    return;
  }
  shower_shape_min_tower_energy_ = min_tower_energy;
}

void Pi0Reconstruction::set_store_shower_shape_tower_patch(bool store_tower_patch)
{
  store_shower_shape_tower_patch_ = store_tower_patch;
}

void Pi0Reconstruction::set_min_cluster_energy(double min_cluster_energy)
{
  min_cluster_energy_ = min_cluster_energy;
}

void Pi0Reconstruction::set_mass_histogram_bins(int nbins, double min, double max)
{
  if (nbins <= 0 || min >= max)
  {
    std::cout << "Pi0Reconstruction::set_mass_histogram_bins - invalid binning requested; keeping current settings" << std::endl;
    return;
  }

  mass_histogram_nbins_ = nbins;
  mass_histogram_min_ = min;
  mass_histogram_max_ = max;
}

int Pi0Reconstruction::Init(PHCompositeNode * /*topNode*/)
{
  create_output_directory();

  ShowerShapeCalculator::Config shower_shape_config;
  shower_shape_config.min_tower_energy = static_cast<float>(shower_shape_min_tower_energy_);
  shower_shape_calculator_ = ShowerShapeCalculator(shower_shape_config);

  output_file_ = TFile::Open(output_file_name_.c_str(), "RECREATE");
  if (!output_file_ || output_file_->IsZombie())
  {
    std::cout << "Pi0Reconstruction::Init - failed to open output file: " << output_file_name_ << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  h_m_gg_ = new TH1D("h_m_gg", "CEMC cluster pair invariant mass;M_{#gamma#gamma} [GeV];Pairs", mass_histogram_nbins_, mass_histogram_min_, mass_histogram_max_);
  h_ncluster_ = new TH1D("h_ncluster", "CEMC clusters per event;N_{cluster};Events", 100, 0.0, 100.0);
  h_cluster_e_ = new TH1D("h_cluster_e", "CEMC cluster energy;E_{cluster} [GeV];Clusters", 200, 0.0, 20.0);
  h_pair_e_asym_ = new TH1D("h_pair_e_asym", "CEMC cluster pair energy asymmetry;(|E_{1}-E_{2}|)/(E_{1}+E_{2});Pairs", 100, -1.0, 1.0);

  event_tree_ = new TTree("event_tree", "Pi0 reconstruction event tree");
  create_tree_branches();

  std::cout << "Pi0Reconstruction::Init - writing output to " << output_file_name_ << std::endl;
  std::cout << "Pi0Reconstruction::Init - cluster node: " << cluster_node_name_ << std::endl;
  std::cout << "Pi0Reconstruction::Init - CEMC tower node: " << cemc_tower_node_name_ << std::endl;
  std::cout << "Pi0Reconstruction::Init - shower-shape tower threshold: " << shower_shape_min_tower_energy_ << " GeV" << std::endl;
  std::cout << "Pi0Reconstruction::Init - process ID: " << process_id_ << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

int Pi0Reconstruction::InitRun(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int Pi0Reconstruction::process_event(PHCompositeNode *topNode)
{
  const unsigned int event_number = event_counter_;
  if (event_number % 200 == 0)
  {
    std::cout << "Pi0Reconstruction::process_event - event " << event_number << std::endl;
  }
  ++event_counter_;

  reset_tree_variables();
  tree_process_id_ = process_id_;
  tree_event_ = event_number;
  tree_event_uid_ = (static_cast<unsigned long long>(process_id_) << 32U) | static_cast<unsigned long long>(event_number); // Fills the upper 32 bits with process_id_ and the lower 32 bits with event_number 
  tree_min_cluster_energy_ = min_cluster_energy_;
  tree_shower_shape_min_tower_energy_ = shower_shape_min_tower_energy_;
  tree_shower_shape_algorithm_version_ = ShowerShapeCalculator::kAlgorithmVersion;
  tree_shower_shape_patch_side_ = ShowerShapeCalculator::kPatchSide;
  tree_store_shower_shape_tower_patch_ = store_shower_shape_tower_patch_;

  RawClusterContainer *cluster_container = findNode::getClass<RawClusterContainer>(topNode, cluster_node_name_);
  if (!cluster_container)
  {
    if (missing_cluster_node_warnings_ < 5)
    {
      std::cout << "Pi0Reconstruction::process_event - missing required cluster node: " << cluster_node_name_ << std::endl;
    }
    ++missing_cluster_node_warnings_;

    if (!abort_on_missing_cluster_node_ && event_tree_)
    {
      event_tree_->Fill();
    }
    return abort_on_missing_cluster_node_ ? Fun4AllReturnCodes::ABORTRUN : Fun4AllReturnCodes::EVENT_OK;
  }

  TowerInfoContainer *cemc_towers = findNode::getClass<TowerInfoContainer>(topNode, cemc_tower_node_name_);
  if (!cemc_towers)
  {
    if (missing_cemc_tower_node_warnings_ < 5)
    {
      std::cout << "Pi0Reconstruction::process_event - missing CEMC tower node required for shower shapes: "
                << cemc_tower_node_name_ << std::endl;
    }
    ++missing_cemc_tower_node_warnings_;
    if (abort_on_missing_cemc_tower_node_)
    {
      return Fun4AllReturnCodes::ABORTRUN;
    }
  }

  std::array<double, 3> vertex = {0.0, 0.0, 0.0};
  if (!get_event_vertex(topNode, vertex))
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }
  tree_vertex_x_ = vertex[0];
  tree_vertex_y_ = vertex[1];
  tree_vertex_z_ = vertex[2];

  std::vector<PhotonCandidate> all_photons;
  std::vector<unsigned int> selected_photon_indices;
  all_photons.reserve(cluster_container->size());
  selected_photon_indices.reserve(cluster_container->size());

  RawClusterContainer::ConstRange cluster_range = cluster_container->getClusters();
  for (RawClusterContainer::ConstIterator iter = cluster_range.first; iter != cluster_range.second; ++iter)
  {
    const RawCluster *cluster = iter->second;
    if (!cluster)
    {
      continue;
    }

    const double energy = cluster->get_energy();
    PhotonCandidate candidate;
    if (!build_photon_candidate(energy, cluster->get_x(), cluster->get_y(), cluster->get_z(), vertex, candidate))
    {
      continue;
    }

    const unsigned int cluster_index = static_cast<unsigned int>(all_photons.size());
    all_photons.push_back(candidate);
    tree_cluster_e_.push_back(candidate.energy);
    tree_cluster_et_.push_back(candidate.transverse_energy);
    tree_cluster_eta_.push_back(candidate.eta);
    tree_cluster_x_.push_back(cluster->get_x());
    tree_cluster_y_.push_back(cluster->get_y());
    tree_cluster_z_.push_back(cluster->get_z());
    tree_cluster_px_.push_back(candidate.momentum[0]);
    tree_cluster_py_.push_back(candidate.momentum[1]);
    tree_cluster_pz_.push_back(candidate.momentum[2]);

    ShowerShapeCalculator::Result shower_shape;
    if (cemc_towers)
    {
      shower_shape = shower_shape_calculator_.calculate(*cluster, *cemc_towers);
    }
    append_shower_shape_result(shower_shape);

    if (energy < min_cluster_energy_)
    {
      continue;
    }

    h_cluster_e_->Fill(energy);
    selected_photon_indices.push_back(cluster_index);
  }

  tree_ncluster_all_ = static_cast<unsigned int>(all_photons.size());
  tree_ncluster_ = static_cast<unsigned int>(selected_photon_indices.size());
  h_ncluster_->Fill(static_cast<double>(tree_ncluster_));
  for (std::size_t i = 0; i < selected_photon_indices.size(); ++i)
  {
    for (std::size_t j = i + 1; j < selected_photon_indices.size(); ++j)
    {
      const unsigned int first_index = selected_photon_indices[i];
      const unsigned int second_index = selected_photon_indices[j];
      const PhotonCandidate &first = all_photons[first_index];
      const PhotonCandidate &second = all_photons[second_index];

      const double total_energy = first.energy + second.energy;
      const double px = first.momentum[0] + second.momentum[0];
      const double py = first.momentum[1] + second.momentum[1];
      const double pz = first.momentum[2] + second.momentum[2];
      const double mass2 = total_energy * total_energy - px * px - py * py - pz * pz;
      const double mass = std::sqrt(std::max(0.0, mass2));

      const double energy_asymmetry = total_energy > 0.0 ? std::abs(first.energy - second.energy) / total_energy : -1.0;

      tree_pair_cluster_i_.push_back(first_index);
      tree_pair_cluster_j_.push_back(second_index);
      tree_pair_m_gg_.push_back(mass);
      tree_pair_e_asym_.push_back(energy_asymmetry);

      h_m_gg_->Fill(mass);
      if (total_energy > 0.0)
      {
        h_pair_e_asym_->Fill(energy_asymmetry);
      }
    }
  }

  if (event_tree_)
  {
    event_tree_->Fill();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int Pi0Reconstruction::ResetEvent(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int Pi0Reconstruction::Reset(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int Pi0Reconstruction::EndRun(const int /*runnumber*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int Pi0Reconstruction::End(PHCompositeNode * /*topNode*/)
{
  if (!output_file_)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  output_file_->cd();
  h_m_gg_->Write();
  h_ncluster_->Write();
  h_cluster_e_->Write();
  h_pair_e_asym_->Write();
  if (event_tree_)
  {
    event_tree_->Write();
  }
  output_file_->Close();
  delete output_file_;
  output_file_ = nullptr;

  std::cout << "Pi0Reconstruction::End - processed " << event_counter_ << " events" << std::endl;
  std::cout << "Pi0Reconstruction::End - wrote " << output_file_name_ << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

bool Pi0Reconstruction::get_event_vertex(PHCompositeNode *topNode, std::array<double, 3> &vertex)
{
  vertex = {0.0, 0.0, 0.0};

  if (vertex_mode_ == VertexMode::Origin)
  {
    return true;
  }

  GlobalVertexMap *vertex_map = findNode::getClass<GlobalVertexMap>(topNode, vertex_node_name_);
  if (!vertex_map || vertex_map->empty())
  {
    if (missing_vertex_node_warnings_ < 5)
    {
      std::cout << "Pi0Reconstruction::get_event_vertex - missing or empty vertex node: " << vertex_node_name_ << "; using origin" << std::endl;
    }
    ++missing_vertex_node_warnings_;
    return !abort_on_missing_vertex_node_;
  }

  const GlobalVertex *global_vertex = vertex_map->begin()->second;
  if (!global_vertex)
  {
    return !abort_on_missing_vertex_node_;
  }

  vertex[0] = global_vertex->get_x();
  vertex[1] = global_vertex->get_y();
  vertex[2] = global_vertex->get_z();

  if (!std::isfinite(vertex[0]) || !std::isfinite(vertex[1]) || !std::isfinite(vertex[2]))
  {
    vertex = {0.0, 0.0, 0.0};
    return !abort_on_missing_vertex_node_;
  }

  return true;
}

bool Pi0Reconstruction::build_photon_candidate(double energy, double x, double y, double z, const std::array<double, 3> &vertex, PhotonCandidate &candidate) const
{
  if (!std::isfinite(energy) || !std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z))
  {
    return false;
  }

  const double dx = x - vertex[0];
  const double dy = y - vertex[1];
  const double dz = z - vertex[2];
  const double distance = std::sqrt(dx * dx + dy * dy + dz * dz);
  if (distance <= std::numeric_limits<double>::epsilon())
  {
    return false;
  }

  candidate.energy = energy;
  candidate.momentum = {energy * dx / distance, energy * dy / distance, energy * dz / distance};
  candidate.transverse_energy = std::hypot(candidate.momentum[0], candidate.momentum[1]);
  candidate.eta = candidate.transverse_energy > std::numeric_limits<double>::epsilon()
      ? std::asinh(candidate.momentum[2] / candidate.transverse_energy)
      : std::numeric_limits<double>::quiet_NaN();
  return true;
}

void Pi0Reconstruction::append_shower_shape_result(const ShowerShapeCalculator::Result &result)
{
  tree_cluster_shower_valid_.push_back(result.valid ? 1U : 0U);
  tree_cluster_shower_full_containment_.push_back(result.full_containment ? 1U : 0U);
  tree_cluster_shower_edge_padded_.push_back(result.edge_padded ? 1U : 0U);
  tree_cluster_shower_tower_data_complete_.push_back(result.tower_data_complete ? 1U : 0U);
  tree_cluster_shower_cog_ieta_.push_back(result.cog_ieta);
  tree_cluster_shower_cog_iphi_.push_back(result.cog_iphi);
  tree_cluster_shower_cluster_e_thresholded_.push_back(result.cluster_energy_above_threshold);
  tree_cluster_shower_owned_patch_e_.push_back(result.owned_patch_energy);
  tree_cluster_shower_w_eta_cogx_.push_back(result.w_eta_cogx);
  tree_cluster_shower_w_phi_cogx_.push_back(result.w_phi_cogx);
  tree_cluster_shower_e11_.push_back(result.e11);
  tree_cluster_shower_e33_.push_back(result.e33);
  tree_cluster_shower_e32_.push_back(result.e32);
  tree_cluster_shower_e35_.push_back(result.e35);
  tree_cluster_shower_e11_over_e33_.push_back(result.e11_over_e33);
  tree_cluster_shower_e32_over_e35_.push_back(result.e32_over_e35);
  tree_cluster_shower_et1_.push_back(result.et1);
  tree_cluster_shower_et2_.push_back(result.et2);
  tree_cluster_shower_et3_.push_back(result.et3);
  tree_cluster_shower_et4_.push_back(result.et4);

  if (store_shower_shape_tower_patch_)
  {
    tree_cluster_shower_patch_e_.insert(
        tree_cluster_shower_patch_e_.end(), result.patch_energy.begin(), result.patch_energy.end());
    tree_cluster_shower_patch_good_.insert(
        tree_cluster_shower_patch_good_.end(), result.patch_good.begin(), result.patch_good.end());
    tree_cluster_shower_patch_owned_.insert(
        tree_cluster_shower_patch_owned_.end(), result.patch_owned.begin(), result.patch_owned.end());
  }
}

void Pi0Reconstruction::reset_tree_variables()
{
  tree_process_id_ = 0;
  tree_event_ = 0;
  tree_event_uid_ = 0;
  tree_ncluster_ = 0;
  tree_ncluster_all_ = 0;
  tree_min_cluster_energy_ = 0.0;
  tree_vertex_x_ = 0.0;
  tree_vertex_y_ = 0.0;
  tree_vertex_z_ = 0.0;
  tree_shower_shape_min_tower_energy_ = 0.0;
  tree_shower_shape_algorithm_version_ = 0;
  tree_shower_shape_patch_side_ = 0;
  tree_store_shower_shape_tower_patch_ = false;

  tree_cluster_e_.clear();
  tree_cluster_et_.clear();
  tree_cluster_eta_.clear();
  tree_cluster_x_.clear();
  tree_cluster_y_.clear();
  tree_cluster_z_.clear();
  tree_cluster_px_.clear();
  tree_cluster_py_.clear();
  tree_cluster_pz_.clear();
  tree_cluster_shower_valid_.clear();
  tree_cluster_shower_full_containment_.clear();
  tree_cluster_shower_edge_padded_.clear();
  tree_cluster_shower_tower_data_complete_.clear();
  tree_cluster_shower_cog_ieta_.clear();
  tree_cluster_shower_cog_iphi_.clear();
  tree_cluster_shower_cluster_e_thresholded_.clear();
  tree_cluster_shower_owned_patch_e_.clear();
  tree_cluster_shower_w_eta_cogx_.clear();
  tree_cluster_shower_w_phi_cogx_.clear();
  tree_cluster_shower_e11_.clear();
  tree_cluster_shower_e33_.clear();
  tree_cluster_shower_e32_.clear();
  tree_cluster_shower_e35_.clear();
  tree_cluster_shower_e11_over_e33_.clear();
  tree_cluster_shower_e32_over_e35_.clear();
  tree_cluster_shower_et1_.clear();
  tree_cluster_shower_et2_.clear();
  tree_cluster_shower_et3_.clear();
  tree_cluster_shower_et4_.clear();
  tree_cluster_shower_patch_e_.clear();
  tree_cluster_shower_patch_good_.clear();
  tree_cluster_shower_patch_owned_.clear();
  tree_pair_cluster_i_.clear();
  tree_pair_cluster_j_.clear();
  tree_pair_m_gg_.clear();
  tree_pair_e_asym_.clear();
}

void Pi0Reconstruction::create_tree_branches()
{
  if (!event_tree_)
  {
    return;
  }

  event_tree_->Branch("process_id", &tree_process_id_);
  event_tree_->Branch("event", &tree_event_);
  event_tree_->Branch("event_uid", &tree_event_uid_);
  event_tree_->Branch("ncluster", &tree_ncluster_);
  event_tree_->Branch("ncluster_all", &tree_ncluster_all_);
  event_tree_->Branch("min_cluster_energy", &tree_min_cluster_energy_);
  event_tree_->Branch("vertex_x", &tree_vertex_x_);
  event_tree_->Branch("vertex_y", &tree_vertex_y_);
  event_tree_->Branch("vertex_z", &tree_vertex_z_);
  event_tree_->Branch("shower_shape_min_tower_energy", &tree_shower_shape_min_tower_energy_);
  event_tree_->Branch("shower_shape_algorithm_version", &tree_shower_shape_algorithm_version_);
  event_tree_->Branch("shower_shape_patch_side", &tree_shower_shape_patch_side_);
  event_tree_->Branch("store_shower_shape_tower_patch", &tree_store_shower_shape_tower_patch_);

  event_tree_->Branch("cluster_e", &tree_cluster_e_);
  event_tree_->Branch("cluster_et", &tree_cluster_et_);
  event_tree_->Branch("cluster_eta", &tree_cluster_eta_);
  event_tree_->Branch("cluster_x", &tree_cluster_x_);
  event_tree_->Branch("cluster_y", &tree_cluster_y_);
  event_tree_->Branch("cluster_z", &tree_cluster_z_);
  event_tree_->Branch("cluster_px", &tree_cluster_px_);
  event_tree_->Branch("cluster_py", &tree_cluster_py_);
  event_tree_->Branch("cluster_pz", &tree_cluster_pz_);
  event_tree_->Branch("cluster_shower_valid", &tree_cluster_shower_valid_);
  event_tree_->Branch("cluster_shower_full_containment", &tree_cluster_shower_full_containment_);
  event_tree_->Branch("cluster_shower_edge_padded", &tree_cluster_shower_edge_padded_);
  event_tree_->Branch("cluster_shower_tower_data_complete", &tree_cluster_shower_tower_data_complete_);
  event_tree_->Branch("cluster_shower_cog_ieta", &tree_cluster_shower_cog_ieta_);
  event_tree_->Branch("cluster_shower_cog_iphi", &tree_cluster_shower_cog_iphi_);
  event_tree_->Branch("cluster_shower_cluster_e_thresholded", &tree_cluster_shower_cluster_e_thresholded_);
  event_tree_->Branch("cluster_shower_owned_patch_e", &tree_cluster_shower_owned_patch_e_);
  event_tree_->Branch("cluster_shower_w_eta_cogx", &tree_cluster_shower_w_eta_cogx_);
  event_tree_->Branch("cluster_shower_w_phi_cogx", &tree_cluster_shower_w_phi_cogx_);
  event_tree_->Branch("cluster_shower_e11", &tree_cluster_shower_e11_);
  event_tree_->Branch("cluster_shower_e33", &tree_cluster_shower_e33_);
  event_tree_->Branch("cluster_shower_e32", &tree_cluster_shower_e32_);
  event_tree_->Branch("cluster_shower_e35", &tree_cluster_shower_e35_);
  event_tree_->Branch("cluster_shower_e11_over_e33", &tree_cluster_shower_e11_over_e33_);
  event_tree_->Branch("cluster_shower_e32_over_e35", &tree_cluster_shower_e32_over_e35_);
  event_tree_->Branch("cluster_shower_et1", &tree_cluster_shower_et1_);
  event_tree_->Branch("cluster_shower_et2", &tree_cluster_shower_et2_);
  event_tree_->Branch("cluster_shower_et3", &tree_cluster_shower_et3_);
  event_tree_->Branch("cluster_shower_et4", &tree_cluster_shower_et4_);
  event_tree_->Branch("cluster_shower_patch_e", &tree_cluster_shower_patch_e_);
  event_tree_->Branch("cluster_shower_patch_good", &tree_cluster_shower_patch_good_);
  event_tree_->Branch("cluster_shower_patch_owned", &tree_cluster_shower_patch_owned_);
  event_tree_->Branch("pair_cluster_i", &tree_pair_cluster_i_);
  event_tree_->Branch("pair_cluster_j", &tree_pair_cluster_j_);
  event_tree_->Branch("pair_m_gg", &tree_pair_m_gg_);
  event_tree_->Branch("pair_e_asym", &tree_pair_e_asym_);
}

void Pi0Reconstruction::create_output_directory() const
{
  const std::string::size_type slash_position = output_file_name_.find_last_of('/');
  if (slash_position == std::string::npos)
  {
    return;
  }

  const std::string directory = output_file_name_.substr(0, slash_position);
  if (!directory.empty())
  {
    gSystem->mkdir(directory.c_str(), true);
  }
}
