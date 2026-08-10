#include "PythiaPhotonAnalysisTree.h"

#include <calobase/RawClusterContainer.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/TowerInfoContainer.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <HepMC/GenEvent.h>
#include <HepMC/SimpleVector.h>

#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>

#include <cmath>
#include <iostream>

PythiaPhotonAnalysisTree::PythiaPhotonAnalysisTree(const std::string& name)
  : SubsysReco(name)
{
}

PythiaPhotonAnalysisTree::~PythiaPhotonAnalysisTree()
{
  close_output_file();
}

int PythiaPhotonAnalysisTree::Init(PHCompositeNode* /*topNode*/)
{
  if (signal_embedding_id_ <= 0 || min_cluster_energy_ < 0.0 ||
      shower_shape_min_tower_energy_ < 0.0)
  {
    std::cout << "PythiaPhotonAnalysisTree::Init - invalid numeric configuration" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  common_.set_min_cluster_energy(min_cluster_energy_);
  common_.set_shower_shape_min_tower_energy(shower_shape_min_tower_energy_);
  common_.set_store_shower_shape_tower_patch(store_shower_shape_tower_patch_);
  if (!common_.initialize())
  {
    std::cout << "PythiaPhotonAnalysisTree::Init - failed to initialize common reco filler" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  truth_matcher_.set_verbosity(verbosity_);

  create_output_directory();
  output_file_ = TFile::Open(output_file_name_.c_str(), "RECREATE");
  if (!output_file_ || output_file_->IsZombie())
  {
    std::cout << "PythiaPhotonAnalysisTree::Init - failed to create " << output_file_name_ << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  create_trees();

  std::cout << "PythiaPhotonAnalysisTree::Init - input: " << input_file_name_ << '\n'
            << "PythiaPhotonAnalysisTree::Init - output: " << output_file_name_ << '\n'
            << "PythiaPhotonAnalysisTree::Init - signal embedding ID: " << signal_embedding_id_ << '\n'
            << "PythiaPhotonAnalysisTree::Init - split/no-split nodes: "
            << split_cluster_node_name_ << "/" << nosplit_cluster_node_name_ << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

int PythiaPhotonAnalysisTree::process_event(PHCompositeNode* topNode)
{
  reset_event();
  b_source_file_id_ = source_file_id_;
  b_event_in_file_ = static_cast<unsigned int>(n_events_processed_);
  b_event_uid_ = (static_cast<unsigned long long>(source_file_id_) << 32U) |
                 static_cast<unsigned long long>(b_event_in_file_);
  ++n_events_processed_;

  auto* truth = findNode::getClass<PHG4TruthInfoContainer>(topNode, truth_node_name_);
  auto* event_map = findNode::getClass<PHHepMCGenEventMap>(topNode, hepmc_event_map_node_name_);
  auto* towers = findNode::getClass<TowerInfoContainer>(topNode, tower_node_name_);
  auto* raw_truth_towers = findNode::getClass<RawTowerContainer>(topNode, raw_truth_tower_node_name_);
  auto* geometry = findNode::getClass<RawTowerGeomContainer>(topNode, tower_geom_node_name_);
  auto* split_clusters = findNode::getClass<RawClusterContainer>(topNode, split_cluster_node_name_);
  auto* nosplit_clusters = findNode::getClass<RawClusterContainer>(topNode, nosplit_cluster_node_name_);

  const bool missing_required_node =
      !truth || !event_map || !towers || !geometry || !split_clusters ||
      (require_nosplit_cluster_node_ && !nosplit_clusters);
  if (missing_required_node)
  {
    ++n_events_invalid_detector_;
    std::cout << "PythiaPhotonAnalysisTree::process_event - missing node in event "
              << b_event_in_file_ << ": truth=" << static_cast<bool>(truth)
              << " hepmc=" << static_cast<bool>(event_map)
              << " towers=" << static_cast<bool>(towers)
              << " raw_truth_towers=" << static_cast<bool>(raw_truth_towers)
              << " geometry=" << static_cast<bool>(geometry)
              << " split=" << static_cast<bool>(split_clusters)
              << " nosplit=" << static_cast<bool>(nosplit_clusters) << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  photon_tree::EventVertex event_vertex;
  if (!fill_hepmc_event_info(event_map, event_vertex))
  {
    ++n_events_invalid_truth_;
    if (verbosity_ > 0)
    {
      std::cout << "PythiaPhotonAnalysisTree::process_event - missing/invalid signal HepMC event "
                << signal_embedding_id_ << " in event " << b_event_in_file_ << std::endl;
    }
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  if (!truth_matcher_.begin_event(topNode))
  {
    ++n_events_invalid_detector_;
    std::cout << "PythiaPhotonAnalysisTree::process_event - CaloEvalStack lacks reduced-sim nodes"
              << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  const bool split_valid = common_.fill_collection(
      split_clusters, towers, geometry, event_vertex, true, false, common_.split());
  const bool nosplit_valid = !nosplit_clusters || common_.fill_collection(
      nosplit_clusters, towers, geometry, event_vertex, true, true, common_.nosplit());
  const bool split_truth_valid = split_valid && fill_cluster_truth(
      common_.split(), true, towers, raw_truth_towers, truth, event_map, split_truth_);
  const bool nosplit_truth_valid = nosplit_valid && (!nosplit_clusters || fill_cluster_truth(
      common_.nosplit(), false, towers, raw_truth_towers, truth, event_map, nosplit_truth_));
  if (!split_valid || !nosplit_valid)
  {
    ++n_events_invalid_detector_;
    std::cout << "PythiaPhotonAnalysisTree::process_event - invalid cluster/tower content in event "
              << b_event_in_file_ << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  if (!split_truth_valid || !nosplit_truth_valid)
  {
    ++n_events_invalid_truth_;
  }

  event_tree_->Fill();
  ++n_events_written_;
  if (verbosity_ > 0 && b_event_in_file_ < 5)
  {
    std::cout << "PythiaPhotonAnalysisTree::process_event - event " << b_event_in_file_
              << " subevents=" << b_hepmc_n_subevent_
              << " split/nosplit=" << common_.split().data.ncluster << "/"
              << common_.nosplit().data.ncluster << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int PythiaPhotonAnalysisTree::End(PHCompositeNode* /*topNode*/)
{
  close_output_file();
  std::cout << "PythiaPhotonAnalysisTree::End - processed/written/invalid_truth/invalid_detector = "
            << n_events_processed_ << "/" << n_events_written_ << "/"
            << n_events_invalid_truth_ << "/" << n_events_invalid_detector_ << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

bool PythiaPhotonAnalysisTree::fill_hepmc_event_info(
    const PHHepMCGenEventMap* event_map,
    photon_tree::EventVertex& event_vertex)
{
  if (!event_map)
  {
    return false;
  }
  b_hepmc_n_subevent_ = static_cast<unsigned int>(event_map->size());
  for (const auto& entry : event_map->get_map())
  {
    b_hepmc_embedding_id_.push_back(entry.first);
  }

  const PHHepMCGenEvent* signal_event = event_map->get(signal_embedding_id_);
  const HepMC::GenEvent* hepmc_event = signal_event ? signal_event->getEvent() : nullptr;
  if (!signal_event || !hepmc_event || !signal_event->is_simulated())
  {
    return false;
  }
  const HepMC::FourVector& vertex = signal_event->get_collision_vertex();
  if (!std::isfinite(vertex.x()) || !std::isfinite(vertex.y()) ||
      !std::isfinite(vertex.z()) || !std::isfinite(vertex.t()))
  {
    return false;
  }

  b_vertex_valid_ = 1U;
  b_vertex_embedding_id_ = signal_embedding_id_;
  b_vertex_x_ = vertex.x();
  b_vertex_y_ = vertex.y();
  b_vertex_z_ = vertex.z();
  b_vertex_t_ = vertex.t();
  b_hepmc_signal_event_number_ = hepmc_event->event_number();
  event_vertex.valid = true;
  event_vertex.x = b_vertex_x_;
  event_vertex.y = b_vertex_y_;
  event_vertex.z = b_vertex_z_;
  event_vertex.source = signal_embedding_id_;
  return true;
}

bool PythiaPhotonAnalysisTree::fill_cluster_truth(
    const photon_tree::FilledCollection& reco,
    bool allocate_split_tower_energy,
    TowerInfoContainer* towers,
    RawTowerContainer* raw_truth_towers,
    PHG4TruthInfoContainer* truth,
    const PHHepMCGenEventMap* event_map,
    photon_tree::ClusterTruthCollection& output) const
{
  output.clear();
  bool all_valid = true;
  for (const RawCluster* cluster : reco.ordered_clusters)
  {
    const photon_tree::ClusterTruthMatch match = truth_matcher_.match(
        cluster, towers, raw_truth_towers, truth, event_map, allocate_split_tower_energy);
    output.append(match);
    all_valid = all_valid && match.valid;
  }
  return all_valid;
}

void PythiaPhotonAnalysisTree::reset_event()
{
  b_source_file_id_ = 0;
  b_event_in_file_ = 0;
  b_event_uid_ = 0;
  b_vertex_valid_ = 0;
  b_vertex_embedding_id_ = -999;
  b_vertex_x_ = b_vertex_y_ = b_vertex_z_ = b_vertex_t_ = photon_tree::kInvalidDouble;
  b_hepmc_n_subevent_ = 0;
  b_hepmc_signal_event_number_ = -999;
  b_hepmc_embedding_id_.clear();
  common_.clear_event();
  split_truth_.clear();
  nosplit_truth_.clear();
}

void PythiaPhotonAnalysisTree::create_trees()
{
  output_file_->cd();
  event_tree_ = new TTree("event_tree", "Pythia cluster and contributor truth analysis tree");
  event_tree_->Branch("source_file_id", &b_source_file_id_);
  event_tree_->Branch("event_in_file", &b_event_in_file_);
  event_tree_->Branch("event_uid", &b_event_uid_);
  event_tree_->Branch("vertex_valid", &b_vertex_valid_);
  event_tree_->Branch("vertex_embedding_id", &b_vertex_embedding_id_);
  event_tree_->Branch("vertex_x", &b_vertex_x_);
  event_tree_->Branch("vertex_y", &b_vertex_y_);
  event_tree_->Branch("vertex_z", &b_vertex_z_);
  event_tree_->Branch("vertex_t", &b_vertex_t_);
  event_tree_->Branch("hepmc_n_subevent", &b_hepmc_n_subevent_);
  event_tree_->Branch("hepmc_signal_event_number", &b_hepmc_signal_event_number_);
  event_tree_->Branch("hepmc_embedding_id", &b_hepmc_embedding_id_);
  event_tree_->Branch("min_cluster_energy", &min_cluster_energy_);
  event_tree_->Branch("shower_shape_min_tower_energy", &shower_shape_min_tower_energy_);
  event_tree_->Branch("shower_shape_algorithm_version", &shower_shape_algorithm_version_);
  event_tree_->Branch("shower_shape_patch_side", &shower_shape_patch_side_);
  event_tree_->Branch("store_shower_shape_tower_patch", &store_shower_shape_tower_patch_);
  common_.create_collection_branches(event_tree_, "split", common_.split().data, true);
  common_.create_collection_branches(event_tree_, "nosplit", common_.nosplit().data, true);
  split_truth_.create_branches(event_tree_, "split");
  nosplit_truth_.create_branches(event_tree_, "nosplit");

  metadata_tree_ = new TTree("metadata", "One entry per Pythia source DST set");
  metadata_tree_->Branch("schema_version", &metadata_schema_version_);
  metadata_tree_->Branch("sample_type", &sample_type_);
  metadata_tree_->Branch("input_file", &input_file_name_);
  metadata_tree_->Branch("output_file", &output_file_name_);
  metadata_tree_->Branch("source_file_id", &source_file_id_);
  metadata_tree_->Branch("signal_embedding_id", &signal_embedding_id_);
  metadata_tree_->Branch("truth_scheme", &truth_scheme_);
  metadata_tree_->Branch("truth_matcher_version", &metadata_matcher_version_);
  metadata_tree_->Branch("photon_classifier_version", &metadata_classifier_version_);
  metadata_tree_->Branch("split_truth_allocation", &split_truth_allocation_);
  metadata_tree_->Branch("nosplit_truth_allocation", &nosplit_truth_allocation_);
  metadata_tree_->Branch("truth_node", &truth_node_name_);
  metadata_tree_->Branch("hepmc_event_map_node", &hepmc_event_map_node_name_);
  metadata_tree_->Branch("tower_node", &tower_node_name_);
  metadata_tree_->Branch("tower_geom_node", &tower_geom_node_name_);
  metadata_tree_->Branch("raw_truth_tower_node", &raw_truth_tower_node_name_);
  metadata_tree_->Branch("split_cluster_node", &split_cluster_node_name_);
  metadata_tree_->Branch("nosplit_cluster_node", &nosplit_cluster_node_name_);
  metadata_tree_->Branch("require_nosplit_cluster_node", &require_nosplit_cluster_node_);
  metadata_tree_->Branch("cluster_ordering", &cluster_ordering_);
  metadata_tree_->Branch("n_events_processed", &n_events_processed_);
  metadata_tree_->Branch("n_events_written", &n_events_written_);
  metadata_tree_->Branch("n_events_invalid_truth", &n_events_invalid_truth_);
  metadata_tree_->Branch("n_events_invalid_detector", &n_events_invalid_detector_);
}

void PythiaPhotonAnalysisTree::create_output_directory() const
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

void PythiaPhotonAnalysisTree::close_output_file()
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
