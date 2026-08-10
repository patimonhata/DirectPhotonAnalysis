#ifndef RYOTARO_PYTHIAPHOTONANALYSISTREE_H_20260810
#define RYOTARO_PYTHIAPHOTONANALYSISTREE_H_20260810

#include "HepMCPhotonClassifier.h"
#include "PhotonTreeCommon.h"
#include "PythiaClusterTruthMatcher.h"

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>

class PHCompositeNode;
class PHG4TruthInfoContainer;
class PHHepMCGenEventMap;
class RawTowerContainer;
class TowerInfoContainer;
class TFile;
class TTree;

class PythiaPhotonAnalysisTree : public SubsysReco
{
 public:
  explicit PythiaPhotonAnalysisTree(const std::string& name = "PythiaPhotonAnalysisTree");
  ~PythiaPhotonAnalysisTree() override;

  int Init(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int End(PHCompositeNode* topNode) override;

  void set_input_file_name(const std::string& value) { input_file_name_ = value; }
  void set_output_file_name(const std::string& value) { output_file_name_ = value; }
  void set_source_file_id(unsigned int value) { source_file_id_ = value; }
  void set_signal_embedding_id(int value) { signal_embedding_id_ = value; }
  void set_truth_node_name(const std::string& value) { truth_node_name_ = value; }
  void set_hepmc_event_map_node_name(const std::string& value) { hepmc_event_map_node_name_ = value; }
  void set_tower_node_name(const std::string& value) { tower_node_name_ = value; }
  void set_raw_truth_tower_node_name(const std::string& value) { raw_truth_tower_node_name_ = value; }
  void set_tower_geom_node_name(const std::string& value) { tower_geom_node_name_ = value; }
  void set_split_cluster_node_name(const std::string& value) { split_cluster_node_name_ = value; }
  void set_nosplit_cluster_node_name(const std::string& value) { nosplit_cluster_node_name_ = value; }
  void set_require_nosplit_cluster_node(bool value) { require_nosplit_cluster_node_ = value; }
  void set_min_cluster_energy(double value) { min_cluster_energy_ = value; }
  void set_shower_shape_min_tower_energy(double value) { shower_shape_min_tower_energy_ = value; }
  void set_store_shower_shape_tower_patch(bool value) { store_shower_shape_tower_patch_ = value; }
  void set_verbosity(int value) { verbosity_ = value; }

 private:
  static constexpr int schema_version_ = 4;

  void create_output_directory() const;
  void create_trees();
  void close_output_file();
  void reset_event();
  bool fill_hepmc_event_info(const PHHepMCGenEventMap* event_map,
                             photon_tree::EventVertex& event_vertex);
  bool fill_cluster_truth(const photon_tree::FilledCollection& reco,
                          bool allocate_split_tower_energy,
                          TowerInfoContainer* towers,
                          RawTowerContainer* raw_truth_towers,
                          PHG4TruthInfoContainer* truth,
                          const PHHepMCGenEventMap* event_map,
                          photon_tree::ClusterTruthCollection& output) const;

  std::string input_file_name_;
  std::string output_file_name_ = "pythia_photon_analysis_tree.root";
  std::string truth_node_name_ = "G4TruthInfo";
  std::string hepmc_event_map_node_name_ = "PHHepMCGenEventMap";
  std::string tower_node_name_ = "TOWERINFO_CALIB_CEMC";
  std::string raw_truth_tower_node_name_ = "TOWER_CALIB_CEMC";
  std::string tower_geom_node_name_ = "TOWERGEOM_CEMC";
  std::string split_cluster_node_name_ = "CLUSTERINFO_CEMC";
  std::string nosplit_cluster_node_name_ = "CLUSTERINFO_CEMC_NO_SPLIT";
  std::string sample_type_ = "pythia";
  std::string truth_scheme_ = "towerinfo_or_rawtower_g4_shower_edep_to_primary_then_hepmc_barcode";
  std::string split_truth_allocation_ = "cluster_tower_energy_over_tower_energy_clamped_0_1";
  std::string nosplit_truth_allocation_ = "full_tower_shower_edep";
  std::string cluster_ordering_ = "energy_descending_then_cluster_id_ascending";
  unsigned int source_file_id_ = 0;
  int signal_embedding_id_ = 1;
  bool require_nosplit_cluster_node_ = false;
  double min_cluster_energy_ = 0.0;
  double shower_shape_min_tower_energy_ = 0.070;
  bool store_shower_shape_tower_patch_ = true;
  int verbosity_ = 0;

  photon_tree::PhotonTreeCommon common_;
  photon_tree::PythiaClusterTruthMatcher truth_matcher_;
  photon_tree::ClusterTruthCollection split_truth_;
  photon_tree::ClusterTruthCollection nosplit_truth_;

  TFile* output_file_ = nullptr;
  TTree* event_tree_ = nullptr;
  TTree* metadata_tree_ = nullptr;
  bool metadata_filled_ = false;

  int metadata_schema_version_ = schema_version_;
  int metadata_matcher_version_ = photon_tree::PythiaClusterTruthMatcher::kAlgorithmVersion;
  int metadata_classifier_version_ = photon_tree::HepMCPhotonClassifier::kAlgorithmVersion;
  int shower_shape_algorithm_version_ = ShowerShapeCalculator::kAlgorithmVersion;
  int shower_shape_patch_side_ = ShowerShapeCalculator::kPatchSide;
  unsigned long long n_events_processed_ = 0;
  unsigned long long n_events_written_ = 0;
  unsigned long long n_events_invalid_truth_ = 0;
  unsigned long long n_events_invalid_detector_ = 0;

  unsigned int b_source_file_id_ = 0;
  unsigned int b_event_in_file_ = 0;
  unsigned long long b_event_uid_ = 0;
  unsigned char b_vertex_valid_ = 0;
  int b_vertex_embedding_id_ = -999;
  double b_vertex_x_ = photon_tree::kInvalidDouble;
  double b_vertex_y_ = photon_tree::kInvalidDouble;
  double b_vertex_z_ = photon_tree::kInvalidDouble;
  double b_vertex_t_ = photon_tree::kInvalidDouble;
  unsigned int b_hepmc_n_subevent_ = 0;
  int b_hepmc_signal_event_number_ = -999;
  std::vector<int> b_hepmc_embedding_id_;
};

#endif
