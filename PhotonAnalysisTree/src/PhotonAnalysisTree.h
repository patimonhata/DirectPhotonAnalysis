#ifndef RYOTARO_PHOTONANALYSISTREE_H_20260721
#define RYOTARO_PHOTONANALYSISTREE_H_20260721

#include <fun4all/SubsysReco.h>

#include "PhotonTreeCommon.h"

#include <string>
#include <vector>

class PHCompositeNode;
class PHG4Particle;
class PHG4TruthInfoContainer;
class RawTowerGeomContainer;
class TFile;
class TTree;

class PhotonAnalysisTree : public SubsysReco
{
 public:
  explicit PhotonAnalysisTree(const std::string& name = "PhotonAnalysisTree");
  ~PhotonAnalysisTree() override;

  int Init(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int End(PHCompositeNode* topNode) override;

  void set_input_file_name(const std::string& value) { input_file_name_ = value; }
  void set_output_file_name(const std::string& value) { output_file_name_ = value; }
  void set_source_file_id(unsigned int value) { source_file_id_ = value; }
  void set_expected_primary_pdg(int value) { expected_primary_pdg_ = value; }
  void set_truth_node_name(const std::string& value) { truth_node_name_ = value; }
  void set_tower_node_name(const std::string& value) { tower_node_name_ = value; }
  void set_tower_geom_node_name(const std::string& value) { tower_geom_node_name_ = value; }
  void set_split_cluster_node_name(const std::string& value) { split_cluster_node_name_ = value; }
  void set_nosplit_cluster_node_name(const std::string& value) { nosplit_cluster_node_name_ = value; }
  void set_require_truth_node(bool value) { require_truth_node_ = value; }
  void set_require_nosplit_cluster_node(bool value) { require_nosplit_cluster_node_ = value; }
  void set_acceptance_eta_max(double value) { acceptance_eta_max_ = value; }
  void set_min_cluster_energy(double value) { min_cluster_energy_ = value; }
  void set_shower_shape_min_tower_energy(double value) { shower_shape_min_tower_energy_ = value; }
  void set_store_shower_shape_tower_patch(bool value) { store_shower_shape_tower_patch_ = value; }
  void set_verbosity(int value) { verbosity_ = value; }

 private:
  static constexpr int schema_version_ = 3;
  static constexpr double invalid_double_ = photon_tree::kInvalidDouble;
  static constexpr int invalid_int_ = photon_tree::kInvalidInt;
  static constexpr double default_cemc_radius_ = 95.0;

  void create_output_directory() const;
  void create_trees();
  void close_output_file();
  void reset_event();
  bool fill_truth(PHG4TruthInfoContainer* truth, RawTowerGeomContainer* geometry);

  static double radius(double x, double y);
  static double eta_from_xyz(double x, double y, double z);
  static double phi_from_xy(double x, double y);
  static double wrap_delta_phi(double value);
  static bool project_to_radius(double x0, double y0, double z0,
                                double px, double py, double pz,
                                double target_radius,
                                double& x1, double& y1, double& z1);
  static double invariant_mass(const PHG4Particle* first, const PHG4Particle* second);
  static double energy_asymmetry(const PHG4Particle* first, const PHG4Particle* second);
  double cemc_radius(RawTowerGeomContainer* geometry) const;

  std::string input_file_name_;
  std::string output_file_name_ = "photon_analysis_tree.root";
  std::string truth_node_name_ = "G4TruthInfo";
  std::string tower_node_name_ = "TOWERINFO_CALIB_CEMC";
  std::string tower_geom_node_name_ = "TOWERGEOM_CEMC";
  std::string split_cluster_node_name_ = "CLUSTERINFO_CEMC";
  std::string nosplit_cluster_node_name_ = "CLUSTERINFO_CEMC_NO_SPLIT";
  std::string cluster_ordering_ = "energy_descending_then_cluster_id_ascending";
  bool require_truth_node_ = true;
  bool require_nosplit_cluster_node_ = true;
  unsigned int source_file_id_ = 0;
  int expected_primary_pdg_ = 111;
  double acceptance_eta_max_ = 1.1;
  double min_cluster_energy_ = 0.0;
  double shower_shape_min_tower_energy_ = 0.070;
  bool store_shower_shape_tower_patch_ = true;
  int verbosity_ = 0;
  int shower_shape_algorithm_version_ = ShowerShapeCalculator::kAlgorithmVersion;
  int shower_shape_patch_side_ = ShowerShapeCalculator::kPatchSide;
  photon_tree::PhotonTreeCommon common_;

  TFile* output_file_ = nullptr;
  TTree* event_tree_ = nullptr;
  TTree* metadata_tree_ = nullptr;
  bool metadata_filled_ = false;

  int metadata_schema_version_ = schema_version_;
  unsigned long long n_events_processed_ = 0;
  unsigned long long n_events_written_ = 0;
  unsigned long long n_events_invalid_truth_ = 0;
  unsigned long long n_events_invalid_detector_ = 0;

  unsigned int b_source_file_id_ = 0;
  unsigned int b_event_in_file_ = 0;
  unsigned long long b_event_uid_ = 0;
  double b_vertex_x_ = 0.0;
  double b_vertex_y_ = 0.0;
  double b_vertex_z_ = 0.0;
  int b_label_ = invalid_int_;
  unsigned char b_truth_valid_ = 0;
  int b_truth_primary_pdg_id_ = invalid_int_;
  int b_truth_primary_track_id_ = invalid_int_;
  double b_truth_e_ = invalid_double_;
  double b_truth_px_ = invalid_double_;
  double b_truth_py_ = invalid_double_;
  double b_truth_pz_ = invalid_double_;
  double b_truth_pt_ = invalid_double_;
  double b_truth_eta_ = invalid_double_;
  double b_truth_phi_ = invalid_double_;
  double b_truth_vx_ = invalid_double_;
  double b_truth_vy_ = invalid_double_;
  double b_truth_vz_ = invalid_double_;
  unsigned int b_truth_n_direct_daughter_ = 0;
  unsigned char b_truth_is_pi0_to_2gamma_ = 0;
  std::vector<int> b_truth_daughter_track_id_;
  std::vector<int> b_truth_daughter_pdg_id_;
  std::vector<double> b_truth_daughter_e_;
  std::vector<double> b_truth_daughter_px_;
  std::vector<double> b_truth_daughter_py_;
  std::vector<double> b_truth_daughter_pz_;
  std::vector<double> b_truth_daughter_pt_;
  std::vector<double> b_truth_daughter_eta_;
  std::vector<double> b_truth_daughter_phi_;
  std::vector<double> b_truth_daughter_projection_eta_;
  std::vector<double> b_truth_daughter_projection_phi_;
  std::vector<unsigned char> b_truth_daughter_projection_valid_;
  std::vector<unsigned char> b_truth_daughter_in_acceptance_;
  unsigned char b_truth_both_gamma_in_acceptance_ = 0;
  unsigned char b_truth_at_least_one_gamma_out_acceptance_ = 0;
  unsigned char b_truth_missing_gamma_projection_ = 0;
  double b_truth_m_gg_ = invalid_double_;
  double b_truth_pair_e_asym_ = invalid_double_;
};

#endif
