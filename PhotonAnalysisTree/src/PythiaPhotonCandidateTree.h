#ifndef RYOTARO_PYTHIAPHOTONCANDIDATETREE_H_20260901
#define RYOTARO_PYTHIAPHOTONCANDIDATETREE_H_20260901

#include "HepMCPhotonClassifier.h"
#include "PhotonTreeCommon.h"
#include "Pi0AnchorTopologyEvaluator.h"
#include "PythiaClusterTruthMatcher.h"

#include <fun4all/SubsysReco.h>

#include <memory>
#include <string>
#include <vector>

class PHCompositeNode;
class PHHepMCGenEventMap;
class RawClusterContainer;
class Fun4AllInputManager;
class TFile;
class TTree;

namespace TMVA::Experimental
{
class RBDT;
}

namespace photon_tree
{
struct PhotonCandidateSelectionBranches
{
  std::vector<float> bdt_score;
  std::vector<unsigned char> bdt_valid;
  std::vector<float> bdt_tight_boundary;
  std::vector<float> bdt_nontight_lower;
  std::vector<float> bdt_nontight_upper;
  std::vector<unsigned char> pass_kinematics;
  std::vector<unsigned char> pass_preselection;
  std::vector<unsigned char> pass_tight;
  std::vector<unsigned char> pass_nontight;

  std::vector<double> iso_raw_et;
  std::vector<double> iso_corrected_et;
  std::vector<double> iso_boundary;
  std::vector<double> noniso_boundary;
  std::vector<unsigned int> iso_topocluster_count;
  std::vector<unsigned char> pass_isolated;
  std::vector<unsigned char> pass_nonisolated;

  std::vector<unsigned char> pass_region_a;
  std::vector<unsigned char> pass_region_b;
  std::vector<unsigned char> pass_region_c;
  std::vector<unsigned char> pass_region_d;
  std::vector<unsigned char> pass_final_photon;

  std::vector<unsigned char> pi0_tag;
  std::vector<unsigned char> eta_tag;
  std::vector<int> pi0_partner_cluster_id;
  std::vector<int> eta_partner_cluster_id;
  std::vector<float> pi0_partner_mass;
  std::vector<float> eta_partner_mass;

  std::vector<unsigned char> truth_prompt_cluster;
  std::vector<unsigned char> pi0_anchor_valid;
  std::vector<int> pi0_anchor_candidate_index;
  std::vector<float> pi0_anchor_main_fraction;
  std::vector<float> pi0_anchor_second_fraction;
  std::vector<float> pi0_anchor_unmatched_max_fraction;
  std::vector<unsigned char> pi0_anchor_ambiguous_main;
  std::vector<int> pi0_anchor_topology;
  std::vector<int> pi0_anchor_reason;
  std::vector<int> pi0_anchor_missing_category;
  std::vector<int> pi0_anchor_missing_detail;
  std::vector<int> pi0_anchor_partner_photon_index;
  std::vector<float> pi0_anchor_partner_diagnostic_mass;

  void clear();
  void create_branches(TTree* tree);
};

struct Pi0TopologyTreeBranches
{
  unsigned int candidate_count = 0;
  std::vector<int> candidate_pathway;
  std::vector<int> candidate_parent_barcode;
  std::vector<int> candidate_g4_parent_track_id;
  std::vector<float> candidate_e;
  std::vector<float> candidate_pt;
  std::vector<float> candidate_eta;
  std::vector<float> candidate_phi;
  std::vector<unsigned char> candidate_evaluated;

  std::vector<int> photon0_track_id;
  std::vector<int> photon1_track_id;
  std::vector<float> photon0_e;
  std::vector<float> photon1_e;
  std::vector<float> photon0_eta;
  std::vector<float> photon1_eta;
  std::vector<float> photon0_phi;
  std::vector<float> photon1_phi;
  std::vector<unsigned char> photon0_projection_valid;
  std::vector<unsigned char> photon1_projection_valid;
  std::vector<float> photon0_projection_eta;
  std::vector<float> photon1_projection_eta;
  std::vector<float> photon0_projection_phi;
  std::vector<float> photon1_projection_phi;
  std::vector<unsigned char> photon0_in_cemc_acceptance;
  std::vector<unsigned char> photon1_in_cemc_acceptance;
  std::vector<unsigned char> photon0_pre_cemc_interaction;
  std::vector<unsigned char> photon1_pre_cemc_interaction;
  std::vector<float> photon0_first_daughter_radius;
  std::vector<float> photon1_first_daughter_radius;
  std::vector<float> photon0_cemc_edep;
  std::vector<float> photon1_cemc_edep;
  std::vector<int> photon0_best_cluster_id;
  std::vector<int> photon1_best_cluster_id;
  std::vector<float> photon0_maximum_edep;
  std::vector<float> photon1_maximum_edep;
  std::vector<float> photon0_reconstructed_e;
  std::vector<float> photon1_reconstructed_e;
  std::vector<unsigned char> photon0_recovered;
  std::vector<unsigned char> photon1_recovered;

  std::vector<unsigned char> photon0_diagnostic_found;
  std::vector<unsigned char> photon1_diagnostic_found;
  std::vector<int> photon0_diagnostic_cluster_id;
  std::vector<int> photon1_diagnostic_cluster_id;
  std::vector<float> photon0_diagnostic_cluster_e;
  std::vector<float> photon1_diagnostic_cluster_e;
  std::vector<float> photon0_diagnostic_delta_r;
  std::vector<float> photon1_diagnostic_delta_r;
  std::vector<float> photon0_diagnostic_recovery;
  std::vector<float> photon1_diagnostic_recovery;
  std::vector<unsigned char> photon0_diagnostic_below_threshold;
  std::vector<unsigned char> photon1_diagnostic_below_threshold;
  std::vector<unsigned char> photon0_diagnostic_direct_deposit;
  std::vector<unsigned char> photon1_diagnostic_direct_deposit;

  unsigned int anchor_count = 0;
  std::vector<unsigned int> anchor_cluster_id;
  std::vector<unsigned int> anchor_candidate_index;
  std::vector<float> anchor_main_fraction;
  std::vector<float> anchor_second_fraction;
  std::vector<float> anchor_unmatched_max_fraction;
  std::vector<unsigned char> anchor_ambiguous_main;
  std::vector<int> anchor_topology;
  std::vector<int> anchor_reason;
  std::vector<int> anchor_missing_category;
  std::vector<int> anchor_missing_detail;
  std::vector<int> anchor_partner_photon_index;
  std::vector<int> anchor_pre_cemc_photon_index;
  std::vector<float> anchor_partner_diagnostic_mass;

  void clear();
  void fill(const Pi0AnchorTopologyEventResult& result);
  void create_branches(TTree* tree);
};
}

class PythiaPhotonCandidateTree : public SubsysReco
{
 public:
  explicit PythiaPhotonCandidateTree(const std::string& name = "PythiaPhotonCandidateTree");
  ~PythiaPhotonCandidateTree() override;

  int Init(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int End(PHCompositeNode* topNode) override;

  void set_input_file_name(const std::string& value) { input_file_name_ = value; }
  void set_output_file_name(const std::string& value) { output_file_name_ = value; }
  void set_model_file_name(const std::string& value) { model_file_name_ = value; }
  void set_sample_name(const std::string& value) { sample_name_ = value; }
  void set_source_file_id(unsigned int value) { source_file_id_ = value; }
  void set_primary_input_manager(Fun4AllInputManager* value) { primary_input_manager_ = value; }
  void set_manifest_path(const std::string& value) { manifest_path_ = value; }
  void set_manifest_range(long long begin, long long end) { manifest_begin_ = begin; manifest_end_ = end; }
  void set_suffix_range(const std::string& first, const std::string& last) { first_input_suffix_ = first; last_input_suffix_ = last; }
  void set_min_cluster_energy(double value) { min_cluster_energy_ = value; }
  void set_map_chunk_id(unsigned int value) { map_chunk_id_ = value; }
  void set_signal_embedding_id(int value) { signal_embedding_id_ = value; }
  void set_truth_jet_node_name(const std::string& value) { truth_jet_node_name_ = value; }
  void set_verbosity(int value) { verbosity_ = value; }

 private:
  static constexpr int schema_version_ = 3;

  bool configure_sample();
  bool fill_event_truth(const PHHepMCGenEventMap* event_map, PHCompositeNode* topNode);
  bool fill_candidate_selection(RawClusterContainer* split_clusters, RawClusterContainer* topo_clusters);
  bool update_input_provenance();
  void fill_topology_cluster_links(const photon_tree::Pi0AnchorTopologyEventResult& result);
  void create_output();
  void close_output();
  void reset_event();

  std::string input_file_name_;
  std::string manifest_path_;
  std::string first_input_suffix_;
  std::string last_input_suffix_;
  std::string current_input_file_;
  std::string output_file_name_ = "pythia_photon_candidate_tree.root";
  std::string model_file_name_;
  std::string sample_name_;
  std::string model_key_ = "myBDT";
  std::string model_sha256_ = "5494946f40213a1035dc9e883bc71ddfdcebf86015d279300c8f44774e846248";
  std::string analysis_release_ = "ana.565";
  std::string split_cluster_node_name_ = "CLUSTERINFO_CEMC";
  std::string topo_cluster_node_name_ = "TOPOCLUSTER_ALLCALO";
  std::string tower_node_name_ = "TOWERINFO_CALIB_CEMC";
  std::string tower_geom_node_name_ = "TOWERGEOM_CEMC";
  std::string hepmc_event_map_node_name_ = "PHHepMCGenEventMap";
  std::string truth_jet_node_name_ = "AntiKt_Truth_r04";
  std::string truth_node_name_ = "G4TruthInfo";
  std::string raw_truth_tower_node_name_ = "TOWER_SIM_CEMC";
  std::string truth_cell_node_name_ = "G4CELL_CEMC";
  std::string truth_hit_node_name_ = "G4HIT_CEMC";
  std::string isolation_definition_ = "sum_topocluster_et_delta_r_lt_0p4_minus_candidate_et";
  std::string event_weight_definition_ = "cross_section_pb_times_generator_weight_divided_by_full_sample_sum_generator_weight_at_merge";
  std::string photon_stitch_definition_ = "leading_status1_terminal_prompt_hepmc_photon_pt";
  std::string jet_stitch_definition_ = "leading_AntiKt_Truth_r04_jet_pt";

  unsigned int source_file_id_ = 0;
  unsigned int map_chunk_id_ = 0;
  long long manifest_begin_ = -1;
  long long manifest_end_ = -1;
  long long input_file_count_ = 0;
  unsigned int next_event_in_file_ = 0;
  int signal_embedding_id_ = 1;
  int verbosity_ = 0;
  double sample_cross_section_pb_ = 0.0;
  double sample_window_min_ = 0.0;
  double sample_window_max_ = 0.0;
  bool sample_is_photonjet_ = false;
  bool sample_upper_unbounded_ = false;

  double min_cluster_energy_ = 0.1;
  double shower_shape_min_tower_energy_ = 0.070;
  double candidate_et_min_ = 5.0;
  double candidate_et_max_ = 35.0;
  std::string topocluster_configuration_ = "ALLCALO_EMCal_HCal_noise_0p006_0p03_0p09_significance_4_2_1_corner_split_localmax_1_2_0p5_Rshower_0p025";
  double candidate_abs_eta_max_ = 0.7;
  double max_abs_vertex_z_ = 60.0;
  double isolation_radius_ = 0.4;
  double isolation_scale_ = 1.2;
  double isolation_offset_ = 0.1;
  double nonisolation_gap_ = 0.8;
  double meson_partner_min_energy_ = 0.5;
  double pi0_mass_min_ = 0.10;
  double pi0_mass_max_ = 0.20;
  double eta_mass_min_ = 0.45;
  double eta_mass_max_ = 0.65;

  photon_tree::PhotonTreeCommon common_;
  photon_tree::Pi0AnchorTopologyEvaluator topology_evaluator_;
  photon_tree::ClusterTruthCollection cluster_truth_;
  photon_tree::PhotonCandidateSelectionBranches selection_;
  photon_tree::Pi0TopologyTreeBranches topology_;
  std::unique_ptr<TMVA::Experimental::RBDT> bdt_;

  TFile* output_file_ = nullptr;
  TTree* event_tree_ = nullptr;
  TTree* metadata_tree_ = nullptr;
  bool metadata_filled_ = false;
  int metadata_schema_version_ = schema_version_;

  unsigned long long n_events_processed_ = 0;
  unsigned long long n_events_written_ = 0;
  unsigned long long n_events_vertex_rejected_ = 0;
  unsigned long long n_events_invalid_ = 0;
  unsigned long long n_events_stitch_pass_ = 0;
  unsigned long long n_clusters_region_a_ = 0;
  unsigned long long n_clusters_region_b_ = 0;
  unsigned long long n_clusters_region_c_ = 0;
  unsigned long long n_clusters_region_d_ = 0;
  unsigned long long n_clusters_final_photon_ = 0;
  double sum_generator_weight_processed_ = 0.0;
  double sum_generator_weight_stitch_pass_ = 0.0;

  unsigned int b_source_file_id_ = 0;
  unsigned int b_event_in_file_ = 0;
  unsigned long long b_event_uid_ = 0;
  int b_hepmc_event_number_ = -999;
  double b_vertex_x_ = photon_tree::kInvalidDouble;
  double b_vertex_y_ = photon_tree::kInvalidDouble;
  double b_vertex_z_ = photon_tree::kInvalidDouble;
  double b_vertex_t_ = photon_tree::kInvalidDouble;
  unsigned char b_event_weight_valid_ = 0U;
  double b_generator_weight_ = 1.0;
  photon_tree::HepMCPhotonClassifier photon_classifier_;
  double b_weight_numerator_pb_ = 0.0;
  double b_leading_truth_photon_pt_ = photon_tree::kInvalidDouble;
  double b_leading_truth_jet_pt_ = photon_tree::kInvalidDouble;
  unsigned char b_sample_stitching_valid_ = 0U;
  unsigned char b_sample_stitching_pass_ = 0U;
  unsigned int b_topocluster_count_ = 0U;
  unsigned int b_region_a_count_ = 0U;
  unsigned int b_region_b_count_ = 0U;
  unsigned int b_region_c_count_ = 0U;
  unsigned int b_region_d_count_ = 0U;
  unsigned int b_final_photon_count_ = 0U;
  Fun4AllInputManager* primary_input_manager_ = nullptr;
};

#endif
