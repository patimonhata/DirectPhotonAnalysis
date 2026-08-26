#ifndef RYOTARO_TOPOLOGYEVENTDISPLAYDUMP_H_20260824
#define RYOTARO_TOPOLOGYEVENTDISPLAYDUMP_H_20260824

#include <Pi0AnchorTopologyEvaluator.h>

#include <fun4all/SubsysReco.h>

#include <array>
#include <string>

class PHCompositeNode;
class PHG4Particle;
class PHG4TruthInfoContainer;
class TFile;
class TTree;

class TopologyEventDisplayDump : public SubsysReco
{
 public:
  explicit TopologyEventDisplayDump(const std::string& name = "TopologyEventDisplayDump");
  ~TopologyEventDisplayDump() override;

  int Init(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;
  int End(PHCompositeNode*) override;

  void set_output_file_name(const std::string& value) { output_file_name_ = value; }
  void set_source_label(const std::string& value) { source_label_ = value; }
  void set_manifest_path(const std::string& value) { manifest_path_ = value; }
  void set_manifest_range(long long begin, long long end)
  {
    manifest_begin_ = begin;
    manifest_end_ = end;
  }
  void set_sample_mode(photon_tree::Pi0SampleMode value) { config_.sample_mode = value; }
  void set_first_event(int value) { first_event_ = value; }
  void set_signal_embedding_id(int value) { config_.signal_embedding_id = value; }
  void set_truth_node_name(const std::string& value) { config_.truth_node_name = value; }
  void set_hepmc_event_map_node_name(const std::string& value) { config_.hepmc_event_map_node_name = value; }
  void set_tower_node_name(const std::string& value) { config_.tower_node_name = value; }
  void set_raw_truth_tower_node_name(const std::string& value) { config_.raw_truth_tower_node_name = value; }
  void set_truth_cell_node_name(const std::string& value) { config_.truth_cell_node_name = value; }
  void set_truth_hit_node_name(const std::string& value) { config_.truth_hit_node_name = value; }
  void set_cluster_node_name(const std::string& value) { config_.cluster_node_name = value; }
  void set_truth_eta_max(double value) { config_.truth_eta_max = value; }
  void set_anchor_cluster_eta_max(double value) { config_.anchor_cluster_eta_max = value; }
  void set_partner_cluster_eta_max(double value) { config_.partner_cluster_eta_max = value; }
  void set_min_cluster_energy(double value) { config_.min_cluster_energy = value; }
  void set_dominant_fraction_min(double value) { config_.dominant_fraction_min = value; }
  void set_anchor_pi0_fraction_min(double value) { config_.anchor_pi0_fraction_min = value; }
  void set_min_energy_contribution_fraction(double value) { config_.min_energy_contribution_fraction = value; }
  void set_min_photon_energy_recovery(double value) { config_.min_photon_energy_recovery = value; }
  void set_min_direct_match_cluster_energy_coverage(double value) { config_.min_direct_match_cluster_energy_coverage = value; }
  void set_missing_diagnostic_max_delta_r(double value) { config_.missing_diagnostic_max_delta_r = value; }
  void set_write_detail(bool value) { write_detail_ = value; }
  void set_verbosity(int value) { verbosity_ = value; config_.verbosity = value; }

 private:
  static constexpr int schema_version_ = 2;
  static constexpr int invalid_int_ = -999;
  static constexpr double invalid_double_ = -999.0;

  void create_output();
  void create_output_directory() const;
  void close_output();
  void fill_event(PHCompositeNode*, const photon_tree::Pi0AnchorTopologyEventResult&);
  void fill_candidates(const photon_tree::Pi0AnchorTopologyEventResult&);
  void fill_anchors(const photon_tree::Pi0AnchorTopologyEventResult&);
  void fill_candidate_cluster_truth(const photon_tree::Pi0AnchorTopologyEventResult&);
  void fill_clusters(PHCompositeNode*, const photon_tree::Pi0AnchorTopologyEventResult&);
  void fill_truth(PHCompositeNode*, const photon_tree::Pi0AnchorTopologyEventResult&);

  static bool project_to_radius(double x0, double y0, double z0,
                                double px, double py, double pz,
                                double radius, double& x1,
                                double& y1, double& z1);
  static std::array<int, 2> family_for_track(int track_id, PHG4TruthInfoContainer*, const photon_tree::Pi0AnchorTopologyEventResult&);

  std::string output_file_name_ = "topology_event_display.root";
  std::string source_label_;
  std::string manifest_path_;
  long long manifest_begin_ = -1;
  long long manifest_end_ = -1;
  bool write_detail_ = true;
  int verbosity_ = 0;
  int first_event_ = 0;
  photon_tree::Pi0AnchorTopologyConfig config_;
  photon_tree::Pi0AnchorTopologyEvaluator evaluator_;

  TFile* output_file_ = nullptr;
  TTree* metadata_tree_ = nullptr;
  TTree* events_tree_ = nullptr;
  TTree* candidates_tree_ = nullptr;
  TTree* anchors_tree_ = nullptr;
  TTree* candidate_cluster_truth_tree_ = nullptr;
  TTree* clusters_tree_ = nullptr;
  TTree* cluster_contributors_tree_ = nullptr;
  TTree* cluster_towers_tree_ = nullptr;
  TTree* truth_particles_tree_ = nullptr;
  TTree* truth_segments_tree_ = nullptr;

  unsigned long long events_processed_ = 0;
  unsigned long long source_events_seen_ = 0;
  unsigned long long events_written_ = 0;
  unsigned long long events_invalid_ = 0;

  int b_event_ = 0;
  int b_candidate_id_ = invalid_int_;
  int b_pathway_ = invalid_int_;
  std::string b_pathway_name_;
  int b_parent_barcode_ = invalid_int_;
  int b_g4_parent_track_id_ = invalid_int_;
  double b_energy_ = invalid_double_;
  double b_et_ = invalid_double_;
  double b_pt_ = invalid_double_;
  double b_eta_ = invalid_double_;
  double b_phi_ = invalid_double_;
  double b_x_ = invalid_double_;
  double b_y_ = invalid_double_;
  double b_z_ = invalid_double_;
  double b_r_ = invalid_double_;
  int b_photon0_track_id_ = invalid_int_;
  int b_photon1_track_id_ = invalid_int_;
  int metadata_sample_mode_ = 0;
  int metadata_write_detail_ = 1;
  double b_photon0_energy_ = invalid_double_;
  double b_photon1_energy_ = invalid_double_;
  double b_photon0_eta_ = invalid_double_;
  double b_photon1_eta_ = invalid_double_;
  double b_photon0_phi_ = invalid_double_;
  double b_photon1_phi_ = invalid_double_;
  int b_best_cluster0_id_ = invalid_int_;
  int b_best_cluster1_id_ = invalid_int_;
  double b_maximum_edep0_ = invalid_double_;
  double b_maximum_edep1_ = invalid_double_;
  double b_reconstructed_photon0_energy_ = invalid_double_;
  double b_reconstructed_photon1_energy_ = invalid_double_;
  int b_recovered0_ = 0;
  int b_recovered1_ = 0;
  int b_topology_evaluated_ = 0;

  int b_anchor_id_ = invalid_int_;
  unsigned int b_cluster_id_ = 0;
  int b_topology_ = 0;
  std::string b_topology_name_;
  int b_reason_ = 0;
  std::string b_reason_name_;
  int b_missing_detail_ = 0;
  std::string b_missing_detail_name_;
  int b_partner_photon_index_ = invalid_int_;
  double b_main_fraction_ = invalid_double_;
  double b_second_fraction_ = invalid_double_;
  double b_unmatched_max_fraction_ = invalid_double_;
  int b_ambiguous_main_ = 0;

  int b_match_valid_ = 0;
  double b_total_edep_ = 0.0;
  int b_match_usable_ = 0;
  int b_match_status_ = 0;
  std::string b_match_status_name_;
  int b_match_failure_ = 0;
  std::string b_match_failure_name_;
  int b_match_failure_ieta_ = invalid_int_;
  int b_match_failure_iphi_ = invalid_int_;
  unsigned int b_match_tower_count_ = 0;
  unsigned int b_match_matched_tower_count_ = 0;
  double b_match_cluster_member_energy_coverage_ = 0.0;
  double b_gamma0_edep_ = 0.0;
  double b_gamma1_edep_ = 0.0;
  double b_other_edep_ = 0.0;
  double b_gamma0_fraction_ = 0.0;
  double b_gamma1_fraction_ = 0.0;
  double b_other_fraction_ = 0.0;
  double b_gamma0_recovery_estimate_ = 0.0;
  double b_gamma1_recovery_estimate_ = 0.0;

  int b_topology_considered_ = 0;
  int b_partner_diagnostic_found_ = 0;
  int b_partner_diagnostic_below_energy_threshold_ = 0;
  int b_partner_diagnostic_has_direct_deposit_ = 0;
  int b_partner_diagnostic_cluster_id_ = invalid_int_;
  double b_partner_diagnostic_cluster_energy_ = invalid_double_;
  double b_partner_diagnostic_cluster_eta_ = invalid_double_;
  double b_partner_diagnostic_cluster_phi_ = invalid_double_;
  double b_partner_diagnostic_delta_r_ = invalid_double_;
  double b_partner_diagnostic_reconstructed_energy_ = 0.0;
  double b_partner_diagnostic_recovery_ = 0.0;
  int b_partner_diagnostic_match_usable_ = 0;
  int b_partner_diagnostic_match_status_ = 0;
  std::string b_partner_diagnostic_match_status_name_;
  int b_partner_diagnostic_match_failure_ = 0;
  std::string b_partner_diagnostic_match_failure_name_;
  int b_partner_diagnostic_failure_ieta_ = invalid_int_;
  int b_partner_diagnostic_failure_iphi_ = invalid_int_;
  double b_partner_diagnostic_match_coverage_ = 0.0;
  int b_topology_cluster_index_ = invalid_int_;
  int b_anchor_acceptance_ = 0;
  int b_ntowers_ = 0;
  int b_truth_valid_ = 0;
  double b_truth_total_edep_ = 0.0;
  int b_n_contributors_ = 0;
  int b_contributor_index_ = invalid_int_;
  int b_g4_track_id_ = invalid_int_;
  int b_g4_pdg_id_ = invalid_int_;
  int b_embedding_id_ = invalid_int_;
  int b_hepmc_barcode_ = invalid_int_;
  double b_contributor_edep_ = 0.0;
  double b_contributor_fraction_ = 0.0;

  unsigned int b_tower_key_ = 0;
  int b_ieta_ = invalid_int_;
  int b_iphi_ = invalid_int_;
  double b_cluster_tower_energy_ = invalid_double_;
  double b_tower_energy_ = invalid_double_;
  double b_allocation_fraction_ = invalid_double_;

  int b_track_id_ = invalid_int_;
  int b_pid_ = invalid_int_;
  int b_parent_id_ = invalid_int_;
  int b_primary_id_ = invalid_int_;
  int b_is_primary_ = 0;
  int b_is_g4_secondary_pi0_ = 0;
  int b_family_candidate_id_ = invalid_int_;
  int b_family_gamma_index_ = invalid_int_;
  double b_px_ = invalid_double_;
  double b_py_ = invalid_double_;
  double b_pz_ = invalid_double_;
  double b_vx_ = invalid_double_;
  double b_vy_ = invalid_double_;
  double b_vz_ = invalid_double_;
  double b_x0_ = invalid_double_;
  double b_y0_ = invalid_double_;
  double b_z0_ = invalid_double_;
  double b_x1_ = invalid_double_;
  double b_y1_ = invalid_double_;
  double b_z1_ = invalid_double_;

  int b_n_candidates_ = 0;
  int b_n_g4_primary_pi0_ = 0;
  int b_n_generator_pi0_ = 0;
  int b_n_g4_secondary_pi0_ = 0;
  int b_n_selected_family_particles_ = 0;
  int b_n_truth_particles_ = 0;
  int b_n_clusters_ = 0;
  int b_n_topology_clusters_ = 0;
  int b_n_anchors_ = 0;
  int b_n_separated_ = 0;
  int b_n_merged_ = 0;
  int b_n_missing_ = 0;
  int b_n_other_ = 0;
  double b_collision_x_ = 0.0;
  double b_collision_y_ = 0.0;
  double b_collision_z_ = 0.0;
};

#endif
