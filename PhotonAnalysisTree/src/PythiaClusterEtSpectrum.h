#ifndef RYOTARO_PYTHIACLUSTERETSPECTRUM_H_20260811
#define RYOTARO_PYTHIACLUSTERETSPECTRUM_H_20260811

#include "Pi0ClusterTruthMatcher.h"
#include "PythiaClusterTruthMatcher.h"

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class TFile;
class TH1D;
class TTree;

class PythiaClusterEtSpectrum : public SubsysReco
{
 public:
  explicit PythiaClusterEtSpectrum(const std::string& name = "PythiaClusterEtSpectrum");
  ~PythiaClusterEtSpectrum() override;

  int Init(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int End(PHCompositeNode* topNode) override;

  void set_output_file_name(const std::string& value) { output_file_name_ = value; }
  void set_manifest_path(const std::string& value) { manifest_path_ = value; }
  void set_manifest_range(long long begin, long long end)
  {
    manifest_begin_ = begin;
    manifest_end_ = end;
  }
  void set_suffix_range(const std::string& first, const std::string& last)
  {
    first_suffix_ = first;
    last_suffix_ = last;
  }
  void set_signal_embedding_id(int value) { signal_embedding_id_ = value; }
  void set_truth_node_name(const std::string& value) { truth_node_name_ = value; }
  void set_hepmc_event_map_node_name(const std::string& value) { hepmc_event_map_node_name_ = value; }
  void set_tower_node_name(const std::string& value) { tower_node_name_ = value; }
  void set_raw_truth_tower_node_name(const std::string& value) { raw_truth_tower_node_name_ = value; }
  void set_truth_cell_node_name(const std::string& value) { truth_cell_node_name_ = value; }
  void set_truth_hit_node_name(const std::string& value) { truth_hit_node_name_ = value; }
  void set_tower_geom_node_name(const std::string& value) { tower_geom_node_name_ = value; }
  void set_split_cluster_node_name(const std::string& value) { split_cluster_node_name_ = value; }
  void set_binning(int n_bins, double et_max)
  {
    n_bins_ = n_bins;
    et_max_ = et_max;
  }
  void set_truth_eta_max(double value) { truth_eta_max_ = value; }
  void set_cluster_eta_max(double value) { cluster_eta_max_ = value; }
  void set_min_cluster_energy(double value) { min_cluster_energy_ = value; }
  void set_dominant_fraction_min(double value) { dominant_fraction_min_ = value; }
  void set_pi0_contributor_fraction_min(double value) { pi0_contributor_fraction_min_ = value; }
  void set_min_energy_contribution_fraction(double value) { min_energy_contribution_fraction_ = value; }
  void set_separated_delta_r_cut(double value) { separated_delta_r_cut_ = value; }
  void set_merged_delta_r_cut(double value) { merged_delta_r_cut_ = value; }
  void set_response_window(double minimum, double maximum)
  {
    response_min_ = minimum;
    response_max_ = maximum;
  }
  void set_verbosity(int value) { verbosity_ = value; }

 private:
  static constexpr int schema_version_ = 2;

  void create_output_directory() const;
  void create_output();
  void close_output();

  std::string output_file_name_ = "pythia_cluster_et_partial.root";
  std::string manifest_path_;
  std::string first_suffix_;
  std::string last_suffix_;
  std::string truth_node_name_ = "G4TruthInfo";
  std::string hepmc_event_map_node_name_ = "PHHepMCGenEventMap";
  std::string tower_node_name_ = "TOWERINFO_CALIB_CEMC";
  std::string raw_truth_tower_node_name_ = "TOWER_SIM_CEMC";
  std::string truth_cell_node_name_ = "G4CELL_CEMC";
  std::string truth_hit_node_name_ = "G4HIT_CEMC";
  std::string tower_geom_node_name_ = "TOWERGEOM_CEMC";
  std::string split_cluster_node_name_ = "CLUSTERINFO_CEMC";
  std::string cluster_collection_ = "split";
  std::string prompt_selection_ = "dominant_prompt_category_1_or_2";
  std::string pi0_selection_ = "g4_pi0_decay_including_secondary_or_generator_photon_with_pi0_origin";
  std::string topology_priority_ = "separated_then_merged_then_missing_then_none";
  std::string projection_scheme_ = "g4_photon_vertex_and_momentum_to_cemc_cylinder";
  std::string energy_topology_priority_ = "separated_then_merged_then_missing_then_none";
  std::string energy_matching_scheme_ = "g4hit_edep_to_direct_pi0_daughter_maximum_deposit";
  std::string energy_candidate_selection_ = "summed_pi0_daughter_primary_contributor_fraction";
  long long manifest_begin_ = -1;
  long long manifest_end_ = -1;
  int signal_embedding_id_ = 1;
  int n_bins_ = 100;
  double et_max_ = 20.0;
  double truth_eta_max_ = 0.7;
  double cluster_eta_max_ = 0.7;
  double min_cluster_energy_ = 0.2;
  double dominant_fraction_min_ = 0.5;
  double pi0_contributor_fraction_min_ = 0.5;
  double min_energy_contribution_fraction_ = 0.0;
  double separated_delta_r_cut_ = 0.03;
  double merged_delta_r_cut_ = 0.06;
  double response_min_ = 0.5;
  double response_max_ = 1.5;
  int verbosity_ = 0;

  int pi0_truth_matching_algorithm_version_ = photon_tree::Pi0ClusterTruthMatcher::kAlgorithmVersion;
  photon_tree::PythiaClusterTruthMatcher truth_matcher_;
  photon_tree::Pi0ClusterTruthMatcher pi0_truth_matcher_;
  TFile* output_file_ = nullptr;
  TH1D* h_prompt_ = nullptr;
  TH1D* h_pi0_ = nullptr;
  TH1D* h_pi0_separated_ = nullptr;
  TH1D* h_pi0_merged_ = nullptr;
  TH1D* h_pi0_missing_ = nullptr;
  TH1D* h_pi0_energy_separated_ = nullptr;
  TH1D* h_pi0_energy_merged_ = nullptr;
  TH1D* h_pi0_energy_missing_ = nullptr;
  TTree* metadata_tree_ = nullptr;

  unsigned long long n_events_processed_ = 0;
  unsigned long long n_events_written_ = 0;
  unsigned long long n_events_invalid_ = 0;
  unsigned long long n_cluster_considered_ = 0;
  unsigned long long n_cluster_invalid_truth_ = 0;
  unsigned long long n_prompt_cluster_ = 0;
  unsigned long long n_pi0_cluster_ = 0;
  unsigned long long n_pi0_cluster_g4_decay_ = 0;
  unsigned long long n_pi0_cluster_generator_decay_ = 0;
  unsigned long long n_pi0_candidate_g4_decay_ = 0;
  unsigned long long n_pi0_candidate_generator_decay_ = 0;
  unsigned long long n_pi0_malformed_daughters_ = 0;
  unsigned long long n_pi0_projection_failure_ = 0;
  unsigned long long n_pi0_separated_ = 0;
  unsigned long long n_pi0_merged_ = 0;
  unsigned long long n_pi0_missing_ = 0;
  unsigned long long n_pi0_none_ = 0;
  unsigned long long n_pi0_ambiguous_ = 0;
  unsigned long long n_pi0_separated_cluster_fill_ = 0;
  unsigned long long n_pi0_merged_cluster_fill_ = 0;
  unsigned long long n_pi0_missing_cluster_fill_ = 0;
  unsigned long long n_pi0_energy_separated_ = 0;
  unsigned long long n_pi0_energy_merged_ = 0;
  unsigned long long n_pi0_energy_missing_ = 0;
  unsigned long long n_pi0_energy_none_ = 0;
  unsigned long long n_pi0_energy_match_invalid_ = 0;
  unsigned long long n_pi0_energy_separated_cluster_fill_ = 0;
  unsigned long long n_pi0_energy_merged_cluster_fill_ = 0;
  unsigned long long n_pi0_energy_missing_cluster_fill_ = 0;
};

#endif
