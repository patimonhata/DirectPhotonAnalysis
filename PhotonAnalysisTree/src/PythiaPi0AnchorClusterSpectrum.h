#ifndef RYOTARO_PYTHIAPI0ANCHORCLUSTERSPECTRUM_H_20260814
#define RYOTARO_PYTHIAPI0ANCHORCLUSTERSPECTRUM_H_20260814

#include "Pi0ClusterTruthMatcher.h"
#include "PythiaClusterTruthMatcher.h"

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class TFile;
class TH1D;
class TTree;

class PythiaPi0AnchorClusterSpectrum : public SubsysReco
{
 public:
  explicit PythiaPi0AnchorClusterSpectrum(
      const std::string& name = "PythiaPi0AnchorClusterSpectrum");
  ~PythiaPi0AnchorClusterSpectrum() override;

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
  void set_hepmc_event_map_node_name(const std::string& value)
  {
    hepmc_event_map_node_name_ = value;
  }
  void set_tower_node_name(const std::string& value) { tower_node_name_ = value; }
  void set_raw_truth_tower_node_name(const std::string& value)
  {
    raw_truth_tower_node_name_ = value;
  }
  void set_truth_cell_node_name(const std::string& value)
  {
    truth_cell_node_name_ = value;
  }
  void set_truth_hit_node_name(const std::string& value)
  {
    truth_hit_node_name_ = value;
  }
  void set_split_cluster_node_name(const std::string& value)
  {
    split_cluster_node_name_ = value;
  }
  void set_binning(int bins, double et_max)
  {
    n_bins_ = bins;
    et_max_ = et_max;
  }
  void set_truth_eta_max(double value) { truth_eta_max_ = value; }
  void set_anchor_cluster_eta_max(double value)
  {
    anchor_cluster_eta_max_ = value;
  }
  // A non-positive value disables the software eta cut for partner lookup.
  void set_partner_cluster_eta_max(double value)
  {
    partner_cluster_eta_max_ = value;
  }
  void set_min_cluster_energy(double value) { min_cluster_energy_ = value; }
  void set_dominant_fraction_min(double value)
  {
    dominant_fraction_min_ = value;
  }
  void set_anchor_pi0_fraction_min(double value)
  {
    anchor_pi0_fraction_min_ = value;
  }
  void set_min_energy_contribution_fraction(double value)
  {
    min_energy_contribution_fraction_ = value;
  }
  void set_verbosity(int value) { verbosity_ = value; }

 private:
  static constexpr int schema_version_ = 1;

  void create_output_directory() const;
  void create_output();
  void close_output();

  std::string output_file_name_ = "pythia_pi0_anchor_cluster_partial.root";
  std::string manifest_path_;
  std::string first_suffix_;
  std::string last_suffix_;
  std::string truth_node_name_ = "G4TruthInfo";
  std::string hepmc_event_map_node_name_ = "PHHepMCGenEventMap";
  std::string tower_node_name_ = "TOWERINFO_CALIB_CEMC";
  std::string raw_truth_tower_node_name_ = "TOWER_SIM_CEMC";
  std::string truth_cell_node_name_ = "G4CELL_CEMC";
  std::string truth_hit_node_name_ = "G4HIT_CEMC";
  std::string split_cluster_node_name_ = "CLUSTERINFO_CEMC";
  std::string cluster_collection_ = "split";
  std::string classification_unit_ =
      "every_cluster_with_selected_pi0_as_grouped_main_contributor";
  std::string pi0_selection_ =
      "signal_g4_primary_pi0_or_generator_pi0_with_exactly_two_g4_photons";
  std::string partner_selection_ =
      "same_energy_cut_as_anchor_partner_eta_cut_configurable";
  std::string topology_definition_ =
      "anchor_membership_in_direct_daughter_maximum_deposit_clusters";
  std::string topology_priority_ =
      "ambiguous_main_to_other_then_merged_then_separated_then_missing_then_other";
  std::string response_policy_ = "not_used_for_classification";
  long long manifest_begin_ = -1;
  long long manifest_end_ = -1;
  int signal_embedding_id_ = 1;
  int n_bins_ = 100;
  double et_max_ = 20.0;
  double truth_eta_max_ = 0.7;
  double anchor_cluster_eta_max_ = 0.7;
  double partner_cluster_eta_max_ = -1.0;
  double min_cluster_energy_ = 0.2;
  double dominant_fraction_min_ = 0.5;
  double anchor_pi0_fraction_min_ = 0.5;
  double min_energy_contribution_fraction_ = 0.3;
  int verbosity_ = 0;

  int pi0_truth_matching_algorithm_version_ =
      photon_tree::Pi0ClusterTruthMatcher::kAlgorithmVersion;
  photon_tree::PythiaClusterTruthMatcher truth_matcher_;
  photon_tree::Pi0ClusterTruthMatcher pi0_truth_matcher_;
  TFile* output_file_ = nullptr;
  TH1D* h_prompt_ = nullptr;
  TH1D* h_anchor_ = nullptr;
  TH1D* h_separated_ = nullptr;
  TH1D* h_merged_ = nullptr;
  TH1D* h_missing_ = nullptr;
  TH1D* h_other_ = nullptr;
  TTree* metadata_tree_ = nullptr;

  unsigned long long n_events_processed_ = 0;
  unsigned long long n_events_written_ = 0;
  unsigned long long n_events_invalid_ = 0;
  unsigned long long n_cluster_considered_ = 0;
  unsigned long long n_cluster_invalid_truth_ = 0;
  unsigned long long n_prompt_cluster_ = 0;
  unsigned long long n_pi0_candidate_g4_decay_ = 0;
  unsigned long long n_pi0_candidate_generator_decay_ = 0;
  unsigned long long n_pi0_malformed_daughters_ = 0;
  unsigned long long n_anchor_cluster_ = 0;
  unsigned long long n_anchor_g4_decay_ = 0;
  unsigned long long n_anchor_generator_decay_ = 0;
  unsigned long long n_anchor_ambiguous_main_ = 0;
  unsigned long long n_energy_match_invalid_ = 0;
  unsigned long long n_separated_ = 0;
  unsigned long long n_merged_ = 0;
  unsigned long long n_missing_ = 0;
  unsigned long long n_other_ = 0;
};

#endif
