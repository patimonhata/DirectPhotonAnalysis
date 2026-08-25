#ifndef RYOTARO_PYTHIAPI0CLUSTERMULTIPLICITY_H_20260813
#define RYOTARO_PYTHIAPI0CLUSTERMULTIPLICITY_H_20260813

#include "PythiaClusterTruthMatcher.h"

#include <fun4all/SubsysReco.h>

#include <array>
#include <string>

class PHCompositeNode;
class TFile;
class TH1D;
class TH2D;
class TTree;

class PythiaPi0ClusterMultiplicity : public SubsysReco
{
 public:
  explicit PythiaPi0ClusterMultiplicity(const std::string& name = "PythiaPi0ClusterMultiplicity");
  ~PythiaPi0ClusterMultiplicity() override;

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
  void set_split_cluster_node_name(const std::string& value)
  {
    split_cluster_node_name_ = value;
  }
  void set_truth_eta_max(double value) { truth_eta_max_ = value; }
  void set_pt_binning(int bins, double maximum)
  {
    pt_bins_ = bins;
    pt_max_ = maximum;
  }
  void set_multiplicity_max(int value) { multiplicity_max_ = value; }
  void set_cluster_energy_binning(int bins, double maximum)
  {
    cluster_energy_bins_ = bins;
    cluster_energy_max_ = maximum;
  }
  void set_verbosity(int value) { verbosity_ = value; }

 private:
  static constexpr int schema_version_ = 3;
  static constexpr std::size_t threshold_count_ = 4U;
  static constexpr std::size_t pathway_count_ = 2U;

  void create_output_directory() const;
  void create_output();
  void close_output();

  std::string output_file_name_ = "pythia_pi0_cluster_multiplicity_partial.root";
  std::string manifest_path_;
  std::string first_suffix_;
  std::string last_suffix_;
  std::string truth_node_name_ = "G4TruthInfo";
  std::string hepmc_event_map_node_name_ = "PHHepMCGenEventMap";
  std::string tower_node_name_ = "TOWERINFO_CALIB_CEMC";
  std::string raw_truth_tower_node_name_ = "TOWER_CALIB_CEMC";
  std::string split_cluster_node_name_ = "CLUSTERINFO_CEMC";
  std::string cluster_collection_ = "split";
  std::string pi0_selection_ = "g4_primary_pi0_decay_or_generator_photon_pair_with_hepmc_pi0_origin";
  std::string cluster_selection_ = "finite_cluster_kinematics_without_eta_or_energy_threshold";
  std::string fraction_definition_ = "sum_of_ancestry_compatible_primary_shower_edep_over_all_cluster_primary_shower_edep";
  std::string zero_threshold_definition_ = "strictly_positive_fraction";
  long long manifest_begin_ = -1;
  long long manifest_end_ = -1;
  int signal_embedding_id_ = 1;
  int pt_bins_ = 100;
  int multiplicity_max_ = 20;
  int cluster_energy_bins_ = 100;
  double pt_max_ = 20.0;
  double cluster_energy_max_ = 20.0;
  double truth_eta_max_ = 0.7;
  int verbosity_ = 0;
  std::array<double, threshold_count_> thresholds_ = {0.0, 0.1, 0.3, 0.5};

  photon_tree::PythiaClusterTruthMatcher truth_matcher_;
  TFile* output_file_ = nullptr;
  std::array<TH1D*, threshold_count_> h_multiplicity_{};
  std::array<TH2D*, threshold_count_> h_multiplicity_vs_pt_{};
  std::array<std::array<TH1D*, threshold_count_>, pathway_count_>
      h_pathway_multiplicity_{};
  TH1D* h_maximum_fraction_ = nullptr;
  TH1D* h_second_fraction_ = nullptr;
  TH2D* h_fraction_vs_cluster_energy_ = nullptr;
  TTree* metadata_tree_ = nullptr;

  unsigned long long n_events_processed_ = 0;
  unsigned long long n_events_written_ = 0;
  unsigned long long n_events_invalid_ = 0;
  unsigned long long n_cluster_considered_ = 0;
  unsigned long long n_cluster_invalid_truth_ = 0;
  unsigned long long n_pi0_candidate_ = 0;
  unsigned long long n_pi0_candidate_g4_primary_ = 0;
  unsigned long long n_pi0_candidate_generator_ = 0;
  unsigned long long n_pi0_malformed_daughters_ = 0;
  unsigned long long n_pi0_cluster_pair_evaluated_ = 0;
  unsigned long long n_pi0_cluster_pair_positive_ = 0;
};

#endif
