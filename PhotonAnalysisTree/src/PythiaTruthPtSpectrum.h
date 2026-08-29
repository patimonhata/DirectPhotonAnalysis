#ifndef RYOTARO_PYTHIATRUTHPTSPECTRUM_H_20260829
#define RYOTARO_PYTHIATRUTHPTSPECTRUM_H_20260829

#include "HepMCPhotonClassifier.h"

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class TFile;
class TH1D;
class TTree;

class PythiaTruthPtSpectrum : public SubsysReco
{
 public:
  explicit PythiaTruthPtSpectrum(const std::string& name = "PythiaTruthPtSpectrum");
  ~PythiaTruthPtSpectrum() override;

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
  void set_input_file_prefix(const std::string& value) { input_file_prefix_ = value; }
  void set_signal_embedding_id(int value) { signal_embedding_id_ = value; }
  void set_hepmc_event_map_node_name(const std::string& value) { hepmc_event_map_node_name_ = value; }
  void set_truth_node_name(const std::string& value) { truth_node_name_ = value; }
  void set_binning(int bins, double pt_max)
  {
    n_bins_ = bins;
    pt_max_ = pt_max;
  }
  void set_max_abs_eta(double value) { max_abs_eta_ = value; }
  void set_use_event_weight(bool value) { use_event_weight_ = value; }
  void set_verbosity(int value) { verbosity_ = value; }

 private:
  static constexpr int schema_version_ = 3;

  void create_output_directory() const;
  void create_output();
  void close_output();

  std::string output_file_name_ = "pythia_truth_pt_partial.root";
  std::string manifest_path_;
  std::string first_suffix_;
  std::string last_suffix_;
  std::string input_file_prefix_ = "G4Hits_";
  std::string hepmc_event_map_node_name_ = "PHHepMCGenEventMap";
  std::string truth_node_name_ = "G4TruthInfo";
  std::string photon_selection_ = "prompt_category_1_or_2";
  std::string pi0_decay_photon_selection_ = "hepmc_final_photon_with_valid_single_pi0_origin_plus_g4_immediate_photon_daughter_of_signal_primary_pi0";
  long long manifest_begin_ = -1;
  long long manifest_end_ = -1;
  long long files_added_ = 0;
  int signal_embedding_id_ = 1;
  int n_bins_ = 100;
  double pt_max_ = 20.0;
  double max_abs_eta_ = 0.7;
  bool use_event_weight_ = false;
  int verbosity_ = 0;

  photon_tree::HepMCPhotonClassifier photon_classifier_;
  TFile* output_file_ = nullptr;
  TH1D* h_prompt_photon_ = nullptr;
  TH1D* h_pi0_ = nullptr;
  TH1D* h_pi0_decay_photon_ = nullptr;
  TH1D* h_hepmc_pi0_decay_photon_ = nullptr;
  TH1D* h_g4_pi0_decay_photon_ = nullptr;
  TTree* metadata_tree_ = nullptr;

  int metadata_schema_version_ = schema_version_;
  unsigned char metadata_use_event_weight_ = 0U;
  unsigned char metadata_bin_width_normalized_ = 0U;
  long long n_events_processed_ = 0;
  long long n_events_written_ = 0;
  unsigned long long n_malformed_events_ = 0;
  unsigned long long n_invalid_weight_events_ = 0;
  unsigned long long n_prompt_photon_ = 0;
  unsigned long long n_pi0_ = 0;
  unsigned long long n_pi0_decay_photon_ = 0;
  unsigned long long n_hepmc_pi0_decay_photon_ = 0;
  unsigned long long n_g4_pi0_decay_photon_ = 0;
};

#endif
