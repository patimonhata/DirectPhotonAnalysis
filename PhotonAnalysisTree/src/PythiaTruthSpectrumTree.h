#ifndef RYOTARO_PYTHIATRUTHSPECTRUMTREE_H_20260810
#define RYOTARO_PYTHIATRUTHSPECTRUMTREE_H_20260810

#include "HepMCPhotonClassifier.h"

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>

class PHCompositeNode;
class PHG4TruthInfoContainer;
class PHHepMCGenEventMap;
class TFile;
class TTree;

namespace HepMC
{
class GenEvent;
}

class PythiaTruthSpectrumTree : public SubsysReco
{
 public:
  explicit PythiaTruthSpectrumTree(const std::string& name = "PythiaTruthSpectrumTree");
  ~PythiaTruthSpectrumTree() override;

  int Init(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int End(PHCompositeNode* topNode) override;

  void set_input_file_name(const std::string& value) { input_file_name_ = value; }
  void set_output_file_name(const std::string& value) { output_file_name_ = value; }
  void set_source_file_id(unsigned int value) { source_file_id_ = value; }
  void set_signal_embedding_id(int value) { signal_embedding_id_ = value; }
  void set_hepmc_event_map_node_name(const std::string& value) { hepmc_event_map_node_name_ = value; }
  void set_truth_node_name(const std::string& value) { truth_node_name_ = value; }
  void set_verbosity(int value) { verbosity_ = value; }

 private:
  static constexpr int schema_version_ = 1;

  void create_output_directory() const;
  void create_trees();
  void close_output_file();
  void reset_event();
  bool fill_truth(const PHHepMCGenEventMap* event_map,
                  const PHG4TruthInfoContainer* truth);

  std::string input_file_name_;
  std::string output_file_name_ = "pythia_truth_spectrum_tree.root";
  std::string hepmc_event_map_node_name_ = "PHHepMCGenEventMap";
  std::string truth_node_name_ = "G4TruthInfo";
  std::string sample_type_ = "pythia_truth_spectrum";
  std::string photon_selection_ = "status_1_and_no_end_vertex";
  std::string pi0_selection_ = "pdg_111_without_pdg_111_daughter";
  std::string photon_origin_scheme_ = "follow_unique_photon_parents_then_inspect_production_vertex";
  std::string pi0_decay_photon_scheme_ = "g4_immediate_photon_daughter_of_signal_primary_pi0_matched_by_hepmc_barcode";
  std::string event_weight_scheme_ = "first_hepmc_weight_or_unit_if_empty";
  unsigned int source_file_id_ = 0;
  int signal_embedding_id_ = 1;
  int verbosity_ = 0;

  photon_tree::HepMCPhotonClassifier photon_classifier_;

  TFile* output_file_ = nullptr;
  TTree* event_tree_ = nullptr;
  TTree* metadata_tree_ = nullptr;
  bool metadata_filled_ = false;

  int metadata_schema_version_ = schema_version_;
  int metadata_classifier_version_ = photon_tree::HepMCPhotonClassifier::kAlgorithmVersion;
  unsigned long long n_events_processed_ = 0;
  unsigned long long n_events_written_ = 0;
  unsigned long long n_events_invalid_truth_ = 0;
  unsigned long long n_hepmc_particle_record_ = 0; /* TOTAL number of hepmc particle recond in INPUT FILES */
  unsigned long long n_final_photon_ = 0;
  unsigned long long n_terminal_pi0_ = 0;
  unsigned long long n_g4_pi0_decay_photon_ = 0;
  unsigned long long n_photon_copy_edge_ = 0;
  unsigned long long n_nonterminal_pi0_copy_ = 0;
  unsigned long long n_invalid_photon_ancestry_ = 0;
  unsigned long long n_invalid_kinematics_ = 0;
  unsigned long long n_unmatched_g4_pi0_decay_photon_ = 0;

  unsigned int b_source_file_id_ = 0;
  unsigned int b_event_in_file_ = 0;
  unsigned long long b_event_uid_ = 0;
  int b_hepmc_event_number_ = -999;
  unsigned char b_event_weight_valid_ = 0;
  double b_event_weight_ = 1.0;
  std::vector<double> b_event_weights_;
  unsigned int b_hepmc_n_particle_record_ = 0; /* number of hepmc particle recond in a FUN4ALL EVENT */

  unsigned int b_truth_photon_n_ = 0;
  std::vector<int> b_truth_photon_barcode_;
  std::vector<int> b_truth_photon_status_;
  std::vector<unsigned char> b_truth_photon_kinematics_valid_;
  std::vector<float> b_truth_photon_e_;
  std::vector<float> b_truth_photon_pt_;
  std::vector<float> b_truth_photon_eta_;
  std::vector<float> b_truth_photon_phi_;
  std::vector<unsigned char> b_truth_photon_classification_valid_;
  std::vector<int> b_truth_photon_category_;
  std::vector<int> b_truth_photon_immediate_parent_count_;
  std::vector<int> b_truth_photon_immediate_parent_pdg_;
  std::vector<unsigned char> b_truth_photon_copy_chain_valid_;
  std::vector<unsigned int> b_truth_photon_copy_depth_;
  std::vector<int> b_truth_photon_origin_parent_count_;
  std::vector<int> b_truth_photon_origin_parent_pdg_;
  std::vector<int> b_truth_photon_origin_parent_barcode_;

  unsigned int b_truth_pi0_n_ = 0;
  std::vector<int> b_truth_pi0_barcode_;
  std::vector<int> b_truth_pi0_status_;
  std::vector<unsigned char> b_truth_pi0_kinematics_valid_;
  std::vector<float> b_truth_pi0_e_;
  std::vector<float> b_truth_pi0_pt_;
  std::vector<float> b_truth_pi0_eta_;
  std::vector<float> b_truth_pi0_phi_;
  std::vector<unsigned int> b_truth_pi0_hepmc_direct_photon_count_;

  unsigned int b_truth_pi0_decay_photon_n_ = 0;
  std::vector<int> b_truth_pi0_decay_photon_g4_track_id_;
  std::vector<int> b_truth_pi0_decay_photon_parent_g4_track_id_;
  std::vector<int> b_truth_pi0_decay_photon_parent_hepmc_barcode_;
  std::vector<unsigned char> b_truth_pi0_decay_photon_kinematics_valid_;
  std::vector<float> b_truth_pi0_decay_photon_e_;
  std::vector<float> b_truth_pi0_decay_photon_pt_;
  std::vector<float> b_truth_pi0_decay_photon_eta_;
  std::vector<float> b_truth_pi0_decay_photon_phi_;
};

#endif
