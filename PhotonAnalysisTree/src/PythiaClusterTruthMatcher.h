#ifndef RYOTARO_PYTHIACLUSTERTRUTHMATCHER_H_20260810
#define RYOTARO_PYTHIACLUSTERTRUTHMATCHER_H_20260810

#include <memory>
#include <string>
#include <vector>

class CaloEvalStack;
class PHCompositeNode;
class PHG4TruthInfoContainer;
class PHHepMCGenEventMap;
class RawCluster;
class RawTowerContainer;
class TowerInfoContainer;
class TTree;

namespace photon_tree
{
struct TruthContributor
{
  int g4_track_id = -999;
  int g4_pdg_id = -999;
  int embedding_id = -999;
  int hepmc_barcode = -999;
  float edep = 0.0F;
  float fraction = 0.0F;
  unsigned char hepmc_valid = 0;
  int hepmc_pdg_id = -999;
  int photon_category = -999;
  int photon_source = -999;
  int immediate_parent_pdg = -999;
  int classification_parent_pdg = -999;
};

struct ClusterTruthMatch
{
  bool valid = false;
  float total_edep = 0.0F;
  std::vector<TruthContributor> contributors;
};

class PythiaClusterTruthMatcher
{
 public:
  static constexpr int kAlgorithmVersion = 2;

  PythiaClusterTruthMatcher();
  ~PythiaClusterTruthMatcher();
  PythiaClusterTruthMatcher(const PythiaClusterTruthMatcher&) = delete;
  PythiaClusterTruthMatcher& operator=(const PythiaClusterTruthMatcher&) = delete;

  void set_verbosity(int value) { verbosity_ = value; }
  bool begin_event(PHCompositeNode* topNode);
  ClusterTruthMatch match(const RawCluster* cluster,
                          TowerInfoContainer* towers,
                          RawTowerContainer* raw_truth_towers,
                          PHG4TruthInfoContainer* truth,
                          const PHHepMCGenEventMap* hepmc_event_map,
                          bool allocate_split_tower_energy) const;

 private:
  std::unique_ptr<CaloEvalStack> eval_stack_;
  bool initialized_ = false;
  int verbosity_ = 0;
};

struct ClusterTruthCollection
{
  std::vector<unsigned char> valid;
  std::vector<float> total_edep;
  std::vector<unsigned int> n_contributor;
  std::vector<int> dominant_g4_track_id;
  std::vector<int> dominant_g4_pdg_id;
  std::vector<int> dominant_embedding_id;
  std::vector<int> dominant_hepmc_barcode;
  std::vector<float> dominant_edep;
  std::vector<float> dominant_fraction;
  std::vector<unsigned char> dominant_hepmc_valid;
  std::vector<int> dominant_hepmc_pdg_id;
  std::vector<int> dominant_photon_category;
  std::vector<int> dominant_photon_source;
  std::vector<int> dominant_immediate_parent_pdg;
  std::vector<int> dominant_classification_parent_pdg;

  // Flat contributor table. offset has ncluster+1 entries.
  std::vector<unsigned int> contributor_offset;
  std::vector<unsigned int> contributor_cluster_index;
  std::vector<int> contributor_g4_track_id;
  std::vector<int> contributor_g4_pdg_id;
  std::vector<int> contributor_embedding_id;
  std::vector<int> contributor_hepmc_barcode;
  std::vector<float> contributor_edep;
  std::vector<float> contributor_fraction;
  std::vector<unsigned char> contributor_hepmc_valid;
  std::vector<int> contributor_hepmc_pdg_id;
  std::vector<int> contributor_photon_category;
  std::vector<int> contributor_photon_source;
  std::vector<int> contributor_immediate_parent_pdg;
  std::vector<int> contributor_classification_parent_pdg;

  void clear();
  void append(const ClusterTruthMatch& match);
  void create_branches(TTree* tree, const std::string& prefix);
};
}

#endif
