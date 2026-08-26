#ifndef RYOTARO_PI0CLUSTERTRUTHMATCHER_H_20260810
#define RYOTARO_PI0CLUSTERTRUTHMATCHER_H_20260810

#include <array>
#include <string>
#include <vector>

class PHG4TruthInfoContainer;
class PHG4CellContainer;
class PHG4HitContainer;
class RawCluster;
class RawTowerContainer;
class TowerInfoContainer;
class TTree;

namespace photon_tree
{
enum class Pi0ClusterTruthMatchStatus : int
{
  invalid = 0,
  partial = 1,
  complete = 2
};

enum class Pi0ClusterTruthMatchFailure : int
{
  none = 0,
  missing_input = 1,
  missing_tower_info = 2,
  invalid_tower_energy = 3,
  missing_truth_tower = 4,
  missing_g4_cell = 5,
  missing_g4_hit = 6,
  invalid_hit_edep = 7
};

struct Pi0ClusterTruthMatch
{
  bool valid = false;
  bool usable = false;
  Pi0ClusterTruthMatchStatus status = Pi0ClusterTruthMatchStatus::invalid;
  Pi0ClusterTruthMatchFailure failure = Pi0ClusterTruthMatchFailure::none;
  int failure_ieta = -999;
  int failure_iphi = -999;
  unsigned int failure_tower_key = 0;
  unsigned long long failure_cell_id = 0;
  unsigned long long failure_hit_id = 0;
  unsigned int tower_count = 0;
  unsigned int matched_tower_count = 0;
  float cluster_member_energy = 0.0F;
  float matched_cluster_member_energy = 0.0F;
  float cluster_member_energy_coverage = 0.0F;
  float total_edep = 0.0F;
  std::array<float, 2> gamma_edep = {0.0F, 0.0F};
  float other_edep = 0.0F;
};

const char* pi0_cluster_truth_match_status_name(Pi0ClusterTruthMatchStatus value);
const char* pi0_cluster_truth_match_failure_name(Pi0ClusterTruthMatchFailure value);

class Pi0ClusterTruthMatcher
{
 public:
  static constexpr int kAlgorithmVersion = 3;

  Pi0ClusterTruthMatch match(
      const RawCluster* cluster,
      TowerInfoContainer* towers,
      RawTowerContainer* raw_truth_towers,
      PHG4CellContainer* cells,
      PHG4HitContainer* hits,
      PHG4TruthInfoContainer* truth,
      const std::array<int, 2>& direct_gamma_track_ids,
      bool allocate_split_tower_energy) const;

 private:
  static int direct_gamma_index(int track_id, PHG4TruthInfoContainer* truth, const std::array<int, 2>& direct_gamma_track_ids);
};

struct Pi0ClusterTruthCollection
{
  std::vector<unsigned char> match_valid;
  std::vector<float> total_edep;
  std::vector<float> gamma0_edep;
  std::vector<float> gamma1_edep;
  std::vector<float> other_edep;
  std::vector<float> gamma0_fraction;
  std::vector<float> gamma1_fraction;
  std::vector<float> other_fraction;
  std::vector<float> gamma0_recovery;
  std::vector<float> gamma1_recovery;

  void clear();
  void append(const Pi0ClusterTruthMatch& match,
              float truth_gamma0_energy,
              float truth_gamma1_energy);
  void create_branches(TTree* tree, const std::string& prefix);
};
}

#endif
