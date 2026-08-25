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
struct Pi0ClusterTruthMatch
{
  bool valid = false;
  float total_edep = 0.0F;
  std::array<float, 2> gamma_edep = {0.0F, 0.0F};
  float other_edep = 0.0F;
};

class Pi0ClusterTruthMatcher
{
 public:
  static constexpr int kAlgorithmVersion = 2;

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
