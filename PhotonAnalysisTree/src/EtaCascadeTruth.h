#ifndef PHOTON_TREE_ETA_CASCADE_TRUTH_H
#define PHOTON_TREE_ETA_CASCADE_TRUTH_H

#include "PythiaClusterTruthMatcher.h"
#include <map>
#include <vector>
class PHG4Particle;
class PHG4TruthInfoContainer;
class PHHepMCGenEventMap;
class PHG4CellContainer;
class PHG4HitContainer;
class RawCluster;
class TowerInfoContainer;
class RawTowerContainer;
namespace photon_tree
{
// Exact direct-daughter modes; never infer a decay from an arbitrary descendant pion.
enum class EtaDecayMode : int { unknown = 0, two_gamma = 1, three_pi0 = 2, charged_pions_pi0 = 3, other = 4 };
class EtaCascadeTruth
{
 public:
  EtaCascadeTruth(PHG4TruthInfoContainer* truth, const PHHepMCGenEventMap* events);
  void refine(ClusterTruthMatch& match, const RawCluster* cluster, TowerInfoContainer* towers, RawTowerContainer* raw,
              PHG4CellContainer* cells, PHG4HitContainer* hits) const;
  bool is_eta_daughter_pi0(const PHG4Particle* particle) const;
 private:
  void annotate(TruthContributor& contributor) const;
  int g4_mode(int eta_track) const;
  PHG4TruthInfoContainer* truth_;
  const PHHepMCGenEventMap* events_;
  std::map<int, std::vector<const PHG4Particle*>> children_;
};
}
#endif
