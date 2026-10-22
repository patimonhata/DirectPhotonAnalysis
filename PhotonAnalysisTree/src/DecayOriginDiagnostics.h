#ifndef PHOTON_TREE_DECAY_ORIGIN_DIAGNOSTICS_H
#define PHOTON_TREE_DECAY_ORIGIN_DIAGNOSTICS_H
#include "PythiaClusterTruthMatcher.h"
#include <memory>
#include <string>
#include <vector>
class PHG4TruthInfoContainer;
class PHHepMCGenEventMap;
class PHG4CellContainer;
class PHG4HitContainer;
class RawCluster;
class TowerInfoContainer;
class RawTowerContainer;
class TTree;
namespace photon_tree
{
// These are diagnostic contributors, never replacements for ClusterTruthMatch or topology anchors.
struct DecayOriginContribution
{
  int pdg = 0, embedding = 0, source = 0, id = 0; // source: 1 HepMC barcode, 2 G4 track ID, 0 unavailable.
  int root_track_id = 0, parent_pdg = 0, parent_source = 0, parent_id = 0;
  int status = 0, pathway = 0; // status: 0 retained, 1 measured split, -1 incomplete, -2 interaction, -3 unknown flight/lifetime.
  int diphoton = 0, long_lived = 0, daughter_count = 0, daughter_pi0_count = 0, daughter_gamma_count = 0, daughter_charged_count = 0;
  float fraction = 0, gamma0_fraction = -1, gamma1_fraction = -1;
  float flight_cm = -1, parent_flight_cm = -1, proper_length_cm = -1;
};
struct DecayOriginMatch
{
  bool valid = false;
  std::vector<DecayOriginContribution> contributors;
};
class DecayOriginDiagnostics
{
 public:
  static constexpr int kAlgorithmVersion = 1;
  static constexpr double kMaxProperLengthCm = 1e-4; // 1 micrometre; excludes weakly decaying hadrons.
  static constexpr double kMaxFlightCm = 1e-2; // 0.1 mm lab displacement; saved in metadata, not a reco selection.
  DecayOriginDiagnostics(PHG4TruthInfoContainer* truth, const PHHepMCGenEventMap* events);
  ~DecayOriginDiagnostics();
  DecayOriginMatch diagnose(const ClusterTruthMatch& original, const RawCluster* cluster, TowerInfoContainer* towers,
                           RawTowerContainer* raw, PHG4CellContainer* cells, PHG4HitContainer* hits) const;
 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};
struct DecayOriginCollection
{
  std::vector<unsigned char> valid;
  std::vector<unsigned int> offset;
  std::vector<int> pdg, embedding, source, id, root_track_id, parent_pdg, parent_source, parent_id, status, pathway, diphoton, long_lived;
  std::vector<int> daughter_count, daughter_pi0_count, daughter_gamma_count, daughter_charged_count;
  std::vector<float> fraction, gamma0_fraction, gamma1_fraction, flight_cm, parent_flight_cm, proper_length_cm;
  void clear();
  void append(const DecayOriginMatch& match);
  void create_branches(TTree* tree);
};
}
#endif
