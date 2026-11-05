#ifndef PHOTON_TREE_COMMON_DECAY_TRUTH_H
#define PHOTON_TREE_COMMON_DECAY_TRUTH_H
#include "CommonTruthComposition.h"
#include "PythiaClusterTruthMatcher.h"
#include <array>
#include <memory>
class PHG4Particle;
class PHG4TruthInfoContainer;
class PHHepMCGenEventMap;
class PHG4CellContainer;
class PHG4HitContainer;
namespace photon_tree
{
struct CommonDiphoton
{
  composition::Particle parent;
  int g4_parent_track = 0;
  std::array<const PHG4Particle*, 2> photons = {nullptr, nullptr};
  std::array<composition::Identity, 2> photon_ids;
  double energy = 0, pt = 0, eta = 0, phi = 0;
};
class CommonDecayTruth
{
 public:
  CommonDecayTruth(PHG4TruthInfoContainer* truth, const PHHepMCGenEventMap* events);
  ~CommonDecayTruth();
  composition::Result resolve(const ClusterTruthMatch& original, const RawCluster* cluster, TowerInfoContainer* towers,
                              RawTowerContainer* raw, PHG4CellContainer* cells, PHG4HitContainer* hits) const;
  std::vector<CommonDiphoton> diphotons(int pdg, int embedding) const;
 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};
// Flat event-local tables. Every row retains its own ancestry and identity; parent rows are overlapping summaries, never extra deposit.
struct CommonTruthCollection
{
  std::vector<int> main, selected_pdg, selected_embedding, selected_source, selected_id, applicability, topology_reason, legacy_origin, matching_status, anchor_matching_status;
  std::vector<unsigned int> quality;
  std::vector<unsigned char> long_lived_ancestor, anchor_matching_usable;
  std::vector<double> total, attributed, unresolved, fraction, max_photon, other_photons, charged, minimum_coverage, anchor_matching_coverage;
  std::vector<unsigned int> offset, ancestor_offset;
  std::vector<int> table, pdg, embedding, source, id, decay, prompt, unresolved_piece;
  std::vector<unsigned int> piece_quality;
  std::vector<double> deposit, flight_cm, proper_length_cm;
  std::vector<std::vector<double>> production_vertex, decay_vertex;
  std::vector<int> ancestor_pdg, ancestor_embedding, ancestor_source, ancestor_id;
  std::vector<double> ancestor_flight_cm, ancestor_proper_length_cm;
  void clear();
  void append(const composition::Result& result, int legacy = 3);
  void create_branches(TTree* tree);
};
}
#endif
