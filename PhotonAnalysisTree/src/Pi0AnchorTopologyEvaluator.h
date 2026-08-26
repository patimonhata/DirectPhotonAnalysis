#ifndef RYOTARO_PI0ANCHORTOPOLOGYEVALUATOR_H_20260824
#define RYOTARO_PI0ANCHORTOPOLOGYEVALUATOR_H_20260824

#include "Pi0ClusterTruthMatcher.h"
#include "PythiaClusterTruthMatcher.h"

#include <array>
#include <cstddef>
#include <string>
#include <vector>

class PHCompositeNode;
class RawCluster;

namespace photon_tree
{
enum class Pi0SampleMode : int
{
  pythia = 1,
  single_particle = 2
};

enum class Pi0Pathway : int
{
  g4_primary_decay = 1,
  generator_decay = 2,
  single_particle_g4_decay = 3
};

enum class Pi0AnchorTopology : int
{
  other = 0,
  separated = 1,
  merged = 2,
  missing = 3
};

enum class Pi0AnchorReason : int
{
  other_not_daughter_maximum = 0,
  separated_distinct_recovered_clusters = 1,
  merged_shared_recovered_cluster = 2,
  missing_unrecovered_partner = 3,
  ambiguous_main_contributor = 4,
  other_best_cluster_below_recovery = 5
};

enum class Pi0MissingDetail : int
{
  not_missing = 0,
  partner_best_below_recovery = 1,
  partner_cluster_below_energy_threshold_recovered = 2,
  partner_cluster_below_energy_threshold_below_recovery = 3,
  partner_direct_match_incomplete = 4,
  partner_no_direct_deposit = 5
};

enum class Pi0TopologyEventStatus : int
{
  invalid = 0,
  vertex_rejected = 1,
  accepted = 2
};

struct Pi0AnchorTopologyConfig
{
  Pi0SampleMode sample_mode = Pi0SampleMode::pythia;
  std::string truth_node_name = "G4TruthInfo";
  std::string hepmc_event_map_node_name = "PHHepMCGenEventMap";
  std::string tower_node_name = "TOWERINFO_CALIB_CEMC";
  std::string raw_truth_tower_node_name = "TOWER_SIM_CEMC";
  std::string truth_cell_node_name = "G4CELL_CEMC";
  std::string truth_hit_node_name = "G4HIT_CEMC";
  std::string cluster_node_name = "CLUSTERINFO_CEMC";
  int signal_embedding_id = 1;
  double anchor_cluster_eta_max = 0.7;
  double partner_cluster_eta_max = -1.0;
  double min_cluster_energy = 0.2;
  double dominant_fraction_min = 0.5;
  double anchor_pi0_fraction_min = 0.5;
  double min_energy_contribution_fraction = 0.0;
  double min_photon_energy_recovery = 0.5;
  double min_direct_match_cluster_energy_coverage = 0.5;
  double missing_diagnostic_max_delta_r = 0.15;
  // A non-positive value disables the event-level collision-z cut.
  double max_abs_vertex_z = -1.0;
  bool evaluate_all_candidates = false;
  bool enable_missing_diagnostics = false;
  int verbosity = 0;
};

struct Pi0TopologyClusterRecord
{
  const RawCluster* cluster = nullptr;
  unsigned int cluster_id = 0;
  double energy = 0.0;
  double et = 0.0;
  double eta = 0.0;
  double phi = 0.0;
  bool anchor_acceptance = false;
  ClusterTruthMatch truth;
};

struct Pi0PartnerDiagnosticRecord
{
  bool found = false;
  bool below_energy_threshold = false;
  bool has_direct_deposit = false;
  unsigned int cluster_id = 0;
  double cluster_energy = 0.0;
  double cluster_eta = 0.0;
  double cluster_phi = 0.0;
  double delta_r = -1.0;
  double reconstructed_photon_energy = 0.0;
  double recovery = 0.0;
  Pi0ClusterTruthMatch match;
};

struct Pi0TopologyCandidateRecord
{
  Pi0Pathway pathway = Pi0Pathway::g4_primary_decay;
  int parent_barcode = 0;
  int g4_parent_track_id = -999;
  double energy = 0.0;
  double pt = 0.0;
  double eta = 0.0;
  double phi = 0.0;
  std::array<int, 2> photon_track_ids = {-999, -999};
  std::array<double, 2> photon_energy = {0.0, 0.0};
  std::array<double, 2> photon_eta = {0.0, 0.0};
  std::array<double, 2> photon_phi = {0.0, 0.0};
  bool topology_evaluated = false;
  std::array<std::size_t, 2> best_cluster = {static_cast<std::size_t>(-1), static_cast<std::size_t>(-1)};
  std::array<double, 2> maximum_edep = {-1.0, -1.0};
  std::array<double, 2> reconstructed_photon_energy = {0.0, 0.0};
  std::array<bool, 2> recovered = {false, false};
  std::vector<Pi0ClusterTruthMatch> cluster_matches;
  std::array<Pi0PartnerDiagnosticRecord, 2> partner_diagnostics;
};

struct Pi0TopologyAnchorRecord
{
  std::size_t cluster_index = static_cast<std::size_t>(-1);
  std::size_t candidate_index = static_cast<std::size_t>(-1);
  double main_fraction = -1.0;
  double second_fraction = -1.0;
  double unmatched_max_fraction = 0.0;
  bool ambiguous_main = false;
  Pi0AnchorTopology topology = Pi0AnchorTopology::other;
  Pi0AnchorReason reason = Pi0AnchorReason::other_not_daughter_maximum;
  Pi0MissingDetail missing_detail = Pi0MissingDetail::not_missing;
  int partner_photon_index = -1;
};

struct Pi0AnchorTopologyEventResult
{
  Pi0TopologyEventStatus status = Pi0TopologyEventStatus::invalid;
  std::array<double, 3> collision_vertex = {0.0, 0.0, 0.0};
  std::vector<Pi0TopologyClusterRecord> clusters;
  std::vector<unsigned char> prompt_cluster;
  std::vector<Pi0TopologyCandidateRecord> candidates;
  std::vector<Pi0TopologyAnchorRecord> anchors;
  unsigned long long cluster_invalid_truth_count = 0;
  unsigned long long malformed_candidate_count = 0;
  unsigned long long g4_candidate_count = 0;
  unsigned long long generator_candidate_count = 0;
  unsigned long long energy_match_invalid_count = 0;
};

const char* pi0_missing_detail_name(Pi0MissingDetail value);
const char* pi0_pathway_name(Pi0Pathway value);
const char* pi0_anchor_topology_name(Pi0AnchorTopology value);
const char* pi0_anchor_reason_name(Pi0AnchorReason value);

class Pi0AnchorTopologyEvaluator
{
 public:
  static constexpr int kAlgorithmVersion = 4;
  void configure(const Pi0AnchorTopologyConfig& config);
  const Pi0AnchorTopologyConfig& config() const { return config_; }
  Pi0AnchorTopologyEventResult evaluate(PHCompositeNode* topNode);

 private:
  Pi0AnchorTopologyConfig config_;
  PythiaClusterTruthMatcher truth_matcher_;
  Pi0ClusterTruthMatcher pi0_truth_matcher_;
};
}

#endif
