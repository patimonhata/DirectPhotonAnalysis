#ifndef RYOTARO_TOPOLOGYEVENTRENDERER_H_20260824
#define RYOTARO_TOPOLOGYEVENTRENDERER_H_20260824

#include <TBox.h>
#include <TCanvas.h>
#include <TEllipse.h>
#include <TFile.h>
#include <TH2F.h>
#include <TLatex.h>
#include <TLine.h>
#include <TPad.h>
#include <TMarker.h>
#include <TStyle.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace topology_display
{
struct Event
{
  int id = -1;
  double cx = 0.0;
  double cy = 0.0;
  double cz = 0.0;
  int candidates = 0;
  int g4_primary = 0;
  int generator = 0;
  int secondary = 0;
  int family_particles = 0;
  int truth_particles = 0;
  int clusters = 0;
  int anchors = 0;
  int separated = 0;
  int merged = 0;
  int single_contaminated = 0;
  int missing = 0;
  int other = 0;
};

struct Candidate
{
  int id = -1;
  int pathway = 0;
  std::string pathway_name;
  int barcode = -999;
  int g4_parent = -999;
  double energy = 0.0;
  double pt = 0.0;
  double eta = 0.0;
  double phi = 0.0;
  int photon_track[2] = {-999, -999};
  double photon_energy[2] = {0.0, 0.0};
  double photon_eta[2] = {0.0, 0.0};
  double photon_phi[2] = {0.0, 0.0};
  int photon_projection_valid[2] = {0, 0};
  double photon_projection_eta[2] = {-999.0, -999.0};
  double photon_projection_phi[2] = {-999.0, -999.0};
  int photon_in_cemc_acceptance[2] = {0, 0};
  int photon_first_daughter_vertex_valid[2] = {0, 0};
  double photon_first_daughter_radius[2] = {-999.0, -999.0};
  int photon_pre_cemc_interaction[2] = {0, 0};
  int best_cluster[2] = {-999, -999};
  double maximum_edep[2] = {-1.0, -1.0};
  double reconstructed[2] = {0.0, 0.0};
  int recovered[2] = {0, 0};
};

struct Anchor
{
  int id = -1;
  int candidate = -1;
  unsigned int cluster = 0;
  double energy = 0.0;
  double et = 0.0;
  int topology = 0;
  std::string topology_name;
  int reason = 0;
  std::string reason_name;
  int missing_category = 0;
  std::string missing_category_name;
  int missing_detail = 0;
  std::string missing_detail_name;
  int partner_photon = -1;
  int pre_cemc_photon = -1;
  double partner_cemc_edep = 0.0;
  double diagnostic_invariant_mass = -999.0;
  double main_fraction = -1.0;
  double second_fraction = -1.0;
  double unmatched_fraction = 0.0;
  int ambiguous = 0;
  int best_cluster[2] = {-999, -999};
  int recovered[2] = {0, 0};
  double truth_energy[2] = {0.0, 0.0};
  double reconstructed[2] = {0.0, 0.0};
  int match_valid = 0;
  int match_usable = 0;
  int match_status = 0;
  std::string match_status_name;
  int match_failure = 0;
  std::string match_failure_name;
  int match_failure_ieta = -999;
  int match_failure_iphi = -999;
  unsigned int match_towers = 0;
  unsigned int match_matched_towers = 0;
  double match_coverage = 0.0;
  double total_edep = 0.0;
  double gamma_edep[2] = {0.0, 0.0};
  double other_edep = 0.0;
  int diagnostic_found = 0;
  int diagnostic_below_threshold = 0;
  int diagnostic_has_direct = 0;
  int diagnostic_cluster = -999;
  double diagnostic_energy = -999.0;
  double diagnostic_eta = -999.0;
  double diagnostic_phi = -999.0;
  double diagnostic_delta_r = -999.0;
  double diagnostic_reconstructed = 0.0;
  double diagnostic_recovery = 0.0;
  int diagnostic_match_usable = 0;
  std::string diagnostic_match_status_name;
  std::string diagnostic_match_failure_name;
  int diagnostic_failure_ieta = -999;
  int diagnostic_failure_iphi = -999;
  double diagnostic_match_coverage = 0.0;
};

struct Segment
{
  int track = -999;
  int pid = 0;
  int family = -1;
  int gamma = -1;
  double x0 = 0.0;
  double y0 = 0.0;
  double z0 = 0.0;
  double x1 = 0.0;
  double y1 = 0.0;
  double z1 = 0.0;
};

struct Cluster
{
  unsigned int id = 0;
  double energy = 0.0;
  double et = 0.0;
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
  double r = 0.0;
  double eta = 0.0;
  double phi = 0.0;
  int ntowers = 0;
  int considered = 0;
};

struct Tower
{
  unsigned int cluster = 0;
  int ieta = -1;
  int iphi = -1;
  double member_energy = 0.0;
  double tower_energy = 0.0;
  double allocation = 0.0;
};

struct CandidateClusterMatch
{
  int candidate = -1;
  unsigned int cluster = 0;
  int valid = 0;
  double gamma_fraction[2] = {0.0, 0.0};
};

struct DisplayData
{
  Event event;
  std::vector<Candidate> candidates;
  std::vector<Anchor> anchors;
  std::vector<Segment> segments;
  std::vector<Cluster> clusters;
  std::vector<Tower> towers;
  std::vector<CandidateClusterMatch> matches;
};

inline int topology_color(int topology)
{
  if (topology == 1) return kGreen + 2;
  if (topology == 2) return kMagenta + 1;
  if (topology == 3) return kOrange + 7;
  if (topology == 4) return kAzure + 2;
  return kGray + 2;
}

constexpr int photon0_line_style = 11;
constexpr int photon1_line_style = 12;
constexpr double cluster_fill_threshold_gev = 0.1;

inline double cluster_marker_size(double energy, double maximum_size)
{
  const double nonnegative_energy = std::max(0.0, energy);
  if (nonnegative_energy < cluster_fill_threshold_gev) return 0.18 + 1.2 * nonnegative_energy;
  return std::min(maximum_size, 1.0 + 0.18 * (nonnegative_energy - cluster_fill_threshold_gev));
}

inline int family_color(int family)
{
  static const int colors[] = {kAzure + 1, kRed + 1, kViolet + 1,
                               kTeal + 2, kOrange + 1, kPink + 7};
  if (family < 0) return kGray + 1;
  return colors[family % 6];
}

inline int family_line_style(int gamma)
{
  if (gamma == 0) return photon0_line_style;
  if (gamma == 1) return photon1_line_style;
  return 1;
}

inline bool is_core_segment(const DisplayData& data, const Segment& segment)
{
  for (const auto& candidate : data.candidates)
    if (segment.track == candidate.g4_parent || segment.track == candidate.photon_track[0] || segment.track == candidate.photon_track[1]) return true;
  return false;
}

inline bool load_event(TFile* file, int event_id, DisplayData& data)
{
  if (!file || file->IsZombie()) return false;
  data = DisplayData{};
  TTree* events = static_cast<TTree*>(file->Get("events"));
  TTree* candidates = static_cast<TTree*>(file->Get("pi0_candidates"));
  TTree* anchors = static_cast<TTree*>(file->Get("anchor_decisions"));
  TTree* matches = static_cast<TTree*>(file->Get("candidate_cluster_truth"));
  TTree* segments = static_cast<TTree*>(file->Get("truth_segments"));
  TTree* clusters = static_cast<TTree*>(file->Get("clusters"));
  TTree* towers = static_cast<TTree*>(file->Get("cluster_towers"));
  if (!events || !candidates || !anchors) return false;

  int event = -1;
  events->SetBranchAddress("event", &event);
  events->SetBranchAddress("collision_x", &data.event.cx);
  events->SetBranchAddress("collision_y", &data.event.cy);
  events->SetBranchAddress("collision_z", &data.event.cz);
  events->SetBranchAddress("n_candidates", &data.event.candidates);
  events->SetBranchAddress("n_g4_primary_pi0", &data.event.g4_primary);
  events->SetBranchAddress("n_generator_pi0", &data.event.generator);
  events->SetBranchAddress("n_g4_secondary_pi0", &data.event.secondary);
  events->SetBranchAddress("n_selected_family_particles", &data.event.family_particles);
  events->SetBranchAddress("n_truth_particles", &data.event.truth_particles);
  events->SetBranchAddress("n_clusters", &data.event.clusters);
  events->SetBranchAddress("n_anchors", &data.event.anchors);
  events->SetBranchAddress("n_separated", &data.event.separated);
  events->SetBranchAddress("n_merged", &data.event.merged);
  if (events->GetBranch("n_single_contaminated")) events->SetBranchAddress("n_single_contaminated", &data.event.single_contaminated);
  events->SetBranchAddress("n_missing", &data.event.missing);
  events->SetBranchAddress("n_other", &data.event.other);
  bool found = false;
  for (Long64_t entry = 0; entry < events->GetEntries(); ++entry)
  {
    events->GetEntry(entry);
    if (event == event_id)
    {
      data.event.id = event;
      found = true;
      break;
    }
  }
  events->ResetBranchAddresses();
  if (!found) return false;

  Candidate candidate;
  std::string* pathway_name = nullptr;
  candidates->SetBranchAddress("event", &event);
  candidates->SetBranchAddress("candidate_id", &candidate.id);
  candidates->SetBranchAddress("pathway", &candidate.pathway);
  candidates->SetBranchAddress("pathway_name", &pathway_name);
  candidates->SetBranchAddress("parent_barcode", &candidate.barcode);
  candidates->SetBranchAddress("g4_parent_track_id", &candidate.g4_parent);
  candidates->SetBranchAddress("energy", &candidate.energy);
  candidates->SetBranchAddress("pt", &candidate.pt);
  candidates->SetBranchAddress("eta", &candidate.eta);
  candidates->SetBranchAddress("phi", &candidate.phi);
  candidates->SetBranchAddress("photon0_track_id", &candidate.photon_track[0]);
  candidates->SetBranchAddress("photon1_track_id", &candidate.photon_track[1]);
  candidates->SetBranchAddress("photon0_energy", &candidate.photon_energy[0]);
  candidates->SetBranchAddress("photon1_energy", &candidate.photon_energy[1]);
  candidates->SetBranchAddress("photon0_eta", &candidate.photon_eta[0]);
  candidates->SetBranchAddress("photon1_eta", &candidate.photon_eta[1]);
  candidates->SetBranchAddress("photon0_phi", &candidate.photon_phi[0]);
  candidates->SetBranchAddress("photon1_phi", &candidate.photon_phi[1]);
  if (candidates->GetBranch("photon0_projection_valid")) candidates->SetBranchAddress("photon0_projection_valid", &candidate.photon_projection_valid[0]);
  if (candidates->GetBranch("photon1_projection_valid")) candidates->SetBranchAddress("photon1_projection_valid", &candidate.photon_projection_valid[1]);
  if (candidates->GetBranch("photon0_projection_eta")) candidates->SetBranchAddress("photon0_projection_eta", &candidate.photon_projection_eta[0]);
  if (candidates->GetBranch("photon1_projection_eta")) candidates->SetBranchAddress("photon1_projection_eta", &candidate.photon_projection_eta[1]);
  if (candidates->GetBranch("photon0_projection_phi")) candidates->SetBranchAddress("photon0_projection_phi", &candidate.photon_projection_phi[0]);
  if (candidates->GetBranch("photon1_projection_phi")) candidates->SetBranchAddress("photon1_projection_phi", &candidate.photon_projection_phi[1]);
  if (candidates->GetBranch("photon0_in_cemc_acceptance")) candidates->SetBranchAddress("photon0_in_cemc_acceptance", &candidate.photon_in_cemc_acceptance[0]);
  if (candidates->GetBranch("photon1_in_cemc_acceptance")) candidates->SetBranchAddress("photon1_in_cemc_acceptance", &candidate.photon_in_cemc_acceptance[1]);
  if (candidates->GetBranch("photon0_first_daughter_vertex_valid")) candidates->SetBranchAddress("photon0_first_daughter_vertex_valid", &candidate.photon_first_daughter_vertex_valid[0]);
  if (candidates->GetBranch("photon1_first_daughter_vertex_valid")) candidates->SetBranchAddress("photon1_first_daughter_vertex_valid", &candidate.photon_first_daughter_vertex_valid[1]);
  if (candidates->GetBranch("photon0_first_daughter_radius")) candidates->SetBranchAddress("photon0_first_daughter_radius", &candidate.photon_first_daughter_radius[0]);
  if (candidates->GetBranch("photon1_first_daughter_radius")) candidates->SetBranchAddress("photon1_first_daughter_radius", &candidate.photon_first_daughter_radius[1]);
  if (candidates->GetBranch("photon0_pre_cemc_interaction")) candidates->SetBranchAddress("photon0_pre_cemc_interaction", &candidate.photon_pre_cemc_interaction[0]);
  if (candidates->GetBranch("photon1_pre_cemc_interaction")) candidates->SetBranchAddress("photon1_pre_cemc_interaction", &candidate.photon_pre_cemc_interaction[1]);
  candidates->SetBranchAddress("best_cluster0_id", &candidate.best_cluster[0]);
  candidates->SetBranchAddress("best_cluster1_id", &candidate.best_cluster[1]);
  candidates->SetBranchAddress("maximum_edep0", &candidate.maximum_edep[0]);
  candidates->SetBranchAddress("maximum_edep1", &candidate.maximum_edep[1]);
  candidates->SetBranchAddress("reconstructed_photon0_energy", &candidate.reconstructed[0]);
  candidates->SetBranchAddress("reconstructed_photon1_energy", &candidate.reconstructed[1]);
  candidates->SetBranchAddress("recovered0", &candidate.recovered[0]);
  candidates->SetBranchAddress("recovered1", &candidate.recovered[1]);
  for (Long64_t entry = 0; entry < candidates->GetEntries(); ++entry)
  {
    candidates->GetEntry(entry);
    if (event == event_id)
    {
      candidate.pathway_name = pathway_name ? *pathway_name : "unknown";
      data.candidates.push_back(candidate);
    }
  }
  candidates->ResetBranchAddresses();

  Anchor anchor;
  const auto bind_anchor = [&](const char* name, auto* address) {
    if (anchors->GetBranch(name)) anchors->SetBranchAddress(name, address);
  };
  std::string* topology_name = nullptr;
  std::string* reason_name = nullptr;
  std::string* missing_category_name = nullptr;
  std::string* missing_detail_name = nullptr;
  std::string* match_status_name = nullptr;
  std::string* match_failure_name = nullptr;
  std::string* diagnostic_match_status_name = nullptr;
  std::string* diagnostic_match_failure_name = nullptr;
  anchors->SetBranchAddress("event", &event);
  anchors->SetBranchAddress("anchor_id", &anchor.id);
  anchors->SetBranchAddress("candidate_id", &anchor.candidate);
  anchors->SetBranchAddress("cluster_id", &anchor.cluster);
  anchors->SetBranchAddress("energy", &anchor.energy);
  anchors->SetBranchAddress("et", &anchor.et);
  anchors->SetBranchAddress("topology", &anchor.topology);
  anchors->SetBranchAddress("topology_name", &topology_name);
  anchors->SetBranchAddress("reason", &anchor.reason);
  anchors->SetBranchAddress("reason_name", &reason_name);
  bind_anchor("missing_category", &anchor.missing_category);
  bind_anchor("missing_category_name", &missing_category_name);
  bind_anchor("missing_detail", &anchor.missing_detail);
  bind_anchor("missing_detail_name", &missing_detail_name);
  bind_anchor("partner_photon_index", &anchor.partner_photon);
  bind_anchor("pre_cemc_photon_index", &anchor.pre_cemc_photon);
  bind_anchor("partner_cemc_edep", &anchor.partner_cemc_edep);
  bind_anchor("partner_diagnostic_invariant_mass", &anchor.diagnostic_invariant_mass);
  anchors->SetBranchAddress("main_fraction", &anchor.main_fraction);
  anchors->SetBranchAddress("second_fraction", &anchor.second_fraction);
  anchors->SetBranchAddress("unmatched_max_fraction", &anchor.unmatched_fraction);
  anchors->SetBranchAddress("ambiguous_main", &anchor.ambiguous);
  anchors->SetBranchAddress("best_cluster0_id", &anchor.best_cluster[0]);
  anchors->SetBranchAddress("best_cluster1_id", &anchor.best_cluster[1]);
  anchors->SetBranchAddress("recovered0", &anchor.recovered[0]);
  anchors->SetBranchAddress("recovered1", &anchor.recovered[1]);
  anchors->SetBranchAddress("photon0_energy", &anchor.truth_energy[0]);
  anchors->SetBranchAddress("photon1_energy", &anchor.truth_energy[1]);
  anchors->SetBranchAddress("reconstructed_photon0_energy", &anchor.reconstructed[0]);
  anchors->SetBranchAddress("reconstructed_photon1_energy", &anchor.reconstructed[1]);
  anchors->SetBranchAddress("match_valid", &anchor.match_valid);
  bind_anchor("match_usable", &anchor.match_usable);
  bind_anchor("match_status", &anchor.match_status);
  bind_anchor("match_status_name", &match_status_name);
  bind_anchor("match_failure", &anchor.match_failure);
  bind_anchor("match_failure_name", &match_failure_name);
  bind_anchor("match_failure_ieta", &anchor.match_failure_ieta);
  bind_anchor("match_failure_iphi", &anchor.match_failure_iphi);
  bind_anchor("match_tower_count", &anchor.match_towers);
  bind_anchor("match_matched_tower_count", &anchor.match_matched_towers);
  bind_anchor("match_cluster_member_energy_coverage", &anchor.match_coverage);
  anchors->SetBranchAddress("total_edep", &anchor.total_edep);
  anchors->SetBranchAddress("gamma0_edep", &anchor.gamma_edep[0]);
  anchors->SetBranchAddress("gamma1_edep", &anchor.gamma_edep[1]);
  anchors->SetBranchAddress("other_edep", &anchor.other_edep);
  bind_anchor("partner_diagnostic_found", &anchor.diagnostic_found);
  bind_anchor("partner_diagnostic_below_energy_threshold", &anchor.diagnostic_below_threshold);
  bind_anchor("partner_diagnostic_has_direct_deposit", &anchor.diagnostic_has_direct);
  bind_anchor("partner_diagnostic_cluster_id", &anchor.diagnostic_cluster);
  bind_anchor("partner_diagnostic_cluster_energy", &anchor.diagnostic_energy);
  bind_anchor("partner_diagnostic_cluster_eta", &anchor.diagnostic_eta);
  bind_anchor("partner_diagnostic_cluster_phi", &anchor.diagnostic_phi);
  bind_anchor("partner_diagnostic_delta_r", &anchor.diagnostic_delta_r);
  bind_anchor("partner_diagnostic_reconstructed_energy", &anchor.diagnostic_reconstructed);
  bind_anchor("partner_diagnostic_recovery", &anchor.diagnostic_recovery);
  bind_anchor("partner_diagnostic_match_usable", &anchor.diagnostic_match_usable);
  bind_anchor("partner_diagnostic_match_status_name", &diagnostic_match_status_name);
  bind_anchor("partner_diagnostic_match_failure_name", &diagnostic_match_failure_name);
  bind_anchor("partner_diagnostic_failure_ieta", &anchor.diagnostic_failure_ieta);
  bind_anchor("partner_diagnostic_failure_iphi", &anchor.diagnostic_failure_iphi);
  bind_anchor("partner_diagnostic_match_coverage", &anchor.diagnostic_match_coverage);
  const bool has_match_usable = anchors->GetBranch("match_usable");
  for (Long64_t entry = 0; entry < anchors->GetEntries(); ++entry)
  {
    anchors->GetEntry(entry);
    if (event == event_id)
    {
      anchor.topology_name = topology_name ? *topology_name : "unknown";
      anchor.reason_name = reason_name ? *reason_name : "unknown";
      anchor.missing_category_name = missing_category_name ? *missing_category_name : "not_recorded";
      anchor.missing_detail_name = missing_detail_name ? *missing_detail_name : "not_recorded";
      anchor.match_status_name = match_status_name ? *match_status_name : (anchor.match_valid ? "complete" : "invalid");
      anchor.match_failure_name = match_failure_name ? *match_failure_name : "not_recorded";
      anchor.diagnostic_match_status_name = diagnostic_match_status_name ? *diagnostic_match_status_name : "not_recorded";
      anchor.diagnostic_match_failure_name = diagnostic_match_failure_name ? *diagnostic_match_failure_name : "not_recorded";
      if (!has_match_usable) anchor.match_usable = anchor.match_valid;
      data.anchors.push_back(anchor);
    }
  }
  anchors->ResetBranchAddresses();
  if (matches)
  {
    CandidateClusterMatch match;
    matches->SetBranchAddress("event", &event);
    matches->SetBranchAddress("candidate_id", &match.candidate);
    matches->SetBranchAddress("cluster_id", &match.cluster);
    if (matches->GetBranch("match_usable"))
      matches->SetBranchAddress("match_usable", &match.valid);
    else
      matches->SetBranchAddress("match_valid", &match.valid);
    matches->SetBranchAddress("gamma0_fraction", &match.gamma_fraction[0]);
    matches->SetBranchAddress("gamma1_fraction", &match.gamma_fraction[1]);
    for (Long64_t entry = 0; entry < matches->GetEntries(); ++entry)
    {
      matches->GetEntry(entry);
      if (event == event_id) data.matches.push_back(match);
    }
    matches->ResetBranchAddresses();
  }

  if (segments)
  {
    Segment segment;
    segments->SetBranchAddress("event", &event);
    segments->SetBranchAddress("track_id", &segment.track);
    segments->SetBranchAddress("pid", &segment.pid);
    segments->SetBranchAddress("family_candidate_id", &segment.family);
    segments->SetBranchAddress("family_gamma_index", &segment.gamma);
    segments->SetBranchAddress("x0", &segment.x0);
    segments->SetBranchAddress("y0", &segment.y0);
    segments->SetBranchAddress("z0", &segment.z0);
    segments->SetBranchAddress("x1", &segment.x1);
    segments->SetBranchAddress("y1", &segment.y1);
    segments->SetBranchAddress("z1", &segment.z1);
    for (Long64_t entry = 0; entry < segments->GetEntries(); ++entry)
    {
      segments->GetEntry(entry);
      if (event == event_id) data.segments.push_back(segment);
    }
    segments->ResetBranchAddresses();
  }

  if (clusters)
  {
    Cluster cluster;
    clusters->SetBranchAddress("event", &event);
    clusters->SetBranchAddress("cluster_id", &cluster.id);
    clusters->SetBranchAddress("energy", &cluster.energy);
    clusters->SetBranchAddress("et", &cluster.et);
    clusters->SetBranchAddress("x", &cluster.x);
    clusters->SetBranchAddress("y", &cluster.y);
    clusters->SetBranchAddress("z", &cluster.z);
    clusters->SetBranchAddress("r", &cluster.r);
    clusters->SetBranchAddress("eta", &cluster.eta);
    clusters->SetBranchAddress("phi", &cluster.phi);
    clusters->SetBranchAddress("ntowers", &cluster.ntowers);
    clusters->SetBranchAddress("topology_considered", &cluster.considered);
    for (Long64_t entry = 0; entry < clusters->GetEntries(); ++entry)
    {
      clusters->GetEntry(entry);
      if (event == event_id) data.clusters.push_back(cluster);
    }
    clusters->ResetBranchAddresses();
  }

  if (towers)
  {
    Tower tower;
    towers->SetBranchAddress("event", &event);
    towers->SetBranchAddress("cluster_id", &tower.cluster);
    towers->SetBranchAddress("ieta", &tower.ieta);
    towers->SetBranchAddress("iphi", &tower.iphi);
    towers->SetBranchAddress("cluster_tower_energy", &tower.member_energy);
    towers->SetBranchAddress("tower_energy", &tower.tower_energy);
    towers->SetBranchAddress("allocation_fraction", &tower.allocation);
    for (Long64_t entry = 0; entry < towers->GetEntries(); ++entry)
    {
      towers->GetEntry(entry);
      if (event == event_id) data.towers.push_back(tower);
    }
    towers->ResetBranchAddresses();
  }
  return true;
}

inline std::vector<int> event_ids(
    TFile* file, int topology_filter = -1, int pathway_filter = -1,
    double vertex_z_min = -std::numeric_limits<double>::infinity(),
    double vertex_z_max = std::numeric_limits<double>::infinity(),
    double truth_pi0_pt_min = -std::numeric_limits<double>::infinity(),
    double truth_pi0_pt_max = std::numeric_limits<double>::infinity())
{
  std::vector<int> result;
  if (!file || file->IsZombie() || vertex_z_min > vertex_z_max || truth_pi0_pt_min > truth_pi0_pt_max) return result;
  TTree* events = static_cast<TTree*>(file->Get("events"));
  TTree* anchors = static_cast<TTree*>(file->Get("anchor_decisions"));
  TTree* candidates = static_cast<TTree*>(file->Get("pi0_candidates"));
  if (!events || !anchors || !candidates) return result;
  struct CandidateFilterValues
  {
    int pathway = -1;
    double pt = -1.0;
  };
  std::map<std::pair<int, int>, CandidateFilterValues> candidate_values;
  std::set<int> candidate_accepted;
  int event = -1;
  int candidate = -1;
  int path = -1;
  double pt = -1.0;
  candidates->SetBranchAddress("event", &event);
  candidates->SetBranchAddress("candidate_id", &candidate);
  candidates->SetBranchAddress("pathway", &path);
  candidates->SetBranchAddress("pt", &pt);
  for (Long64_t entry = 0; entry < candidates->GetEntries(); ++entry)
  {
    candidates->GetEntry(entry);
    candidate_values[{event, candidate}] = {path, pt};
    const bool pathway_ok = pathway_filter < 0 || path == pathway_filter;
    const bool pt_ok = pt >= truth_pi0_pt_min && pt <= truth_pi0_pt_max;
    if (pathway_ok && pt_ok) candidate_accepted.insert(event);
  }
  candidates->ResetBranchAddresses();
  std::set<int> topology_accepted;
  int topology = -1;
  anchors->SetBranchAddress("event", &event);
  anchors->SetBranchAddress("candidate_id", &candidate);
  anchors->SetBranchAddress("topology", &topology);
  for (Long64_t entry = 0; entry < anchors->GetEntries(); ++entry)
  {
    anchors->GetEntry(entry);
    const bool topology_ok = topology_filter < 0 || topology == topology_filter;
    const auto found = candidate_values.find({event, candidate});
    const bool pathway_ok = pathway_filter < 0 ||
        (found != candidate_values.end() && found->second.pathway == pathway_filter);
    const bool pt_ok = found != candidate_values.end() &&
        found->second.pt >= truth_pi0_pt_min && found->second.pt <= truth_pi0_pt_max;
    if (topology_ok && pathway_ok && pt_ok) topology_accepted.insert(event);
  }
  anchors->ResetBranchAddresses();
  double collision_z = 0.0;
  events->SetBranchAddress("event", &event);
  events->SetBranchAddress("collision_z", &collision_z);
  const bool candidate_filter_enabled = pathway_filter >= 0 ||
      truth_pi0_pt_min != -std::numeric_limits<double>::infinity() ||
      truth_pi0_pt_max != std::numeric_limits<double>::infinity();
  for (Long64_t entry = 0; entry < events->GetEntries(); ++entry)
  {
    events->GetEntry(entry);
    const bool vertex_ok = collision_z >= vertex_z_min && collision_z <= vertex_z_max;
    const bool candidate_topology_ok = topology_filter >= 0
        ? topology_accepted.count(event) != 0U
        : (!candidate_filter_enabled || candidate_accepted.count(event) != 0U);
    if (vertex_ok && candidate_topology_ok) result.push_back(event);
  }
  events->ResetBranchAddresses();
  return result;
}

inline const Candidate* find_candidate(const DisplayData& data, int id)
{
  for (const auto& candidate : data.candidates)
    if (candidate.id == id) return &candidate;
  return nullptr;
}

inline int cluster_topology(const DisplayData& data, unsigned int id)
{
  for (const auto& anchor : data.anchors)
    if (anchor.cluster == id) return anchor.topology;
  return -1;
}

inline int cluster_family(const DisplayData& data, unsigned int id,
                          int selected_family = -1)
{
  int family = -1;
  double maximum_fraction = 0.5;
  for (const auto& match : data.matches)
  {
    if (!match.valid || match.cluster != id ||
        (selected_family >= 0 && match.candidate != selected_family))
      continue;
    const double fraction = match.gamma_fraction[0] + match.gamma_fraction[1];
    if (fraction >= maximum_fraction)
    {
      family = match.candidate;
      maximum_fraction = fraction;
    }
  }
  return family;
}

inline void draw_cluster_marker(const DisplayData& data, const Cluster& cluster,
                                double first, double second,
                                double marker_size,
                                int selected_family = -1)
{
  const int topology = cluster_topology(data, cluster.id);
  auto* marker = new TMarker(first, second, 20);
  marker->SetMarkerColor(topology >= 0 ? topology_color(topology) : kGray + 2);
  marker->SetMarkerSize(marker_size);
  marker->Draw();
  const int family = cluster_family(data, cluster.id, selected_family);
  if (family >= 0)
  {
    auto* family_ring = new TMarker(first, second, 24);
    family_ring->SetMarkerColor(family_color(family));
    family_ring->SetMarkerSize(marker_size + 0.35);
    family_ring->Draw();
  }
}

inline void draw_xy(const DisplayData& data, int selected_family = -1)
{
  auto* frame = new TH2F(Form("xy_%d_%d", data.event.id, selected_family),
      ";x [cm];y [cm]", 100, -130.0, 130.0, 100, -130.0, 130.0);
  frame->SetDirectory(nullptr);
  frame->Draw();
  auto* outer = new TEllipse(0.0, 0.0, 113.0, 113.0);
  outer->SetFillStyle(0); outer->SetLineColor(kGray + 2); outer->SetLineStyle(2); outer->Draw();
  auto* inner = new TEllipse(0.0, 0.0, 93.0, 93.0);
  inner->SetFillStyle(0); inner->SetLineColor(kGray + 2); inner->SetLineStyle(2); inner->Draw();
  for (const auto& segment : data.segments)
  {
    if (selected_family >= 0 && segment.family != selected_family) continue;
    auto* line = new TLine(segment.x0, segment.y0, segment.x1, segment.y1);
    line->SetLineColor(family_color(segment.family));
    line->SetLineWidth(segment.family >= 0 && is_core_segment(data, segment) ? 2 : 1);
    line->SetLineStyle(segment.family >= 0 ? family_line_style(segment.gamma) : 3);
    line->Draw();
  }
  for (const auto& candidate : data.candidates)
  {
    if (selected_family >= 0 && candidate.id != selected_family) continue;
    auto* direction = new TLine(data.event.cx, data.event.cy,
        data.event.cx + 100.0 * std::cos(candidate.phi),
        data.event.cy + 100.0 * std::sin(candidate.phi));
    direction->SetLineColor(family_color(candidate.id));
    direction->SetLineStyle(1); direction->SetLineWidth(2); direction->Draw();
  }
  for (const auto& cluster : data.clusters)
  {
    draw_cluster_marker(data, cluster, cluster.x, cluster.y,
        cluster_marker_size(cluster.energy, 2.2),
        selected_family);
  }
}

inline void draw_zr(const DisplayData& data, int selected_family = -1)
{
  auto* frame = new TH2F(Form("zr_%d_%d", data.event.id, selected_family),
      ";z [cm];r [cm]", 120, -180.0, 180.0, 100, 0.0, 140.0);
  frame->SetDirectory(nullptr); frame->Draw();
  for (const double radius : {93.0, 113.0})
  {
    auto* surface = new TLine(-180.0, radius, 180.0, radius);
    surface->SetLineColor(kGray + 2); surface->SetLineStyle(2); surface->Draw();
  }
  for (const auto& segment : data.segments)
  {
    if (selected_family >= 0 && segment.family != selected_family) continue;
    auto* line = new TLine(segment.z0, std::hypot(segment.x0, segment.y0),
                           segment.z1, std::hypot(segment.x1, segment.y1));
    line->SetLineColor(family_color(segment.family));
    line->SetLineWidth(segment.family >= 0 && is_core_segment(data, segment) ? 2 : 1);
    line->SetLineStyle(segment.family >= 0 ? family_line_style(segment.gamma) : 3); line->Draw();
  }
  for (const auto& candidate : data.candidates)
  {
    if (selected_family >= 0 && candidate.id != selected_family) continue;
    auto* direction = new TLine(data.event.cz, std::hypot(data.event.cx, data.event.cy),
                                data.event.cz + 100.0 * std::sinh(candidate.eta), 100.0);
    direction->SetLineColor(family_color(candidate.id)); direction->SetLineStyle(1);
    direction->SetLineWidth(2); direction->Draw();
  }
  for (const auto& cluster : data.clusters)
  {
    draw_cluster_marker(data, cluster, cluster.z, cluster.r,
        cluster_marker_size(cluster.energy, 2.2),
        selected_family);
  }
}

inline void draw_eta_phi(const DisplayData& data, int selected_family = -1)
{
  auto* frame = new TH2F(Form("etaphi_%d_%d", data.event.id, selected_family),
      ";#eta;#phi", 100, -1.2, 1.2, 128, -3.2, 3.2);
  frame->SetDirectory(nullptr); frame->Draw();
  for (const auto& cluster : data.clusters)
  {
    if (!std::isfinite(cluster.eta) || !std::isfinite(cluster.phi)) continue;
    draw_cluster_marker(data, cluster, cluster.eta, cluster.phi,
        cluster_marker_size(cluster.energy, 2.0),
        selected_family);
  }
  for (const auto& candidate : data.candidates)
  {
    if (selected_family >= 0 && candidate.id != selected_family) continue;
    for (int gamma = 0; gamma < 2; ++gamma)
    {
      const double eta = candidate.photon_eta[gamma];
      const double phi = candidate.photon_phi[gamma];
      auto* marker = new TMarker(eta, phi, 29);
      marker->SetMarkerColor(family_color(candidate.id));
      marker->SetMarkerSize(1.6); marker->Draw();
    }
  }
}

inline int unwrap_iphi(int iphi, int reference)
{
  int result = iphi;
  while (result - reference > 128) result -= 256;
  while (result - reference < -128) result += 256;
  return result;
}

inline void draw_anchor_towers(const DisplayData& data, const Anchor& anchor)
{
  std::set<unsigned int> roles = {anchor.cluster};
  if (anchor.best_cluster[0] >= 0) roles.insert(static_cast<unsigned int>(anchor.best_cluster[0]));
  if (anchor.best_cluster[1] >= 0) roles.insert(static_cast<unsigned int>(anchor.best_cluster[1]));
  int reference = 0;
  bool have_reference = false;
  for (const auto& tower : data.towers)
    if (tower.cluster == anchor.cluster) { reference = tower.iphi; have_reference = true; break; }
  if (!have_reference)
  {
    auto* frame = new TH2F(Form("empty_towers_%d_%d", data.event.id, anchor.id),
        "No detailed tower data;tower i#eta;tower i#phi", 10, 0, 10, 10, 0, 10);
    frame->SetDirectory(nullptr); frame->Draw(); return;
  }
  int min_eta = 96, max_eta = -1, min_phi = reference, max_phi = reference;
  std::map<std::pair<int, int>, double> tower_energy;
  for (const auto& tower : data.towers)
  {
    if (!roles.count(tower.cluster)) continue;
    const int phi = unwrap_iphi(tower.iphi, reference);
    min_eta = std::min(min_eta, tower.ieta); max_eta = std::max(max_eta, tower.ieta);
    min_phi = std::min(min_phi, phi); max_phi = std::max(max_phi, phi);
    auto& stored = tower_energy[{tower.ieta, phi}];
    stored = std::max(stored, tower.tower_energy);
  }
  min_eta -= 3; max_eta += 3; min_phi -= 3; max_phi += 3;
  auto* heat = new TH2F(Form("tower_heat_%d_%d", data.event.id, anchor.id),
      ";tower i#eta;tower i#phi;calibrated tower E [GeV]",
      max_eta - min_eta + 1, min_eta - 0.5, max_eta + 0.5,
      max_phi - min_phi + 1, min_phi - 0.5, max_phi + 0.5);
  heat->SetDirectory(nullptr);
  for (const auto& [index, energy] : tower_energy)
    heat->Fill(index.first, index.second, energy);
  heat->Draw("colz");

  const auto draw_outline = [&](int cluster_id, double half_width, int color, int style, int width)
  {
    if (cluster_id < 0) return;
    for (const auto& tower : data.towers)
    {
      if (static_cast<int>(tower.cluster) != cluster_id) continue;
      const int phi = unwrap_iphi(tower.iphi, reference);
      auto* box = new TBox(tower.ieta - half_width, phi - half_width,
                           tower.ieta + half_width, phi + half_width);
      box->SetFillStyle(0); box->SetLineColor(color); box->SetLineStyle(style); box->SetLineWidth(width); box->Draw();
    }
  };
  draw_outline(static_cast<int>(anchor.cluster), 0.47, topology_color(anchor.topology), 1, 3);
  draw_outline(anchor.best_cluster[0], 0.36, family_color(anchor.candidate), photon0_line_style, 2);
  draw_outline(anchor.best_cluster[1], 0.25, family_color(anchor.candidate), photon1_line_style, 2);
}

inline int anchor_display_role(const DisplayData& data, const Anchor& anchor)
{
  const Candidate* candidate = find_candidate(data, anchor.candidate);
  if (!candidate) return 2;
  if (static_cast<int>(anchor.cluster) == candidate->best_cluster[0]) return 0;
  if (static_cast<int>(anchor.cluster) == candidate->best_cluster[1]) return 1;
  return 2;
}

inline std::vector<const Anchor*> ordered_anchors(const DisplayData& data)
{
  std::vector<const Anchor*> result;
  result.reserve(data.anchors.size());
  for (const auto& anchor : data.anchors) result.push_back(&anchor);
  std::stable_sort(result.begin(), result.end(), [&](const Anchor* left, const Anchor* right)
  {
    return std::make_tuple(left->candidate, anchor_display_role(data, *left), left->id) <
        std::make_tuple(right->candidate, anchor_display_role(data, *right), right->id);
  });
  return result;
}

inline void draw_legend_line(double y, int color, int style, int width, const char* label, double text_size)
{
  auto* line = new TLine(0.06, y, 0.25, y);
  line->SetLineColor(color); line->SetLineStyle(style); line->SetLineWidth(width); line->Draw();
  TLatex text; text.SetTextSize(text_size); text.SetTextAlign(12); text.DrawLatex(0.30, y, label);
}

inline void draw_display_legend(int family = -1, bool tower_detail = false,
                                double start_y = 0.92, double step = 0.075,
                                double text_size = 0.038)
{
  gPad->Range(0.0, 0.0, 1.0, 1.0);
  const int color = family >= 0 ? family_color(family) : kAzure + 1;
  TLatex title; title.SetTextSize(text_size + 0.012); title.SetTextFont(62); title.DrawLatex(0.05, start_y, "Legend");
  double y = start_y - step;
  draw_legend_line(y, color, 1, 3, "#pi^{0} direction / segment", text_size); y -= step;
  draw_legend_line(y, color, photon0_line_style, 2, "photon: fine dash", text_size); y -= step;
  draw_legend_line(y, color, photon1_line_style, 2, "partner photon: coarse dash", text_size); y -= step;
  draw_legend_line(y, color, photon0_line_style, 1, "pre-CEMC descendants: thin", text_size); y -= step;
  if (family < 0)
  {
    draw_legend_line(y, kGray + 1, 3, 1, "other G4 truth segment", text_size);
  }
  else
  {
    TLatex note; note.SetTextSize(text_size); note.SetTextAlign(12); note.DrawLatex(0.06, y, "other G4 truth segments hidden");
  }
  y -= step;

  auto* star = new TMarker(0.15, y, 29);
  star->SetMarkerColor(color); star->SetMarkerSize(1.4); star->Draw();
  TLatex text; text.SetTextSize(text_size); text.SetTextAlign(12); text.DrawLatex(0.30, y, "truth photon momentum in #eta-#phi"); y -= step;

  auto* low_energy_cluster = new TMarker(0.09, y, 20);
  low_energy_cluster->SetMarkerColor(kGray + 2); low_energy_cluster->SetMarkerSize(0.3); low_energy_cluster->Draw();
  auto* cluster = new TMarker(0.18, y, 20);
  cluster->SetMarkerColor(kGray + 2); cluster->SetMarkerSize(1.2); cluster->Draw();
  auto* ring = new TMarker(0.18, y, 24);
  ring->SetMarkerColor(color); ring->SetMarkerSize(1.55); ring->Draw();
  text.DrawLatex(0.28, y, "cluster size: E<0.1 small, E#geq0.1 large"); y -= step;

  draw_legend_line(y, kGray + 2, 2, 1, "CEMC inner / outer surface", text_size); y -= step;
  text.SetTextColor(kGreen + 2); text.DrawLatex(0.05, y, "S separated");
  text.SetTextColor(kMagenta + 1); text.DrawLatex(0.24, y, "M merged");
  text.SetTextColor(kAzure + 2); text.DrawLatex(0.43, y, "C contaminated");
  text.SetTextColor(kOrange + 7); text.DrawLatex(0.67, y, "X missing");
  text.SetTextColor(kGray + 2); text.DrawLatex(0.84, y, "O other");
  text.SetTextColor(kBlack);

  if (!tower_detail) return;
  y -= step;
  text.DrawLatex(0.05, y, "Tower map: heat = calibrated tower E"); y -= step;
  draw_legend_line(y, topology_color(1), 1, 3, "anchor cluster outline", text_size); y -= step;
  draw_legend_line(y, color, photon0_line_style, 2, "best-cluster outline: fine dash", text_size); y -= step;
  draw_legend_line(y, color, photon1_line_style, 2, "best-cluster outline: coarse dash", text_size);
}

inline void draw_overview_text(const DisplayData& data)
{
  gPad->Range(0.0, 0.0, 1.0, 1.0);
  TLatex text; text.SetNDC(); text.SetTextSize(0.043);
  double y = 0.95;
  text.DrawLatex(0.04, y, Form("vertex = (%.2f, %.2f, %.2f) cm", data.event.cx, data.event.cy, data.event.cz)); y -= 0.07;
  text.DrawLatex(0.04, y, Form("selected #pi^{0}: %d   G4-primary: %d   generator: %d",
      data.event.candidates, data.event.g4_primary, data.event.generator)); y -= 0.065;
  text.DrawLatex(0.04, y, Form("clusters: %d   anchors: %d   S/M/C/X/O = %d/%d/%d/%d/%d",
      data.event.clusters, data.event.anchors, data.event.separated,
      data.event.merged, data.event.single_contaminated, data.event.missing, data.event.other)); y -= 0.075;

  const auto anchors = ordered_anchors(data);
  text.SetTextSize(0.038);
  std::size_t shown = 0;
  for (const Anchor* anchor : anchors)
  {
    if (y < 0.57) break;
    text.SetTextColor(topology_color(anchor->topology));
    text.DrawLatex(0.04, y, Form("P%d  A%d  C%u  %-9s  E_{T}=%.2f  f=%.3f",
        anchor->candidate, anchor->id, anchor->cluster, anchor->topology_name.c_str(), anchor->et, anchor->main_fraction));
    text.SetTextColor(kBlack); y -= 0.052; ++shown;
  }
  if (shown < anchors.size()) text.DrawLatex(0.04, y, Form("... %zu more anchors", anchors.size() - shown));
  draw_display_legend(-1, false, 0.48, 0.052, 0.033);
}

inline void draw_anchor_text(const DisplayData& data, const Anchor& anchor)
{
  gPad->Range(0.0, 0.0, 1.0, 1.0);
  const Candidate* candidate = find_candidate(data, anchor.candidate);
  TLatex text; text.SetNDC(); text.SetTextSize(0.040);
  double y = 0.96;
  text.SetTextFont(62); text.DrawLatex(0.04, y, "Anchor classification"); text.SetTextFont(42); y -= 0.055;
  text.SetTextColor(topology_color(anchor.topology));
  text.DrawLatex(0.04, y, Form("%s", anchor.topology_name.c_str()));
  text.SetTextColor(kBlack); y -= 0.050;
  text.SetTextSize(0.030); text.DrawLatex(0.04, y, Form("reason: %s", anchor.reason_name.c_str())); y -= 0.045;
  if (anchor.topology == 4)
  {
    text.DrawLatex(0.04, y, Form("pre-CEMC interacting photon: #gamma%d", anchor.pre_cemc_photon)); y -= 0.045;
  }
  if (anchor.topology == 3)
  {
    text.DrawLatex(0.04, y, Form("missing: %s / %s  (partner #gamma%d)", anchor.missing_category_name.c_str(), anchor.missing_detail_name.c_str(), anchor.partner_photon)); y -= 0.045;
    text.DrawLatex(0.04, y, Form("partner CEMC truth Edep = %.4g GeV", anchor.partner_cemc_edep)); y -= 0.045;
    if (candidate && anchor.partner_photon >= 0 && anchor.partner_photon < 2)
    {
      const int partner = anchor.partner_photon;
      if (candidate->photon_projection_valid[partner])
      {
        text.DrawLatex(0.04, y, Form("partner projection: #eta=%.4f  #phi=%.4f  in acceptance=%s", candidate->photon_projection_eta[partner],
            candidate->photon_projection_phi[partner], candidate->photon_in_cemc_acceptance[partner] ? "yes" : "no")); y -= 0.045;
      }
      else
      {
        text.DrawLatex(0.04, y, "partner projection: invalid"); y -= 0.045;
      }
    }
  }
  text.SetTextSize(0.036);
  text.DrawLatex(0.04, y, Form("E = %.3f GeV    E_{T} = %.3f GeV", anchor.energy, anchor.et)); y -= 0.048;
  text.DrawLatex(0.04, y, Form("main/second/unmatched f = %.4f / %.4f / %.4f",
      anchor.main_fraction, anchor.second_fraction, anchor.unmatched_fraction)); y -= 0.044;
  if (candidate)
  {
    text.DrawLatex(0.04, y, Form("candidate P%d    p_{T} = %.3f GeV", candidate->id, candidate->pt)); y -= 0.044;
    text.SetTextSize(0.030); text.DrawLatex(0.04, y, candidate->pathway_name.c_str()); y -= 0.050;
  }
  text.SetTextSize(0.032);
  for (int photon = 0; photon < 2; ++photon)
  {
    const double recovery = anchor.truth_energy[photon] > 0.0 ? anchor.reconstructed[photon] / anchor.truth_energy[photon] : 0.0;
    const int style = photon == 0 ? photon0_line_style : photon1_line_style;
    auto* sample = new TLine(0.04, y + 0.008, 0.16, y + 0.008);
    sample->SetLineColor(family_color(anchor.candidate)); sample->SetLineStyle(style); sample->SetLineWidth(2); sample->Draw();
    text.DrawLatex(0.19, y, Form("#gamma%d Cbest=%d  recovered=%s", photon, anchor.best_cluster[photon],
        anchor.recovered[photon] ? "yes" : "no")); y -= 0.038;
    text.DrawLatex(0.19, y, Form("Etruth=%.3f  Erec=%.3f  ratio=%.3f",
        anchor.truth_energy[photon], anchor.reconstructed[photon], recovery)); y -= 0.038;
    if (candidate->photon_first_daughter_vertex_valid[photon])
      text.DrawLatex(0.19, y, Form("first daughter r=%.3f cm  pre-CEMC=%s", candidate->photon_first_daughter_radius[photon],
          candidate->photon_pre_cemc_interaction[photon] ? "yes" : "no"));
    else
      text.DrawLatex(0.19, y, "first daughter vertex: not stored");
    y -= 0.043;
  }
  text.DrawLatex(0.04, y, Form("anchor match: strict=%d usable=%d  %s  coverage=%.3f",
      anchor.match_valid, anchor.match_usable, anchor.match_status_name.c_str(), anchor.match_coverage)); y -= 0.038;
  text.DrawLatex(0.04, y, Form("towers=%d/%d  failure=%s at (%d,%d)",
      anchor.match_matched_towers, anchor.match_towers, anchor.match_failure_name.c_str(),
      anchor.match_failure_ieta, anchor.match_failure_iphi)); y -= 0.038;
  text.DrawLatex(0.04, y, Form("anchor edep: total=%.4g  other=%.4g", anchor.total_edep, anchor.other_edep)); y -= 0.038;
  text.DrawLatex(0.04, y, Form("fine/coarse photon edep = %.4g / %.4g", anchor.gamma_edep[0], anchor.gamma_edep[1])); y -= 0.038;
  if (anchor.diagnostic_found)
  {
    text.DrawLatex(0.04, y, Form("partner diagnostic C%d: E=%.3f  #DeltaR=%.4f  below=%d direct=%d",
        anchor.diagnostic_cluster, anchor.diagnostic_energy, anchor.diagnostic_delta_r,
        anchor.diagnostic_below_threshold, anchor.diagnostic_has_direct)); y -= 0.038;
    text.DrawLatex(0.04, y, Form("Erec=%.3f  recovery=%.3f  match=%s usable=%d coverage=%.3f",
        anchor.diagnostic_reconstructed, anchor.diagnostic_recovery,
        anchor.diagnostic_match_status_name.c_str(), anchor.diagnostic_match_usable,
        anchor.diagnostic_match_coverage)); y -= 0.038;
    if (anchor.diagnostic_invariant_mass >= 0.0)
    {
      text.DrawLatex(0.04, y, Form("anchor+diagnostic m_{#gamma#gamma}=%.4f GeV", anchor.diagnostic_invariant_mass)); y -= 0.038;
    }
    text.DrawLatex(0.04, y, Form("diagnostic failure=%s at (%d,%d)",
        anchor.diagnostic_match_failure_name.c_str(), anchor.diagnostic_failure_ieta, anchor.diagnostic_failure_iphi));
  }
}

inline TPad* make_display_pad(TCanvas* canvas, const std::string& name,
                              double x0, double y0, double x1, double y1,
                              int fill_color = kWhite)
{
  canvas->cd();
  auto* pad = new TPad(name.c_str(), "", x0, y0, x1, y1);
  pad->SetFillColor(fill_color); pad->SetBorderMode(0); pad->Draw();
  return pad;
}

inline void draw_page_header(TCanvas* canvas, const DisplayData& data,
                             const Anchor* anchor = nullptr,
                             std::size_t anchor_position = 0,
                             std::size_t anchor_count = 0)
{
  auto* header = make_display_pad(canvas, Form("header_%d_%d", data.event.id, anchor ? anchor->id : -1),
      0.0, 0.93, 1.0, 1.0, anchor ? kOrange - 9 : kAzure - 9);
  header->cd(); header->Range(0.0, 0.0, 1.0, 1.0);
  TLatex text; text.SetTextAlign(12); text.SetTextFont(62); text.SetTextSize(anchor ? 0.30 : 0.34);
  if (!anchor)
  {
    text.DrawLatex(0.02, 0.50, Form("EVENT OVERVIEW   |   Event %d", data.event.id));
    return;
  }
  text.DrawLatex(0.02, 0.50, Form("ANCHOR DETAIL %zu/%zu   |   Event %d   |   P%d   A%d   C%u",
      anchor_position, anchor_count, data.event.id, anchor->candidate, anchor->id, anchor->cluster));
  auto* badge = new TBox(0.83, 0.13, 0.98, 0.87);
  badge->SetFillColor(topology_color(anchor->topology)); badge->SetLineColor(topology_color(anchor->topology)); badge->Draw();
  text.SetTextAlign(22); text.SetTextColor(kWhite); text.SetTextSize(0.25);
  text.DrawLatex(0.905, 0.50, anchor->topology_name.c_str());
}

inline void configure_plot_pad(TPad* pad, bool color_bar = false)
{
  pad->SetLeftMargin(0.13); pad->SetBottomMargin(0.13); pad->SetTopMargin(0.08);
  pad->SetRightMargin(color_bar ? 0.17 : 0.05);
}

inline std::vector<std::unique_ptr<TCanvas>> make_event_pages(TFile* file, int event_id)
{
  DisplayData data;
  std::vector<std::unique_ptr<TCanvas>> pages;
  if (!load_event(file, event_id, data)) return pages;
  gStyle->SetOptStat(0);
  gStyle->SetLineStyleString(photon0_line_style, "4 4");
  gStyle->SetLineStyleString(photon1_line_style, "16 8");

  auto overview = std::make_unique<TCanvas>(Form("overview_%d", event_id), "topology event overview", 1800, 1100);
  draw_page_header(overview.get(), data);
  auto* overview_xy = make_display_pad(overview.get(), Form("overview_xy_%d", event_id), 0.00, 0.49, 0.50, 0.93);
  auto* overview_zr = make_display_pad(overview.get(), Form("overview_zr_%d", event_id), 0.50, 0.49, 1.00, 0.93);
  auto* overview_eta_phi = make_display_pad(overview.get(), Form("overview_eta_phi_%d", event_id), 0.00, 0.00, 0.50, 0.49);
  auto* overview_text = make_display_pad(overview.get(), Form("overview_text_%d", event_id), 0.50, 0.00, 1.00, 0.49);
  configure_plot_pad(overview_xy); configure_plot_pad(overview_zr); configure_plot_pad(overview_eta_phi);
  overview_xy->cd(); draw_xy(data);
  overview_zr->cd(); draw_zr(data);
  overview_eta_phi->cd(); draw_eta_phi(data);
  overview_text->cd(); draw_overview_text(data);
  pages.push_back(std::move(overview));

  const auto anchors = ordered_anchors(data);
  for (std::size_t position = 0; position < anchors.size(); ++position)
  {
    const Anchor& anchor = *anchors[position];
    auto page = std::make_unique<TCanvas>(Form("anchor_%d_%d", event_id, anchor.id), "topology anchor detail", 1800, 1100);
    draw_page_header(page.get(), data, &anchor, position + 1, anchors.size());
    auto* xy = make_display_pad(page.get(), Form("anchor_xy_%d_%d", event_id, anchor.id), 0.00, 0.49, 0.333, 0.93);
    auto* zr = make_display_pad(page.get(), Form("anchor_zr_%d_%d", event_id, anchor.id), 0.333, 0.49, 0.666, 0.93);
    auto* eta_phi = make_display_pad(page.get(), Form("anchor_eta_phi_%d_%d", event_id, anchor.id), 0.666, 0.49, 1.00, 0.93);
    auto* towers = make_display_pad(page.get(), Form("anchor_towers_%d_%d", event_id, anchor.id), 0.00, 0.00, 0.42, 0.49);
    auto* info = make_display_pad(page.get(), Form("anchor_info_%d_%d", event_id, anchor.id), 0.42, 0.00, 0.72, 0.49);
    auto* legend = make_display_pad(page.get(), Form("anchor_legend_%d_%d", event_id, anchor.id), 0.72, 0.00, 1.00, 0.49);
    configure_plot_pad(xy); configure_plot_pad(zr); configure_plot_pad(eta_phi); configure_plot_pad(towers, true);
    xy->cd(); draw_xy(data, anchor.candidate);
    zr->cd(); draw_zr(data, anchor.candidate);
    eta_phi->cd(); draw_eta_phi(data, anchor.candidate);
    towers->cd(); draw_anchor_towers(data, anchor);
    info->cd(); draw_anchor_text(data, anchor);
    legend->cd(); draw_display_legend(anchor.candidate, true, 0.94, 0.068, 0.035);
    pages.push_back(std::move(page));
  }
  return pages;
}

inline bool print_event_pages(TFile* file, int event_id,
                              const std::string& output_pdf,
                              bool open_pdf, bool close_pdf)
{
  auto pages = make_event_pages(file, event_id);
  if (pages.empty()) return false;
  if (open_pdf) pages.front()->Print((output_pdf + "[").c_str());
  for (auto& page : pages) page->Print(output_pdf.c_str());
  if (close_pdf) pages.back()->Print((output_pdf + "]").c_str());
  return true;
}
}

#endif
