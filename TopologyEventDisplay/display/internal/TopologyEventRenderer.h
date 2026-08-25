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
  double main_fraction = -1.0;
  double second_fraction = -1.0;
  double unmatched_fraction = 0.0;
  int ambiguous = 0;
  int best_cluster[2] = {-999, -999};
  int recovered[2] = {0, 0};
  double truth_energy[2] = {0.0, 0.0};
  double reconstructed[2] = {0.0, 0.0};
  int match_valid = 0;
  double total_edep = 0.0;
  double gamma_edep[2] = {0.0, 0.0};
  double other_edep = 0.0;
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
  return kGray + 2;
}

inline int family_color(int family, int gamma = -1)
{
  static const int colors[] = {kAzure + 1, kRed + 1, kViolet + 1,
                               kTeal + 2, kOrange + 1, kPink + 7};
  if (family < 0) return kGray + 1;
  int color = colors[family % 6];
  if (gamma == 1) color += 2;
  return color;
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
  std::string* topology_name = nullptr;
  std::string* reason_name = nullptr;
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
  anchors->SetBranchAddress("total_edep", &anchor.total_edep);
  anchors->SetBranchAddress("gamma0_edep", &anchor.gamma_edep[0]);
  anchors->SetBranchAddress("gamma1_edep", &anchor.gamma_edep[1]);
  anchors->SetBranchAddress("other_edep", &anchor.other_edep);
  for (Long64_t entry = 0; entry < anchors->GetEntries(); ++entry)
  {
    anchors->GetEntry(entry);
    if (event == event_id)
    {
      anchor.topology_name = topology_name ? *topology_name : "unknown";
      anchor.reason_name = reason_name ? *reason_name : "unknown";
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
    if (selected_family >= 0 && segment.family >= 0 && segment.family != selected_family) continue;
    auto* line = new TLine(segment.x0, segment.y0, segment.x1, segment.y1);
    line->SetLineColor(family_color(segment.family, segment.gamma));
    line->SetLineWidth(segment.family >= 0 ? 2 : 1);
    line->SetLineStyle(segment.family >= 0 ? 1 : 3);
    line->Draw();
  }
  for (const auto& candidate : data.candidates)
  {
    if (selected_family >= 0 && candidate.id != selected_family) continue;
    auto* direction = new TLine(data.event.cx, data.event.cy,
        data.event.cx + 100.0 * std::cos(candidate.phi),
        data.event.cy + 100.0 * std::sin(candidate.phi));
    direction->SetLineColor(family_color(candidate.id));
    direction->SetLineStyle(2); direction->SetLineWidth(2); direction->Draw();
  }
  for (const auto& cluster : data.clusters)
  {
    draw_cluster_marker(data, cluster, cluster.x, cluster.y,
        std::min(2.2, 0.7 + 0.12 * std::max(0.0, cluster.energy)),
        selected_family);
  }
}

inline void draw_zr(const DisplayData& data, int selected_family = -1)
{
  auto* frame = new TH2F(Form("zr_%d_%d", data.event.id, selected_family),
      ";z [cm];r [cm]", 120, -180.0, 180.0, 100, 0.0, 140.0);
  frame->SetDirectory(nullptr); frame->Draw();
  for (const auto& segment : data.segments)
  {
    if (selected_family >= 0 && segment.family >= 0 && segment.family != selected_family) continue;
    auto* line = new TLine(segment.z0, std::hypot(segment.x0, segment.y0),
                           segment.z1, std::hypot(segment.x1, segment.y1));
    line->SetLineColor(family_color(segment.family, segment.gamma));
    line->SetLineWidth(segment.family >= 0 ? 2 : 1);
    line->SetLineStyle(segment.family >= 0 ? 1 : 3); line->Draw();
  }
  for (const auto& candidate : data.candidates)
  {
    if (selected_family >= 0 && candidate.id != selected_family) continue;
    auto* direction = new TLine(data.event.cz, std::hypot(data.event.cx, data.event.cy),
                                data.event.cz + 100.0 * std::sinh(candidate.eta), 100.0);
    direction->SetLineColor(family_color(candidate.id)); direction->SetLineStyle(2);
    direction->SetLineWidth(2); direction->Draw();
  }
  for (const auto& cluster : data.clusters)
  {
    draw_cluster_marker(data, cluster, cluster.z, cluster.r,
        std::min(2.2, 0.7 + 0.12 * std::max(0.0, cluster.energy)),
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
        std::min(2.0, 0.7 + 0.12 * std::max(0.0, cluster.energy)),
        selected_family);
  }
  for (const auto& candidate : data.candidates)
  {
    if (selected_family >= 0 && candidate.id != selected_family) continue;
    for (int gamma = 0; gamma < 2; ++gamma)
    {
      auto* marker = new TMarker(candidate.photon_eta[gamma], candidate.photon_phi[gamma], 29);
      marker->SetMarkerColor(family_color(candidate.id, gamma));
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
  for (const auto& tower : data.towers)
  {
    if (!roles.count(tower.cluster)) continue;
    const int phi = unwrap_iphi(tower.iphi, reference);
    double half_width = 0.47;
    if (tower.cluster != anchor.cluster &&
        static_cast<int>(tower.cluster) == anchor.best_cluster[0])
      half_width = 0.36;
    else if (tower.cluster != anchor.cluster &&
             static_cast<int>(tower.cluster) == anchor.best_cluster[1])
      half_width = 0.25;
    auto* box = new TBox(tower.ieta - half_width, phi - half_width,
                         tower.ieta + half_width, phi + half_width);
    box->SetFillStyle(0);
    int color = kBlack;
    if (tower.cluster == anchor.cluster) color = topology_color(anchor.topology);
    else if (static_cast<int>(tower.cluster) == anchor.best_cluster[0]) color = kAzure + 1;
    else if (static_cast<int>(tower.cluster) == anchor.best_cluster[1]) color = kRed + 1;
    box->SetLineColor(color); box->SetLineWidth(2); box->Draw();
  }
}

inline void draw_overview_text(const DisplayData& data)
{
  TLatex text; text.SetNDC(); text.SetTextSize(0.035);
  double y = 0.94;
  text.DrawLatex(0.04, y, Form("Event %d  vertex=(%.2f, %.2f, %.2f) cm",
      data.event.id, data.event.cx, data.event.cy, data.event.cz)); y -= 0.07;
  text.DrawLatex(0.04, y, Form("selected #pi^{0}: %d  G4-primary: %d  generator: %d",
      data.event.candidates, data.event.g4_primary, data.event.generator)); y -= 0.06;
  text.DrawLatex(0.04, y, Form("G4-secondary #pi^{0}: %d  selected-family particles: %d / %d",
      data.event.secondary, data.event.family_particles, data.event.truth_particles)); y -= 0.06;
  text.DrawLatex(0.04, y, Form("clusters: %d  anchors: %d  S/M/X/O=%d/%d/%d/%d",
      data.event.clusters, data.event.anchors, data.event.separated,
      data.event.merged, data.event.missing, data.event.other)); y -= 0.08;
  for (const auto& anchor : data.anchors)
  {
    text.SetTextColor(topology_color(anchor.topology));
    text.DrawLatex(0.04, y, Form("A%d C%u P%d %-9s E_{T}=%.2f f=%.3f",
        anchor.id, anchor.cluster, anchor.candidate,
        anchor.topology_name.c_str(), anchor.et, anchor.main_fraction));
    text.SetTextColor(kBlack); y -= 0.05;
    if (y < 0.08) break;
  }
  text.SetTextSize(0.025);
  text.DrawLatex(0.04, 0.055, "Cluster: open ring = selected #pi^{0} family (#geq 50% direct-daughter edep); fill = anchor topology.");
  text.DrawLatex(0.04, 0.025, "Dashed directions are projections; charged-particle lines are not propagated trajectories.");
}

inline void draw_anchor_text(const DisplayData& data, const Anchor& anchor)
{
  const Candidate* candidate = find_candidate(data, anchor.candidate);
  TLatex text; text.SetNDC(); text.SetTextSize(0.033);
  double y = 0.94;
  text.SetTextColor(topology_color(anchor.topology));
  text.DrawLatex(0.04, y, Form("Anchor %d  cluster %u  %s",
      anchor.id, anchor.cluster, anchor.topology_name.c_str()));
  text.SetTextColor(kBlack); y -= 0.07;
  text.DrawLatex(0.04, y, Form("reason: %s", anchor.reason_name.c_str())); y -= 0.06;
  text.DrawLatex(0.04, y, Form("E=%.3f GeV  E_{T}=%.3f GeV", anchor.energy, anchor.et)); y -= 0.06;
  text.DrawLatex(0.04, y, Form("main/second/unmatched f = %.4f / %.4f / %.4f",
      anchor.main_fraction, anchor.second_fraction, anchor.unmatched_fraction)); y -= 0.07;
  if (candidate)
  {
    text.DrawLatex(0.04, y, Form("candidate %d  %s  p_{T}=%.3f GeV",
        candidate->id, candidate->pathway_name.c_str(), candidate->pt)); y -= 0.06;
  }
  for (int gamma = 0; gamma < 2; ++gamma)
  {
    const double recovery = anchor.truth_energy[gamma] > 0.0
        ? anchor.reconstructed[gamma] / anchor.truth_energy[gamma] : 0.0;
    text.DrawLatex(0.04, y, Form("#gamma%d: Cbest=%d recovered=%d Etruth=%.3f Erec=%.3f ratio=%.3f",
        gamma, anchor.best_cluster[gamma], anchor.recovered[gamma],
        anchor.truth_energy[gamma], anchor.reconstructed[gamma], recovery)); y -= 0.06;
  }
  y -= 0.02;
  text.DrawLatex(0.04, y, Form("anchor edep: total=%.4g #gamma0=%.4g #gamma1=%.4g other=%.4g valid=%d",
      anchor.total_edep, anchor.gamma_edep[0], anchor.gamma_edep[1],
      anchor.other_edep, anchor.match_valid));
}

inline std::vector<std::unique_ptr<TCanvas>> make_event_pages(
    TFile* file, int event_id)
{
  DisplayData data;
  std::vector<std::unique_ptr<TCanvas>> pages;
  if (!load_event(file, event_id, data)) return pages;
  gStyle->SetOptStat(0);
  auto overview = std::make_unique<TCanvas>(Form("overview_%d", event_id),
      "topology event overview", 1600, 1000);
  overview->Divide(2, 2);
  overview->cd(1); draw_xy(data);
  overview->cd(2); draw_zr(data);
  overview->cd(3); draw_eta_phi(data);
  overview->cd(4); draw_overview_text(data);
  pages.push_back(std::move(overview));
  for (const auto& anchor : data.anchors)
  {
    auto page = std::make_unique<TCanvas>(Form("anchor_%d_%d", event_id, anchor.id),
        "topology anchor detail", 1600, 1000);
    page->Divide(2, 2);
    page->cd(1); draw_xy(data, anchor.candidate);
    page->cd(2); draw_eta_phi(data, anchor.candidate);
    page->cd(3); gPad->SetRightMargin(0.14); draw_anchor_towers(data, anchor);
    page->cd(4); draw_anchor_text(data, anchor);
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
