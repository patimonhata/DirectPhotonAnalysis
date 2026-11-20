#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <tuple>

namespace
{
bool require_branch(TTree* tree, const char* name)
{
  if (tree && tree->GetBranch(name)) return true;
  std::cerr << "CheckTopologyEventDump - missing branch "
            << (tree ? tree->GetName() : "<null>") << "/" << name
            << std::endl;
  return false;
}

struct CandidateProjection
{
  std::array<int, 2> valid = {0, 0};
  std::array<int, 2> no_intersection = {0, 0};
  int interaction_info_valid = 1;
  int parent_pdg = 111;
  std::array<int, 3> identity{};
  std::array<double, 2> eta = {-999.0, -999.0};
  std::array<double, 2> phi = {-999.0, -999.0};
  std::array<int, 2> in_acceptance = {0, 0};
  std::array<int, 2> first_daughter_valid = {0, 0};
  std::array<double, 2> first_daughter_radius = {-999.0, -999.0};
  std::array<int, 2> pre_cemc = {0, 0};
  std::array<double, 2> cemc_edep = {0.0, 0.0};
};

struct ClusterOrigin
{
  int main = 0;
  int pdg = 0;
  std::array<int, 3> identity{};
  int applicability = 0;
  int reason = 0;
  int match_status = -1;
  double coverage = -1.0;
  bool usable = false;
};
}

int CheckTopologyEventDump(const char* input_file)
{
  TFile* file = TFile::Open(input_file, "READ");
  if (!file || file->IsZombie())
  {
    std::cerr << "CheckTopologyEventDump - cannot open " << input_file
              << std::endl;
    return 1;
  }
  TTree* metadata = static_cast<TTree*>(file->Get("metadata"));
  TTree* events = static_cast<TTree*>(file->Get("events"));
  auto* candidates = file->Get<TTree>("meson_candidates");
  if (!candidates) candidates = file->Get<TTree>("pi0_candidates");
  auto* origins = file->Get<TTree>("cluster_origins");
  TTree* anchors = static_cast<TTree*>(file->Get("anchor_decisions"));
  TTree* clusters = static_cast<TTree*>(file->Get("clusters"));
  TTree* matches = static_cast<TTree*>(file->Get("candidate_cluster_truth"));
  bool ok = metadata && events && candidates && anchors && clusters && matches;
  ok &= metadata && metadata->GetEntries() == 1;
  for (const char* branch : {
           "schema_version", "sample_mode", "write_detail", "first_event", "tower_geom_node",
           "cemc_acceptance_eta_max", "pre_cemc_interaction_radius", "min_direct_match_cluster_energy_coverage",
           "missing_diagnostic_max_delta_r", "enable_missing_diagnostics",
           "events_processed", "events_written", "events_invalid"})
    ok &= require_branch(metadata, branch);
  for (const char* branch : {"event", "n_anchors", "n_separated",
                             "n_merged", "n_single_contaminated", "n_missing", "n_other"})
    ok &= require_branch(events, branch);
  for (const char* branch : {"event", "candidate_id", "pathway",
                             "best_cluster0_id", "best_cluster1_id",
                             "photon0_projection_valid", "photon1_projection_valid",
                             "photon0_projection_eta", "photon1_projection_eta",
                             "photon0_projection_phi", "photon1_projection_phi",
                             "photon0_in_cemc_acceptance", "photon1_in_cemc_acceptance",
                             "photon0_first_daughter_vertex_valid", "photon1_first_daughter_vertex_valid",
                             "photon0_first_daughter_radius", "photon1_first_daughter_radius",
                             "photon0_pre_cemc_interaction", "photon1_pre_cemc_interaction",
                             "photon0_cemc_edep", "photon1_cemc_edep"})
    ok &= require_branch(candidates, branch);
  for (const char* branch : {
           "event", "candidate_id", "cluster_id", "topology", "reason", "main_fraction",
           "missing_category", "missing_detail", "partner_photon_index", "pre_cemc_photon_index",
           "partner_cemc_edep", "partner_diagnostic_invariant_mass",
           "best_cluster0_id", "best_cluster1_id", "recovered0", "recovered1",
           "match_valid", "match_usable", "match_status", "match_failure",
           "match_failure_ieta", "match_failure_iphi", "match_tower_count",
           "match_matched_tower_count", "match_cluster_member_energy_coverage",
           "partner_diagnostic_found", "partner_diagnostic_below_energy_threshold",
           "partner_diagnostic_has_direct_deposit", "partner_diagnostic_cluster_id",
           "partner_diagnostic_cluster_energy", "partner_diagnostic_delta_r",
           "partner_diagnostic_recovery", "partner_diagnostic_match_usable",
           "partner_diagnostic_match_status", "partner_diagnostic_match_failure",
           "partner_diagnostic_failure_ieta", "partner_diagnostic_failure_iphi",
           "partner_diagnostic_match_coverage"})
    ok &= require_branch(anchors, branch);
  for (const char* branch : {
           "match_valid", "match_usable", "match_status", "match_failure",
           "match_failure_ieta", "match_failure_iphi", "match_tower_count",
           "match_matched_tower_count", "match_cluster_member_energy_coverage",
           "total_edep", "gamma0_edep", "gamma1_edep", "other_edep"})
    ok &= require_branch(matches, branch);
  if (!ok)
  {
    file->Close();
    return 2;
  }
  int schema_version = 0;
  double cemc_acceptance_eta_max = 0.0;
  double pre_cemc_interaction_radius = 0.0;
  double missing_diagnostic_max_delta_r = 0.0;
  bool enable_missing_diagnostics = false;
  metadata->SetBranchAddress("schema_version", &schema_version);
  metadata->SetBranchAddress("cemc_acceptance_eta_max", &cemc_acceptance_eta_max);
  metadata->SetBranchAddress("pre_cemc_interaction_radius", &pre_cemc_interaction_radius);
  metadata->SetBranchAddress("missing_diagnostic_max_delta_r", &missing_diagnostic_max_delta_r);
  metadata->SetBranchAddress("enable_missing_diagnostics", &enable_missing_diagnostics);
  metadata->GetEntry(0);
  metadata->ResetBranchAddresses();
  if (schema_version != 6 && schema_version != 7)
  {
    std::cerr << "CheckTopologyEventDump - expected schema version 6 or 7, got " << schema_version << std::endl;
    ok = false;
  }
  const bool current = schema_version == 7;
  int common_version = 0;
  double min_match_coverage = 0.0, missing_min = 0.0, missing_max = 0.0;
  std::array<double, 2> mass_min{}, mass_max{};
  if (current)
  {
    for (const char* branch : {"common_composition_algorithm_version", "missing_energy_min", "missing_energy_max",
                              "pi0_mass_min", "pi0_mass_max", "eta_mass_min", "eta_mass_max"})
      ok &= require_branch(metadata, branch);
    for (const char* branch : {"n_topology_clusters", "n_unevaluable_anchors"}) ok &= require_branch(events, branch);
    for (const char* branch : {"parent_pdg_id", "parent_embedding", "parent_source", "parent_id",
                              "photon0_projection_no_intersection", "photon1_projection_no_intersection", "interaction_info_valid"})
      ok &= require_branch(candidates, branch);
    for (const char* branch : {"applicability", "topology_reason", "match_ancestry_complete",
                              "truth_partner_cluster_id", "truth_partner_cluster_energy", "truth_partner_invariant_mass"})
      ok &= require_branch(anchors, branch);
    for (const char* branch : {"event", "cluster_id", "common_main", "selected_pdg", "selected_embedding", "selected_source", "selected_id",
                              "applicability", "topology_reason", "anchor_matching_status", "anchor_matching_coverage", "anchor_matching_usable"})
      ok &= require_branch(origins, branch);
    ok &= require_branch(matches, "match_ancestry_complete");
    if (!ok) return 2;
    metadata->SetBranchAddress("common_composition_algorithm_version", &common_version);
    metadata->SetBranchAddress("min_direct_match_cluster_energy_coverage", &min_match_coverage);
    metadata->SetBranchAddress("missing_energy_min", &missing_min);
    metadata->SetBranchAddress("missing_energy_max", &missing_max);
    metadata->SetBranchAddress("pi0_mass_min", &mass_min[0]);
    metadata->SetBranchAddress("pi0_mass_max", &mass_max[0]);
    metadata->SetBranchAddress("eta_mass_min", &mass_min[1]);
    metadata->SetBranchAddress("eta_mass_max", &mass_max[1]);
    metadata->GetEntry(0);
    metadata->ResetBranchAddresses();
  }
  const bool common = common_version > 0;
  if (!std::isfinite(cemc_acceptance_eta_max) || !(cemc_acceptance_eta_max > 0.0) ||
      !std::isfinite(pre_cemc_interaction_radius) || !(pre_cemc_interaction_radius > 0.0) ||
      !std::isfinite(missing_diagnostic_max_delta_r) || !(missing_diagnostic_max_delta_r > 0.0))
  {
    std::cerr << "CheckTopologyEventDump - invalid CEMC acceptance or pre-CEMC radius" << std::endl;
    ok = false;
  }

  int event = -1;
  int n_anchor = 0, n_separated = 0, n_merged = 0, n_single_contaminated = 0, n_missing = 0, n_other = 0;
  int n_unevaluable = 0, n_topology_clusters = 0;
  std::map<int, std::array<int, 7>> remaining_anchor_counts;
  std::map<int, int> remaining_origin_counts;
  events->SetBranchAddress("event", &event);
  if (current)
  {
    events->SetBranchAddress("n_unevaluable_anchors", &n_unevaluable);
    events->SetBranchAddress("n_topology_clusters", &n_topology_clusters);
  }
  events->SetBranchAddress("n_anchors", &n_anchor);
  events->SetBranchAddress("n_separated", &n_separated);
  events->SetBranchAddress("n_merged", &n_merged);
  events->SetBranchAddress("n_single_contaminated", &n_single_contaminated);
  events->SetBranchAddress("n_missing", &n_missing);
  events->SetBranchAddress("n_other", &n_other);
  for (Long64_t entry = 0; entry < events->GetEntries(); ++entry)
  {
    events->GetEntry(entry);
    if (current)
    {
      if (!remaining_anchor_counts.emplace(event, std::array<int, 7>{
          n_other, n_separated, n_merged, n_missing, n_single_contaminated, n_unevaluable, n_anchor}).second) ok = false;
      remaining_origin_counts.emplace(event, common ? n_topology_clusters : 0);
    }
    if (n_anchor != n_separated + n_merged + n_single_contaminated + n_missing + n_other + n_unevaluable)
    {
      std::cerr << "CheckTopologyEventDump - category closure failed at event "
                << event << std::endl;
      ok = false;
    }
  }
  events->ResetBranchAddresses();

  std::set<std::pair<int, int>> candidate_keys;
  std::map<std::pair<int, int>, CandidateProjection> candidate_projection;
  int candidate_id = -1;
  CandidateProjection projection;
  candidates->SetBranchAddress("event", &event);
  candidates->SetBranchAddress("candidate_id", &candidate_id);
  if (current)
  {
    candidates->SetBranchAddress("parent_pdg_id", &projection.parent_pdg);
    candidates->SetBranchAddress("parent_embedding", &projection.identity[0]);
    candidates->SetBranchAddress("parent_source", &projection.identity[1]);
    candidates->SetBranchAddress("parent_id", &projection.identity[2]);
    candidates->SetBranchAddress("photon0_projection_no_intersection", &projection.no_intersection[0]);
    candidates->SetBranchAddress("photon1_projection_no_intersection", &projection.no_intersection[1]);
    candidates->SetBranchAddress("interaction_info_valid", &projection.interaction_info_valid);
  }
  candidates->SetBranchAddress("photon0_projection_valid", &projection.valid[0]);
  candidates->SetBranchAddress("photon1_projection_valid", &projection.valid[1]);
  candidates->SetBranchAddress("photon0_projection_eta", &projection.eta[0]);
  candidates->SetBranchAddress("photon1_projection_eta", &projection.eta[1]);
  candidates->SetBranchAddress("photon0_projection_phi", &projection.phi[0]);
  candidates->SetBranchAddress("photon1_projection_phi", &projection.phi[1]);
  candidates->SetBranchAddress("photon0_in_cemc_acceptance", &projection.in_acceptance[0]);
  candidates->SetBranchAddress("photon1_in_cemc_acceptance", &projection.in_acceptance[1]);
  candidates->SetBranchAddress("photon0_first_daughter_vertex_valid", &projection.first_daughter_valid[0]);
  candidates->SetBranchAddress("photon1_first_daughter_vertex_valid", &projection.first_daughter_valid[1]);
  candidates->SetBranchAddress("photon0_first_daughter_radius", &projection.first_daughter_radius[0]);
  candidates->SetBranchAddress("photon1_first_daughter_radius", &projection.first_daughter_radius[1]);
  candidates->SetBranchAddress("photon0_pre_cemc_interaction", &projection.pre_cemc[0]);
  candidates->SetBranchAddress("photon1_pre_cemc_interaction", &projection.pre_cemc[1]);
  candidates->SetBranchAddress("photon0_cemc_edep", &projection.cemc_edep[0]);
  candidates->SetBranchAddress("photon1_cemc_edep", &projection.cemc_edep[1]);
  for (Long64_t entry = 0; entry < candidates->GetEntries(); ++entry)
  {
    candidates->GetEntry(entry);
    const std::pair<int, int> key = {event, candidate_id};
    bool invalid_projection = !candidate_keys.insert(key).second ||
        (projection.parent_pdg != 111 && projection.parent_pdg != 221) ||
        (common && ((projection.identity[1] != 1 && projection.identity[1] != 2) || projection.identity[2] == 0));
    for (int photon = 0; photon < 2; ++photon)
    {
      invalid_projection |= (projection.valid[photon] != 0 && projection.valid[photon] != 1) ||
          (projection.in_acceptance[photon] != 0 && projection.in_acceptance[photon] != 1) ||
          (projection.no_intersection[photon] != 0 && projection.no_intersection[photon] != 1) ||
          (projection.valid[photon] && projection.no_intersection[photon]);
      if (projection.valid[photon])
      {
        const bool expected_in_acceptance = std::abs(projection.eta[photon]) < cemc_acceptance_eta_max;
        invalid_projection |= !std::isfinite(projection.eta[photon]) || !std::isfinite(projection.phi[photon]) ||
            projection.in_acceptance[photon] != static_cast<int>(expected_in_acceptance);
      }
      else
      {
        invalid_projection |= projection.in_acceptance[photon] != 0;
      }
      invalid_projection |= !std::isfinite(projection.cemc_edep[photon]) || projection.cemc_edep[photon] < 0.0;
      invalid_projection |= (projection.first_daughter_valid[photon] != 0 && projection.first_daughter_valid[photon] != 1) ||
          (projection.pre_cemc[photon] != 0 && projection.pre_cemc[photon] != 1);
      if (projection.first_daughter_valid[photon])
      {
        invalid_projection |= !std::isfinite(projection.first_daughter_radius[photon]) || projection.first_daughter_radius[photon] < 0.0 ||
            projection.pre_cemc[photon] != static_cast<int>(projection.first_daughter_radius[photon] < pre_cemc_interaction_radius);
      }
      else
      {
        invalid_projection |= projection.pre_cemc[photon] != 0;
      }
    }
    if (invalid_projection)
    {
      std::cerr << "CheckTopologyEventDump - invalid candidate projection at entry " << entry << std::endl;
      ok = false;
    }
    candidate_projection[key] = projection;
  }
  candidates->ResetBranchAddresses();

  std::set<std::pair<int, unsigned int>> cluster_keys;
  unsigned int cluster_id = 0;
  clusters->SetBranchAddress("event", &event);
  clusters->SetBranchAddress("cluster_id", &cluster_id);
  for (Long64_t entry = 0; entry < clusters->GetEntries(); ++entry)
  {
    clusters->GetEntry(entry);
    if (!cluster_keys.insert({event, cluster_id}).second) ok = false;
  }
  clusters->ResetBranchAddresses();

  std::map<std::pair<int, unsigned int>, ClusterOrigin> cluster_origins;
  if (current)
  {
    ClusterOrigin origin;
    origins->SetBranchAddress("event", &event);
    origins->SetBranchAddress("cluster_id", &cluster_id);
    origins->SetBranchAddress("common_main", &origin.main);
    origins->SetBranchAddress("selected_pdg", &origin.pdg);
    origins->SetBranchAddress("selected_embedding", &origin.identity[0]);
    origins->SetBranchAddress("selected_source", &origin.identity[1]);
    origins->SetBranchAddress("selected_id", &origin.identity[2]);
    origins->SetBranchAddress("applicability", &origin.applicability);
    origins->SetBranchAddress("topology_reason", &origin.reason);
    origins->SetBranchAddress("anchor_matching_status", &origin.match_status);
    origins->SetBranchAddress("anchor_matching_coverage", &origin.coverage);
    origins->SetBranchAddress("anchor_matching_usable", &origin.usable);
    for (Long64_t entry = 0; entry < origins->GetEntries(); ++entry)
    {
      origins->GetEntry(entry);
      const auto counts = remaining_origin_counts.find(event);
      if (counts != remaining_origin_counts.end()) --counts->second;
      if (!common || counts == remaining_origin_counts.end() || !cluster_origins.emplace(std::make_pair(event, cluster_id), origin).second ||
          (!cluster_keys.empty() && !cluster_keys.count({event, cluster_id})) ||
          origin.main < 0 || origin.main > 11 || origin.applicability < 0 || origin.applicability > 2 || origin.reason < 0 || origin.reason > 6 ||
          origin.match_status < -1 || origin.match_status > 2 || !std::isfinite(origin.coverage) ||
          (origin.match_status < 0 ? origin.coverage != -1 || origin.usable : origin.coverage < 0.0 || origin.coverage > 1.0 + 1e-6) ||
          (origin.applicability == 1 && (!origin.usable || origin.reason != 0 || (origin.main != 2 && origin.main != 3))) ||
          (origin.applicability == 2 && origin.reason != 1))
      {
        std::cerr << "CheckTopologyEventDump - invalid cluster origin at entry " << entry << std::endl;
        ok = false;
      }
    }
    origins->ResetBranchAddresses();
  }

  int topology = -1, reason = -1, missing_category = 0, missing_detail = 0, partner_photon = -1, pre_cemc_photon = -1;
  int applicability = -1, topology_reason = 0, truth_partner = -999;
  double truth_partner_energy = -999.0, truth_partner_mass = -999.0;
  bool ancestry_complete = true;
  int best_cluster0 = -999, best_cluster1 = -999, recovered0 = 0, recovered1 = 0;
  int match_valid = 0, match_usable = 0, match_status = 0, match_failure = 0;
  int match_failure_ieta = -999, match_failure_iphi = -999;
  unsigned int match_towers = 0, match_matched_towers = 0;
  double match_coverage = 0.0;
  double partner_cemc_edep = 0.0, diagnostic_invariant_mass = -999.0;
  int diagnostic_found = 0, diagnostic_below_threshold = 0, diagnostic_has_direct = 0, diagnostic_match_usable = 0;
  int diagnostic_match_status = 0, diagnostic_match_failure = 0;
  double diagnostic_delta_r = -999.0, diagnostic_match_coverage = 0.0;
  anchors->SetBranchAddress("event", &event);
  anchors->SetBranchAddress("candidate_id", &candidate_id);
  anchors->SetBranchAddress("cluster_id", &cluster_id);
  anchors->SetBranchAddress("topology", &topology);
  if (current)
  {
    anchors->SetBranchAddress("applicability", &applicability);
    anchors->SetBranchAddress("topology_reason", &topology_reason);
    anchors->SetBranchAddress("truth_partner_cluster_id", &truth_partner);
    anchors->SetBranchAddress("truth_partner_cluster_energy", &truth_partner_energy);
    anchors->SetBranchAddress("truth_partner_invariant_mass", &truth_partner_mass);
    anchors->SetBranchAddress("match_ancestry_complete", &ancestry_complete);
  }
  anchors->SetBranchAddress("reason", &reason);
  anchors->SetBranchAddress("missing_category", &missing_category);
  anchors->SetBranchAddress("missing_detail", &missing_detail);
  anchors->SetBranchAddress("partner_photon_index", &partner_photon);
  anchors->SetBranchAddress("pre_cemc_photon_index", &pre_cemc_photon);
  anchors->SetBranchAddress("partner_cemc_edep", &partner_cemc_edep);
  anchors->SetBranchAddress("partner_diagnostic_invariant_mass", &diagnostic_invariant_mass);
  anchors->SetBranchAddress("best_cluster0_id", &best_cluster0);
  anchors->SetBranchAddress("best_cluster1_id", &best_cluster1);
  anchors->SetBranchAddress("recovered0", &recovered0);
  anchors->SetBranchAddress("recovered1", &recovered1);
  anchors->SetBranchAddress("match_valid", &match_valid);
  anchors->SetBranchAddress("match_usable", &match_usable);
  anchors->SetBranchAddress("match_status", &match_status);
  anchors->SetBranchAddress("match_failure", &match_failure);
  anchors->SetBranchAddress("match_failure_ieta", &match_failure_ieta);
  anchors->SetBranchAddress("match_failure_iphi", &match_failure_iphi);
  anchors->SetBranchAddress("match_tower_count", &match_towers);
  anchors->SetBranchAddress("match_matched_tower_count", &match_matched_towers);
  anchors->SetBranchAddress("match_cluster_member_energy_coverage", &match_coverage);
  anchors->SetBranchAddress("partner_diagnostic_found", &diagnostic_found);
  anchors->SetBranchAddress("partner_diagnostic_below_energy_threshold", &diagnostic_below_threshold);
  anchors->SetBranchAddress("partner_diagnostic_has_direct_deposit", &diagnostic_has_direct);
  anchors->SetBranchAddress("partner_diagnostic_delta_r", &diagnostic_delta_r);
  anchors->SetBranchAddress("partner_diagnostic_match_usable", &diagnostic_match_usable);
  anchors->SetBranchAddress("partner_diagnostic_match_status", &diagnostic_match_status);
  anchors->SetBranchAddress("partner_diagnostic_match_failure", &diagnostic_match_failure);
  anchors->SetBranchAddress("partner_diagnostic_match_coverage", &diagnostic_match_coverage);
  std::set<std::pair<int, unsigned int>> anchor_keys;
  const auto expected_usable = [&](int valid, int status, double coverage) {
    return ancestry_complete && (valid || (status == 1 && coverage >= min_match_coverage));
  };
  for (Long64_t entry = 0; entry < anchors->GetEntries(); ++entry)
  {
    anchors->GetEntry(entry);
    const auto projection_it = candidate_projection.find({event, candidate_id});
    bool invalid_reference = projection_it == candidate_projection.end() ||
        (!cluster_keys.empty() && !cluster_keys.count({event, cluster_id})) ||
        topology < 0 || topology > 4;
    if (current && !invalid_reference)
    {
      const auto counts = remaining_anchor_counts.find(event);
      invalid_reference = counts == remaining_anchor_counts.end() || !anchor_keys.emplace(event, cluster_id).second;
      if (counts != remaining_anchor_counts.end())
      {
        --counts->second[6];
        --counts->second[common && applicability != 1 ? 5 : topology];
      }
      if (common)
      {
        const auto found = cluster_origins.find({event, cluster_id});
        const auto& particle = projection_it->second;
        const bool projected = (particle.valid[0] || particle.no_intersection[0]) && (particle.valid[1] || particle.no_intersection[1]);
        const int expected_reason = !match_usable ? 4 : !projected ? 5 : !particle.interaction_info_valid ? 6 : 0;
        invalid_reference |= found == cluster_origins.end() || topology_reason != expected_reason || applicability != (expected_reason == 0);
        if (found != cluster_origins.end())
        {
          const auto& origin = found->second;
          invalid_reference |= origin.pdg != particle.parent_pdg || origin.main != (particle.parent_pdg == 111 ? 2 : 3) ||
              origin.identity != particle.identity || origin.applicability != applicability || origin.reason != topology_reason ||
              origin.match_status != match_status || origin.coverage != match_coverage || origin.usable != static_cast<bool>(match_usable);
        }
      }
      else invalid_reference |= applicability != -1;
    }
    bool invalid_missing_detail = false;
    if (topology == 3)
    {
      invalid_missing_detail = missing_category < 1 || missing_category > 7 ||
          missing_detail < 1 || missing_detail > 11 || partner_photon < 0 || partner_photon > 1;
      invalid_missing_detail |= projection_it == candidate_projection.end();
      if (!invalid_missing_detail)
      {
        const bool projection_valid = projection_it->second.valid[partner_photon];
        const bool in_acceptance = projection_it->second.in_acceptance[partner_photon];
        const double expected_cemc_edep = projection_it->second.cemc_edep[partner_photon];
        invalid_missing_detail |= !std::isfinite(partner_cemc_edep) || partner_cemc_edep < 0.0 ||
            std::abs(partner_cemc_edep - expected_cemc_edep) > 1e-9 * std::max(1.0, expected_cemc_edep);
        if (current)
        {
          int category = 3, detail = 1;
          if (projection_it->second.no_intersection[partner_photon] || (projection_valid && !in_acceptance)) { category = 2; detail = 6; }
          else if (!projection_valid) detail = 7;
          else if (!enable_missing_diagnostics) detail = 11;
          else if (truth_partner >= 0)
          {
            if (truth_partner != static_cast<int>(cluster_id))
            {
              const int meson = projection_it->second.parent_pdg == 221 ? 1 : 0;
              if (truth_partner_energy <= missing_min) category = 5;
              else if (truth_partner_energy <= missing_max && std::isfinite(truth_partner_mass) && truth_partner_mass >= 0.0)
                category = truth_partner_mass > mass_min[meson] && truth_partner_mass < mass_max[meson] ? 1 : 4;
            }
          }
          else if (diagnostic_found && !diagnostic_match_usable) detail = 4;
          else
          {
            const int best = partner_photon == 0 ? best_cluster0 : best_cluster1;
            category = 6;
            detail = expected_cemc_edep > 0.0 || best >= 0 || diagnostic_has_direct ? 10 : 5;
          }
          invalid_missing_detail |= missing_category != category || missing_detail != detail;
        }
        else if (missing_category == 2)
        {
          invalid_missing_detail |= missing_detail != 6 || !projection_valid || in_acceptance;
        }
        else if (missing_category == 1)
        {
          invalid_missing_detail |= (missing_detail != 2 && missing_detail != 3) || !projection_valid || !in_acceptance ||
              !diagnostic_found || !diagnostic_below_threshold || !diagnostic_has_direct ||
              diagnostic_delta_r > missing_diagnostic_max_delta_r ||
              !std::isfinite(diagnostic_invariant_mass) || diagnostic_invariant_mass < 0.0;
        }
        else if (missing_category == 4)
        {
          invalid_missing_detail |= (missing_detail != 8 && missing_detail != 9) || !projection_valid || !in_acceptance ||
              !diagnostic_found || !diagnostic_below_threshold || !diagnostic_has_direct ||
              diagnostic_delta_r <= missing_diagnostic_max_delta_r ||
              !std::isfinite(diagnostic_invariant_mass) || diagnostic_invariant_mass < 0.0;
        }
        else if (missing_category == 5)
        {
          invalid_missing_detail |= missing_detail != 5 || !projection_valid || !in_acceptance || partner_cemc_edep > 0.0;
        }
        else if (missing_category == 6)
        {
          invalid_missing_detail |= missing_detail != 10 || !projection_valid || !in_acceptance || !(partner_cemc_edep > 0.0);
        }
        else if (missing_category == 7)
        {
          invalid_missing_detail |= missing_detail != 4 || !projection_valid || !in_acceptance ||
              !diagnostic_found || diagnostic_match_usable;
        }
        else
        {
          invalid_missing_detail |= (missing_detail != 1 && missing_detail != 7 && missing_detail != 11) ||
              (missing_detail == 7 && projection_valid) ||
              (missing_detail == 1 && (!projection_valid || !in_acceptance)) ||
              (missing_detail == 11 && enable_missing_diagnostics);
        }
      }
    }
    else
    {
      invalid_missing_detail = missing_category != 0 || missing_detail != 0 ||
          (topology == 1 ? (partner_photon < 0 || partner_photon > 1) : partner_photon != -1);
    }
    bool invalid_pre_cemc = false;
    if (topology == 4)
    {
      invalid_pre_cemc = reason != 6 || pre_cemc_photon < 0 || pre_cemc_photon > 1 || projection_it == candidate_projection.end() ||
          !recovered0 || !recovered1 || best_cluster0 != static_cast<int>(cluster_id) || best_cluster1 != static_cast<int>(cluster_id);
      if (!invalid_pre_cemc)
      {
        const int other_photon = 1 - pre_cemc_photon;
        invalid_pre_cemc = !projection_it->second.pre_cemc[pre_cemc_photon] || projection_it->second.pre_cemc[other_photon];
      }
    }
    else
    {
      invalid_pre_cemc = pre_cemc_photon != -1;
    }
    const bool invalid_match = match_status < 0 || match_status > 2 ||
        match_failure < 0 || match_failure > 7 ||
        !std::isfinite(match_coverage) || match_coverage < 0.0 || match_coverage > 1.0 + 1e-6 ||
        match_matched_towers > match_towers ||
        (match_valid != (match_status == 2)) ||
        (common ? static_cast<bool>(match_usable) != expected_usable(match_valid, match_status, match_coverage) : match_valid && !match_usable) ||
        (match_usable && match_status == 0) || (match_status == 2 && match_failure != 0);
    const bool invalid_diagnostic = diagnostic_found &&
        (!std::isfinite(diagnostic_delta_r) || diagnostic_delta_r < 0.0 ||
         diagnostic_match_status < 0 || diagnostic_match_status > 2 ||
         diagnostic_match_failure < 0 || diagnostic_match_failure > 7 ||
         !std::isfinite(diagnostic_match_coverage) ||
         diagnostic_match_coverage < 0.0 || diagnostic_match_coverage > 1.0 + 1e-6 ||
         (diagnostic_match_usable && diagnostic_match_status == 0));
    if (invalid_reference || invalid_missing_detail || invalid_pre_cemc || invalid_match || invalid_diagnostic)
    {
      std::cerr << "CheckTopologyEventDump - invalid anchor record at entry "
                << entry << std::endl;
      ok = false;
    }
  }
  anchors->ResetBranchAddresses();
  for (const auto& counts : remaining_anchor_counts)
    if (std::any_of(counts.second.begin(), counts.second.end(), [](int remaining) { return remaining != 0; })) ok = false;
  for (const auto& counts : remaining_origin_counts) if (counts.second != 0) ok = false;
  for (const auto& entry : cluster_origins)
    if ((entry.second.match_status >= 0) != (anchor_keys.count(entry.first) != 0)) ok = false;

  match_valid = 0;
  match_usable = 0;
  match_status = 0;
  match_failure = 0;
  match_towers = match_matched_towers = 0;
  match_coverage = 0.0;
  double total = 0.0, gamma0 = 0.0, gamma1 = 0.0, other = 0.0;
  matches->SetBranchAddress("match_valid", &match_valid);
  matches->SetBranchAddress("match_usable", &match_usable);
  if (current) matches->SetBranchAddress("match_ancestry_complete", &ancestry_complete);
  matches->SetBranchAddress("match_status", &match_status);
  matches->SetBranchAddress("match_failure", &match_failure);
  matches->SetBranchAddress("match_tower_count", &match_towers);
  matches->SetBranchAddress("match_matched_tower_count", &match_matched_towers);
  matches->SetBranchAddress("match_cluster_member_energy_coverage", &match_coverage);
  matches->SetBranchAddress("total_edep", &total);
  matches->SetBranchAddress("gamma0_edep", &gamma0);
  matches->SetBranchAddress("gamma1_edep", &gamma1);
  matches->SetBranchAddress("other_edep", &other);
  for (Long64_t entry = 0; entry < matches->GetEntries(); ++entry)
  {
    matches->GetEntry(entry);
    const double scale = std::max(1.0, std::abs(total));
    const bool invalid_match = match_status < 0 || match_status > 2 ||
        match_failure < 0 || match_failure > 7 ||
        !std::isfinite(match_coverage) || match_coverage < 0.0 || match_coverage > 1.0 + 1e-6 ||
        match_matched_towers > match_towers ||
        (match_valid != (match_status == 2)) ||
        (common ? static_cast<bool>(match_usable) != expected_usable(match_valid, match_status, match_coverage) : match_valid && !match_usable) ||
        (match_usable && match_status == 0) || (match_status == 2 && match_failure != 0);
    const bool invalid_closure = match_usable &&
        std::abs(total - gamma0 - gamma1 - other) > 1e-5 * scale;
    if (invalid_match || invalid_closure)
    {
      std::cerr << "CheckTopologyEventDump - invalid direct match at entry "
                << entry << std::endl;
      ok = false;
    }
  }
  matches->ResetBranchAddresses();
  file->Close();
  std::cout << "CheckTopologyEventDump - " << (ok ? "OK" : "FAILED")
            << std::endl;
  return ok ? 0 : 3;
}
