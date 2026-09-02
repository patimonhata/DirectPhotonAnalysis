#include <TFile.h>
#include <TTree.h>

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
  std::array<double, 2> eta = {-999.0, -999.0};
  std::array<double, 2> phi = {-999.0, -999.0};
  std::array<int, 2> in_acceptance = {0, 0};
  std::array<int, 2> first_daughter_valid = {0, 0};
  std::array<double, 2> first_daughter_radius = {-999.0, -999.0};
  std::array<int, 2> pre_cemc = {0, 0};
  std::array<double, 2> cemc_edep = {0.0, 0.0};
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
  TTree* candidates = static_cast<TTree*>(file->Get("pi0_candidates"));
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
  if (schema_version != 6)
  {
    std::cerr << "CheckTopologyEventDump - expected schema version 6, got " << schema_version << std::endl;
    ok = false;
  }
  if (!std::isfinite(cemc_acceptance_eta_max) || !(cemc_acceptance_eta_max > 0.0) ||
      !std::isfinite(pre_cemc_interaction_radius) || !(pre_cemc_interaction_radius > 0.0) ||
      !std::isfinite(missing_diagnostic_max_delta_r) || !(missing_diagnostic_max_delta_r > 0.0))
  {
    std::cerr << "CheckTopologyEventDump - invalid CEMC acceptance or pre-CEMC radius" << std::endl;
    ok = false;
  }

  int event = -1;
  int n_anchor = 0, n_separated = 0, n_merged = 0, n_single_contaminated = 0, n_missing = 0, n_other = 0;
  events->SetBranchAddress("event", &event);
  events->SetBranchAddress("n_anchors", &n_anchor);
  events->SetBranchAddress("n_separated", &n_separated);
  events->SetBranchAddress("n_merged", &n_merged);
  events->SetBranchAddress("n_single_contaminated", &n_single_contaminated);
  events->SetBranchAddress("n_missing", &n_missing);
  events->SetBranchAddress("n_other", &n_other);
  for (Long64_t entry = 0; entry < events->GetEntries(); ++entry)
  {
    events->GetEntry(entry);
    if (n_anchor != n_separated + n_merged + n_single_contaminated + n_missing + n_other)
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
    bool invalid_projection = !candidate_keys.insert(key).second;
    for (int photon = 0; photon < 2; ++photon)
    {
      invalid_projection |= (projection.valid[photon] != 0 && projection.valid[photon] != 1) ||
          (projection.in_acceptance[photon] != 0 && projection.in_acceptance[photon] != 1);
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

  int topology = -1, reason = -1, missing_category = 0, missing_detail = 0, partner_photon = -1, pre_cemc_photon = -1;
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
  for (Long64_t entry = 0; entry < anchors->GetEntries(); ++entry)
  {
    anchors->GetEntry(entry);
    const bool invalid_reference = !candidate_keys.count({event, candidate_id}) ||
        (!cluster_keys.empty() && !cluster_keys.count({event, cluster_id})) ||
        topology < 0 || topology > 4;
    bool invalid_missing_detail = false;
    if (topology == 3)
    {
      invalid_missing_detail = missing_category < 1 || missing_category > 7 ||
          missing_detail < 1 || missing_detail > 11 || partner_photon < 0 || partner_photon > 1;
      const auto projection_it = candidate_projection.find({event, candidate_id});
      invalid_missing_detail |= projection_it == candidate_projection.end();
      if (!invalid_missing_detail)
      {
        const bool projection_valid = projection_it->second.valid[partner_photon];
        const bool in_acceptance = projection_it->second.in_acceptance[partner_photon];
        const double expected_cemc_edep = projection_it->second.cemc_edep[partner_photon];
        invalid_missing_detail |= !std::isfinite(partner_cemc_edep) || partner_cemc_edep < 0.0 ||
            std::abs(partner_cemc_edep - expected_cemc_edep) > 1e-9 * std::max(1.0, expected_cemc_edep);
        if (missing_category == 2)
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
    const auto projection_it = candidate_projection.find({event, candidate_id});
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
        (match_valid != (match_status == 2)) || (match_valid && !match_usable) ||
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

  match_valid = 0;
  match_usable = 0;
  match_status = 0;
  match_failure = 0;
  match_towers = match_matched_towers = 0;
  match_coverage = 0.0;
  double total = 0.0, gamma0 = 0.0, gamma1 = 0.0, other = 0.0;
  matches->SetBranchAddress("match_valid", &match_valid);
  matches->SetBranchAddress("match_usable", &match_usable);
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
        (match_valid != (match_status == 2)) || (match_valid && !match_usable) ||
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
