#include <TFile.h>
#include <TTree.h>

#include <cmath>
#include <iostream>
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
           "schema_version", "sample_mode", "write_detail", "first_event", "min_direct_match_cluster_energy_coverage",
           "missing_diagnostic_max_delta_r", "enable_missing_diagnostics",
           "events_processed", "events_written", "events_invalid"})
    ok &= require_branch(metadata, branch);
  for (const char* branch : {"event", "n_anchors", "n_separated",
                             "n_merged", "n_missing", "n_other"})
    ok &= require_branch(events, branch);
  for (const char* branch : {"event", "candidate_id", "pathway",
                             "best_cluster0_id", "best_cluster1_id"})
    ok &= require_branch(candidates, branch);
  for (const char* branch : {
           "event", "candidate_id", "cluster_id", "topology", "reason", "main_fraction",
           "missing_detail", "partner_photon_index",
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
  int schema_version = 0;
  metadata->SetBranchAddress("schema_version", &schema_version);
  metadata->GetEntry(0);
  metadata->ResetBranchAddresses();
  if (schema_version != 2)
  {
    std::cerr << "CheckTopologyEventDump - expected schema version 2, got "
              << schema_version << std::endl;
    ok = false;
  }

  }

  int event = -1;
  int n_anchor = 0, n_separated = 0, n_merged = 0, n_missing = 0, n_other = 0;
  events->SetBranchAddress("event", &event);
  events->SetBranchAddress("n_anchors", &n_anchor);
  events->SetBranchAddress("n_separated", &n_separated);
  events->SetBranchAddress("n_merged", &n_merged);
  events->SetBranchAddress("n_missing", &n_missing);
  events->SetBranchAddress("n_other", &n_other);
  for (Long64_t entry = 0; entry < events->GetEntries(); ++entry)
  {
    events->GetEntry(entry);
    if (n_anchor != n_separated + n_merged + n_missing + n_other)
    {
      std::cerr << "CheckTopologyEventDump - category closure failed at event "
                << event << std::endl;
      ok = false;
    }
  }
  events->ResetBranchAddresses();

  std::set<std::pair<int, int>> candidate_keys;
  int candidate_id = -1;
  candidates->SetBranchAddress("event", &event);
  candidates->SetBranchAddress("candidate_id", &candidate_id);
  for (Long64_t entry = 0; entry < candidates->GetEntries(); ++entry)
  {
    candidates->GetEntry(entry);
    if (!candidate_keys.insert({event, candidate_id}).second) ok = false;
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

  int topology = -1, missing_detail = 0, partner_photon = -1;
  int match_valid = 0, match_usable = 0, match_status = 0, match_failure = 0;
  int match_failure_ieta = -999, match_failure_iphi = -999;
  unsigned int match_towers = 0, match_matched_towers = 0;
  double match_coverage = 0.0;
  int diagnostic_found = 0, diagnostic_match_usable = 0;
  int diagnostic_match_status = 0, diagnostic_match_failure = 0;
  double diagnostic_delta_r = -999.0, diagnostic_match_coverage = 0.0;
  anchors->SetBranchAddress("event", &event);
  anchors->SetBranchAddress("candidate_id", &candidate_id);
  anchors->SetBranchAddress("cluster_id", &cluster_id);
  anchors->SetBranchAddress("topology", &topology);
  anchors->SetBranchAddress("missing_detail", &missing_detail);
  anchors->SetBranchAddress("partner_photon_index", &partner_photon);
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
        topology < 0 || topology > 3;
    const bool invalid_missing_detail =
        (topology == 3 && (missing_detail < 1 || missing_detail > 5 || partner_photon < 0 || partner_photon > 1)) ||
        (topology != 3 && (missing_detail != 0 || partner_photon != -1));
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
    if (invalid_reference || invalid_missing_detail || invalid_match || invalid_diagnostic)
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
