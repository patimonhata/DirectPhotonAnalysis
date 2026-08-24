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
  for (const char* branch : {"schema_version", "sample_mode", "write_detail",
                             "events_processed", "events_written", "events_invalid"})
    ok &= require_branch(metadata, branch);
  for (const char* branch : {"event", "n_anchors", "n_separated",
                             "n_merged", "n_missing", "n_other"})
    ok &= require_branch(events, branch);
  for (const char* branch : {"event", "candidate_id", "pathway",
                             "best_cluster0_id", "best_cluster1_id"})
    ok &= require_branch(candidates, branch);
  for (const char* branch : {"event", "candidate_id", "cluster_id",
                             "topology", "reason", "main_fraction"})
    ok &= require_branch(anchors, branch);
  if (!ok)
  {
    file->Close();
    return 2;
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

  int topology = -1;
  anchors->SetBranchAddress("event", &event);
  anchors->SetBranchAddress("candidate_id", &candidate_id);
  anchors->SetBranchAddress("cluster_id", &cluster_id);
  anchors->SetBranchAddress("topology", &topology);
  for (Long64_t entry = 0; entry < anchors->GetEntries(); ++entry)
  {
    anchors->GetEntry(entry);
    if (!candidate_keys.count({event, candidate_id}) ||
        (!cluster_keys.empty() && !cluster_keys.count({event, cluster_id})) ||
        topology < 0 || topology > 3)
    {
      std::cerr << "CheckTopologyEventDump - invalid anchor reference at entry "
                << entry << std::endl;
      ok = false;
    }
  }
  anchors->ResetBranchAddresses();

  int match_valid = 0;
  double total = 0.0, gamma0 = 0.0, gamma1 = 0.0, other = 0.0;
  matches->SetBranchAddress("match_valid", &match_valid);
  matches->SetBranchAddress("total_edep", &total);
  matches->SetBranchAddress("gamma0_edep", &gamma0);
  matches->SetBranchAddress("gamma1_edep", &gamma1);
  matches->SetBranchAddress("other_edep", &other);
  for (Long64_t entry = 0; entry < matches->GetEntries(); ++entry)
  {
    matches->GetEntry(entry);
    const double scale = std::max(1.0, std::abs(total));
    if (match_valid && std::abs(total - gamma0 - gamma1 - other) > 1e-5 * scale)
    {
      std::cerr << "CheckTopologyEventDump - edep closure failed at entry "
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
