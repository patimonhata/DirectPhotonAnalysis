#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <set>
#include <vector>

namespace
{
template <class T>
bool bind(TTree* tree, const char* name, T* address)
{
  return tree->GetBranch(name) && tree->SetBranchAddress(name, address) >= 0;
}

bool close_enough(double first, double second)
{
  return std::abs(first - second) <= 2e-5 * std::max({1.0, std::abs(first), std::abs(second)});
}
}

int check_pythia_photon_candidate_tree(const char* path)
{
  TFile file(path, "READ");
  auto* events = file.Get<TTree>("event_tree");
  auto* metadata = file.Get<TTree>("metadata");
  if (file.IsZombie() || !events || !metadata || metadata->GetEntries() != 1)
  {
    std::cerr << "Missing file, event_tree, or one-entry metadata tree" << std::endl;
    return 1;
  }

  unsigned int ncluster = 0;
  double vertex_z = 0.0;
  std::vector<unsigned int>* cluster_id = nullptr;
  std::vector<double>* cluster_e = nullptr;
  std::vector<float>* score = nullptr;
  std::vector<unsigned char>* valid = nullptr;
  std::vector<unsigned char>* kinematics = nullptr;
  std::vector<unsigned char>* preselection = nullptr;
  std::vector<unsigned char>* tight = nullptr;
  std::vector<unsigned char>* nontight = nullptr;
  std::vector<float>* iso_raw = nullptr;
  std::vector<float>* iso_corrected = nullptr;
  std::vector<float>* iso_boundary = nullptr;
  std::vector<float>* noniso_boundary = nullptr;
  std::vector<unsigned char>* isolated = nullptr;
  std::vector<unsigned char>* nonisolated = nullptr;
  std::vector<unsigned char>* region_a = nullptr;
  std::vector<unsigned char>* region_b = nullptr;
  std::vector<unsigned char>* region_c = nullptr;
  std::vector<unsigned char>* region_d = nullptr;
  std::vector<unsigned char>* final_photon = nullptr;
  std::vector<unsigned char>* pi0_tag = nullptr;
  std::vector<unsigned char>* eta_tag = nullptr;
  std::vector<float>* patch = nullptr;
  std::vector<unsigned int>* pair_i = nullptr;
  unsigned int topology_anchor_count = 0;
  std::vector<unsigned int>* topology_anchor_cluster_id = nullptr;

  bool ok = true;
  ok &= bind(events, "split_ncluster", &ncluster);
  ok &= bind(events, "vertex_z", &vertex_z);
  ok &= bind(events, "split_cluster_id", &cluster_id);
  ok &= bind(events, "split_cluster_e", &cluster_e);
  ok &= bind(events, "split_cluster_bdt_score", &score);
  ok &= bind(events, "split_cluster_bdt_valid", &valid);
  ok &= bind(events, "split_cluster_pass_kinematics", &kinematics);
  ok &= bind(events, "split_cluster_pass_preselection", &preselection);
  ok &= bind(events, "split_cluster_pass_tight", &tight);
  ok &= bind(events, "split_cluster_pass_nontight", &nontight);
  ok &= bind(events, "split_cluster_iso_raw_et", &iso_raw);
  ok &= bind(events, "split_cluster_iso_corrected_et", &iso_corrected);
  ok &= bind(events, "split_cluster_iso_boundary", &iso_boundary);
  ok &= bind(events, "split_cluster_noniso_boundary", &noniso_boundary);
  ok &= bind(events, "split_cluster_pass_isolated", &isolated);
  ok &= bind(events, "split_cluster_pass_nonisolated", &nonisolated);
  ok &= bind(events, "split_cluster_pass_region_a", &region_a);
  ok &= bind(events, "split_cluster_pass_region_b", &region_b);
  ok &= bind(events, "split_cluster_pass_region_c", &region_c);
  ok &= bind(events, "split_cluster_pass_region_d", &region_d);
  ok &= bind(events, "split_cluster_pass_final_photon", &final_photon);
  ok &= bind(events, "split_cluster_pi0_tag", &pi0_tag);
  ok &= bind(events, "split_cluster_eta_tag", &eta_tag);
  ok &= bind(events, "split_cluster_shower_patch_e", &patch);
  ok &= bind(events, "split_pair_cluster_i", &pair_i);
  ok &= bind(events, "pi0_topology_anchor_count", &topology_anchor_count);
  ok &= bind(events, "pi0_topology_anchor_cluster_id", &topology_anchor_cluster_id);
  if (!ok)
  {
    std::cerr << "Missing required branch" << std::endl;
    return 2;
  }

  unsigned long long malformed = 0;
  unsigned long long checked_clusters = 0;
  for (Long64_t entry = 0; entry < events->GetEntries(); ++entry)
  {
    events->GetEntry(entry);
    const std::vector<std::size_t> sizes = {
        cluster_id->size(), cluster_e->size(), score->size(), valid->size(), kinematics->size(),
        preselection->size(), tight->size(), nontight->size(), iso_raw->size(), iso_corrected->size(),
        iso_boundary->size(), noniso_boundary->size(), isolated->size(), nonisolated->size(),
        region_a->size(), region_b->size(), region_c->size(), region_d->size(), final_photon->size(),
        pi0_tag->size(), eta_tag->size()};
    bool event_ok = std::all_of(sizes.begin(), sizes.end(), [ncluster](std::size_t size) { return size == ncluster; });
    event_ok &= std::abs(vertex_z) < 60.0 && patch->empty() && pair_i->empty();
    event_ok &= topology_anchor_cluster_id->size() == topology_anchor_count;
    const std::set<unsigned int> stored_ids(cluster_id->begin(), cluster_id->end());
    for (unsigned int id : *topology_anchor_cluster_id)
    {
      event_ok &= stored_ids.count(id) == 1U;
    }

    for (std::size_t cluster = 0; event_ok && cluster < ncluster; ++cluster)
    {
      event_ok &= (*cluster_e)[cluster] > 0.1;
      event_ok &= close_enough((*iso_corrected)[cluster], 1.2 * (*iso_raw)[cluster] + 0.1);
      event_ok &= close_enough((*noniso_boundary)[cluster], (*iso_boundary)[cluster] + 0.8);
      event_ok &= static_cast<bool>((*isolated)[cluster]) == ((*iso_corrected)[cluster] < (*iso_boundary)[cluster]);
      event_ok &= static_cast<bool>((*nonisolated)[cluster]) == ((*iso_corrected)[cluster] > (*noniso_boundary)[cluster]);
      const bool common = (*kinematics)[cluster] && (*preselection)[cluster];
      event_ok &= static_cast<bool>((*region_a)[cluster]) == (common && (*isolated)[cluster] && (*tight)[cluster]);
      event_ok &= static_cast<bool>((*region_b)[cluster]) == (common && (*nonisolated)[cluster] && (*tight)[cluster]);
      event_ok &= static_cast<bool>((*region_c)[cluster]) == (common && (*isolated)[cluster] && (*nontight)[cluster]);
      event_ok &= static_cast<bool>((*region_d)[cluster]) == (common && (*nonisolated)[cluster] && (*nontight)[cluster]);
      event_ok &= static_cast<bool>((*final_photon)[cluster]) ==
          ((*region_a)[cluster] && !(*pi0_tag)[cluster] && !(*eta_tag)[cluster]);
      event_ok &= !(*valid)[cluster] || std::isfinite((*score)[cluster]);
    }
    malformed += event_ok ? 0ULL : 1ULL;
    checked_clusters += ncluster;
  }

  std::cout << "check_pythia_photon_candidate_tree - events/clusters/malformed = "
            << events->GetEntries() << "/" << checked_clusters << "/" << malformed << std::endl;
  return malformed == 0 ? 0 : 3;
}
