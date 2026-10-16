#include "PythiaPhotonCandidateTree.h"
#include <TFile.h>
#include <TTree.h>
#include <cassert>
#include <cmath>
#include <iostream>

int main(int argc, char** argv)
{
  assert(argc == 2);
  using namespace photon_tree;
  Pi0AnchorTopologyConfig config;
  assert(config.parent_pdg_id == 111 && !config.require_hepmc_two_photon_decay);
  config.parent_pdg_id = 221; config.require_hepmc_two_photon_decay = true;
  config.enable_missing_diagnostics = true; config.tagging_partner_min_cluster_energy = .5;
  config.tagging_pi0_mass_min = .45; config.tagging_pi0_mass_max = .65;
  Pi0TopologyCandidateRecord candidate;
  candidate.photon_track_ids = {101, 102}; candidate.photon_energy = {10, 1};
  candidate.best_cluster[0] = 0; candidate.recovered[0] = true;
  candidate.photon_projection_valid = {true, true}; candidate.photon_in_cemc_acceptance = {true, true};
  auto& partner = candidate.truth_partner_clusters[1];
  partner.found = true; partner.cluster_id = 20; partner.cluster_energy = .3;
  partner.cluster_phi = std::acos(1 - .55 * .55 / (2 * 10 * .3));
  Pi0TopologyClusterRecord cluster; cluster.cluster_id = 10; cluster.energy = 10;
  Pi0TopologyAnchorRecord anchor; anchor.cluster_index = 0;
  classify_pi0_anchor(config, candidate, cluster, anchor);
  assert(anchor.topology == Pi0AnchorTopology::missing);
  assert(anchor.missing_category == Pi0MissingCategory::energy_band_taggable);
  config.tagging_pi0_mass_max = anchor.truth_partner_invariant_mass;
  classify_pi0_anchor(config, candidate, cluster, anchor);
  assert(anchor.missing_category == Pi0MissingCategory::energy_band_not_taggable);
  config.tagging_pi0_mass_max = .65;
  partner.cluster_energy = 1;
  partner.cluster_phi = std::acos(1 - .55 * .55 / (2 * 10));
  classify_pi0_anchor(config, candidate, cluster, anchor);
  assert(anchor.topology == Pi0AnchorTopology::separated && anchor.truth_partner_tag_status == Pi0TruthPartnerTagStatus::taggable);
  PhotonCandidateSelectionBranches pi0_links, eta_links;
  eta_links.pi0_anchor_valid = {1}; eta_links.pi0_anchor_topology = {1};
  Pi0AnchorTopologyEventResult result;
  result.candidates.push_back(candidate); result.clusters.push_back(cluster); result.anchors.push_back(anchor);
  Pi0TopologyTreeBranches pi0, eta;
  eta.fill(result);
  {
    TFile output(argv[1], "RECREATE"); TTree tree("event_tree", "Branch naming regression");
    pi0_links.create_branches(&tree); eta_links.create_branches(&tree, "eta", true);
    pi0.create_branches(&tree); eta.create_branches(&tree, "eta");
    assert(tree.GetBranch("split_cluster_pi0_anchor_topology") && tree.GetBranch("split_cluster_eta_anchor_topology"));
    assert(tree.GetBranch("pi0_topology_candidate_count") && tree.GetBranch("eta_topology_candidate_count"));
    assert(!tree.GetBranch("split_cluster_eta_bdt_score"));
    tree.Fill(); tree.Write();
  }
  TFile input(argv[1]); auto* tree = input.Get<TTree>("event_tree");
  unsigned int count = 0; std::vector<int>* track = nullptr; std::vector<int>* topology = nullptr;
  tree->SetBranchAddress("eta_topology_candidate_count", &count);
  tree->SetBranchAddress("eta_topology_photon0_track_id", &track);
  tree->SetBranchAddress("split_cluster_eta_anchor_topology", &topology);
  tree->GetEntry(0);
  assert(count == 1 && track && track->at(0) == 101 && topology && topology->at(0) == 1);
  std::cout << "Eta mass boundaries, partner tagging and map branch roundtrip: passed" << std::endl;
}
