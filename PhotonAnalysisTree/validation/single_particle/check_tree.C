#include <TFile.h>
#include <TKey.h>
#include <TTree.h>

#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace
{
template <class T>
bool bind_branch(TTree* tree, const std::string& name, T* address)
{
  return tree->GetBranch(name.c_str()) &&
      tree->SetBranchAddress(name.c_str(), address) >= 0;
}

template <class T>
bool aligned(const std::vector<T>* values, const std::size_t expected)
{
  return values && values->size() == expected;
}

struct CollectionBranches
{
  UInt_t ncluster = 0;
  UInt_t ntower = 0;
  std::vector<unsigned int>* cluster_id = nullptr;
  std::vector<int>* cluster_ntower = nullptr;
  std::vector<double>* cluster_e = nullptr;
  std::vector<unsigned char>* shower_valid = nullptr;
  std::vector<float>* shower_patch_e = nullptr;
  std::vector<unsigned char>* shower_patch_good = nullptr;
  std::vector<unsigned char>* shower_patch_owned = nullptr;
  std::vector<unsigned int>* pair_i = nullptr;
  std::vector<unsigned int>* pair_j = nullptr;
  std::vector<double>* pair_mass = nullptr;
  std::vector<double>* pair_asymmetry = nullptr;
  std::vector<int>* tower_cluster_index = nullptr;
  std::vector<unsigned int>* tower_key = nullptr;
  std::vector<double>* tower_energy = nullptr;
  std::vector<double>* tower_cluster_value = nullptr;
  std::vector<unsigned char>* truth_match_valid = nullptr;
  std::vector<float>* truth_total_edep = nullptr;
  std::vector<float>* truth_gamma0_edep = nullptr;
  std::vector<float>* truth_gamma1_edep = nullptr;
  std::vector<float>* truth_other_edep = nullptr;
  std::vector<float>* truth_gamma0_fraction = nullptr;
  std::vector<float>* truth_gamma1_fraction = nullptr;
  std::vector<float>* truth_other_fraction = nullptr;
  std::vector<float>* truth_gamma0_recovery = nullptr;
  std::vector<float>* truth_gamma1_recovery = nullptr;
};

bool bind_collection(TTree* tree, const std::string& prefix, CollectionBranches& c)
{
  const auto name = [&prefix](const char* suffix) {
    return prefix + "_" + suffix;
  };
  bool ok = true;
  ok &= bind_branch(tree, name("ncluster"), &c.ncluster);
  ok &= bind_branch(tree, name("ntower"), &c.ntower);
  ok &= bind_branch(tree, name("cluster_id"), &c.cluster_id);
  ok &= bind_branch(tree, name("cluster_ntower"), &c.cluster_ntower);
  ok &= bind_branch(tree, name("cluster_e"), &c.cluster_e);
  ok &= bind_branch(tree, name("cluster_shower_valid"), &c.shower_valid);
  ok &= bind_branch(tree, name("cluster_shower_patch_e"), &c.shower_patch_e);
  ok &= bind_branch(tree, name("cluster_shower_patch_good"), &c.shower_patch_good);
  ok &= bind_branch(tree, name("cluster_shower_patch_owned"), &c.shower_patch_owned);
  ok &= bind_branch(tree, name("pair_cluster_i"), &c.pair_i);
  ok &= bind_branch(tree, name("pair_cluster_j"), &c.pair_j);
  ok &= bind_branch(tree, name("pair_m_gg"), &c.pair_mass);
  ok &= bind_branch(tree, name("pair_e_asym"), &c.pair_asymmetry);
  ok &= bind_branch(tree, name("tower_cluster_index"), &c.tower_cluster_index);
  ok &= bind_branch(tree, name("tower_key"), &c.tower_key);
  ok &= bind_branch(tree, name("tower_energy"), &c.tower_energy);
  ok &= bind_branch(tree, name("tower_cluster_value"), &c.tower_cluster_value);
  ok &= bind_branch(tree, prefix + "_cluster_truth_match_valid", &c.truth_match_valid);
  ok &= bind_branch(tree, prefix + "_cluster_truth_total_edep", &c.truth_total_edep);
  ok &= bind_branch(tree, prefix + "_cluster_truth_gamma0_edep", &c.truth_gamma0_edep);
  ok &= bind_branch(tree, prefix + "_cluster_truth_gamma1_edep", &c.truth_gamma1_edep);
  ok &= bind_branch(tree, prefix + "_cluster_truth_other_edep", &c.truth_other_edep);
  ok &= bind_branch(tree, prefix + "_cluster_truth_gamma0_fraction", &c.truth_gamma0_fraction);
  ok &= bind_branch(tree, prefix + "_cluster_truth_gamma1_fraction", &c.truth_gamma1_fraction);
  ok &= bind_branch(tree, prefix + "_cluster_truth_other_fraction", &c.truth_other_fraction);
  ok &= bind_branch(tree, prefix + "_cluster_truth_gamma0_recovery", &c.truth_gamma0_recovery);
  ok &= bind_branch(tree, prefix + "_cluster_truth_gamma1_recovery", &c.truth_gamma1_recovery);
  return ok;
}

bool collection_aligned(const CollectionBranches& c,
                        const bool store_patch,
                        const int patch_side)
{
  const std::size_t ncluster = c.ncluster;
  const std::size_t ntower = c.ntower;
  const std::size_t npair =
      ncluster * (ncluster > 0U ? ncluster - 1U : 0U) / 2U;
  const std::size_t patch_size = store_patch
      ? ncluster * static_cast<std::size_t>(patch_side) *
            static_cast<std::size_t>(patch_side)
      : 0U;

  bool ok = aligned(c.cluster_id, ncluster) &&
      aligned(c.cluster_ntower, ncluster) &&
      aligned(c.cluster_e, ncluster) &&
      aligned(c.shower_valid, ncluster) &&
      aligned(c.shower_patch_e, patch_size) &&
      aligned(c.shower_patch_good, patch_size) &&
      aligned(c.shower_patch_owned, patch_size) &&
      aligned(c.pair_i, npair) &&
      aligned(c.pair_j, npair) &&
      aligned(c.pair_mass, npair) &&
      aligned(c.pair_asymmetry, npair) &&
      aligned(c.tower_cluster_index, ntower) &&
      aligned(c.tower_key, ntower) &&
      aligned(c.tower_energy, ntower) &&
      aligned(c.tower_cluster_value, ntower) &&
      aligned(c.truth_match_valid, ncluster) &&
      aligned(c.truth_total_edep, ncluster) &&
      aligned(c.truth_gamma0_edep, ncluster) &&
      aligned(c.truth_gamma1_edep, ncluster) &&
      aligned(c.truth_other_edep, ncluster) &&
      aligned(c.truth_gamma0_fraction, ncluster) &&
      aligned(c.truth_gamma1_fraction, ncluster) &&
      aligned(c.truth_other_fraction, ncluster) &&
      aligned(c.truth_gamma0_recovery, ncluster) &&
      aligned(c.truth_gamma1_recovery, ncluster);
  if (!ok)
  {
    return false;
  }

  for (std::size_t pair = 0; pair < npair; ++pair)
  {
    if ((*c.pair_i)[pair] >= ncluster || (*c.pair_j)[pair] >= ncluster ||
        (*c.pair_i)[pair] >= (*c.pair_j)[pair])
    {
      return false;
    }
  }
  for (const int cluster_index : *c.tower_cluster_index)
  {
    if (cluster_index < 0 ||
        static_cast<std::size_t>(cluster_index) >= ncluster)
    {
      return false;
    }
  }
  return true;
}

bool clean_layout(TFile& file)
{
  bool has_event_tree = false;
  bool has_metadata = false;
  bool ok = file.GetListOfKeys()->GetSize() == 2;
  TIter next(file.GetListOfKeys());
  while (auto* key = dynamic_cast<TKey*>(next()))
  {
    const std::string name = key->GetName();
    const bool is_tree = std::string(key->GetClassName()) == "TTree";
    has_event_tree |= name == "event_tree" && is_tree;
    has_metadata |= name == "metadata" && is_tree;
    ok &= (name == "event_tree" || name == "metadata") && is_tree;
  }
  return ok && has_event_tree && has_metadata;
}
}

int check_tree(const char* input_path, const bool require_unscored = true)
{
  if (!input_path)
  {
    return 1;
  }
  std::unique_ptr<TFile> input(TFile::Open(input_path, "READ"));
  if (!input || input->IsZombie() || !clean_layout(*input))
  {
    std::cerr << "check_tree - expected only event_tree and metadata TTrees"
              << std::endl;
    return 2;
  }
  TTree* tree = input->Get<TTree>("event_tree");
  TTree* metadata = input->Get<TTree>("metadata");
  if (!tree || !metadata)
  {
    return 2;
  }

  Int_t schema_version = 0;
  UInt_t metadata_source_file_id = 0;
  Int_t expected_primary_pdg = 0;
  ULong64_t n_events_processed = 0;
  ULong64_t n_events_written = 0;
  ULong64_t n_events_invalid_truth = 0;
  ULong64_t n_events_invalid_detector = 0;
  const bool metadata_ok = metadata->GetEntries() == 1 &&
      bind_branch(metadata, "schema_version", &schema_version) &&
      bind_branch(metadata, "source_file_id", &metadata_source_file_id) &&
      bind_branch(metadata, "expected_primary_pdg", &expected_primary_pdg) &&
      bind_branch(metadata, "n_events_processed", &n_events_processed) &&
      bind_branch(metadata, "n_events_written", &n_events_written) &&
      bind_branch(metadata, "n_events_invalid_truth", &n_events_invalid_truth) &&
      bind_branch(metadata, "n_events_invalid_detector", &n_events_invalid_detector) &&
      metadata->GetEntry(0) > 0 &&
      schema_version == 4 &&
      (expected_primary_pdg == 22 || expected_primary_pdg == 111 ||
       expected_primary_pdg == 221) &&
      n_events_written == static_cast<ULong64_t>(tree->GetEntries()) &&
      n_events_processed >= n_events_written &&
      n_events_invalid_truth <= n_events_written &&
      n_events_invalid_detector <= n_events_processed;
  if (!metadata_ok)
  {
    std::cerr << "check_tree - unreadable or inconsistent metadata" << std::endl;
    return 3;
  }

  const std::vector<std::string> score_branches = {
      "split_cluster_bdt_base_v3E_score",
      "split_cluster_bdt_base_v3E_valid",
      "split_cluster_bdt_ppg15v1_score",
      "split_cluster_bdt_ppg15v1_valid",
      "nosplit_cluster_bdt_base_v3E_score",
      "nosplit_cluster_bdt_base_v3E_valid",
      "split_cluster_p_gamma",
      "split_cluster_p_gamma_valid",
      "nosplit_cluster_p_gamma",
      "nosplit_cluster_p_gamma_valid"};
  if (require_unscored)
  {
    for (const std::string& branch : score_branches)
    {
      if (tree->GetBranch(branch.c_str()))
      {
        std::cerr << "check_tree - score branch already exists: " << branch
                  << std::endl;
        return 4;
      }
    }
  }

  UInt_t source_file_id = 0;
  UInt_t event_in_file = 0;
  ULong64_t event_uid = 0;
  UInt_t truth_n_direct_daughter = 0;
  Bool_t store_patch = false;
  Int_t patch_side = 0;
  std::vector<int>* truth_daughter_track_id = nullptr;
  std::vector<int>* truth_daughter_pdg_id = nullptr;
  std::vector<double>* truth_daughter_e = nullptr;
  std::vector<double>* truth_daughter_projection_eta = nullptr;
  std::vector<double>* truth_daughter_projection_phi = nullptr;
  std::vector<unsigned char>* truth_daughter_projection_valid = nullptr;
  std::vector<unsigned char>* truth_daughter_in_acceptance = nullptr;
  CollectionBranches split;
  CollectionBranches nosplit;

  bool branches_ok = true;
  branches_ok &= bind_branch(tree, "source_file_id", &source_file_id);
  branches_ok &= bind_branch(tree, "event_in_file", &event_in_file);
  branches_ok &= bind_branch(tree, "event_uid", &event_uid);
  branches_ok &= bind_branch(tree, "truth_n_direct_daughter", &truth_n_direct_daughter);
  branches_ok &= bind_branch(tree, "store_shower_shape_tower_patch", &store_patch);
  branches_ok &= bind_branch(tree, "shower_shape_patch_side", &patch_side);
  branches_ok &= bind_branch(tree, "truth_daughter_track_id", &truth_daughter_track_id);
  branches_ok &= bind_branch(tree, "truth_daughter_pdg_id", &truth_daughter_pdg_id);
  branches_ok &= bind_branch(tree, "truth_daughter_e", &truth_daughter_e);
  branches_ok &= bind_branch(tree, "truth_daughter_projection_eta",
                      &truth_daughter_projection_eta);
  branches_ok &= bind_branch(tree, "truth_daughter_projection_phi",
                      &truth_daughter_projection_phi);
  branches_ok &= bind_branch(tree, "truth_daughter_projection_valid",
                      &truth_daughter_projection_valid);
  branches_ok &= bind_branch(tree, "truth_daughter_in_acceptance",
                      &truth_daughter_in_acceptance);
  branches_ok &= bind_collection(tree, "split", split);
  branches_ok &= bind_collection(tree, "nosplit", nosplit);
  if (!branches_ok)
  {
    std::cerr << "check_tree - missing required base branch" << std::endl;
    return 5;
  }

  Long64_t malformed = 0;
  for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry)
  {
    if (tree->GetEntry(entry) <= 0)
    {
      ++malformed;
      continue;
    }
    const std::size_t ndaughter = truth_n_direct_daughter;
    const ULong64_t expected_uid =
        (static_cast<ULong64_t>(source_file_id) << 32U) |
        static_cast<ULong64_t>(event_in_file);
    const bool truth_aligned =
        aligned(truth_daughter_track_id, ndaughter) &&
        aligned(truth_daughter_pdg_id, ndaughter) &&
        aligned(truth_daughter_e, ndaughter) &&
        aligned(truth_daughter_projection_eta, ndaughter) &&
        aligned(truth_daughter_projection_phi, ndaughter) &&
        aligned(truth_daughter_projection_valid, ndaughter) &&
        aligned(truth_daughter_in_acceptance, ndaughter);
    const bool event_ok = source_file_id == metadata_source_file_id &&
        event_uid == expected_uid && patch_side > 0 &&
        truth_aligned &&
        collection_aligned(split, store_patch, patch_side) &&
        collection_aligned(nosplit, store_patch, patch_side);
    malformed += event_ok ? 0 : 1;
  }

  std::cout << "check_tree - events/malformed = " << tree->GetEntries() << "/"
            << malformed << std::endl;
  std::cout << "check_tree - metadata schema/source/processed/written/"
               "invalid_truth/invalid_detector = "
            << schema_version << "/" << metadata_source_file_id << "/"
            << n_events_processed << "/" << n_events_written << "/"
            << n_events_invalid_truth << "/" << n_events_invalid_detector
            << std::endl;
  return malformed == 0 ? 0 : 6;
}
