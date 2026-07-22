#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace
{
template <class T>
bool bind(TTree* tree, const char* name, T* address)
{
  return tree->GetBranch(name) && tree->SetBranchAddress(name, address) >= 0;
}
}

int check_scored_tree(const char* input_path)
{
  std::unique_ptr<TFile> input(TFile::Open(input_path, "READ"));
  if (!input || input->IsZombie())
  {
    return 1;
  }
  TTree* tree = input->Get<TTree>("event_tree");
  TTree* metadata = input->Get<TTree>("metadata");
  if (!tree || !metadata)
  {
    std::cerr << "check_scored_tree - missing event_tree or metadata" << std::endl;
    return 2;
  }

  UInt_t metadata_source_file_id = 0;
  ULong64_t n_events_processed = 0;
  ULong64_t n_events_written = 0;
  const bool metadata_ok = metadata->GetEntries() == 1 &&
      bind(metadata, "source_file_id", &metadata_source_file_id) &&
      bind(metadata, "n_events_processed", &n_events_processed) &&
      bind(metadata, "n_events_written", &n_events_written) &&
      metadata->GetEntry(0) > 0 &&
      n_events_written == static_cast<ULong64_t>(tree->GetEntries());
  if (!metadata_ok)
  {
    std::cerr << "check_scored_tree - unreadable or inconsistent metadata" << std::endl;
    return 3;
  }

  UInt_t source_file_id = 0;
  UInt_t split_ncluster = 0;
  UInt_t nosplit_ncluster = 0;
  UInt_t nosplit_ntower = 0;
  std::vector<int>* truth_daughter_id = nullptr;
  std::vector<double>* truth_daughter_e = nullptr;
  std::vector<double>* split_e = nullptr;
  std::vector<float>* split_patch = nullptr;
  std::vector<unsigned int>* split_pair_i = nullptr;
  std::vector<double>* split_pair_mass = nullptr;
  std::vector<double>* nosplit_e = nullptr;
  std::vector<float>* nosplit_patch = nullptr;
  std::vector<unsigned int>* nosplit_pair_i = nullptr;
  std::vector<double>* nosplit_pair_mass = nullptr;
  std::vector<int>* tower_cluster_index = nullptr;
  std::vector<float>* bdt_score = nullptr;
  std::vector<unsigned char>* bdt_valid = nullptr;
  std::vector<float>* gamma_score = nullptr;
  std::vector<unsigned char>* gamma_valid = nullptr;

  bool ok = true;
  ok &= bind(tree, "source_file_id", &source_file_id);
  ok &= bind(tree, "split_ncluster", &split_ncluster);
  ok &= bind(tree, "nosplit_ncluster", &nosplit_ncluster);
  ok &= bind(tree, "nosplit_ntower", &nosplit_ntower);
  ok &= bind(tree, "truth_daughter_track_id", &truth_daughter_id);
  ok &= bind(tree, "truth_daughter_e", &truth_daughter_e);
  ok &= bind(tree, "split_cluster_e", &split_e);
  ok &= bind(tree, "split_cluster_shower_patch_e", &split_patch);
  ok &= bind(tree, "split_pair_cluster_i", &split_pair_i);
  ok &= bind(tree, "split_pair_m_gg", &split_pair_mass);
  ok &= bind(tree, "nosplit_cluster_e", &nosplit_e);
  ok &= bind(tree, "nosplit_cluster_shower_patch_e", &nosplit_patch);
  ok &= bind(tree, "nosplit_pair_cluster_i", &nosplit_pair_i);
  ok &= bind(tree, "nosplit_pair_m_gg", &nosplit_pair_mass);
  ok &= bind(tree, "nosplit_tower_cluster_index", &tower_cluster_index);
  ok &= bind(tree, "split_cluster_bdt_base_v3E_score", &bdt_score);
  ok &= bind(tree, "split_cluster_bdt_base_v3E_valid", &bdt_valid);
  ok &= bind(tree, "nosplit_cluster_p_gamma", &gamma_score);
  ok &= bind(tree, "nosplit_cluster_p_gamma_valid", &gamma_valid);
  if (!ok)
  {
    std::cerr << "check_scored_tree - missing required branch" << std::endl;
    return 4;
  }

  Long64_t malformed = 0;
  Long64_t valid_bdt = 0;
  Long64_t valid_gamma = 0;
  for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry)
  {
    tree->GetEntry(entry);
    if (source_file_id != metadata_source_file_id)
    {
      ++malformed;
    }
    const std::size_t split_pairs = static_cast<std::size_t>(split_ncluster) *
                                    (split_ncluster > 0 ? split_ncluster - 1U : 0U) / 2U;
    const std::size_t nosplit_pairs = static_cast<std::size_t>(nosplit_ncluster) *
                                      (nosplit_ncluster > 0 ? nosplit_ncluster - 1U : 0U) / 2U;
    bool event_ok = truth_daughter_id && truth_daughter_e &&
        truth_daughter_id->size() == truth_daughter_e->size() &&
        split_e && split_e->size() == split_ncluster &&
        split_patch && split_patch->size() == 49U * split_ncluster &&
        split_pair_i && split_pair_i->size() == split_pairs &&
        split_pair_mass && split_pair_mass->size() == split_pairs &&
        nosplit_e && nosplit_e->size() == nosplit_ncluster &&
        nosplit_patch && nosplit_patch->size() == 49U * nosplit_ncluster &&
        nosplit_pair_i && nosplit_pair_i->size() == nosplit_pairs &&
        nosplit_pair_mass && nosplit_pair_mass->size() == nosplit_pairs &&
        tower_cluster_index && tower_cluster_index->size() == nosplit_ntower &&
        bdt_score && bdt_score->size() == split_ncluster &&
        bdt_valid && bdt_valid->size() == split_ncluster &&
        gamma_score && gamma_score->size() == nosplit_ncluster &&
        gamma_valid && gamma_valid->size() == nosplit_ncluster;
    if (!event_ok)
    {
      ++malformed;
    }
    valid_bdt += bdt_valid ? std::count(bdt_valid->begin(), bdt_valid->end(), 1U) : 0;
    valid_gamma += gamma_valid ? std::count(gamma_valid->begin(), gamma_valid->end(), 1U) : 0;
  }

  std::cout << "check_scored_tree - events/valid_bdt/valid_gamma/malformed = "
            << tree->GetEntries() << "/" << valid_bdt << "/" << valid_gamma
            << "/" << malformed << std::endl;
  std::cout << "check_scored_tree - metadata source/processed/written = "
            << metadata_source_file_id << "/" << n_events_processed << "/"
            << n_events_written << std::endl;
  return malformed == 0 ? 0 : 5;
}
