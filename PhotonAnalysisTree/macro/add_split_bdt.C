#include <TMVA/RBDT.hxx>

#include <TBranch.h>
#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

namespace
{
constexpr float invalid_score = -999.0F;

template <class T>
bool bind(TTree* tree, const char* name, T* address)
{
  if (!tree->GetBranch(name))
  {
    std::cerr << "Missing branch: " << name << std::endl;
    return false;
  }
  return tree->SetBranchAddress(name, address) >= 0;
}

}

int add_split_bdt(
    const char* input_path,
    const char* model_path =
        "/sphenix/user/shuhangli/ppg12/FunWithxgboost/binned_models/model_base_v3E_split_single_tmva.root")
{
  constexpr const char* tree_name = "event_tree";
  constexpr const char* score_branch = "split_cluster_bdt_base_v3E_score";
  constexpr const char* valid_branch = "split_cluster_bdt_base_v3E_valid";
  constexpr const char* model_key = "myBDT";
  if (!input_path || !model_path)
  {
    std::cerr << "Input and model paths must be non-empty" << std::endl;
    return 1;
  }

  std::unique_ptr<TFile> file(TFile::Open(input_path, "UPDATE"));
  if (!file || file->IsZombie() || !file->IsWritable())
  {
    std::cerr << "Cannot open input for update: " << input_path << std::endl;
    return 2;
  }
  TTree* tree = file->Get<TTree>(tree_name);
  if (!tree)
  {
    std::cerr << "Missing event_tree" << std::endl;
    return 3;
  }
  if (tree->GetBranch(score_branch) || tree->GetBranch(valid_branch))
  {
    std::cerr << "SPLIT BDT score or valid branch already exists" << std::endl;
    return 4;
  }

  UInt_t ncluster = 0;
  Double_t vertex_z = 0.0;
  std::vector<double>* cluster_et = nullptr;
  std::vector<double>* cluster_eta = nullptr;
  std::vector<unsigned char>* shower_valid = nullptr;
  std::vector<float>* weta = nullptr;
  std::vector<float>* wphi = nullptr;
  std::vector<float>* e11_over_e33 = nullptr;
  std::vector<float>* e32_over_e35 = nullptr;
  std::vector<float>* et1 = nullptr;
  std::vector<float>* et2 = nullptr;
  std::vector<float>* et3 = nullptr;
  std::vector<float>* et4 = nullptr;

  bool ok = true;
  ok &= bind(tree, "split_ncluster", &ncluster);
  ok &= bind(tree, "vertex_z", &vertex_z);
  ok &= bind(tree, "split_cluster_et", &cluster_et);
  ok &= bind(tree, "split_cluster_eta", &cluster_eta);
  ok &= bind(tree, "split_cluster_shower_valid", &shower_valid);
  ok &= bind(tree, "split_cluster_shower_w_eta_cogx", &weta);
  ok &= bind(tree, "split_cluster_shower_w_phi_cogx", &wphi);
  ok &= bind(tree, "split_cluster_shower_e11_over_e33", &e11_over_e33);
  ok &= bind(tree, "split_cluster_shower_e32_over_e35", &e32_over_e35);
  ok &= bind(tree, "split_cluster_shower_et1", &et1);
  ok &= bind(tree, "split_cluster_shower_et2", &et2);
  ok &= bind(tree, "split_cluster_shower_et3", &et3);
  ok &= bind(tree, "split_cluster_shower_et4", &et4);
  if (!ok)
  {
    return 5;
  }

  std::unique_ptr<TMVA::Experimental::RBDT> bdt;
  try
  {
    bdt = std::make_unique<TMVA::Experimental::RBDT>(model_key, model_path);
  }
  catch (const std::exception& error)
  {
    std::cerr << "Cannot load BDT: " << error.what() << std::endl;
    return 6;
  }

  file->cd();
  std::vector<float> scores;
  std::vector<unsigned char> valid;
  TBranch* score_output = tree->Branch(score_branch, &scores);
  TBranch* valid_output = tree->Branch(valid_branch, &valid);
  if (!score_output || !valid_output)
  {
    std::cerr << "Cannot create SPLIT BDT output branches" << std::endl;
    return 7;
  }

  Long64_t total_clusters = 0;
  Long64_t valid_scores = 0;
  Long64_t malformed_events = 0;
  for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry)
  {
    tree->GetEntry(entry);
    scores.assign(ncluster, invalid_score);
    valid.assign(ncluster, 0U);
    total_clusters += ncluster;
    const std::vector<std::size_t> sizes = {
        cluster_et ? cluster_et->size() : 0U,
        cluster_eta ? cluster_eta->size() : 0U,
        shower_valid ? shower_valid->size() : 0U,
        weta ? weta->size() : 0U,
        wphi ? wphi->size() : 0U,
        e11_over_e33 ? e11_over_e33->size() : 0U,
        e32_over_e35 ? e32_over_e35->size() : 0U,
        et1 ? et1->size() : 0U,
        et2 ? et2->size() : 0U,
        et3 ? et3->size() : 0U,
        et4 ? et4->size() : 0U};
    const std::size_t available = *std::min_element(sizes.begin(), sizes.end());
    if (available != ncluster)
    {
      ++malformed_events;
    }

    for (std::size_t cluster = 0; cluster < std::min<std::size_t>(ncluster, available); ++cluster)
    {
      std::vector<float> features = {
          static_cast<float>((*cluster_et)[cluster]),
          (*weta)[cluster], (*wphi)[cluster], static_cast<float>(vertex_z),
          static_cast<float>((*cluster_eta)[cluster]), (*e11_over_e33)[cluster],
          (*et1)[cluster], (*et2)[cluster], (*et3)[cluster], (*et4)[cluster],
          (*e32_over_e35)[cluster]};
      const bool finite = std::all_of(features.begin(), features.end(),
                                      [](float value) { return std::isfinite(value); });
      if (!finite)
      {
        continue;
      }
      const auto result = bdt->Compute(features);
      if (result.empty() || !std::isfinite(result[0]))
      {
        continue;
      }
      scores[cluster] = result[0];
      valid[cluster] = (*shower_valid)[cluster] ? 1U : 0U;
      valid_scores += valid[cluster] ? 1 : 0;
    }
    if (score_output->Fill() < 0 || valid_output->Fill() < 0)
    {
      std::cerr << "Failed to fill SPLIT BDT output branches at entry " << entry << std::endl;
      return 9;
    }
  }

  const Long64_t entries = tree->GetEntries();
  if (score_output->GetEntries() != entries || valid_output->GetEntries() != entries)
  {
    std::cerr << "SPLIT BDT output branch entry count mismatch" << std::endl;
    return 9;
  }
  file->cd();
  if (tree->Write(tree_name, TObject::kOverwrite) <= 0)
  {
    std::cerr << "Failed to write updated event_tree" << std::endl;
    return 9;
  }
  file->Flush();
  file->Close();

  std::cout << "add_split_bdt - events/clusters/valid/malformed = "
            << entries << "/" << total_clusters << "/"
            << valid_scores << "/" << malformed_events << std::endl;
  return malformed_events == 0 ? 0 : 8;
}
