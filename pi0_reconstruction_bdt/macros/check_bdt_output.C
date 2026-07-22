#include <TFile.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

int check_bdt_output(const char *file_name)
{
  TFile file(file_name, "READ");
  if (file.IsZombie())
  {
    std::cerr << "Cannot open file: " << file_name << std::endl;
    return 1;
  }

  auto *tree = dynamic_cast<TTree *>(file.Get("event_tree"));
  if (!tree)
  {
    std::cerr << "Missing event_tree" << std::endl;
    return 2;
  }

  constexpr const char *kScoreBranch = "cluster_bdt_base_v3E_split";
  constexpr const char *kValidBranch = "cluster_bdt_base_v3E_split_valid";
  if (!tree->GetBranch(kScoreBranch) || !tree->GetBranch(kValidBranch))
  {
    std::cerr << "Missing BDT score or validity branch" << std::endl;
    return 3;
  }

  UInt_t ncluster = 0;
  std::vector<float> *scores = nullptr;
  std::vector<unsigned char> *valid = nullptr;
  tree->SetBranchAddress("ncluster", &ncluster);
  tree->SetBranchAddress(kScoreBranch, &scores);
  tree->SetBranchAddress(kValidBranch, &valid);

  Long64_t total_clusters = 0;
  Long64_t valid_clusters = 0;
  Long64_t malformed_events = 0;
  double score_min = std::numeric_limits<double>::infinity();
  double score_max = -std::numeric_limits<double>::infinity();

  for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry)
  {
    tree->GetEntry(entry);
    if (!scores || !valid || scores->size() != ncluster || valid->size() != ncluster)
    {
      ++malformed_events;
      continue;
    }

    total_clusters += ncluster;
    for (size_t cluster = 0; cluster < ncluster; ++cluster)
    {
      if ((*valid)[cluster] == 0U)
      {
        continue;
      }
      ++valid_clusters;
      score_min = std::min(score_min, static_cast<double>((*scores)[cluster]));
      score_max = std::max(score_max, static_cast<double>((*scores)[cluster]));
    }
  }

  auto *model = dynamic_cast<TNamed *>(file.Get("bdt_model_file"));
  auto *features = dynamic_cast<TNamed *>(file.Get("bdt_feature_order"));
  auto *stored_total = dynamic_cast<TParameter<Long64_t> *>(file.Get("bdt_total_clusters"));

  std::cout << "File: " << file_name << '\n'
            << "Events: " << tree->GetEntries() << '\n'
            << "Clusters: " << total_clusters << '\n'
            << "Valid scores: " << valid_clusters << '\n'
            << "Malformed events: " << malformed_events << '\n';
  if (valid_clusters > 0)
  {
    std::cout << "Valid score range: " << score_min << " to " << score_max << '\n';
  }
  std::cout << "Model: " << (model ? model->GetTitle() : "metadata missing") << '\n'
            << "Features: " << (features ? features->GetTitle() : "metadata missing") << std::endl;

  if (malformed_events != 0 || (stored_total && stored_total->GetVal() != total_clusters))
  {
    return 4;
  }
  return 0;
}
