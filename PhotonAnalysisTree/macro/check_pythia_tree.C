#include <TFile.h>
#include <TTree.h>

#include <algorithm>
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
  return tree->GetBranch(name.c_str()) && tree->SetBranchAddress(name.c_str(), address) >= 0;
}

bool check_collection(TTree* tree, const std::string& prefix)
{
  unsigned int ncluster = 0;
  std::vector<unsigned char>* valid = nullptr;
  std::vector<float>* total_edep = nullptr;
  std::vector<unsigned int>* n_contributor = nullptr;
  std::vector<int>* dominant_track = nullptr;
  std::vector<float>* dominant_fraction = nullptr;
  std::vector<unsigned int>* offset = nullptr;
  std::vector<unsigned int>* contributor_cluster = nullptr;
  std::vector<int>* contributor_track = nullptr;
  std::vector<float>* contributor_edep = nullptr;
  std::vector<float>* contributor_fraction = nullptr;

  const std::string truth = prefix + "_cluster_truth_";
  bool ok = true;
  ok &= bind_branch(tree, prefix + "_ncluster", &ncluster);
  ok &= bind_branch(tree, truth + "valid", &valid);
  ok &= bind_branch(tree, truth + "total_edep", &total_edep);
  ok &= bind_branch(tree, truth + "n_contributor", &n_contributor);
  ok &= bind_branch(tree, truth + "dominant_g4_track_id", &dominant_track);
  ok &= bind_branch(tree, truth + "dominant_fraction", &dominant_fraction);
  ok &= bind_branch(tree, truth + "contributor_offset", &offset);
  ok &= bind_branch(tree, truth + "contributor_cluster_index", &contributor_cluster);
  ok &= bind_branch(tree, truth + "contributor_g4_track_id", &contributor_track);
  ok &= bind_branch(tree, truth + "contributor_edep", &contributor_edep);
  ok &= bind_branch(tree, truth + "contributor_fraction", &contributor_fraction);
  if (!ok)
  {
    std::cerr << "check_pythia_tree - missing " << prefix << " truth branch" << std::endl;
    return false;
  }

  for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry)
  {
    tree->GetEntry(entry);
    const bool aligned = valid && total_edep && n_contributor && dominant_track &&
        dominant_fraction && offset && contributor_cluster && contributor_track &&
        contributor_edep && contributor_fraction &&
        valid->size() == ncluster && total_edep->size() == ncluster &&
        n_contributor->size() == ncluster && dominant_track->size() == ncluster &&
        dominant_fraction->size() == ncluster && offset->size() == ncluster + 1U &&
        contributor_cluster->size() == contributor_track->size() &&
        contributor_track->size() == contributor_edep->size() &&
        contributor_edep->size() == contributor_fraction->size();
    if (!aligned || offset->front() != 0U || offset->back() != contributor_track->size())
    {
      std::cerr << "check_pythia_tree - " << prefix << " alignment failed at entry "
                << entry << std::endl;
      return false;
    }
    for (unsigned int cluster = 0; cluster < ncluster; ++cluster)
    {
      if ((*offset)[cluster] > (*offset)[cluster + 1U] ||
          (*offset)[cluster + 1U] - (*offset)[cluster] != (*n_contributor)[cluster])
      {
        std::cerr << "check_pythia_tree - " << prefix << " offset failed at entry "
                  << entry << ", cluster " << cluster << std::endl;
        return false;
      }
      float edep_sum = 0.0F;
      float fraction_sum = 0.0F;
      for (unsigned int index = (*offset)[cluster]; index < (*offset)[cluster + 1U]; ++index)
      {
        if ((*contributor_cluster)[index] != cluster || !std::isfinite((*contributor_edep)[index]) ||
            !std::isfinite((*contributor_fraction)[index]) || (*contributor_edep)[index] < 0.0F ||
            (*contributor_fraction)[index] < 0.0F)
        {
          return false;
        }
        edep_sum += (*contributor_edep)[index];
        fraction_sum += (*contributor_fraction)[index];
      }
      const float tolerance = 1.0e-4F * std::max(1.0F, (*total_edep)[cluster]);
      if (std::abs(edep_sum - (*total_edep)[cluster]) > tolerance ||
          (edep_sum > 0.0F && std::abs(fraction_sum - 1.0F) > 1.0e-4F))
      {
        return false;
      }
    }
  }
  tree->ResetBranchAddresses();
  return true;
}
}

int check_pythia_tree(const std::string& file_name)
{
  std::unique_ptr<TFile> file(TFile::Open(file_name.c_str(), "READ"));
  if (!file || file->IsZombie())
  {
    std::cerr << "check_pythia_tree - cannot open " << file_name << std::endl;
    return 1;
  }
  auto* tree = dynamic_cast<TTree*>(file->Get("event_tree"));
  auto* metadata = dynamic_cast<TTree*>(file->Get("metadata"));
  if (!tree || !metadata || metadata->GetEntries() != 1)
  {
    std::cerr << "check_pythia_tree - missing event_tree or metadata" << std::endl;
    return 2;
  }
  int schema_version = 0;
  std::string* sample_type = nullptr;
  if (!bind_branch(metadata, "schema_version", &schema_version) ||
      !bind_branch(metadata, "sample_type", &sample_type))
  {
    return 3;
  }
  metadata->GetEntry(0);
  if (schema_version != 4 || !sample_type || *sample_type != "pythia")
  {
    std::cerr << "check_pythia_tree - metadata identifies a different schema" << std::endl;
    return 4;
  }
  metadata->ResetBranchAddresses();
  if (!check_collection(tree, "split") || !check_collection(tree, "nosplit"))
  {
    return 5;
  }
  std::cout << "check_pythia_tree - validated " << tree->GetEntries() << " events" << std::endl;
  return 0;
}
