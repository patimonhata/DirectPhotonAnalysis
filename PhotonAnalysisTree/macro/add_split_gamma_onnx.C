#if defined(__CLING__)
R__LOAD_LIBRARY(libPi0GammaOnnx.so)
#endif

#include "../src/Pi0GammaOnnx.h"

#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <memory>
#include <string>
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

double delta_phi(double lhs, double rhs)
{
  constexpr double pi = 3.14159265358979323846;
  double value = lhs - rhs;
  while (value > pi) value -= 2.0 * pi;
  while (value <= -pi) value += 2.0 * pi;
  return value;
}

bool copy_metadata_tree(TFile& input, TFile& output)
{
  TTree* source_tree = input.Get<TTree>("metadata");
  if (!source_tree)
  {
    std::cerr << "Missing metadata TTree" << std::endl;
    return false;
  }
  output.cd();
  std::unique_ptr<TTree> cloned_tree(source_tree->CloneTree(-1));
  if (!cloned_tree)
  {
    std::cerr << "Cannot clone metadata TTree" << std::endl;
    return false;
  }
  cloned_tree->SetDirectory(&output);
  const bool written = cloned_tree->Write("metadata", TObject::kOverwrite) > 0;
  cloned_tree->SetDirectory(nullptr);
  return written;
}
}

int add_split_gamma_onnx(
    const char* input_path,
    const char* output_path,
    const char* model_path = nullptr)
{
  constexpr const char* tree_name = "event_tree";
  constexpr const char* score_branch = "split_cluster_p_gamma";
  constexpr const char* valid_branch = "split_cluster_p_gamma_valid";
  if (!input_path || !output_path || !model_path)
  {
    std::cerr << "Input, output, and split ONNX model paths must be provided" << std::endl;
    return 1;
  }
  std::error_code input_error;
  std::error_code output_error;
  const auto input_canonical = std::filesystem::weakly_canonical(input_path, input_error);
  const auto output_canonical = std::filesystem::weakly_canonical(output_path, output_error);
  if (!input_error && !output_error && input_canonical == output_canonical)
  {
    std::cerr << "Input and output must differ" << std::endl;
    return 2;
  }
  if (std::filesystem::exists(output_path))
  {
    std::cerr << "Refusing to overwrite output: " << output_path << std::endl;
    return 3;
  }
  const std::filesystem::path output_fs(output_path);
  if (!output_fs.parent_path().empty())
  {
    std::filesystem::create_directories(output_fs.parent_path());
  }

  std::unique_ptr<TFile> input(TFile::Open(input_path, "READ"));
  if (!input || input->IsZombie())
  {
    return 4;
  }
  TTree* tree = input->Get<TTree>(tree_name);
  if (!tree || tree->GetBranch(score_branch))
  {
    std::cerr << "Missing event_tree or gamma branch already exists" << std::endl;
    return 5;
  }

  UInt_t ncluster = 0;
  UInt_t ntower = 0;
  std::vector<int>* cluster_ntower = nullptr;
  std::vector<double>* cluster_e = nullptr;
  std::vector<double>* cluster_eta = nullptr;
  std::vector<double>* cluster_phi = nullptr;
  std::vector<int>* tower_cluster_index = nullptr;
  std::vector<double>* tower_eta = nullptr;
  std::vector<double>* tower_phi = nullptr;
  // For split clusters, use the energy contribution stored in RawCluster rather than
  // the full calibrated tower energy, since one tower may be shared between clusters.
  std::vector<double>* tower_assigned_energy = nullptr;

  bool ok = true;
  ok &= bind(tree, "split_ncluster", &ncluster);
  ok &= bind(tree, "split_ntower", &ntower);
  ok &= bind(tree, "split_cluster_ntower", &cluster_ntower);
  ok &= bind(tree, "split_cluster_e", &cluster_e);
  ok &= bind(tree, "split_cluster_eta", &cluster_eta);
  ok &= bind(tree, "split_cluster_phi", &cluster_phi);
  ok &= bind(tree, "split_tower_cluster_index", &tower_cluster_index);
  ok &= bind(tree, "split_tower_eta", &tower_eta);
  ok &= bind(tree, "split_tower_phi", &tower_phi);
  ok &= bind(tree, "split_tower_cluster_value", &tower_assigned_energy);
  if (!ok)
  {
    return 6;
  }

  std::unique_ptr<pi0gamma::Pi0GammaOnnx> classifier;
  try
  {
    classifier = std::make_unique<pi0gamma::Pi0GammaOnnx>(model_path);
  }
  catch (const std::exception& error)
  {
    std::cerr << "Cannot load ONNX model: " << error.what() << std::endl;
    return 7;
  }

  std::unique_ptr<TFile> output(TFile::Open(output_path, "RECREATE"));
  if (!output || output->IsZombie())
  {
    return 8;
  }
  output->cd();
  TTree* output_tree = tree->CloneTree(0);
  std::vector<float> scores;
  std::vector<unsigned char> valid;
  output_tree->Branch(score_branch, &scores);
  output_tree->Branch(valid_branch, &valid);

  Long64_t total_clusters = 0;
  Long64_t valid_scores = 0;
  Long64_t invalid_clusters = 0;
  Long64_t multi_cluster_events = 0;
  Long64_t malformed_events = 0;
  const Long64_t entries = tree->GetEntries();
  for (Long64_t entry = 0; entry < entries; ++entry)
  {
    tree->GetEntry(entry);
    scores.assign(ncluster, invalid_score);
    valid.assign(ncluster, 0U);
    total_clusters += ncluster;
    multi_cluster_events += ncluster > 1 ? 1 : 0;

    const bool cluster_lengths_ok = cluster_ntower && cluster_e && cluster_eta && cluster_phi &&
        cluster_ntower->size() == ncluster && cluster_e->size() == ncluster &&
        cluster_eta->size() == ncluster && cluster_phi->size() == ncluster;
    const bool tower_lengths_ok = tower_cluster_index && tower_eta && tower_phi && tower_assigned_energy &&
        tower_cluster_index->size() == ntower && tower_eta->size() == ntower &&
        tower_phi->size() == ntower && tower_assigned_energy->size() == ntower;
    if (!cluster_lengths_ok || !tower_lengths_ok)
    {
      ++malformed_events;
      invalid_clusters += ncluster;
      output_tree->Fill();
      continue;
    }

    std::vector<std::vector<std::size_t>> tower_indices(ncluster);
    bool indices_ok = true;
    for (std::size_t tower = 0; tower < ntower; ++tower)
    {
      const int index = (*tower_cluster_index)[tower];
      if (index < 0 || static_cast<unsigned int>(index) >= ncluster)
      {
        indices_ok = false;
        break;
      }
      tower_indices[static_cast<std::size_t>(index)].push_back(tower);
    }
    if (!indices_ok)
    {
      ++malformed_events;
      invalid_clusters += ncluster;
      output_tree->Fill();
      continue;
    }

    for (std::size_t cluster = 0; cluster < ncluster; ++cluster)
    {
      const double energy = (*cluster_e)[cluster];
      const double eta = (*cluster_eta)[cluster];
      const double phi = (*cluster_phi)[cluster];
      const int declared_ntower = (*cluster_ntower)[cluster];
      const auto& indices = tower_indices[cluster];
      bool features_valid = std::isfinite(energy) && energy > 0.0 &&
          std::isfinite(eta) && std::isfinite(phi) && declared_ntower > 0 &&
          static_cast<std::size_t>(declared_ntower) == indices.size();
      double tower_assigned_energy_sum = 0.0;
      std::vector<pi0gamma::Pi0GammaOnnx::TowerFeatures> points;
      points.reserve(indices.size());
      for (const std::size_t tower : indices)
      {
        const double tower_e = (*tower_assigned_energy)[tower];
        const double tower_eta_value = (*tower_eta)[tower];
        const double tower_phi_value = (*tower_phi)[tower];
        features_valid = features_valid && std::isfinite(tower_e) && tower_e > 0.0 &&
                         std::isfinite(tower_eta_value) && std::isfinite(tower_phi_value);
        if (!features_valid)
        {
          break;
        }
        tower_assigned_energy_sum += tower_e;
        points.push_back({
            static_cast<float>(tower_eta_value - eta),
            static_cast<float>(delta_phi(tower_phi_value, phi)),
            static_cast<float>(std::log1p(tower_e)),
            static_cast<float>(tower_e / energy)});
      }
      features_valid = features_valid &&
          std::abs(tower_assigned_energy_sum - energy) <= 1.0e-6 + 1.0e-5 * std::abs(energy);
      if (!features_valid)
      {
        ++invalid_clusters;
        continue;
      }

      const double transverse_energy = energy / std::cosh(eta);
      try
      {
        scores[cluster] = classifier->predict(
            {static_cast<float>(eta),
             static_cast<float>(std::log1p(transverse_energy)),
             static_cast<float>(std::log1p(declared_ntower))},
            points);
        valid[cluster] = 1U;
        ++valid_scores;
      }
      catch (const std::exception& error)
      {
        if (invalid_clusters < 5)
        {
          std::cerr << "Inference failed at entry " << entry << ", cluster " << cluster
                    << ": " << error.what() << std::endl;
        }
        ++invalid_clusters;
      }
    }
    output_tree->Fill();
  }

  output->cd();
  output_tree->Write();
  if (!copy_metadata_tree(*input, *output))
  {
    output->Close();
    return 10;
  }
  output->Close();

  std::cout << "add_split_gamma_onnx - events/clusters/valid/invalid/multicluster/malformed = "
            << entries << "/" << total_clusters << "/" << valid_scores << "/"
            << invalid_clusters << "/" << multi_cluster_events << "/" << malformed_events
            << std::endl;
  return malformed_events == 0 ? 0 : 9;
}
