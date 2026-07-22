#if defined(__CLING__)
R__LOAD_LIBRARY(libPi0GammaOnnx.so)
#endif

#include "../src/Pi0GammaOnnx.h"

#include <TFile.h>
#include <TKey.h>
#include <TNamed.h>
#include <TObject.h>
#include <TParameter.h>
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

bool copy_other_keys(TFile& input, TFile& output, const std::string& tree_name)
{
  TIter next(input.GetListOfKeys());
  while (auto* key = dynamic_cast<TKey*>(next()))
  {
    if (key->GetName() == tree_name)
    {
      continue;
    }
    std::unique_ptr<TObject> object(key->ReadObj());
    if (!object)
    {
      std::cerr << "Cannot read input key: " << key->GetName() << std::endl;
      return false;
    }
    output.cd();
    if (auto* source_tree = dynamic_cast<TTree*>(object.get()))
    {
      std::unique_ptr<TTree> cloned_tree(source_tree->CloneTree(-1));
      if (!cloned_tree)
      {
        std::cerr << "Cannot clone input TTree: " << key->GetName() << std::endl;
        return false;
      }
      cloned_tree->SetName(key->GetName());
      cloned_tree->SetDirectory(&output);
      const bool written = cloned_tree->Write(key->GetName(), TObject::kOverwrite) > 0;
      cloned_tree->SetDirectory(nullptr);
      if (!written)
      {
        std::cerr << "Cannot write cloned TTree: " << key->GetName() << std::endl;
        return false;
      }
    }
    else if (object->Write(key->GetName(), TObject::kOverwrite) <= 0)
    {
      std::cerr << "Cannot copy input key: " << key->GetName() << std::endl;
      return false;
    }
  }
  return true;
}
}

int add_nosplit_gamma_onnx(
    const char* input_path,
    const char* output_path,
    const char* model_path =
        "/sphenix/user/ryotaro/DirectPhotonAnalysis/PhotonAnalysisTree/models/best_model.onnx")
{
  constexpr const char* tree_name = "event_tree";
  constexpr const char* score_branch = "nosplit_cluster_p_gamma";
  constexpr const char* valid_branch = "nosplit_cluster_p_gamma_valid";
  constexpr const char* global_features =
      "cluster_eta,log1p(cluster_energy/cosh(cluster_eta)),log1p(cluster_ntower)";
  constexpr const char* point_features =
      "tower_eta-cluster_eta,wrapped(tower_phi-cluster_phi),log1p(tower_energy),tower_energy/cluster_energy";

  if (!input_path || !output_path || !model_path)
  {
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
  std::vector<double>* tower_energy = nullptr;

  bool ok = true;
  ok &= bind(tree, "nosplit_ncluster", &ncluster);
  ok &= bind(tree, "nosplit_ntower", &ntower);
  ok &= bind(tree, "nosplit_cluster_ntower", &cluster_ntower);
  ok &= bind(tree, "nosplit_cluster_e", &cluster_e);
  ok &= bind(tree, "nosplit_cluster_eta", &cluster_eta);
  ok &= bind(tree, "nosplit_cluster_phi", &cluster_phi);
  ok &= bind(tree, "nosplit_tower_cluster_index", &tower_cluster_index);
  ok &= bind(tree, "nosplit_tower_eta", &tower_eta);
  ok &= bind(tree, "nosplit_tower_phi", &tower_phi);
  ok &= bind(tree, "nosplit_tower_energy", &tower_energy);
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
    const bool tower_lengths_ok = tower_cluster_index && tower_eta && tower_phi && tower_energy &&
        tower_cluster_index->size() == ntower && tower_eta->size() == ntower &&
        tower_phi->size() == ntower && tower_energy->size() == ntower;
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
      double tower_energy_sum = 0.0;
      std::vector<pi0gamma::Pi0GammaOnnx::TowerFeatures> points;
      points.reserve(indices.size());
      for (const std::size_t tower : indices)
      {
        const double tower_e = (*tower_energy)[tower];
        const double tower_eta_value = (*tower_eta)[tower];
        const double tower_phi_value = (*tower_phi)[tower];
        features_valid = features_valid && std::isfinite(tower_e) && tower_e > 0.0 &&
                         std::isfinite(tower_eta_value) && std::isfinite(tower_phi_value);
        if (!features_valid)
        {
          break;
        }
        tower_energy_sum += tower_e;
        points.push_back({
            static_cast<float>(tower_eta_value - eta),
            static_cast<float>(delta_phi(tower_phi_value, phi)),
            static_cast<float>(std::log1p(tower_e)),
            static_cast<float>(tower_e / energy)});
      }
      features_valid = features_valid &&
          std::abs(tower_energy_sum - energy) <= 1.0e-6 + 1.0e-5 * std::abs(energy);
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
  if (!copy_other_keys(*input, *output, tree_name))
  {
    output->Close();
    return 10;
  }
  TNamed("nosplit_gamma_onnx_model_file", model_path).Write();
  TNamed("nosplit_gamma_global_features", global_features).Write();
  TNamed("nosplit_gamma_point_features", point_features).Write();
  TNamed("nosplit_gamma_training_domain_warning",
         "The model was trained only on events with exactly one no-split cluster; multi-cluster scores require separate validation.").Write();
  TParameter<Long64_t>("nosplit_gamma_total_clusters", total_clusters).Write();
  TParameter<Long64_t>("nosplit_gamma_valid_scores", valid_scores).Write();
  TParameter<Long64_t>("nosplit_gamma_invalid_clusters", invalid_clusters).Write();
  TParameter<Long64_t>("nosplit_gamma_multi_cluster_events", multi_cluster_events).Write();
  TParameter<Long64_t>("nosplit_gamma_malformed_events", malformed_events).Write();
  output->Close();

  std::cout << "add_nosplit_gamma_onnx - events/clusters/valid/invalid/multicluster/malformed = "
            << entries << "/" << total_clusters << "/" << valid_scores << "/"
            << invalid_clusters << "/" << multi_cluster_events << "/" << malformed_events
            << std::endl;
  return malformed_events == 0 ? 0 : 9;
}
