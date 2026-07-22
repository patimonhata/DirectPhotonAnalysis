#include <ROOT/RVec.hxx>
#include <TMVA/RBDT.hxx>

#include <TBranch.h>
#include <TFile.h>
#include <TH1F.h>
#include <TKey.h>
#include <TNamed.h>
#include <TObject.h>
#include <TParameter.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace
{
constexpr float kInvalidScore = -999.F;

template <class T>
bool bind_branch(TTree *tree, const char *name, T *address)
{
  if (!tree->GetBranch(name))
  {
    std::cerr << "Missing branch: " << name << std::endl;
    return false;
  }
  return tree->SetBranchAddress(name, address) >= 0;
}

template <class T>
bool bind_vector_branch(TTree *tree, const char *name, std::vector<T> **address)
{
  return bind_branch(tree, name, address);
}
}  // namespace

int add_ppg15_bdt_to_pi0_tree(
    const char *input_file =
        "/sphenix/user/jaein213/photon/BDT/PPG15PhotonAN/mc_analysis/pi0_reconstruction_bdt/inputs/ForJaein_pi0_reconstruction_SPLIT_000000.root",
    const char *output_file =
        "/sphenix/user/jaein213/photon/BDT/PPG15PhotonAN/mc_analysis/pi0_reconstruction_bdt/output/ForJaein_pi0_reconstruction_SPLIT_000000_with_bdt.root",
    const char *model_file =
        "/sphenix/user/shuhangli/ppg12/FunWithxgboost/binned_models/model_base_v3E_split_single_tmva.root",
    bool overwrite_output = false)
{
  constexpr const char *kTreeName = "event_tree";
  constexpr const char *kModelKey = "myBDT";
  constexpr const char *kScoreBranch = "cluster_bdt_base_v3E_split";
  constexpr const char *kValidBranch = "cluster_bdt_base_v3E_split_valid";
  constexpr const char *kFeatureOrder =
      "ET,weta_cogx,wphi_cogx,vertex_z,cluster_eta,e11_over_e33,"
      "et1,et2,et3,et4,e32_over_e35";

  std::error_code input_path_error;
  std::error_code output_path_error;
  const auto normalized_input = std::filesystem::weakly_canonical(input_file, input_path_error);
  const auto normalized_output = std::filesystem::weakly_canonical(output_file, output_path_error);
  if (!input_path_error && !output_path_error && normalized_input == normalized_output)
  {
    std::cerr << "Input and output paths must differ; refusing to overwrite the input file."
              << std::endl;
    return 7;
  }
  if (!overwrite_output && std::filesystem::exists(output_file))
  {
    std::cerr << "Output file already exists; refusing to overwrite it: " << output_file << std::endl;
    return 9;
  }

  const std::filesystem::path output_path(output_file);
  if (!output_path.parent_path().empty())
  {
    std::error_code directory_error;
    std::filesystem::create_directories(output_path.parent_path(), directory_error);
    if (directory_error)
    {
      std::cerr << "Cannot create output directory: " << output_path.parent_path()
                << " (" << directory_error.message() << ")" << std::endl;
      return 8;
    }
  }

  std::unique_ptr<TFile> input(TFile::Open(input_file, "READ"));
  if (!input || input->IsZombie())
  {
    std::cerr << "Cannot open input file: " << input_file << std::endl;
    return 1;
  }

  auto *tree = dynamic_cast<TTree *>(input->Get(kTreeName));
  if (!tree)
  {
    std::cerr << "Cannot find tree '" << kTreeName << "'" << std::endl;
    return 2;
  }
  if (tree->GetBranch(kScoreBranch))
  {
    std::cerr << "Input already contains branch '" << kScoreBranch << "'" << std::endl;
    return 3;
  }

  UInt_t ncluster = 0;
  Double_t vertex_z = 0.;
  std::vector<double> *cluster_et = nullptr;
  std::vector<double> *cluster_eta = nullptr;
  std::vector<unsigned char> *shower_valid = nullptr;
  std::vector<float> *weta_cogx = nullptr;
  std::vector<float> *wphi_cogx = nullptr;
  std::vector<float> *e11_over_e33 = nullptr;
  std::vector<float> *e32_over_e35 = nullptr;
  std::vector<float> *et1 = nullptr;
  std::vector<float> *et2 = nullptr;
  std::vector<float> *et3 = nullptr;
  std::vector<float> *et4 = nullptr;

  bool branches_ok = true;
  branches_ok &= bind_branch(tree, "ncluster", &ncluster);
  branches_ok &= bind_branch(tree, "vertex_z", &vertex_z);
  branches_ok &= bind_vector_branch(tree, "cluster_et", &cluster_et);
  branches_ok &= bind_vector_branch(tree, "cluster_eta", &cluster_eta);
  branches_ok &= bind_vector_branch(tree, "cluster_shower_valid", &shower_valid);
  branches_ok &= bind_vector_branch(tree, "cluster_shower_w_eta_cogx", &weta_cogx);
  branches_ok &= bind_vector_branch(tree, "cluster_shower_w_phi_cogx", &wphi_cogx);
  branches_ok &= bind_vector_branch(tree, "cluster_shower_e11_over_e33", &e11_over_e33);
  branches_ok &= bind_vector_branch(tree, "cluster_shower_e32_over_e35", &e32_over_e35);
  branches_ok &= bind_vector_branch(tree, "cluster_shower_et1", &et1);
  branches_ok &= bind_vector_branch(tree, "cluster_shower_et2", &et2);
  branches_ok &= bind_vector_branch(tree, "cluster_shower_et3", &et3);
  branches_ok &= bind_vector_branch(tree, "cluster_shower_et4", &et4);
  if (!branches_ok)
  {
    return 4;
  }

  using RBDT = TMVA::Experimental::RBDT;
  std::unique_ptr<RBDT> bdt;
  try
  {
    bdt = std::make_unique<RBDT>(kModelKey, model_file);
  }
  catch (const std::exception &error)
  {
    std::cerr << "Failed to load BDT model: " << error.what() << std::endl;
    return 5;
  }

  std::unique_ptr<TFile> output(TFile::Open(output_file, "RECREATE"));
  if (!output || output->IsZombie())
  {
    std::cerr << "Cannot create output file: " << output_file << std::endl;
    return 6;
  }
  output->cd();

  auto *output_tree = tree->CloneTree(0);
  std::vector<float> scores;
  std::vector<unsigned char> scores_valid;
  output_tree->Branch(kScoreBranch, &scores);
  output_tree->Branch(kValidBranch, &scores_valid);

  TH1F h_score("h_cluster_bdt_base_v3E_split",
               "PPG15 base_v3E split BDT score;BDT score;Clusters", 100, 0., 1.);

  Long64_t total_clusters = 0;
  Long64_t scored_clusters = 0;
  Long64_t invalid_shape_clusters = 0;
  Long64_t substituted_features = 0;
  Long64_t malformed_events = 0;

  const Long64_t entries = tree->GetEntries();
  for (Long64_t entry = 0; entry < entries; ++entry)
  {
    tree->GetEntry(entry);
    scores.assign(ncluster, kInvalidScore);
    scores_valid.assign(ncluster, 0U);
    total_clusters += ncluster;

    const std::vector<size_t> sizes = {
        cluster_et ? cluster_et->size() : 0U,
        cluster_eta ? cluster_eta->size() : 0U,
        shower_valid ? shower_valid->size() : 0U,
        weta_cogx ? weta_cogx->size() : 0U,
        wphi_cogx ? wphi_cogx->size() : 0U,
        e11_over_e33 ? e11_over_e33->size() : 0U,
        e32_over_e35 ? e32_over_e35->size() : 0U,
        et1 ? et1->size() : 0U,
        et2 ? et2->size() : 0U,
        et3 ? et3->size() : 0U,
        et4 ? et4->size() : 0U};
    const size_t available = *std::min_element(sizes.begin(), sizes.end());
    if (available < ncluster)
    {
      ++malformed_events;
    }

    for (size_t cluster = 0; cluster < std::min<size_t>(ncluster, available); ++cluster)
    {
      std::vector<float> features = {
          static_cast<float>((*cluster_et)[cluster]),
          (*weta_cogx)[cluster],
          (*wphi_cogx)[cluster],
          static_cast<float>(vertex_z),
          static_cast<float>((*cluster_eta)[cluster]),
          (*e11_over_e33)[cluster],
          (*et1)[cluster],
          (*et2)[cluster],
          (*et3)[cluster],
          (*et4)[cluster],
          (*e32_over_e35)[cluster]};

      bool all_features_finite = true;
      for (float &feature : features)
      {
        if (!std::isfinite(feature))
        {
          feature = 0.F;
          all_features_finite = false;
          ++substituted_features;
        }
      }

      const auto result = bdt->Compute(features);
      if (result.empty() || !std::isfinite(result[0]))
      {
        continue;
      }

      scores[cluster] = result[0];
      scores_valid[cluster] = ((*shower_valid)[cluster] != 0U && all_features_finite) ? 1U : 0U;
      if ((*shower_valid)[cluster] == 0U)
      {
        ++invalid_shape_clusters;
      }
      ++scored_clusters;
      h_score.Fill(result[0]);
    }

    output_tree->Fill();
  }

  output->cd();
  output_tree->Write();
  h_score.Write();

  TIter next_key(input->GetListOfKeys());
  while (auto *key = dynamic_cast<TKey *>(next_key()))
  {
    if (std::string(key->GetName()) == kTreeName)
    {
      continue;
    }
    std::unique_ptr<TObject> object(key->ReadObj());
    if (object)
    {
      object->Write(key->GetName());
    }
  }

  TNamed("bdt_model_file", model_file).Write();
  TNamed("bdt_model_key", kModelKey).Write();
  TNamed("bdt_feature_order", kFeatureOrder).Write();
  TNamed("bdt_score_branch", kScoreBranch).Write();
  TNamed("bdt_input_compatibility",
         "Uses stored ShowerShapeCalculator algorithm_version=1 variables; exact source-formula "
         "comparison with PhotonAna was unavailable because ShowerShapeCalculator.cc was not readable.")
      .Write();
  TNamed("bdt_kinematic_warning",
         "The base_v3E model has documented performance bins beginning at ET=6 GeV; scores below "
         "that range are extrapolations and are not validated for physics use.")
      .Write();
  TParameter<Long64_t>("bdt_total_clusters", total_clusters).Write();
  TParameter<Long64_t>("bdt_scored_clusters", scored_clusters).Write();
  TParameter<Long64_t>("bdt_invalid_shape_clusters", invalid_shape_clusters).Write();
  TParameter<Long64_t>("bdt_substituted_features", substituted_features).Write();
  TParameter<Long64_t>("bdt_malformed_events", malformed_events).Write();

  output->Close();

  std::cout << "Wrote: " << output_file << '\n'
            << "Events: " << entries << '\n'
            << "Clusters: " << total_clusters << '\n'
            << "Scored clusters: " << scored_clusters << '\n'
            << "Invalid shower shapes: " << invalid_shape_clusters << '\n'
            << "Non-finite features replaced with zero: " << substituted_features << '\n'
            << "Malformed events: " << malformed_events << std::endl;
  return 0;
}
