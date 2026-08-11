#include <TFile.h>
#include <TKey.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <set>
#include <string>
#include <vector>

namespace
{
template <class T>
bool bind_branch(TTree* tree, const char* name, T* address)
{
  if (!tree->GetBranch(name))
  {
    std::cerr << "check_pythia_truth_spectrum_tree - missing branch: " << name << std::endl;
    return false;
  }
  return tree->SetBranchAddress(name, address) >= 0;
}

bool has_only_expected_trees(TFile* file)
{
  bool has_event_tree = false;
  bool has_metadata = false;
  bool valid = file != nullptr;
  TIter next(file ? file->GetListOfKeys() : nullptr);
  while (auto* key = dynamic_cast<TKey*>(next()))
  {
    const std::string name = key->GetName();
    const bool is_tree = std::string(key->GetClassName()) == "TTree";
    has_event_tree |= name == "event_tree" && is_tree;
    has_metadata |= name == "metadata" && is_tree;
    valid &= (name == "event_tree" || name == "metadata") && is_tree;
  }
  return valid && has_event_tree && has_metadata;
}

template <class T>
bool aligned(const std::vector<T>* values, const std::size_t expected)
{
  return values && values->size() == expected;
}

bool finite_kinematics(const float e, const float pt, const float eta, const float phi)
{
  return std::isfinite(e) && std::isfinite(pt) && std::isfinite(eta) &&
      std::isfinite(phi) && e >= 0.0F && pt >= 0.0F;
}
}

int check_pythia_truth_spectrum_tree(const char* input_path)
{
  if (!input_path)
  {
    return 1;
  }
  std::unique_ptr<TFile> input(TFile::Open(input_path, "READ"));
  if (!input || input->IsZombie() || !has_only_expected_trees(input.get()))
  {
    std::cerr << "check_pythia_truth_spectrum_tree - expected only event_tree and metadata TTree keys"
              << std::endl;
    return 2;
  }
  TTree* tree = input->Get<TTree>("event_tree");
  TTree* metadata = input->Get<TTree>("metadata");
  if (!tree || !metadata || metadata->GetEntries() != 1)
  {
    return 2;
  }

  int schema_version = 0;
  std::string* sample_type = nullptr;
  UInt_t metadata_source_file_id = 0U;
  ULong64_t n_events_processed = 0ULL;
  ULong64_t n_events_written = 0ULL;
  ULong64_t metadata_n_particle = 0ULL;
  ULong64_t metadata_n_photon = 0ULL;
  ULong64_t metadata_n_pi0 = 0ULL;
  ULong64_t metadata_n_pi0_decay_photon = 0ULL;
  bool metadata_ok = true;
  metadata_ok &= bind_branch(metadata, "schema_version", &schema_version);
  metadata_ok &= bind_branch(metadata, "sample_type", &sample_type);
  metadata_ok &= bind_branch(metadata, "source_file_id", &metadata_source_file_id);
  metadata_ok &= bind_branch(metadata, "n_events_processed", &n_events_processed);
  metadata_ok &= bind_branch(metadata, "n_events_written", &n_events_written);
  metadata_ok &= bind_branch(metadata, "n_hepmc_particle_record", &metadata_n_particle);
  metadata_ok &= bind_branch(metadata, "n_final_photon", &metadata_n_photon);
  metadata_ok &= bind_branch(metadata, "n_terminal_pi0", &metadata_n_pi0);
  metadata_ok &= bind_branch(metadata, "n_g4_pi0_decay_photon", &metadata_n_pi0_decay_photon);
  metadata_ok &= metadata->GetEntry(0) > 0;
  metadata_ok &= schema_version == 1 && sample_type &&
      *sample_type == "pythia_truth_spectrum";
  metadata_ok &= n_events_written == static_cast<ULong64_t>(tree->GetEntries());
  metadata_ok &= n_events_processed >= n_events_written;
  if (!metadata_ok)
  {
    std::cerr << "check_pythia_truth_spectrum_tree - invalid metadata" << std::endl;
    return 3;
  }
  metadata->ResetBranchAddresses();

  UInt_t source_file_id = 0U;
  UInt_t event_in_file = 0U;
  ULong64_t event_uid = 0ULL;
  UChar_t event_weight_valid = 0U;
  double event_weight = 1.0;
  std::vector<double>* event_weights = nullptr;
  UInt_t hepmc_n_particle_record = 0U;

  UInt_t truth_photon_n = 0U;
  std::vector<int>* photon_barcode = nullptr;
  std::vector<int>* photon_status = nullptr;
  std::vector<unsigned char>* photon_kinematics_valid = nullptr;
  std::vector<float>* photon_e = nullptr;
  std::vector<float>* photon_pt = nullptr;
  std::vector<float>* photon_eta = nullptr;
  std::vector<float>* photon_phi = nullptr;
  std::vector<unsigned char>* photon_classification_valid = nullptr;
  std::vector<int>* photon_category = nullptr;
  std::vector<int>* photon_immediate_parent_count = nullptr;
  std::vector<int>* photon_immediate_parent_pdg = nullptr;
  std::vector<unsigned char>* photon_copy_chain_valid = nullptr;
  std::vector<unsigned int>* photon_copy_depth = nullptr;
  std::vector<int>* photon_origin_parent_count = nullptr;
  std::vector<int>* photon_origin_parent_pdg = nullptr;
  std::vector<int>* photon_origin_parent_barcode = nullptr;

  UInt_t truth_pi0_n = 0U;
  std::vector<int>* pi0_barcode = nullptr;
  std::vector<int>* pi0_status = nullptr;
  std::vector<unsigned char>* pi0_kinematics_valid = nullptr;
  std::vector<float>* pi0_e = nullptr;
  std::vector<float>* pi0_pt = nullptr;
  std::vector<float>* pi0_eta = nullptr;
  std::vector<float>* pi0_phi = nullptr;
  std::vector<unsigned int>* pi0_hepmc_direct_photon_count = nullptr;

  UInt_t truth_pi0_decay_photon_n = 0U;
  std::vector<int>* pi0_decay_photon_g4_track_id = nullptr;
  std::vector<int>* pi0_decay_photon_parent_g4_track_id = nullptr;
  std::vector<int>* pi0_decay_photon_parent_hepmc_barcode = nullptr;
  std::vector<unsigned char>* pi0_decay_photon_kinematics_valid = nullptr;
  std::vector<float>* pi0_decay_photon_e = nullptr;
  std::vector<float>* pi0_decay_photon_pt = nullptr;
  std::vector<float>* pi0_decay_photon_eta = nullptr;
  std::vector<float>* pi0_decay_photon_phi = nullptr;

  bool ok = true;
  ok &= bind_branch(tree, "source_file_id", &source_file_id);
  ok &= bind_branch(tree, "event_in_file", &event_in_file);
  ok &= bind_branch(tree, "event_uid", &event_uid);
  ok &= bind_branch(tree, "event_weight_valid", &event_weight_valid);
  ok &= bind_branch(tree, "event_weight", &event_weight);
  ok &= bind_branch(tree, "event_weights", &event_weights);
  ok &= bind_branch(tree, "hepmc_n_particle_record", &hepmc_n_particle_record);
  ok &= bind_branch(tree, "truth_photon_n", &truth_photon_n);
  ok &= bind_branch(tree, "truth_photon_barcode", &photon_barcode);
  ok &= bind_branch(tree, "truth_photon_status", &photon_status);
  ok &= bind_branch(tree, "truth_photon_kinematics_valid", &photon_kinematics_valid);
  ok &= bind_branch(tree, "truth_photon_e", &photon_e);
  ok &= bind_branch(tree, "truth_photon_pt", &photon_pt);
  ok &= bind_branch(tree, "truth_photon_eta", &photon_eta);
  ok &= bind_branch(tree, "truth_photon_phi", &photon_phi);
  ok &= bind_branch(tree, "truth_photon_classification_valid", &photon_classification_valid);
  ok &= bind_branch(tree, "truth_photon_category", &photon_category);
  ok &= bind_branch(tree, "truth_photon_immediate_parent_count", &photon_immediate_parent_count);
  ok &= bind_branch(tree, "truth_photon_immediate_parent_pdg", &photon_immediate_parent_pdg);
  ok &= bind_branch(tree, "truth_photon_copy_chain_valid", &photon_copy_chain_valid);
  ok &= bind_branch(tree, "truth_photon_copy_depth", &photon_copy_depth);
  ok &= bind_branch(tree, "truth_photon_origin_parent_count", &photon_origin_parent_count);
  ok &= bind_branch(tree, "truth_photon_origin_parent_pdg", &photon_origin_parent_pdg);
  ok &= bind_branch(tree, "truth_photon_origin_parent_barcode", &photon_origin_parent_barcode);
  ok &= bind_branch(tree, "truth_pi0_n", &truth_pi0_n);
  ok &= bind_branch(tree, "truth_pi0_barcode", &pi0_barcode);
  ok &= bind_branch(tree, "truth_pi0_status", &pi0_status);
  ok &= bind_branch(tree, "truth_pi0_kinematics_valid", &pi0_kinematics_valid);
  ok &= bind_branch(tree, "truth_pi0_e", &pi0_e);
  ok &= bind_branch(tree, "truth_pi0_pt", &pi0_pt);
  ok &= bind_branch(tree, "truth_pi0_eta", &pi0_eta);
  ok &= bind_branch(tree, "truth_pi0_phi", &pi0_phi);
  ok &= bind_branch(tree, "truth_pi0_hepmc_direct_photon_count", &pi0_hepmc_direct_photon_count);
  ok &= bind_branch(tree, "truth_pi0_decay_photon_n", &truth_pi0_decay_photon_n);
  ok &= bind_branch(tree, "truth_pi0_decay_photon_g4_track_id", &pi0_decay_photon_g4_track_id);
  ok &= bind_branch(tree, "truth_pi0_decay_photon_parent_g4_track_id", &pi0_decay_photon_parent_g4_track_id);
  ok &= bind_branch(tree, "truth_pi0_decay_photon_parent_hepmc_barcode", &pi0_decay_photon_parent_hepmc_barcode);
  ok &= bind_branch(tree, "truth_pi0_decay_photon_kinematics_valid", &pi0_decay_photon_kinematics_valid);
  ok &= bind_branch(tree, "truth_pi0_decay_photon_e", &pi0_decay_photon_e);
  ok &= bind_branch(tree, "truth_pi0_decay_photon_pt", &pi0_decay_photon_pt);
  ok &= bind_branch(tree, "truth_pi0_decay_photon_eta", &pi0_decay_photon_eta);
  ok &= bind_branch(tree, "truth_pi0_decay_photon_phi", &pi0_decay_photon_phi);
  if (!ok)
  {
    return 4;
  }

  ULong64_t observed_n_particle = 0ULL;
  ULong64_t observed_n_photon = 0ULL;
  ULong64_t observed_n_pi0 = 0ULL;
  ULong64_t observed_n_pi0_decay_photon = 0ULL;
  for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry)
  {
    if (tree->GetEntry(entry) <= 0)
    {
      return 5;
    }
    const std::size_t n_photon = truth_photon_n;
    const std::size_t n_pi0 = truth_pi0_n;
    const std::size_t n_pi0_decay_photon = truth_pi0_decay_photon_n;
    const bool photon_aligned =
        aligned(photon_barcode, n_photon) && aligned(photon_status, n_photon) &&
        aligned(photon_kinematics_valid, n_photon) && aligned(photon_e, n_photon) &&
        aligned(photon_pt, n_photon) && aligned(photon_eta, n_photon) &&
        aligned(photon_phi, n_photon) && aligned(photon_classification_valid, n_photon) &&
        aligned(photon_category, n_photon) && aligned(photon_immediate_parent_count, n_photon) &&
        aligned(photon_immediate_parent_pdg, n_photon) && aligned(photon_copy_chain_valid, n_photon) &&
        aligned(photon_copy_depth, n_photon) && aligned(photon_origin_parent_count, n_photon) &&
        aligned(photon_origin_parent_pdg, n_photon) && aligned(photon_origin_parent_barcode, n_photon);
    const bool pi0_aligned =
        aligned(pi0_barcode, n_pi0) && aligned(pi0_status, n_pi0) &&
        aligned(pi0_kinematics_valid, n_pi0) && aligned(pi0_e, n_pi0) &&
        aligned(pi0_pt, n_pi0) && aligned(pi0_eta, n_pi0) &&
        aligned(pi0_phi, n_pi0) && aligned(pi0_hepmc_direct_photon_count, n_pi0);
    const bool pi0_decay_photon_aligned =
        aligned(pi0_decay_photon_g4_track_id, n_pi0_decay_photon) &&
        aligned(pi0_decay_photon_parent_g4_track_id, n_pi0_decay_photon) &&
        aligned(pi0_decay_photon_parent_hepmc_barcode, n_pi0_decay_photon) &&
        aligned(pi0_decay_photon_kinematics_valid, n_pi0_decay_photon) &&
        aligned(pi0_decay_photon_e, n_pi0_decay_photon) && aligned(pi0_decay_photon_pt, n_pi0_decay_photon) &&
        aligned(pi0_decay_photon_eta, n_pi0_decay_photon) && aligned(pi0_decay_photon_phi, n_pi0_decay_photon);
    if (!photon_aligned || !pi0_aligned || !pi0_decay_photon_aligned || !event_weights ||
        source_file_id != metadata_source_file_id ||
        event_uid != ((static_cast<ULong64_t>(source_file_id) << 32U) | event_in_file) ||
        event_weight_valid > 1U || !std::isfinite(event_weight) ||
        (event_weight_valid && !std::all_of(event_weights->begin(), event_weights->end(),
            [](const double value) { return std::isfinite(value); })))
    {
      std::cerr << "check_pythia_truth_spectrum_tree - malformed event " << entry << std::endl;
      return 5;
    }

    std::set<int> photon_barcodes;
    for (std::size_t particle = 0; particle < n_photon; ++particle)
    {
      if ((*photon_status)[particle] != 1 ||
          (*photon_kinematics_valid)[particle] > 1U ||
          (*photon_classification_valid)[particle] > 1U ||
          (*photon_copy_chain_valid)[particle] > 1U ||
          ((*photon_kinematics_valid)[particle] &&
           !finite_kinematics((*photon_e)[particle], (*photon_pt)[particle],
                              (*photon_eta)[particle], (*photon_phi)[particle])) ||
          ((*photon_classification_valid)[particle] &&
           ((*photon_category)[particle] < 0 || (*photon_category)[particle] > 3)) ||
          (*photon_immediate_parent_count)[particle] < 0 ||
          (*photon_origin_parent_count)[particle] < 0 ||
          !photon_barcodes.insert((*photon_barcode)[particle]).second)
      {
        std::cerr << "check_pythia_truth_spectrum_tree - invalid photon at event/index "
                  << entry << "/" << particle << std::endl;
        return 6;
      }
    }

    std::set<int> pi0_barcodes;
    for (std::size_t particle = 0; particle < n_pi0; ++particle)
    {
      if ((*pi0_kinematics_valid)[particle] > 1U ||
          ((*pi0_kinematics_valid)[particle] &&
           !finite_kinematics((*pi0_e)[particle], (*pi0_pt)[particle],
                              (*pi0_eta)[particle], (*pi0_phi)[particle])) ||
          !pi0_barcodes.insert((*pi0_barcode)[particle]).second)
      {
        std::cerr << "check_pythia_truth_spectrum_tree - invalid pi0 at event/index "
                  << entry << "/" << particle << std::endl;
        return 7;
      }
    }
    std::set<int> g4_decay_photon_tracks;
    for (std::size_t particle = 0; particle < n_pi0_decay_photon; ++particle)
    {
      if ((*pi0_decay_photon_kinematics_valid)[particle] > 1U ||
          ((*pi0_decay_photon_kinematics_valid)[particle] &&
           !finite_kinematics((*pi0_decay_photon_e)[particle], (*pi0_decay_photon_pt)[particle],
                              (*pi0_decay_photon_eta)[particle], (*pi0_decay_photon_phi)[particle])) ||
          !g4_decay_photon_tracks.insert((*pi0_decay_photon_g4_track_id)[particle]).second ||
          pi0_barcodes.count((*pi0_decay_photon_parent_hepmc_barcode)[particle]) != 1U)
      {
        std::cerr << "check_pythia_truth_spectrum_tree - invalid pi0 decay photon at event/index "
                  << entry << "/" << particle << std::endl;
        return 8;
      }
    }
    observed_n_particle += hepmc_n_particle_record;
    observed_n_photon += truth_photon_n;
    observed_n_pi0 += truth_pi0_n;
    observed_n_pi0_decay_photon += truth_pi0_decay_photon_n;
  }
  if (observed_n_particle != metadata_n_particle ||
      observed_n_photon != metadata_n_photon || observed_n_pi0 != metadata_n_pi0 ||
      observed_n_pi0_decay_photon != metadata_n_pi0_decay_photon)
  {
    std::cerr << "check_pythia_truth_spectrum_tree - metadata aggregate mismatch" << std::endl;
    return 8;
  }

  std::cout << "check_pythia_truth_spectrum_tree - events/particles/final photons/terminal pi0/"
               "G4 pi0 decay photons = "
            << tree->GetEntries() << "/" << observed_n_particle << "/"
            << observed_n_photon << "/" << observed_n_pi0 << "/"
            << observed_n_pi0_decay_photon << std::endl;
  return 0;
}
