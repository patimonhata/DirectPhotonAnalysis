#include <TChain.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace
{
constexpr int kMapSchema = 4;
constexpr int kTopologyVersion = 8;
constexpr double kDiagnosticFloor = 0.1;
constexpr int kEtBins = 40;
constexpr double kEtMax = 40.0;
constexpr Long64_t kCacheSize = 64LL * 1024LL * 1024LL;

struct Sample
{
  const char* name;
  long long expected_end;
};

const std::array<Sample, 7> kSamples = {{{"jet3", 10000}, {"jet5", 10000}, {"jet8", 10000}, {"jet12", 100000},
                                         {"jet20", 10000}, {"jet30", 10000}, {"jet40", 10000}}};

constexpr std::array<const char*, 41> kFlowKeys = {
    "all_before", "all_pi0_only_veto", "all_eta_only_veto", "all_both_veto", "all_after",
    "prompt_before", "prompt_pi0_only_veto", "prompt_eta_only_veto", "prompt_both_veto", "prompt_after",
    "background_before", "background_pi0_only_veto", "background_eta_only_veto", "background_both_veto", "background_after",
    "truth_pi0_before", "truth_pi0_pi0_veto", "truth_pi0_after",
    "pi0_anchor_before", "pi0_anchor_pi0_veto", "pi0_anchor_truth_taggable_veto", "pi0_anchor_combinatorial_only_veto",
    "pi0_anchor_selected_truth_partner", "pi0_anchor_selected_other_partner", "pi0_anchor_displaced_selected_truth_partner", "pi0_anchor_after",
    "missing_energy_threshold_before", "missing_energy_threshold_after", "missing_acceptance_before", "missing_acceptance_after",
    "missing_other_before", "missing_other_after", "missing_displaced_partner_before", "missing_displaced_partner_after",
    "missing_no_cemc_deposit_before", "missing_no_cemc_deposit_after", "missing_unclustered_deposit_before", "missing_unclustered_deposit_after",
    "missing_match_incomplete_before", "missing_match_incomplete_after", "selected_partner_join_failed"};

constexpr std::array<const char*, 4> kFeatureGroups = {
    "prompt_combinatorial", "truth_pi0_selected_truth_partner", "truth_pi0_selected_other_partner", "truth_pi0_combinatorial_only"};
constexpr std::array<const char*, 14> kFeatureKeys = {
    "partner_e", "partner_et", "partner_eta", "delta_r", "energy_asymmetry", "mass", "event_cluster_multiplicity", "partner_ntower",
    "partner_shower_valid", "partner_full_containment", "partner_tower_data_complete", "partner_e11_over_e33", "partner_e32_over_e35", "partner_bdt_score"};

template <class T>
bool bind_active(TTree& tree, const char* name, T* address)
{
  if (!tree.GetBranch(name)) return false;
  tree.SetBranchStatus(name, true);
  tree.AddBranchToCache(name, true);
  return tree.SetBranchAddress(name, address) >= 0;
}

bool close(double a, double b)
{
  return std::abs(a - b) <= 1e-11 * std::max({1.0, std::abs(a), std::abs(b)});
}

struct Input
{
  std::string pattern;
  unsigned long long map_count = 0;
  double sumw = 0;
  std::string release;
  std::string model;
};

bool inspect(const Sample& sample, const std::string& pattern, double threshold, bool require_complete, Input& input)
{
  TChain tree("metadata");
  if (tree.Add(pattern.c_str()) <= 0) return false;
  int schema = -1, topology_version = -1;
  unsigned int chunk = 0;
  long long begin = -1, end = -1;
  double min_energy = -1, diagnostic_floor = -1, tag_energy = -1, sumw = 0;
  std::string* sample_name = nullptr;
  std::string* release = nullptr;
  std::string* model = nullptr;
  const bool ok = tree.SetBranchAddress("schema_version", &schema) >= 0 && tree.SetBranchAddress("sample_name", &sample_name) >= 0 &&
      tree.SetBranchAddress("analysis_release", &release) >= 0 && tree.SetBranchAddress("model_sha256", &model) >= 0 &&
      tree.SetBranchAddress("manifest_begin", &begin) >= 0 && tree.SetBranchAddress("manifest_end", &end) >= 0 &&
      tree.SetBranchAddress("map_chunk_id", &chunk) >= 0 && tree.SetBranchAddress("sum_generator_weight_processed", &sumw) >= 0 &&
      tree.SetBranchAddress("min_cluster_energy", &min_energy) >= 0 &&
      tree.SetBranchAddress("partner_diagnostic_min_cluster_energy", &diagnostic_floor) >= 0 &&
      tree.SetBranchAddress("meson_partner_min_energy", &tag_energy) >= 0 &&
      tree.SetBranchAddress("pi0_topology_algorithm_version", &topology_version) >= 0;
  if (!ok) return false;
  long long expected_begin = 0;
  for (Long64_t entry = 0; entry < tree.GetEntries(); ++entry)
  {
    if (tree.GetEntry(entry) <= 0 || !sample_name || !release || !model || schema != kMapSchema || *sample_name != sample.name ||
        begin != expected_begin || end <= begin || chunk != static_cast<unsigned int>(entry) || !close(min_energy, threshold) ||
        !close(tag_energy, threshold) || !close(diagnostic_floor, kDiagnosticFloor) || topology_version != kTopologyVersion || !std::isfinite(sumw)) return false;
    if (input.release.empty())
    {
      input.release = *release;
      input.model = *model;
    }
    else if (input.release != *release || input.model != *model) return false;
    input.sumw += sumw;
    expected_begin = end;
  }
  if (expected_begin > sample.expected_end || (require_complete && expected_begin != sample.expected_end) || !(input.sumw > 0)) return false;
  input.pattern = pattern;
  input.map_count = tree.GetEntries();
  return true;
}

struct FlowHistograms
{
  std::array<std::unique_ptr<TH1D>, kFlowKeys.size()> count;
  std::array<std::unique_ptr<TH1D>, kFlowKeys.size()> pb;

  FlowHistograms()
  {
    for (std::size_t i = 0; i < kFlowKeys.size(); ++i)
    {
      count[i] = std::make_unique<TH1D>((std::string("h_") + kFlowKeys[i] + "_et_count").c_str(), "", kEtBins, 0, kEtMax);
      pb[i] = std::make_unique<TH1D>((std::string("h_") + kFlowKeys[i] + "_et_pb").c_str(), "", kEtBins, 0, kEtMax);
      count[i]->Sumw2();
      pb[i]->Sumw2();
    }
  }

  void fill(std::size_t key, double et, double weight)
  {
    count[key]->Fill(et);
    pb[key]->Fill(et, weight);
  }

  void write()
  {
    for (std::size_t i = 0; i < kFlowKeys.size(); ++i)
    {
      count[i]->Write();
      pb[i]->Write();
    }
  }
};

struct FeatureHistograms
{
  std::array<std::unique_ptr<TH1D>, kFeatureKeys.size()> count;
  std::array<std::unique_ptr<TH1D>, kFeatureKeys.size()> pb;

  explicit FeatureHistograms(const std::string& group)
  {
    const std::array<int, kFeatureKeys.size()> bins = {100, 100, 120, 100, 100, 120, 100, 30, 2, 2, 2, 120, 120, 120};
    const std::array<double, kFeatureKeys.size()> low = {0, 0, -1.5, 0, 0, 0, 0, 0, -0.5, -0.5, -0.5, 0, 0, -1.2};
    const std::array<double, kFeatureKeys.size()> high = {2, 2, 1.5, 1, 1, 0.3, 100, 30, 1.5, 1.5, 1.5, 1.2, 1.2, 1.2};
    for (std::size_t i = 0; i < kFeatureKeys.size(); ++i)
    {
      const std::string base = "h_" + group + "_" + kFeatureKeys[i];
      count[i] = std::make_unique<TH1D>((base + "_count").c_str(), "", bins[i], low[i], high[i]);
      pb[i] = std::make_unique<TH1D>((base + "_pb").c_str(), "", bins[i], low[i], high[i]);
      count[i]->Sumw2();
      pb[i]->Sumw2();
    }
  }

  void fill(const std::array<double, kFeatureKeys.size()>& values, double weight)
  {
    for (std::size_t i = 0; i < kFeatureKeys.size(); ++i)
    {
      if (!std::isfinite(values[i])) continue;
      count[i]->Fill(values[i]);
      pb[i]->Fill(values[i], weight);
    }
  }

  void write()
  {
    for (std::size_t i = 0; i < kFeatureKeys.size(); ++i)
    {
      count[i]->Write();
      pb[i]->Write();
    }
  }
};

constexpr std::size_t missing_before_key(int category)
{
  return category == 1 ? 26 : category == 2 ? 28 : category == 3 ? 30 : category == 4 ? 32 : category == 5 ? 34 : category == 6 ? 36 : category == 7 ? 38 : kFlowKeys.size();
}
}

int ReduceTaggingPurityDiagnostic(const std::string map_root, const std::string output_file, const double expected_threshold,
    const std::string sample_name, const int shard_index, const int shard_count, const bool require_complete = true, const Long64_t max_events = 0)
{
  const auto sample = std::find_if(kSamples.begin(), kSamples.end(), [&](const auto& value) { return sample_name == value.name; });
  if (sample == kSamples.end() || map_root.empty() || output_file.empty() || !std::isfinite(expected_threshold) || expected_threshold < 0.2 ||
      shard_index < 0 || shard_count <= 0 || shard_index >= shard_count || max_events < 0) return 1;
  Input input;
  std::string root = map_root;
  while (root.size() > 1 && root.back() == '/') root.pop_back();
  if (!inspect(*sample, root + "/" + sample_name + "/map_*.root", expected_threshold, require_complete, input)) return 2;

  TChain tree("event_tree");
  if (tree.Add(input.pattern.c_str()) <= 0) return 3;
  tree.SetBranchStatus("*", false);
  tree.SetCacheSize(kCacheSize);
  unsigned char weight_valid = 0, stitch_valid = 0, stitch_pass = 0;
  double weight_numerator = 0;
  unsigned int ncluster = 0;
  std::vector<unsigned int>* id = nullptr;
  std::vector<int>* ntower = nullptr;
  std::vector<double>* energy = nullptr;
  std::vector<double>* et = nullptr;
  std::vector<double>* eta = nullptr;
  std::vector<double>* phi = nullptr;
  std::vector<unsigned char>* shower_valid = nullptr;
  std::vector<unsigned char>* full_containment = nullptr;
  std::vector<unsigned char>* tower_complete = nullptr;
  std::vector<float>* e11_over_e33 = nullptr;
  std::vector<float>* e32_over_e35 = nullptr;
  std::vector<float>* bdt_score = nullptr;
  std::vector<unsigned char>* region_a = nullptr;
  std::vector<unsigned char>* pi0_tag = nullptr;
  std::vector<unsigned char>* eta_tag = nullptr;
  std::vector<int>* pi0_partner_id = nullptr;
  std::vector<float>* pi0_partner_mass = nullptr;
  std::vector<unsigned char>* truth_prompt = nullptr;
  std::vector<float>* dominant_fraction = nullptr;
  std::vector<unsigned char>* anchor_valid = nullptr;
  std::vector<float>* anchor_fraction = nullptr;
  std::vector<int>* topology = nullptr;
  std::vector<int>* missing_category = nullptr;
  std::vector<int>* alignment = nullptr;
  std::vector<int>* truth_tag_status = nullptr;
  std::vector<int>* tag_result = nullptr;
  std::vector<unsigned char>* selected_matches = nullptr;
  const bool ok = bind_active(tree, "event_weight_valid", &weight_valid) && bind_active(tree, "sample_stitching_valid", &stitch_valid) &&
      bind_active(tree, "sample_stitching_pass", &stitch_pass) && bind_active(tree, "weight_numerator_pb", &weight_numerator) && bind_active(tree, "split_ncluster", &ncluster) &&
      bind_active(tree, "split_cluster_id", &id) && bind_active(tree, "split_cluster_ntower", &ntower) && bind_active(tree, "split_cluster_e", &energy) &&
      bind_active(tree, "split_cluster_et", &et) && bind_active(tree, "split_cluster_eta", &eta) && bind_active(tree, "split_cluster_phi", &phi) &&
      bind_active(tree, "split_cluster_shower_valid", &shower_valid) && bind_active(tree, "split_cluster_shower_full_containment", &full_containment) &&
      bind_active(tree, "split_cluster_shower_tower_data_complete", &tower_complete) && bind_active(tree, "split_cluster_shower_e11_over_e33", &e11_over_e33) &&
      bind_active(tree, "split_cluster_shower_e32_over_e35", &e32_over_e35) && bind_active(tree, "split_cluster_bdt_score", &bdt_score) &&
      bind_active(tree, "split_cluster_pass_region_a", &region_a) && bind_active(tree, "split_cluster_pi0_tag", &pi0_tag) && bind_active(tree, "split_cluster_eta_tag", &eta_tag) &&
      bind_active(tree, "split_cluster_pi0_partner_cluster_id", &pi0_partner_id) && bind_active(tree, "split_cluster_pi0_partner_mass", &pi0_partner_mass) &&
      bind_active(tree, "split_cluster_truth_prompt_cluster", &truth_prompt) && bind_active(tree, "split_cluster_truth_dominant_fraction", &dominant_fraction) &&
      bind_active(tree, "split_cluster_pi0_anchor_valid", &anchor_valid) && bind_active(tree, "split_cluster_pi0_anchor_main_fraction", &anchor_fraction) &&
      bind_active(tree, "split_cluster_pi0_anchor_topology", &topology) && bind_active(tree, "split_cluster_pi0_anchor_missing_category", &missing_category) &&
      bind_active(tree, "split_cluster_pi0_anchor_partner_alignment", &alignment) &&
      bind_active(tree, "split_cluster_pi0_anchor_truth_partner_tag_status", &truth_tag_status) && bind_active(tree, "split_cluster_pi0_anchor_tag_result", &tag_result) &&
      bind_active(tree, "split_cluster_pi0_anchor_selected_tag_partner_matches_truth_partner", &selected_matches);
  if (!ok) return 3;
  tree.StopCacheLearningPhase();

  FlowHistograms flow;
  std::array<std::unique_ptr<FeatureHistograms>, kFeatureGroups.size()> features;
  for (std::size_t i = 0; i < features.size(); ++i) features[i] = std::make_unique<FeatureHistograms>(kFeatureGroups[i]);
  const Long64_t total_entries = tree.GetEntries();
  const Long64_t entry_begin = total_entries * shard_index / shard_count;
  Long64_t entry_end = total_entries * (shard_index + 1) / shard_count;
  if (max_events > 0) entry_end = std::min(entry_end, entry_begin + max_events);
  unsigned long long stitched_events = 0, selected_anchors = 0, partner_joins_failed = 0;
  for (Long64_t entry = entry_begin; entry < entry_end; ++entry)
  {
    if (tree.GetEntry(entry) <= 0) return 4;
    const std::array<std::size_t, 24> sizes = {id ? id->size() : 0, ntower ? ntower->size() : 0, energy ? energy->size() : 0, et ? et->size() : 0,
        eta ? eta->size() : 0, phi ? phi->size() : 0, shower_valid ? shower_valid->size() : 0, full_containment ? full_containment->size() : 0,
        tower_complete ? tower_complete->size() : 0, e11_over_e33 ? e11_over_e33->size() : 0, e32_over_e35 ? e32_over_e35->size() : 0,
        bdt_score ? bdt_score->size() : 0, region_a ? region_a->size() : 0, pi0_tag ? pi0_tag->size() : 0, eta_tag ? eta_tag->size() : 0,
        pi0_partner_id ? pi0_partner_id->size() : 0, pi0_partner_mass ? pi0_partner_mass->size() : 0, truth_prompt ? truth_prompt->size() : 0,
        dominant_fraction ? dominant_fraction->size() : 0, anchor_valid ? anchor_valid->size() : 0, anchor_fraction ? anchor_fraction->size() : 0,
        topology ? topology->size() : 0, missing_category ? missing_category->size() : 0, alignment ? alignment->size() : 0};
    if (std::any_of(sizes.begin(), sizes.end(), [&](std::size_t size) { return size != ncluster; }) || !truth_tag_status || !tag_result || !selected_matches ||
        truth_tag_status->size() != ncluster || tag_result->size() != ncluster || selected_matches->size() != ncluster) return 4;
    if (!weight_valid || !stitch_valid || !stitch_pass) continue;
    ++stitched_events;
    const double weight = weight_numerator / input.sumw;
    if (!std::isfinite(weight)) return 4;
    std::unordered_map<unsigned int, std::size_t> by_id;
    by_id.reserve(ncluster);
    for (std::size_t cluster = 0; cluster < ncluster; ++cluster) by_id.emplace((*id)[cluster], cluster);
    for (std::size_t cluster = 0; cluster < ncluster; ++cluster)
    {
      if (!(*region_a)[cluster]) continue;
      ++selected_anchors;
      const bool ptag = (*pi0_tag)[cluster] != 0;
      const bool etag = (*eta_tag)[cluster] != 0;
      const bool survived = !ptag && !etag;
      const bool prompt = (*truth_prompt)[cluster] && (*dominant_fraction)[cluster] > 0.5F;
      const bool truth_pi0 = (*anchor_valid)[cluster] && (*anchor_fraction)[cluster] > 0.5F;
      flow.fill(0, (*et)[cluster], weight);
      if (ptag && !etag) flow.fill(1, (*et)[cluster], weight);
      if (!ptag && etag) flow.fill(2, (*et)[cluster], weight);
      if (ptag && etag) flow.fill(3, (*et)[cluster], weight);
      if (survived) flow.fill(4, (*et)[cluster], weight);
      const std::size_t base = prompt ? 5 : 10;
      flow.fill(base, (*et)[cluster], weight);
      if (ptag && !etag) flow.fill(base + 1, (*et)[cluster], weight);
      if (!ptag && etag) flow.fill(base + 2, (*et)[cluster], weight);
      if (ptag && etag) flow.fill(base + 3, (*et)[cluster], weight);
      if (survived) flow.fill(base + 4, (*et)[cluster], weight);
      if (truth_pi0)
      {
        flow.fill(15, (*et)[cluster], weight);
        if (ptag) flow.fill(16, (*et)[cluster], weight);
        if (survived) flow.fill(17, (*et)[cluster], weight);
      }
      if ((*anchor_valid)[cluster])
      {
        flow.fill(18, (*et)[cluster], weight);
        if (ptag) flow.fill(19, (*et)[cluster], weight);
        if (ptag && (*tag_result)[cluster] == 2) flow.fill(20, (*et)[cluster], weight);
        if (ptag && (*tag_result)[cluster] == 3) flow.fill(21, (*et)[cluster], weight);
        if (ptag && (*selected_matches)[cluster]) flow.fill(22, (*et)[cluster], weight);
        if (ptag && (*truth_tag_status)[cluster] == 1 && !(*selected_matches)[cluster]) flow.fill(23, (*et)[cluster], weight);
        if (ptag && (*selected_matches)[cluster] && (*alignment)[cluster] == 2) flow.fill(24, (*et)[cluster], weight);
        if (survived) flow.fill(25, (*et)[cluster], weight);
        if ((*topology)[cluster] == 3)
        {
          const std::size_t before = missing_before_key((*missing_category)[cluster]);
          if (before < kFlowKeys.size())
          {
            flow.fill(before, (*et)[cluster], weight);
            if (survived) flow.fill(before + 1, (*et)[cluster], weight);
          }
        }
      }
      if (!ptag) continue;
      const int selected_partner_id = (*pi0_partner_id)[cluster];
      const auto partner = selected_partner_id < 0 ? by_id.end() : by_id.find(static_cast<unsigned int>(selected_partner_id));
      if (partner == by_id.end())
      {
        ++partner_joins_failed;
        flow.fill(40, (*et)[cluster], weight);
        continue;
      }
      const std::size_t j = partner->second;
      const double dphi = std::remainder((*phi)[j] - (*phi)[cluster], 2.0 * std::acos(-1.0));
      const double sum_e = (*energy)[cluster] + (*energy)[j];
      const double asym = sum_e > 0 ? std::abs((*energy)[cluster] - (*energy)[j]) / sum_e : NAN;
      const std::array<double, kFeatureKeys.size()> values = {(*energy)[j], (*et)[j], (*eta)[j], std::hypot((*eta)[j] - (*eta)[cluster], dphi),
          asym, (*pi0_partner_mass)[cluster], static_cast<double>(ncluster), static_cast<double>((*ntower)[j]), static_cast<double>((*shower_valid)[j]),
          static_cast<double>((*full_containment)[j]), static_cast<double>((*tower_complete)[j]), (*e11_over_e33)[j], (*e32_over_e35)[j], (*bdt_score)[j]};
      if (prompt) features[0]->fill(values, weight);
      if (truth_pi0 && (*truth_tag_status)[cluster] == 1 && (*selected_matches)[cluster]) features[1]->fill(values, weight);
      if (truth_pi0 && (*truth_tag_status)[cluster] == 1 && !(*selected_matches)[cluster]) features[2]->fill(values, weight);
      if (truth_pi0 && (*truth_tag_status)[cluster] != 1) features[3]->fill(values, weight);
    }
  }

  const std::size_t slash = output_file.find_last_of('/');
  if (slash != std::string::npos)
  {
    const std::string directory = output_file.substr(0, slash);
    if (gSystem->mkdir(directory.c_str(), true) != 0 && gSystem->AccessPathName(directory.c_str())) return 5;
  }
  TFile output(output_file.c_str(), "RECREATE");
  if (output.IsZombie()) return 5;
  output.mkdir("flow")->cd();
  flow.write();
  TDirectory* feature_root = output.mkdir("selected_partner_features");
  for (std::size_t i = 0; i < features.size(); ++i)
  {
    feature_root->mkdir(kFeatureGroups[i])->cd();
    features[i]->write();
  }
  output.cd();
  int schema_version = 1, source_schema = kMapSchema, topology_version = kTopologyVersion;
  int metadata_shard_index = shard_index, metadata_shard_count = shard_count;
  bool complete = require_complete;
  Long64_t metadata_total_entries = total_entries, metadata_entry_begin = entry_begin, metadata_entry_end = entry_end;
  Long64_t metadata_max_events = max_events;
  double production_threshold = expected_threshold, diagnostic_floor = kDiagnosticFloor;
  std::string metadata_sample = sample_name, metadata_map_root = root;
  TTree metadata("metadata", "Tagging-purity diagnostic partial metadata");
  metadata.Branch("schema_version", &schema_version);
  metadata.Branch("source_map_schema_version", &source_schema);
  metadata.Branch("map_root", &metadata_map_root);
  metadata.Branch("sample_name", &metadata_sample);
  metadata.Branch("require_complete", &complete);
  metadata.Branch("production_threshold", &production_threshold);
  metadata.Branch("diagnostic_floor", &diagnostic_floor);
  metadata.Branch("pi0_topology_algorithm_version", &topology_version);
  metadata.Branch("analysis_release", &input.release);
  metadata.Branch("model_sha256", &input.model);
  metadata.Branch("map_count", &input.map_count);
  metadata.Branch("sample_sum_generator_weights", &input.sumw);
  metadata.Branch("shard_index", &metadata_shard_index);
  metadata.Branch("shard_count", &metadata_shard_count);
  metadata.Branch("total_entries", &metadata_total_entries);
  metadata.Branch("entry_begin", &metadata_entry_begin);
  metadata.Branch("entry_end", &metadata_entry_end);
  metadata.Branch("max_events", &metadata_max_events);
  metadata.Branch("stitched_events", &stitched_events);
  metadata.Branch("region_a_anchors", &selected_anchors);
  metadata.Branch("selected_partner_joins_failed", &partner_joins_failed);
  metadata.Fill();
  metadata.Write();
  output.Close();
  std::cout << "ReduceTaggingPurityDiagnostic - threshold/sample/shard/events/anchors/join_failures/output = " << expected_threshold << "/" << sample_name << "/"
            << shard_index << "/" << (entry_end - entry_begin) << "/" << selected_anchors << "/" << partner_joins_failed << "/" << output_file << std::endl;
  return partner_joins_failed == 0 ? 0 : 6;
}
