#include "ReduceTaggingPurityDiagnostic.C"

namespace
{
constexpr std::array<const char*, 4> kWorkingPoints = {
    "both_0p5", "both_0p2", "pi0_0p2_eta_0p5", "pi0_0p5_eta_0p2"};
constexpr std::array<const char*, 18> kHybridFlowKeys = {
    "all_before", "all_pi0_only_veto", "all_eta_only_veto", "all_both_veto", "all_after",
    "prompt_before", "prompt_pi0_only_veto", "prompt_eta_only_veto", "prompt_both_veto", "prompt_after",
    "background_before", "background_pi0_only_veto", "background_eta_only_veto", "background_both_veto", "background_after",
    "truth_pi0_before", "truth_pi0_pi0_veto", "truth_pi0_after"};

struct HybridFlow
{
  std::array<std::unique_ptr<TH1D>, kHybridFlowKeys.size()> count;
  std::array<std::unique_ptr<TH1D>, kHybridFlowKeys.size()> pb;

  explicit HybridFlow(const std::string& working_point)
  {
    for (std::size_t i = 0; i < kHybridFlowKeys.size(); ++i)
    {
      const std::string base = "h_" + working_point + "_" + kHybridFlowKeys[i] + "_et";
      count[i] = std::make_unique<TH1D>((base + "_count").c_str(), "", kEtBins, 0, kEtMax);
      pb[i] = std::make_unique<TH1D>((base + "_pb").c_str(), "", kEtBins, 0, kEtMax);
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
    for (std::size_t i = 0; i < kHybridFlowKeys.size(); ++i)
    {
      count[i]->Write();
      pb[i]->Write();
    }
  }
};

struct OverlapHistograms
{
  static constexpr std::size_t kPairs = 6;
  std::array<std::unique_ptr<TH1D>, kPairs> all;
  std::array<std::unique_ptr<TH1D>, kPairs> prompt;
  std::array<std::unique_ptr<TH1D>, kPairs> background;

  OverlapHistograms()
  {
    std::size_t pair = 0;
    for (std::size_t first = 0; first < kWorkingPoints.size(); ++first)
    {
      for (std::size_t second = first + 1; second < kWorkingPoints.size(); ++second, ++pair)
      {
        const std::string suffix = std::string(kWorkingPoints[first]) + "__" + kWorkingPoints[second] + "_et_pb";
        all[pair] = std::make_unique<TH1D>(("h_all_after_overlap_" + suffix).c_str(), "", kEtBins, 0, kEtMax);
        prompt[pair] = std::make_unique<TH1D>(("h_prompt_after_overlap_" + suffix).c_str(), "", kEtBins, 0, kEtMax);
        background[pair] = std::make_unique<TH1D>(("h_background_after_overlap_" + suffix).c_str(), "", kEtBins, 0, kEtMax);
        all[pair]->Sumw2();
        prompt[pair]->Sumw2();
        background[pair]->Sumw2();
      }
    }
  }

  void fill(const std::array<bool, kWorkingPoints.size()>& survived, bool is_prompt, double et, double weight)
  {
    std::size_t pair = 0;
    for (std::size_t first = 0; first < kWorkingPoints.size(); ++first)
    {
      for (std::size_t second = first + 1; second < kWorkingPoints.size(); ++second, ++pair)
      {
        if (!survived[first] || !survived[second]) continue;
        all[pair]->Fill(et, weight);
        (is_prompt ? prompt[pair] : background[pair])->Fill(et, weight);
      }
    }
  }

  void write()
  {
    for (std::size_t pair = 0; pair < kPairs; ++pair)
    {
      all[pair]->Write();
      prompt[pair]->Write();
      background[pair]->Write();
    }
  }
};
}

int ReduceHybridThresholdDiagnostic(const std::string low_map_root, const std::string high_map_root, const std::string output_file,
    const std::string sample_name, const int shard_index, const int shard_count, const bool require_complete = true, const Long64_t max_events = 0)
{
  const auto sample = std::find_if(kSamples.begin(), kSamples.end(), [&](const auto& value) { return sample_name == value.name; });
  if (sample == kSamples.end() || low_map_root.empty() || high_map_root.empty() || output_file.empty() || shard_index < 0 || shard_count <= 0 ||
      shard_index >= shard_count || max_events < 0) return 1;
  std::string low_root = low_map_root, high_root = high_map_root;
  while (low_root.size() > 1 && low_root.back() == '/') low_root.pop_back();
  while (high_root.size() > 1 && high_root.back() == '/') high_root.pop_back();
  Input low_input, high_input;
  if (!inspect(*sample, low_root + "/" + sample_name + "/map_*.root", 0.2, require_complete, low_input) ||
      !inspect(*sample, high_root + "/" + sample_name + "/map_*.root", 0.5, require_complete, high_input)) return 2;
  if (low_input.map_count != high_input.map_count || !close(low_input.sumw, high_input.sumw) || low_input.release != high_input.release ||
      low_input.model != high_input.model) return 2;

  TChain low("event_tree"), high("event_tree");
  if (low.Add(low_input.pattern.c_str()) <= 0 || high.Add(high_input.pattern.c_str()) <= 0 || low.GetEntries() != high.GetEntries()) return 3;
  low.SetBranchStatus("*", false);
  high.SetBranchStatus("*", false);
  low.SetCacheSize(kCacheSize);
  high.SetCacheSize(kCacheSize);

  ULong64_t low_uid = 0, high_uid = 0;
  unsigned char weight_valid = 0, stitch_valid = 0, stitch_pass = 0;
  double weight_numerator = 0;
  unsigned int low_ncluster = 0, high_ncluster = 0;
  std::vector<unsigned int>* low_id = nullptr;
  std::vector<unsigned int>* high_id = nullptr;
  std::vector<double>* et = nullptr;
  std::vector<unsigned char>* low_region = nullptr;
  std::vector<unsigned char>* high_region = nullptr;
  std::vector<unsigned char>* low_pi0 = nullptr;
  std::vector<unsigned char>* high_pi0 = nullptr;
  std::vector<unsigned char>* low_eta = nullptr;
  std::vector<unsigned char>* high_eta = nullptr;
  std::vector<unsigned char>* truth_prompt = nullptr;
  std::vector<float>* dominant_fraction = nullptr;
  std::vector<unsigned char>* anchor_valid = nullptr;
  std::vector<float>* anchor_fraction = nullptr;

  const bool low_ok = bind_active(low, "event_uid", &low_uid) && bind_active(low, "event_weight_valid", &weight_valid) &&
      bind_active(low, "sample_stitching_valid", &stitch_valid) && bind_active(low, "sample_stitching_pass", &stitch_pass) &&
      bind_active(low, "weight_numerator_pb", &weight_numerator) && bind_active(low, "split_ncluster", &low_ncluster) &&
      bind_active(low, "split_cluster_id", &low_id) && bind_active(low, "split_cluster_et", &et) &&
      bind_active(low, "split_cluster_pass_region_a", &low_region) && bind_active(low, "split_cluster_pi0_tag", &low_pi0) &&
      bind_active(low, "split_cluster_eta_tag", &low_eta) && bind_active(low, "split_cluster_truth_prompt_cluster", &truth_prompt) &&
      bind_active(low, "split_cluster_truth_dominant_fraction", &dominant_fraction) && bind_active(low, "split_cluster_pi0_anchor_valid", &anchor_valid) &&
      bind_active(low, "split_cluster_pi0_anchor_main_fraction", &anchor_fraction);
  const bool high_ok = bind_active(high, "event_uid", &high_uid) && bind_active(high, "split_ncluster", &high_ncluster) &&
      bind_active(high, "split_cluster_id", &high_id) && bind_active(high, "split_cluster_pass_region_a", &high_region) &&
      bind_active(high, "split_cluster_pi0_tag", &high_pi0) && bind_active(high, "split_cluster_eta_tag", &high_eta);
  if (!low_ok || !high_ok) return 3;
  low.StopCacheLearningPhase();
  high.StopCacheLearningPhase();

  std::array<std::unique_ptr<HybridFlow>, kWorkingPoints.size()> flow;
  for (std::size_t wp = 0; wp < flow.size(); ++wp) flow[wp] = std::make_unique<HybridFlow>(kWorkingPoints[wp]);
  OverlapHistograms overlap;
  const Long64_t total_entries = low.GetEntries();
  const Long64_t entry_begin = total_entries * shard_index / shard_count;
  Long64_t entry_end = total_entries * (shard_index + 1) / shard_count;
  if (max_events > 0) entry_end = std::min(entry_end, entry_begin + max_events);
  unsigned long long stitched_events = 0, region_a_anchors = 0, uid_mismatches = 0, missing_anchors = 0, region_mismatches = 0;

  for (Long64_t entry = entry_begin; entry < entry_end; ++entry)
  {
    if (low.GetEntry(entry) <= 0 || high.GetEntry(entry) <= 0) return 4;
    if (low_uid != high_uid)
    {
      ++uid_mismatches;
      continue;
    }
    const std::array<std::size_t, 9> low_sizes = {low_id ? low_id->size() : 0, et ? et->size() : 0, low_region ? low_region->size() : 0,
        low_pi0 ? low_pi0->size() : 0, low_eta ? low_eta->size() : 0, truth_prompt ? truth_prompt->size() : 0,
        dominant_fraction ? dominant_fraction->size() : 0, anchor_valid ? anchor_valid->size() : 0, anchor_fraction ? anchor_fraction->size() : 0};
    const std::array<std::size_t, 4> high_sizes = {high_id ? high_id->size() : 0, high_region ? high_region->size() : 0,
        high_pi0 ? high_pi0->size() : 0, high_eta ? high_eta->size() : 0};
    if (std::any_of(low_sizes.begin(), low_sizes.end(), [&](std::size_t size) { return size != low_ncluster; }) ||
        std::any_of(high_sizes.begin(), high_sizes.end(), [&](std::size_t size) { return size != high_ncluster; })) return 4;
    if (!weight_valid || !stitch_valid || !stitch_pass) continue;
    ++stitched_events;
    const double weight = weight_numerator / low_input.sumw;
    if (!std::isfinite(weight)) return 4;

    std::unordered_map<unsigned int, std::size_t> high_by_id;
    high_by_id.reserve(high_ncluster);
    for (std::size_t cluster = 0; cluster < high_ncluster; ++cluster) high_by_id.emplace((*high_id)[cluster], cluster);
    for (std::size_t cluster = 0; cluster < low_ncluster; ++cluster)
    {
      if (!(*low_region)[cluster]) continue;
      ++region_a_anchors;
      const auto found = high_by_id.find((*low_id)[cluster]);
      if (found == high_by_id.end())
      {
        ++missing_anchors;
        continue;
      }
      const std::size_t high_cluster = found->second;
      if (!(*high_region)[high_cluster])
      {
        ++region_mismatches;
        continue;
      }
      const bool is_prompt = (*truth_prompt)[cluster] && (*dominant_fraction)[cluster] > 0.5F;
      const bool truth_pi0 = (*anchor_valid)[cluster] && (*anchor_fraction)[cluster] > 0.5F;
      const std::array<bool, kWorkingPoints.size()> ptag = {
          (*high_pi0)[high_cluster] != 0, (*low_pi0)[cluster] != 0, (*low_pi0)[cluster] != 0, (*high_pi0)[high_cluster] != 0};
      const std::array<bool, kWorkingPoints.size()> etag = {
          (*high_eta)[high_cluster] != 0, (*low_eta)[cluster] != 0, (*high_eta)[high_cluster] != 0, (*low_eta)[cluster] != 0};
      std::array<bool, kWorkingPoints.size()> survived{};
      for (std::size_t wp = 0; wp < flow.size(); ++wp)
      {
        survived[wp] = !ptag[wp] && !etag[wp];
        flow[wp]->fill(0, (*et)[cluster], weight);
        if (ptag[wp] && !etag[wp]) flow[wp]->fill(1, (*et)[cluster], weight);
        if (!ptag[wp] && etag[wp]) flow[wp]->fill(2, (*et)[cluster], weight);
        if (ptag[wp] && etag[wp]) flow[wp]->fill(3, (*et)[cluster], weight);
        if (survived[wp]) flow[wp]->fill(4, (*et)[cluster], weight);
        const std::size_t base = is_prompt ? 5 : 10;
        flow[wp]->fill(base, (*et)[cluster], weight);
        if (ptag[wp] && !etag[wp]) flow[wp]->fill(base + 1, (*et)[cluster], weight);
        if (!ptag[wp] && etag[wp]) flow[wp]->fill(base + 2, (*et)[cluster], weight);
        if (ptag[wp] && etag[wp]) flow[wp]->fill(base + 3, (*et)[cluster], weight);
        if (survived[wp]) flow[wp]->fill(base + 4, (*et)[cluster], weight);
        if (truth_pi0)
        {
          flow[wp]->fill(15, (*et)[cluster], weight);
          if (ptag[wp]) flow[wp]->fill(16, (*et)[cluster], weight);
          if (survived[wp]) flow[wp]->fill(17, (*et)[cluster], weight);
        }
      }
      overlap.fill(survived, is_prompt, (*et)[cluster], weight);
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
  TDirectory* flow_root = output.mkdir("working_points");
  for (std::size_t wp = 0; wp < flow.size(); ++wp)
  {
    flow_root->mkdir(kWorkingPoints[wp])->cd();
    flow[wp]->write();
  }
  output.mkdir("overlap")->cd();
  overlap.write();
  output.cd();

  int schema_version = 1, source_schema = kMapSchema, topology_version = kTopologyVersion;
  int metadata_shard_index = shard_index, metadata_shard_count = shard_count;
  bool complete = require_complete;
  Long64_t metadata_total_entries = total_entries, metadata_entry_begin = entry_begin, metadata_entry_end = entry_end, metadata_max_events = max_events;
  double low_production_threshold = 0.2, high_production_threshold = 0.5;
  std::string metadata_sample = sample_name;
  TTree metadata("metadata", "Hybrid tagging-threshold diagnostic partial metadata");
  metadata.Branch("schema_version", &schema_version);
  metadata.Branch("source_map_schema_version", &source_schema);
  metadata.Branch("low_map_root", &low_root);
  metadata.Branch("high_map_root", &high_root);
  metadata.Branch("sample_name", &metadata_sample);
  metadata.Branch("require_complete", &complete);
  metadata.Branch("low_production_threshold", &low_production_threshold);
  metadata.Branch("high_production_threshold", &high_production_threshold);
  metadata.Branch("pi0_topology_algorithm_version", &topology_version);
  metadata.Branch("analysis_release", &low_input.release);
  metadata.Branch("model_sha256", &low_input.model);
  metadata.Branch("map_count", &low_input.map_count);
  metadata.Branch("sample_sum_generator_weights", &low_input.sumw);
  metadata.Branch("shard_index", &metadata_shard_index);
  metadata.Branch("shard_count", &metadata_shard_count);
  metadata.Branch("total_entries", &metadata_total_entries);
  metadata.Branch("entry_begin", &metadata_entry_begin);
  metadata.Branch("entry_end", &metadata_entry_end);
  metadata.Branch("max_events", &metadata_max_events);
  metadata.Branch("stitched_events", &stitched_events);
  metadata.Branch("region_a_anchors", &region_a_anchors);
  metadata.Branch("uid_mismatches", &uid_mismatches);
  metadata.Branch("missing_anchors", &missing_anchors);
  metadata.Branch("region_mismatches", &region_mismatches);
  metadata.Fill();
  metadata.Write();
  output.Close();

  std::cout << "ReduceHybridThresholdDiagnostic - sample/shard/events/anchors/uid_missing_anchor_region_mismatch/output = " << sample_name << "/" << shard_index << "/"
            << (entry_end - entry_begin) << "/" << region_a_anchors << "/" << uid_mismatches << "_" << missing_anchors << "_" << region_mismatches << "/"
            << output_file << std::endl;
  return uid_mismatches == 0 && missing_anchors == 0 && region_mismatches == 0 ? 0 : 6;
}
