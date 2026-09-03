#include "ReducePythiaPhotonCandidateSelection.C"

namespace candidate_composition_merge
{
struct Metadata
{
  int schema_version = -1;
  int source_map_schema_version = -1;
  std::string family;
  std::string map_root;
  bool require_complete = false;
  int n_bins = -1;
  double et_max = -1.0;
  double min_cluster_energy = -1.0;
  double partner_diagnostic_min_cluster_energy = -1.0;
  double meson_partner_min_energy = -1.0;
  int pi0_topology_algorithm_version = -1;
  int signal_embedding_id = -1;
  double majority_threshold = -1.0;
  std::string majority_comparison;
  std::string eta_definition;
  std::string other_definition;
  std::string weight_definition;
  std::string analysis_release;
  std::string model_sha256;
  std::vector<std::string> sample_names;
  std::vector<unsigned long long> sample_map_counts;
  std::vector<double> sample_sum_generator_weights;
  double sample_cross_section_pb = -1.0;
  unsigned long long events_written = 0;
  unsigned long long events_stitch_pass = 0;
  unsigned long long region_a_clusters = 0;
  unsigned long long region_a_prompt_clusters = 0;
  unsigned long long region_a_anchor_clusters = 0;
};

struct CompositionSummary
{
  std::array<candidate_composition::Summary, kSelectionCount> selections;
};

struct TopologySummary
{
  std::array<unsigned long long, kSelectionCount> selected_clusters = {};
  std::array<unsigned long long, kSelectionCount> selected_anchor_clusters = {};
};

struct ShardMetadata
{
  std::string sample_name;
  int shard_index = -1;
  int shard_count = -1;
  unsigned long long map_index_begin = 0;
  unsigned long long map_index_end = 0;
  unsigned long long total_map_count = 0;
};

bool read_metadata(TFile& file, Metadata& value)
{
  auto* tree = file.Get<TTree>("metadata");
  if (!tree || tree->GetEntries() != 1) return false;
  std::string* family = nullptr;
  std::string* map_root = nullptr;
  std::string* majority_comparison = nullptr;
  std::string* eta_definition = nullptr;
  std::string* other_definition = nullptr;
  std::string* weight_definition = nullptr;
  std::string* analysis_release = nullptr;
  std::string* model_sha256 = nullptr;
  std::vector<std::string>* sample_names = nullptr;
  std::vector<unsigned long long>* sample_map_counts = nullptr;
  std::vector<double>* sample_sum_generator_weights = nullptr;
  bool ok = true;
  ok &= bind(tree, "schema_version", &value.schema_version);
  ok &= bind(tree, "source_map_schema_version", &value.source_map_schema_version);
  ok &= bind(tree, "family", &family);
  ok &= bind(tree, "map_root", &map_root);
  ok &= bind(tree, "require_complete", &value.require_complete);
  ok &= bind(tree, "n_bins", &value.n_bins);
  ok &= bind(tree, "et_max", &value.et_max);
  ok &= bind(tree, "min_cluster_energy", &value.min_cluster_energy);
  ok &= bind(tree, "partner_diagnostic_min_cluster_energy", &value.partner_diagnostic_min_cluster_energy);
  ok &= bind(tree, "meson_partner_min_energy", &value.meson_partner_min_energy);
  ok &= bind(tree, "pi0_topology_algorithm_version", &value.pi0_topology_algorithm_version);
  ok &= bind(tree, "signal_embedding_id", &value.signal_embedding_id);
  ok &= bind(tree, "majority_threshold", &value.majority_threshold);
  ok &= bind(tree, "majority_comparison", &majority_comparison);
  ok &= bind(tree, "eta_definition", &eta_definition);
  ok &= bind(tree, "other_definition", &other_definition);
  ok &= bind(tree, "weight_definition", &weight_definition);
  ok &= bind(tree, "analysis_release", &analysis_release);
  ok &= bind(tree, "model_sha256", &model_sha256);
  ok &= bind(tree, "sample_names", &sample_names);
  ok &= bind(tree, "sample_map_counts", &sample_map_counts);
  ok &= bind(tree, "sample_sum_generator_weights", &sample_sum_generator_weights);
  ok &= bind(tree, "sample_cross_section_pb", &value.sample_cross_section_pb);
  ok &= bind(tree, "events_written", &value.events_written);
  ok &= bind(tree, "events_stitch_pass", &value.events_stitch_pass);
  ok &= bind(tree, "region_a_clusters", &value.region_a_clusters);
  ok &= bind(tree, "region_a_prompt_clusters", &value.region_a_prompt_clusters);
  ok &= bind(tree, "region_a_anchor_clusters", &value.region_a_anchor_clusters);
  if (!ok || tree->GetEntry(0) <= 0 || !family || !map_root || !majority_comparison || !eta_definition ||
      !other_definition || !weight_definition || !analysis_release || !model_sha256 || !sample_names || !sample_map_counts ||
      !sample_sum_generator_weights) return false;
  value.family = *family;
  value.map_root = *map_root;
  value.majority_comparison = *majority_comparison;
  value.eta_definition = *eta_definition;
  value.other_definition = *other_definition;
  value.weight_definition = *weight_definition;
  value.analysis_release = *analysis_release;
  value.model_sha256 = *model_sha256;
  value.sample_names = *sample_names;
  value.sample_map_counts = *sample_map_counts;
  value.sample_sum_generator_weights = *sample_sum_generator_weights;
  return true;
}

bool read_composition_summary(TFile& file, CompositionSummary& value)
{
  auto* tree = file.Get<TTree>("composition_summary");
  if (!tree || tree->GetEntries() != static_cast<Long64_t>(kSelectionCount)) return false;
  unsigned int selection_index = 0;
  std::string* selection_key = nullptr;
  candidate_composition::Summary summary;
  const bool ok = bind(tree, "selection_index", &selection_index) && bind(tree, "selection_key", &selection_key) &&
      bind(tree, "selected_clusters", &summary.selected) && bind(tree, "prompt_clusters", &summary.prompt) &&
      bind(tree, "pi0_clusters", &summary.pi0) && bind(tree, "eta_clusters", &summary.eta) && bind(tree, "other_clusters", &summary.other) &&
      bind(tree, "overlap_clusters", &summary.overlap) && bind(tree, "half_boundary_clusters", &summary.half_boundary) &&
      bind(tree, "invalid_truth_clusters", &summary.invalid_truth);
  if (!ok) return false;
  for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry)
  {
    if (tree->GetEntry(entry) <= 0 || !selection_key || selection_index != static_cast<unsigned int>(entry) || selection_index >= kSelectionCount ||
        *selection_key != kSelectionKeys[selection_index] || summary.selected != summary.prompt + summary.pi0 + summary.eta + summary.other) return false;
    value.selections[selection_index] = summary;
  }
  return true;
}

bool read_topology_summary(TFile& file, TopologySummary& value)
{
  auto* tree = file.Get<TTree>("topology_summary");
  if (!tree || tree->GetEntries() != static_cast<Long64_t>(kSelectionCount)) return false;
  unsigned int selection_index = 0;
  std::string* selection_key = nullptr;
  unsigned long long selected = 0, anchor = 0;
  if (!bind(tree, "selection_index", &selection_index) || !bind(tree, "selection_key", &selection_key) ||
      !bind(tree, "selected_clusters", &selected) || !bind(tree, "selected_anchor_clusters", &anchor)) return false;
  for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry)
  {
    if (tree->GetEntry(entry) <= 0 || !selection_key || selection_index != static_cast<unsigned int>(entry) ||
        selection_index >= kSelectionCount || *selection_key != kSelectionKeys[selection_index] || anchor > selected) return false;
    value.selected_clusters[selection_index] = selected;
    value.selected_anchor_clusters[selection_index] = anchor;
  }
  return true;
}

bool read_shard_metadata(TFile& file, ShardMetadata& value)
{
  auto* tree = file.Get<TTree>("shard_metadata");
  if (!tree || tree->GetEntries() != 1) return false;
  std::string* sample_name = nullptr;
  const bool ok = bind(tree, "sample_name", &sample_name) && bind(tree, "shard_index", &value.shard_index) &&
      bind(tree, "shard_count", &value.shard_count) && bind(tree, "map_index_begin", &value.map_index_begin) &&
      bind(tree, "map_index_end", &value.map_index_end) && bind(tree, "total_map_count", &value.total_map_count);
  if (!ok || tree->GetEntry(0) <= 0 || !sample_name) return false;
  value.sample_name = *sample_name;
  return true;
}

bool compatible(const Metadata& value, const Metadata& reference)
{
  return value.schema_version == reference.schema_version && value.source_map_schema_version == reference.source_map_schema_version &&
      value.family == reference.family && value.map_root == reference.map_root &&
      value.require_complete == reference.require_complete && value.n_bins == reference.n_bins && same_double(value.et_max, reference.et_max) &&
      same_double(value.min_cluster_energy, reference.min_cluster_energy) &&
      same_double(value.partner_diagnostic_min_cluster_energy, reference.partner_diagnostic_min_cluster_energy) &&
      same_double(value.meson_partner_min_energy, reference.meson_partner_min_energy) &&
      value.pi0_topology_algorithm_version == reference.pi0_topology_algorithm_version &&
      value.signal_embedding_id == reference.signal_embedding_id && same_double(value.majority_threshold, reference.majority_threshold) &&
      value.majority_comparison == reference.majority_comparison && value.eta_definition == reference.eta_definition &&
      value.other_definition == reference.other_definition && value.weight_definition == reference.weight_definition &&
      value.analysis_release == reference.analysis_release && value.model_sha256 == reference.model_sha256;
}

bool valid_histogram(const TH1D* histogram, const Metadata& metadata)
{
  return histogram && histogram->GetNbinsX() == metadata.n_bins && same_double(histogram->GetXaxis()->GetXmin(), 0.0) &&
      same_double(histogram->GetXaxis()->GetXmax(), metadata.et_max);
}

void add_summary(candidate_composition::Summary& total, const candidate_composition::Summary& value)
{
  total.selected += value.selected;
  total.prompt += value.prompt;
  total.pi0 += value.pi0;
  total.eta += value.eta;
  total.other += value.other;
  total.overlap += value.overlap;
  total.half_boundary += value.half_boundary;
  total.invalid_truth += value.invalid_truth;
}

void add_metadata_counters(Metadata& total, const Metadata& value)
{
  total.events_written += value.events_written;
  total.events_stitch_pass += value.events_stitch_pass;
  total.region_a_clusters += value.region_a_clusters;
  total.region_a_prompt_clusters += value.region_a_prompt_clusters;
  total.region_a_anchor_clusters += value.region_a_anchor_clusters;
}
}

int MergePythiaPhotonCandidateSelection(
    const std::string family,
    const std::string partial_root,
    const std::string composition_output_base,
    const std::string topology_output_base)
{
  using namespace candidate_composition;
  using namespace candidate_composition_merge;
  const std::vector<SampleDefinition> samples = sample_definitions(family);
  if (samples.empty() || partial_root.empty() || composition_output_base.empty() || topology_output_base.empty()) return 1;
  std::string normalized_partial_root = partial_root;
  while (normalized_partial_root.size() > 1U && normalized_partial_root.back() == '/') normalized_partial_root.pop_back();
  if (gSystem->AccessPathName(composition_output_base.c_str()) && gSystem->mkdir(composition_output_base.c_str(), true) != 0) return 2;

  std::array<std::unique_ptr<Histograms>, kSelectionCount> composition_histograms;
  std::array<Summary, kSelectionCount> composition_summaries;
  std::array<std::unique_ptr<Spectra>, kSelectionCount> topology_histograms;
  std::vector<unsigned long long> sample_events_written, sample_events_stitch_pass, sample_region_a_clusters;
  std::vector<unsigned long long> sample_region_a_prompt_clusters, sample_region_a_anchor_clusters;
  std::array<std::vector<unsigned long long>, kSelectionCount> sample_selected_clusters, sample_selected_anchor_clusters;
  std::vector<double> sample_cross_sections_pb;
  Metadata combined;
  bool have_reference = false;
  for (const SampleDefinition& sample : samples)
  {
    const int shard_count = required_shard_count(sample.name);
    unsigned long long next_map_index = 0, total_map_count = 0;
    double sample_sumw = 0.0;
    unsigned long long sample_event_count = 0, sample_stitch_count = 0, sample_region_a_count = 0;
    unsigned long long sample_region_a_prompt_count = 0, sample_region_a_anchor_count = 0;
    std::array<unsigned long long, kSelectionCount> sample_selected = {}, sample_selected_anchor = {};
    for (int shard_index = 0; shard_index < shard_count; ++shard_index)
    {
      const std::string path = normalized_partial_root + "/" + sample.name + "/shard_" + std::to_string(shard_index) + "/photon_candidate_selection.root";
      TFile file(path.c_str(), "READ");
      Metadata metadata;
      ShardMetadata shard;
      CompositionSummary composition_summary;
      TopologySummary topology_summary;
      if (file.IsZombie() || !read_metadata(file, metadata) || !read_shard_metadata(file, shard) || !read_composition_summary(file, composition_summary) || !read_topology_summary(file, topology_summary))
      {
        std::cerr << "Missing or invalid composition partial: " << path << std::endl;
        return 3;
      }
      if (metadata.schema_version != 4 || metadata.source_map_schema_version != 4 || metadata.family != family ||
          metadata.sample_names.size() != 1U || metadata.sample_names.front() != sample.name ||
          metadata.sample_map_counts.size() != 1U || metadata.sample_sum_generator_weights.size() != 1U || !same_double(metadata.sample_cross_section_pb, sample.cross_section_pb) ||
          metadata.events_stitch_pass > metadata.events_written || metadata.region_a_prompt_clusters > metadata.region_a_clusters ||
          metadata.region_a_anchor_clusters > metadata.region_a_clusters ||
          shard.sample_name != sample.name || shard.shard_index != shard_index || shard.shard_count != shard_count ||
          shard.total_map_count != metadata.sample_map_counts.front() || shard.total_map_count == 0 ||
          shard.map_index_begin != shard.total_map_count * static_cast<unsigned long long>(shard_index) / static_cast<unsigned long long>(shard_count) ||
          shard.map_index_end != shard.total_map_count * static_cast<unsigned long long>(shard_index + 1) / static_cast<unsigned long long>(shard_count) ||
          shard.map_index_begin != next_map_index || shard.map_index_begin >= shard.map_index_end)
      {
        std::cerr << "Unexpected partial metadata: " << path << std::endl;
        return 3;
      }
      if (!have_reference)
      {
        combined = metadata;
        combined.sample_names.clear();
        combined.sample_map_counts.clear();
        combined.sample_sum_generator_weights.clear();
        combined.events_written = 0;
        combined.events_stitch_pass = 0;
        combined.region_a_clusters = 0;
        combined.region_a_prompt_clusters = 0;
        combined.region_a_anchor_clusters = 0;
        have_reference = true;
      }
      else if (!compatible(metadata, combined))
      {
        std::cerr << "Incompatible composition partial metadata: " << path << std::endl;
        return 4;
      }
      if (shard_index == 0)
      {
        total_map_count = shard.total_map_count;
        sample_sumw = metadata.sample_sum_generator_weights.front();
        if (!std::isfinite(sample_sumw) || sample_sumw <= 0.0)
        {
          std::cerr << "Invalid full-sample normalization metadata: " << path << std::endl;
          return 4;
        }
      }
      else if (shard.total_map_count != total_map_count || !same_double(metadata.sample_sum_generator_weights.front(), sample_sumw))
      {
        std::cerr << "Inconsistent shard normalization metadata: " << path << std::endl;
        return 4;
      }

      for (std::size_t composition_selection = 0; composition_selection < kSelectionCount; ++composition_selection)
      {
        if (!composition_histograms[composition_selection])
        {
          composition_histograms[composition_selection] = std::make_unique<Histograms>(metadata.n_bins, metadata.et_max);
          for (std::size_t index = 0; index < category_count; ++index)
          {
            composition_histograms[composition_selection]->counts[index]->Reset("ICES");
            composition_histograms[composition_selection]->weighted[index]->Reset("ICES");
          }
        }
        for (std::size_t index = 0; index < category_count; ++index)
        {
          const std::string prefix = std::string("composition/") + kSelectionKeys[composition_selection] + "/h_candidate_" + candidate_composition::kKeys[index] + "_et_";
          TH1D* count = file.Get<TH1D>((prefix + "count").c_str());
          TH1D* weighted_pb = file.Get<TH1D>((prefix + "pb").c_str());
          if (!valid_histogram(count, metadata) || !valid_histogram(weighted_pb, metadata) ||
              !composition_histograms[composition_selection]->counts[index]->Add(count) ||
              !composition_histograms[composition_selection]->weighted[index]->Add(weighted_pb))
          {
            std::cerr << "Missing or incompatible candidate-composition histogram in " << path << std::endl;
            return 5;
          }
        }
        add_summary(composition_summaries[composition_selection], composition_summary.selections[composition_selection]);
      }
      for (std::size_t topology_selection = 0; topology_selection < kSelectionCount; ++topology_selection)
      {
        if (!topology_histograms[topology_selection])
        {
          topology_histograms[topology_selection] = std::make_unique<Spectra>(metadata.n_bins, metadata.et_max);
          for (std::size_t index = 0; index < kSpectrumCount; ++index)
          {
            topology_histograms[topology_selection]->counts[index]->Reset("ICES");
            topology_histograms[topology_selection]->weighted_pb[index]->Reset("ICES");
            topology_histograms[topology_selection]->counts[index]->SetDirectory(nullptr);
            topology_histograms[topology_selection]->weighted_pb[index]->SetDirectory(nullptr);
          }
        }
        for (std::size_t index = 0; index < kSpectrumCount; ++index)
        {
          const std::string prefix = std::string("anchor_topology/") + kSelectionKeys[topology_selection] + "/h_region_a_" + ::kKeys[index] + "_et_";
          TH1D* count = file.Get<TH1D>((prefix + "count").c_str());
          TH1D* weighted_pb = file.Get<TH1D>((prefix + "pb").c_str());
          if (!valid_histogram(count, metadata) || !valid_histogram(weighted_pb, metadata) ||
              !topology_histograms[topology_selection]->counts[index]->Add(count) ||
              !topology_histograms[topology_selection]->weighted_pb[index]->Add(weighted_pb))
          {
            std::cerr << "Missing or incompatible anchor-topology histogram in " << path << std::endl;
            return 5;
          }
        }
        sample_selected[topology_selection] += topology_summary.selected_clusters[topology_selection];
        sample_selected_anchor[topology_selection] += topology_summary.selected_anchor_clusters[topology_selection];
      }
      sample_event_count += metadata.events_written;
      sample_stitch_count += metadata.events_stitch_pass;
      sample_region_a_count += metadata.region_a_clusters;
      sample_region_a_prompt_count += metadata.region_a_prompt_clusters;
      sample_region_a_anchor_count += metadata.region_a_anchor_clusters;
      add_metadata_counters(combined, metadata);
      next_map_index = shard.map_index_end;
    }
    if (next_map_index != total_map_count) return 6;
    combined.sample_names.push_back(sample.name);
    combined.sample_map_counts.push_back(total_map_count);
    combined.sample_sum_generator_weights.push_back(sample_sumw);
    sample_cross_sections_pb.push_back(sample.cross_section_pb);
    sample_events_written.push_back(sample_event_count);
    sample_events_stitch_pass.push_back(sample_stitch_count);
    sample_region_a_clusters.push_back(sample_region_a_count);
    sample_region_a_prompt_clusters.push_back(sample_region_a_prompt_count);
    sample_region_a_anchor_clusters.push_back(sample_region_a_anchor_count);
    for (std::size_t topology_selection = 0; topology_selection < kSelectionCount; ++topology_selection)
    {
      sample_selected_clusters[topology_selection].push_back(sample_selected[topology_selection]);
      sample_selected_anchor_clusters[topology_selection].push_back(sample_selected_anchor[topology_selection]);
    }
  }

  if (!have_reference) return 6;
  for (std::size_t selection_index = 0; selection_index < kSelectionCount; ++selection_index)
  {
    const Summary& summary = composition_summaries[selection_index];
    if (!composition_histograms[selection_index] || summary.selected != summary.prompt + summary.pi0 + summary.eta + summary.other ||
        !candidate_composition::valid_partition(composition_histograms[selection_index]->counts) ||
        !candidate_composition::valid_partition(composition_histograms[selection_index]->weighted) ||
        !topology_histograms[selection_index] || !::valid_partition(topology_histograms[selection_index]->counts) ||
        !::valid_partition(topology_histograms[selection_index]->weighted_pb)) return 6;
  }

  TFile composition_output((composition_output_base + "/selection_comparison.root").c_str(), "RECREATE");
  if (composition_output.IsZombie()) return 7;
  for (std::size_t composition_selection = 0; composition_selection < kSelectionCount; ++composition_selection)
  {
    Histograms& histograms = *composition_histograms[composition_selection];
    std::array<std::unique_ptr<TH1D>, category_count> fractions;
    for (std::size_t index = 1; index < category_count; ++index)
    {
      const std::string name = index == prompt ? "h_photon_candidate_purity" : std::string("h_candidate_") + candidate_composition::kKeys[index] + "_fraction";
      fractions[index] = fraction_histogram(*histograms.weighted[index], *histograms.weighted[denominator], name);
    }
    const std::string selection_output_base = composition_output_base + "/" + kSelectionKeys[composition_selection] + "/photon_candidate_composition";
    if (!make_output_directory(selection_output_base)) return 7;
    draw_stack(fractions, selection_output_base + "_category_fraction_stack.pdf", family, kSelectionKeys[composition_selection], combined.min_cluster_energy, false);
    draw_stack(fractions, selection_output_base + "_category_fraction_stack_detailed.pdf", family, kSelectionKeys[composition_selection], combined.min_cluster_energy, true);
    TDirectory* directory = composition_output.mkdir(kSelectionKeys[composition_selection]);
    if (!directory) return 7;
    directory->cd();
    for (std::size_t index = 0; index < category_count; ++index)
    {
      histograms.counts[index]->Write();
      histograms.weighted[index]->Write();
      if (index > 0) fractions[index]->Write();
    }
    composition_output.cd();
  }
  std::vector<std::string> composition_selection_keys(kSelectionKeys.begin(), kSelectionKeys.end());
  std::vector<std::string> composition_selection_labels(kSelectionLabels.begin(), kSelectionLabels.end());
  std::vector<std::string> composition_selection_definitions(kSelectionDefinitions.begin(), kSelectionDefinitions.end());
  TTree metadata("metadata", "Photon-candidate composition merge metadata");
  metadata.Branch("schema_version", &combined.schema_version);
  metadata.Branch("source_map_schema_version", &combined.source_map_schema_version);
  metadata.Branch("family", &combined.family);
  metadata.Branch("map_root", &combined.map_root);
  metadata.Branch("require_complete", &combined.require_complete);
  metadata.Branch("n_bins", &combined.n_bins);
  metadata.Branch("et_max", &combined.et_max);
  metadata.Branch("selection_keys", &composition_selection_keys);
  metadata.Branch("selection_labels", &composition_selection_labels);
  metadata.Branch("selection_definitions", &composition_selection_definitions);
  metadata.Branch("min_cluster_energy", &combined.min_cluster_energy);
  metadata.Branch("partner_diagnostic_min_cluster_energy", &combined.partner_diagnostic_min_cluster_energy);
  metadata.Branch("meson_partner_min_energy", &combined.meson_partner_min_energy);
  metadata.Branch("pi0_topology_algorithm_version", &combined.pi0_topology_algorithm_version);
  metadata.Branch("signal_embedding_id", &combined.signal_embedding_id);
  metadata.Branch("majority_threshold", &combined.majority_threshold);
  metadata.Branch("majority_comparison", &combined.majority_comparison);
  metadata.Branch("eta_definition", &combined.eta_definition);
  metadata.Branch("other_definition", &combined.other_definition);
  metadata.Branch("weight_definition", &combined.weight_definition);
  metadata.Branch("analysis_release", &combined.analysis_release);
  metadata.Branch("model_sha256", &combined.model_sha256);
  metadata.Branch("sample_names", &combined.sample_names);
  metadata.Branch("sample_map_counts", &combined.sample_map_counts);
  metadata.Branch("sample_sum_generator_weights", &combined.sample_sum_generator_weights);
  metadata.Fill();
  metadata.Write();
  TTree composition_summary_output("selection_summary", "Per-selection candidate-composition merge summary");
  unsigned int composition_selection_index = 0;
  std::string composition_selection_key, composition_selection_label;
  unsigned long long composition_selected = 0, composition_prompt = 0, composition_pi0 = 0, composition_eta = 0, composition_other = 0;
  unsigned long long composition_overlap = 0, composition_half_boundary = 0, composition_invalid_truth = 0;
  composition_summary_output.Branch("selection_index", &composition_selection_index);
  composition_summary_output.Branch("selection_key", &composition_selection_key);
  composition_summary_output.Branch("selection_label", &composition_selection_label);
  composition_summary_output.Branch("selected_clusters", &composition_selected);
  composition_summary_output.Branch("prompt_clusters", &composition_prompt);
  composition_summary_output.Branch("pi0_clusters", &composition_pi0);
  composition_summary_output.Branch("eta_clusters", &composition_eta);
  composition_summary_output.Branch("other_clusters", &composition_other);
  composition_summary_output.Branch("overlap_clusters", &composition_overlap);
  composition_summary_output.Branch("half_boundary_clusters", &composition_half_boundary);
  composition_summary_output.Branch("invalid_truth_clusters", &composition_invalid_truth);
  for (composition_selection_index = 0; composition_selection_index < kSelectionCount; ++composition_selection_index)
  {
    composition_selection_key = kSelectionKeys[composition_selection_index];
    composition_selection_label = kSelectionLabels[composition_selection_index];
    const Summary& summary = composition_summaries[composition_selection_index];
    composition_selected = summary.selected;
    composition_prompt = summary.prompt;
    composition_pi0 = summary.pi0;
    composition_eta = summary.eta;
    composition_other = summary.other;
    composition_overlap = summary.overlap;
    composition_half_boundary = summary.half_boundary;
    composition_invalid_truth = summary.invalid_truth;
    composition_summary_output.Fill();
  }
  composition_summary_output.Write();
  composition_output.Close();

  if (gSystem->AccessPathName(topology_output_base.c_str()) && gSystem->mkdir(topology_output_base.c_str(), true) != 0) return 7;
  TFile topology_output((topology_output_base + "/selection_comparison.root").c_str(), "RECREATE");
  if (topology_output.IsZombie()) return 7;
  const std::string family_label = family == "jet" ? "Pythia8 p+p Jet samples" : "Pythia8 p+p PhotonJet samples";
  for (std::size_t topology_selection = 0; topology_selection < kSelectionCount; ++topology_selection)
  {
    Spectra& spectra = *topology_histograms[topology_selection];
    const std::string selection_output_base = topology_output_base + "/" + kSelectionKeys[topology_selection] + "/region_a_pi0_anchor_topology";
    if (!make_output_directory(selection_output_base)) return 7;
    std::array<std::unique_ptr<TH1D>, kSpectrumCount> density;
    for (std::size_t index = 0; index < kSpectrumCount; ++index)
    {
      density[index].reset(static_cast<TH1D*>(spectra.weighted_pb[index]->Clone((std::string("h_") + kSelectionKeys[topology_selection] + "_" + ::kKeys[index] + "_et_pb_per_gev").c_str())));
      density[index]->SetDirectory(nullptr);
      density[index]->Scale(1.0, "width");
      density[index]->SetLineColor(::kColors[index]);
      density[index]->SetLineWidth(index < 2 ? 3 : 2);
      style_axes(density[index].get(), "d#sigma/dE_{T} [pb/GeV]");
      spectra.counts[index]->SetLineColor(::kColors[index]);
      spectra.weighted_pb[index]->SetLineColor(::kColors[index]);
    }
    const std::vector<std::size_t> detailed_spectrum_indices = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
    const std::vector<std::size_t> summary_spectrum_indices(kSummarySpectrum.begin(), kSummarySpectrum.end());
    const std::vector<std::size_t> detailed_category_indices(kDetailedCategories.begin(), kDetailedCategories.end());
    const std::vector<std::size_t> summary_category_indices(kSummaryCategories.begin(), kSummaryCategories.end());
    auto detailed_fractions = make_fractions(spectra.weighted_pb, detailed_category_indices, std::string("h_") + kSelectionKeys[topology_selection] + "_detailed_");
    auto summary_fractions = make_fractions(spectra.weighted_pb, summary_category_indices, std::string("h_") + kSelectionKeys[topology_selection] + "_summary_");
    const std::string selection_label = family_label + ", " + kSelectionLabels[topology_selection] + ", E_{cluster} > " + std::to_string(combined.min_cluster_energy) + " GeV";
    draw_spectrum(density, detailed_spectrum_indices, selection_output_base + "_detailed.pdf", selection_label, true);
    draw_spectrum(density, summary_spectrum_indices, selection_output_base + ".pdf", selection_label, false);
    draw_fraction_lines(detailed_fractions, detailed_category_indices, selection_output_base + "_category_fractions_detailed.pdf", selection_label, true);
    draw_fraction_lines(summary_fractions, summary_category_indices, selection_output_base + "_category_fractions.pdf", selection_label, false);
    draw_fraction_stack(detailed_fractions, detailed_category_indices, selection_output_base + "_category_fraction_stack_detailed.pdf", selection_label, true);
    draw_fraction_stack(summary_fractions, summary_category_indices, selection_output_base + "_category_fraction_stack.pdf", selection_label, false);

    TDirectory* directory = topology_output.mkdir(kSelectionKeys[topology_selection]);
    if (!directory) return 7;
    directory->cd();
    for (std::size_t index = 0; index < kSpectrumCount; ++index)
    {
      spectra.counts[index]->Write();
      spectra.weighted_pb[index]->Write();
      density[index]->Write();
    }
    for (auto& histogram : detailed_fractions) histogram->Write();
    for (auto& histogram : summary_fractions) histogram->Write();
    topology_output.cd();
  }

  int topology_schema_version = 3;
  int topology_source_schema_version = combined.source_map_schema_version;
  std::vector<std::string> topology_selection_keys(kSelectionKeys.begin(), kSelectionKeys.end());
  std::vector<std::string> topology_selection_labels(kSelectionLabels.begin(), kSelectionLabels.end());
  std::vector<std::string> topology_selection_definitions(kSelectionDefinitions.begin(), kSelectionDefinitions.end());
  std::string topology_prompt_definition = "selected_and_truth_prompt_cluster";
  std::string topology_source = "stored_map_topology_evaluated_with_strict_min_cluster_energy_from_metadata";
  unsigned int topology_sample_count = combined.sample_names.size();
  TTree topology_metadata("metadata", "Pi0-anchor topology merge metadata");
  topology_metadata.Branch("schema_version", &topology_schema_version);
  topology_metadata.Branch("source_map_schema_version", &topology_source_schema_version);
  topology_metadata.Branch("family", &combined.family);
  topology_metadata.Branch("map_root", &combined.map_root);
  topology_metadata.Branch("require_complete", &combined.require_complete);
  topology_metadata.Branch("n_bins", &combined.n_bins);
  topology_metadata.Branch("et_max", &combined.et_max);
  topology_metadata.Branch("selection_keys", &topology_selection_keys);
  topology_metadata.Branch("selection_labels", &topology_selection_labels);
  topology_metadata.Branch("selection_definitions", &topology_selection_definitions);
  topology_metadata.Branch("prompt_definition", &topology_prompt_definition);
  topology_metadata.Branch("topology_source", &topology_source);
  topology_metadata.Branch("analysis_release", &combined.analysis_release);
  topology_metadata.Branch("model_sha256", &combined.model_sha256);
  topology_metadata.Branch("weight_definition", &combined.weight_definition);
  topology_metadata.Branch("sample_count", &topology_sample_count);
  topology_metadata.Branch("sample_names", &combined.sample_names);
  topology_metadata.Branch("min_cluster_energy", &combined.min_cluster_energy);
  topology_metadata.Branch("partner_diagnostic_min_cluster_energy", &combined.partner_diagnostic_min_cluster_energy);
  topology_metadata.Branch("meson_partner_min_energy", &combined.meson_partner_min_energy);
  topology_metadata.Branch("pi0_topology_algorithm_version", &combined.pi0_topology_algorithm_version);
  topology_metadata.Branch("sample_map_counts", &combined.sample_map_counts);
  topology_metadata.Branch("sample_events_written", &sample_events_written);
  topology_metadata.Branch("sample_events_stitch_pass", &sample_events_stitch_pass);
  topology_metadata.Branch("sample_region_a_clusters", &sample_region_a_clusters);
  topology_metadata.Branch("sample_region_a_prompt_clusters", &sample_region_a_prompt_clusters);
  topology_metadata.Branch("sample_region_a_anchor_clusters", &sample_region_a_anchor_clusters);
  topology_metadata.Branch("sample_cross_sections_pb", &sample_cross_sections_pb);
  topology_metadata.Branch("sample_sum_generator_weights", &combined.sample_sum_generator_weights);
  topology_metadata.Fill();
  topology_metadata.Write();

  TTree selection_summary("selection_summary", "Per-selection merge summary");
  unsigned int topology_selection_index = 0;
  std::string topology_selection_key, topology_selection_label;
  std::vector<unsigned long long> selected_clusters, selected_anchor_clusters;
  selection_summary.Branch("selection_index", &topology_selection_index);
  selection_summary.Branch("selection_key", &topology_selection_key);
  selection_summary.Branch("selection_label", &topology_selection_label);
  selection_summary.Branch("sample_selected_clusters", &selected_clusters);
  selection_summary.Branch("sample_selected_anchor_clusters", &selected_anchor_clusters);
  for (topology_selection_index = 0; topology_selection_index < kSelectionCount; ++topology_selection_index)
  {
    topology_selection_key = kSelectionKeys[topology_selection_index];
    topology_selection_label = kSelectionLabels[topology_selection_index];
    selected_clusters = sample_selected_clusters[topology_selection_index];
    selected_anchor_clusters = sample_selected_anchor_clusters[topology_selection_index];
    selection_summary.Fill();
  }
  selection_summary.Write();
  topology_output.Close();

  std::cout << "MergePythiaPhotonCandidateSelection - family/selections/samples/output = "
            << family << "/" << kSelectionCount << "/" << combined.sample_names.size() << "/" << composition_output_base << "/" << topology_output_base << std::endl;
  const bool have_overlap = std::any_of(composition_summaries.begin(), composition_summaries.end(), [](const Summary& summary) { return summary.overlap != 0ULL; });
  return have_overlap ? 9 : 0;
}
