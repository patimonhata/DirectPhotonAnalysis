#include "ReducePythiaPhotonCandidateComposition.C"

namespace candidate_composition_merge
{
struct Metadata
{
  int schema_version = -1;
  int source_map_schema_version = -1;
  std::string family;
  std::string selection;
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
  unsigned long long selected_cluster_count = 0;
  unsigned long long prompt_cluster_count = 0;
  unsigned long long pi0_cluster_count = 0;
  unsigned long long eta_cluster_count = 0;
  unsigned long long other_cluster_count = 0;
  unsigned long long overlap_cluster_count = 0;
  unsigned long long half_boundary_cluster_count = 0;
  unsigned long long invalid_truth_cluster_count = 0;
};

bool read_metadata(TFile& file, Metadata& value)
{
  auto* tree = file.Get<TTree>("metadata");
  if (!tree || tree->GetEntries() != 1) return false;
  std::string* family = nullptr;
  std::string* selection = nullptr;
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
  ok &= bind(tree, "selection", &selection);
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
  ok &= bind(tree, "selected_cluster_count", &value.selected_cluster_count);
  ok &= bind(tree, "prompt_cluster_count", &value.prompt_cluster_count);
  ok &= bind(tree, "pi0_cluster_count", &value.pi0_cluster_count);
  ok &= bind(tree, "eta_cluster_count", &value.eta_cluster_count);
  ok &= bind(tree, "other_cluster_count", &value.other_cluster_count);
  ok &= bind(tree, "overlap_cluster_count", &value.overlap_cluster_count);
  ok &= bind(tree, "half_boundary_cluster_count", &value.half_boundary_cluster_count);
  ok &= bind(tree, "invalid_truth_cluster_count", &value.invalid_truth_cluster_count);
  if (!ok || tree->GetEntry(0) <= 0 || !family || !selection || !map_root || !majority_comparison || !eta_definition ||
      !other_definition || !weight_definition || !analysis_release || !model_sha256 || !sample_names || !sample_map_counts ||
      !sample_sum_generator_weights) return false;
  value.family = *family;
  value.selection = *selection;
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

bool compatible(const Metadata& value, const Metadata& reference)
{
  return value.schema_version == reference.schema_version && value.source_map_schema_version == reference.source_map_schema_version &&
      value.family == reference.family && value.selection == reference.selection && value.map_root == reference.map_root &&
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

void add_counters(Metadata& total, const Metadata& value)
{
  total.selected_cluster_count += value.selected_cluster_count;
  total.prompt_cluster_count += value.prompt_cluster_count;
  total.pi0_cluster_count += value.pi0_cluster_count;
  total.eta_cluster_count += value.eta_cluster_count;
  total.other_cluster_count += value.other_cluster_count;
  total.overlap_cluster_count += value.overlap_cluster_count;
  total.half_boundary_cluster_count += value.half_boundary_cluster_count;
  total.invalid_truth_cluster_count += value.invalid_truth_cluster_count;
}
}

int MergePythiaPhotonCandidateComposition(
    const std::string family = "jet",
    const std::string partial_root = "",
    const std::string output_base = "",
    const std::string selection = "region_a")
{
  using namespace candidate_composition;
  using namespace candidate_composition_merge;
  const std::vector<SampleDefinition> samples = sample_definitions(family);
  if (samples.empty() || partial_root.empty() || output_base.empty() || !find_selection(selection)) return 1;
  std::string normalized_partial_root = partial_root;
  while (normalized_partial_root.size() > 1U && normalized_partial_root.back() == '/') normalized_partial_root.pop_back();
  if (!make_output_directory(output_base)) return 2;

  std::array<std::unique_ptr<TH1D>, category_count> counts;
  std::array<std::unique_ptr<TH1D>, category_count> weighted;
  Metadata combined;
  bool have_reference = false;
  for (const SampleDefinition& sample : samples)
  {
    const std::string path = normalized_partial_root + "/" + sample.name + "/photon_candidate_composition.root";
    TFile file(path.c_str(), "READ");
    Metadata metadata;
    if (file.IsZombie() || !read_metadata(file, metadata))
    {
      std::cerr << "Missing or invalid composition partial: " << path << std::endl;
      return 3;
    }
    if (metadata.schema_version != 1 || metadata.source_map_schema_version != 4 || metadata.family != family ||
        metadata.selection != selection || metadata.sample_names.size() != 1U || metadata.sample_names.front() != sample.name ||
        metadata.sample_map_counts.size() != 1U || metadata.sample_sum_generator_weights.size() != 1U ||
        metadata.selected_cluster_count != metadata.prompt_cluster_count + metadata.pi0_cluster_count +
            metadata.eta_cluster_count + metadata.other_cluster_count)
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
      combined.selected_cluster_count = 0;
      combined.prompt_cluster_count = 0;
      combined.pi0_cluster_count = 0;
      combined.eta_cluster_count = 0;
      combined.other_cluster_count = 0;
      combined.overlap_cluster_count = 0;
      combined.half_boundary_cluster_count = 0;
      combined.invalid_truth_cluster_count = 0;
      have_reference = true;
    }
    else if (!compatible(metadata, combined))
    {
      std::cerr << "Incompatible composition partial metadata: " << path << std::endl;
      return 4;
    }

    for (std::size_t index = 0; index < category_count; ++index)
    {
      const std::string count_name = std::string("h_candidate_") + candidate_composition::kKeys[index] + "_et_count";
      const std::string weighted_name = std::string("h_candidate_") + candidate_composition::kKeys[index] + "_et_pb";
      TH1D* count = file.Get<TH1D>(count_name.c_str());
      TH1D* weighted_pb = file.Get<TH1D>(weighted_name.c_str());
      if (!valid_histogram(count, metadata) || !valid_histogram(weighted_pb, metadata))
      {
        std::cerr << "Missing or incompatible histogram in " << path << ": " << count_name << " / " << weighted_name << std::endl;
        return 5;
      }
      if (!counts[index])
      {
        counts[index].reset(static_cast<TH1D*>(count->Clone()));
        weighted[index].reset(static_cast<TH1D*>(weighted_pb->Clone()));
        counts[index]->SetDirectory(nullptr);
        weighted[index]->SetDirectory(nullptr);
      }
      else if (!counts[index]->Add(count) || !weighted[index]->Add(weighted_pb))
      {
        std::cerr << "Could not add histograms from " << path << std::endl;
        return 5;
      }
    }
    combined.sample_names.push_back(metadata.sample_names.front());
    combined.sample_map_counts.push_back(metadata.sample_map_counts.front());
    combined.sample_sum_generator_weights.push_back(metadata.sample_sum_generator_weights.front());
    add_counters(combined, metadata);
  }

  if (!have_reference || combined.selected_cluster_count != combined.prompt_cluster_count + combined.pi0_cluster_count +
          combined.eta_cluster_count + combined.other_cluster_count ||
      !valid_partition(counts) || !valid_partition(weighted)) return 6;
  std::array<std::unique_ptr<TH1D>, category_count> fractions;
  for (std::size_t index = 1; index < category_count; ++index)
  {
    const std::string name = index == prompt ? "h_photon_candidate_purity" : std::string("h_candidate_") + candidate_composition::kKeys[index] + "_fraction";
    fractions[index] = fraction_histogram(*weighted[index], *weighted[denominator], name);
  }
  draw_stack(fractions, output_base + "_category_fraction_stack.pdf", family, selection, combined.min_cluster_energy, false);
  draw_stack(fractions, output_base + "_category_fraction_stack_detailed.pdf", family, selection, combined.min_cluster_energy, true);

  TFile output((output_base + ".root").c_str(), "RECREATE");
  if (output.IsZombie()) return 7;
  for (std::size_t index = 0; index < category_count; ++index)
  {
    counts[index]->Write();
    weighted[index]->Write();
    if (index > 0) fractions[index]->Write();
  }
  TTree metadata("metadata", "Photon-candidate composition reduce metadata");
  metadata.Branch("schema_version", &combined.schema_version);
  metadata.Branch("source_map_schema_version", &combined.source_map_schema_version);
  metadata.Branch("family", &combined.family);
  metadata.Branch("selection", &combined.selection);
  metadata.Branch("map_root", &combined.map_root);
  metadata.Branch("require_complete", &combined.require_complete);
  metadata.Branch("n_bins", &combined.n_bins);
  metadata.Branch("et_max", &combined.et_max);
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
  metadata.Branch("selected_cluster_count", &combined.selected_cluster_count);
  metadata.Branch("prompt_cluster_count", &combined.prompt_cluster_count);
  metadata.Branch("pi0_cluster_count", &combined.pi0_cluster_count);
  metadata.Branch("eta_cluster_count", &combined.eta_cluster_count);
  metadata.Branch("other_cluster_count", &combined.other_cluster_count);
  metadata.Branch("overlap_cluster_count", &combined.overlap_cluster_count);
  metadata.Branch("half_boundary_cluster_count", &combined.half_boundary_cluster_count);
  metadata.Branch("invalid_truth_cluster_count", &combined.invalid_truth_cluster_count);
  metadata.Fill();
  metadata.Write();
  output.Close();
  std::cout << "MergePythiaPhotonCandidateComposition - family/selection/samples/selected/output = "
            << family << "/" << selection << "/" << combined.sample_names.size() << "/" << combined.selected_cluster_count << "/" << output_base << std::endl;
  return combined.overlap_cluster_count == 0ULL ? 0 : 9;
}
