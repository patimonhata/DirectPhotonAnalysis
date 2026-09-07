#include "PlotDisplacedPartnerMassDiagnostic.C"

namespace
{
int required_shards(const std::string& sample)
{
  return sample == "jet12" ? 10 : 1;
}

struct PartialMetadata
{
  std::vector<double> selection_settings;
  int schema = -1;
  int source_schema = -1;
  int shard_index = -1;
  int shard_count = -1;
  Long64_t total_entries = -1;
  Long64_t entry_begin = -1;
  Long64_t entry_end = -1;
  Long64_t max_events = -1;
  bool require_complete = false;
  double topology_threshold = -1;
  double production_tag_threshold = -1;
  double diagnostic_floor = -1;
  double emulated_tag_threshold = -1;
  double mass_min = -1;
  double mass_max = -1;
  std::string family;
  std::string sample;
  std::string release;
  std::string model;
  unsigned long long map_count = 0;
  double sumw = 0;
};

bool read_partial_metadata(TFile& file, PartialMetadata& value)
{
  auto* tree = file.Get<TTree>("metadata");
  if (!tree || tree->GetEntries() != 1) return false;
  std::string* family = nullptr;
  std::string* sample_filter = nullptr;
  std::string* release = nullptr;
  std::string* model = nullptr;
  std::vector<std::string>* sample_names = nullptr;
  std::vector<unsigned long long>* map_counts = nullptr;
  std::vector<double>* sum_weights = nullptr;
  std::vector<double>* settings = nullptr;
  const bool ok = bind_branch(*tree, "selection_settings", &settings) && bind_branch(*tree, "schema_version", &value.schema) && bind_branch(*tree, "source_map_schema_version", &value.source_schema) &&
      bind_branch(*tree, "family", &family) && bind_branch(*tree, "sample_filter", &sample_filter) && bind_branch(*tree, "require_complete", &value.require_complete) &&
      bind_branch(*tree, "max_events_per_sample", &value.max_events) && bind_branch(*tree, "shard_index", &value.shard_index) &&
      bind_branch(*tree, "shard_count", &value.shard_count) && bind_branch(*tree, "total_entries", &value.total_entries) && bind_branch(*tree, "entry_begin", &value.entry_begin) && bind_branch(*tree, "entry_end", &value.entry_end) &&
      bind_branch(*tree, "topology_min_cluster_energy", &value.topology_threshold) &&
      bind_branch(*tree, "production_tagging_partner_min_energy", &value.production_tag_threshold) &&
      bind_branch(*tree, "diagnostic_min_cluster_energy", &value.diagnostic_floor) &&
      bind_branch(*tree, "emulated_tagging_partner_min_energy", &value.emulated_tag_threshold) && bind_branch(*tree, "pi0_mass_min", &value.mass_min) &&
      bind_branch(*tree, "pi0_mass_max", &value.mass_max) && bind_branch(*tree, "analysis_release", &release) && bind_branch(*tree, "model_sha256", &model) &&
      bind_branch(*tree, "sample_names", &sample_names) && bind_branch(*tree, "sample_map_counts", &map_counts) &&
      bind_branch(*tree, "sample_sum_generator_weights", &sum_weights);
  if (!ok || tree->GetEntry(0) <= 0 || !family || !sample_filter || !release || !model || !sample_names || !map_counts || !sum_weights ||
      sample_names->size() != 1 || map_counts->size() != 1 || sum_weights->size() != 1 || *sample_filter != sample_names->front()) return false;
  if (!settings) return false;
  value.selection_settings = *settings;
  value.family = *family;
  value.sample = sample_names->front();
  value.release = *release;
  value.model = *model;
  value.map_count = map_counts->front();
  value.sumw = sum_weights->front();
  return true;
}

bool compatible_partial(const PartialMetadata& value, const PartialMetadata& reference, bool require_complete_partials)
{
  return photon_candidate_settings::same(value.selection_settings, reference.selection_settings) && value.schema == 1 && value.source_schema == kMapSchema && value.family == reference.family && value.release == reference.release &&
      value.model == reference.model && (!require_complete_partials || value.require_complete) && value.max_events == 0 && close(value.topology_threshold, kTopologyThreshold) &&
      close(value.diagnostic_floor, kDiagnosticFloor) &&
      close(value.emulated_tag_threshold, 0.2) && std::isfinite(value.mass_min) && std::isfinite(value.mass_max) && value.mass_min >= 0.0 && value.mass_min < value.mass_max &&
      close(value.topology_threshold, reference.topology_threshold) &&
      close(value.production_tag_threshold, reference.production_tag_threshold) && close(value.diagnostic_floor, reference.diagnostic_floor) &&
      close(value.emulated_tag_threshold, reference.emulated_tag_threshold) && close(value.mass_min, reference.mass_min) && close(value.mass_max, reference.mass_max);
}
}

int MergeDisplacedPartnerMassDiagnostic(const std::string family, const std::string partial_root, const std::string output_base,
    const std::string sample_filter = "", const int shard_count_override = 0, const bool require_complete_partials = true)
{
  auto definitions = samples(family);
  if (definitions.empty() || partial_root.empty() || output_base.empty() || shard_count_override < 0 || (shard_count_override > 0 && sample_filter.empty())) return 1;
  if (!sample_filter.empty())
  {
    definitions.erase(std::remove_if(definitions.begin(), definitions.end(), [&](const auto& sample) { return sample_filter != sample.name; }), definitions.end());
    if (definitions.empty()) return 1;
  }
  if (gSystem->mkdir(output_base.c_str(), true) != 0 && gSystem->AccessPathName(output_base.c_str())) return 2;
  gROOT->cd();
  std::array<std::unique_ptr<Histograms>, kSelectionKeys.size()> histograms;
  for (std::size_t selection = 0; selection < histograms.size(); ++selection) histograms[selection] = std::make_unique<Histograms>(kSelectionKeys[selection]);
  PartialMetadata reference;
  bool have_reference = false;
  std::vector<std::string> sample_names;
  std::vector<unsigned long long> map_counts;
  std::vector<double> sum_weights;
  unsigned long long partial_count = 0;

  for (const auto& definition : definitions)
  {
    const int shard_count = shard_count_override > 0 ? shard_count_override : required_shards(definition.name);
    Long64_t expected_begin = 0, sample_total_entries = -1;
    unsigned long long sample_map_count = 0;
    double sample_sumw = 0;
    for (int shard = 0; shard < shard_count; ++shard)
    {
      const std::string path = partial_root + "/" + definition.name + "/shard_" + std::to_string(shard) + "/displaced_partner_mass_diagnostic.root";
      TFile file(path.c_str(), "READ");
      PartialMetadata metadata;
      if (file.IsZombie() || !read_partial_metadata(file, metadata) || metadata.family != family || metadata.sample != definition.name ||
          metadata.shard_index != shard || metadata.shard_count != shard_count || metadata.total_entries <= 0 || metadata.entry_begin != expected_begin || metadata.entry_end <= metadata.entry_begin)
      {
        std::cerr << "Invalid partial or shard coverage: " << path << std::endl;
        return 3;
      }
      if (!have_reference)
      {
        reference = metadata;
        kPi0MassMin = metadata.mass_min;
        kPi0MassMax = metadata.mass_max;
        have_reference = true;
      }
      if (!compatible_partial(metadata, reference, require_complete_partials)) return 3;
      if (shard == 0)
      {
        sample_total_entries = metadata.total_entries;
        sample_map_count = metadata.map_count;
        sample_sumw = metadata.sumw;
      }
      else if (sample_total_entries != metadata.total_entries || sample_map_count != metadata.map_count || !close(sample_sumw, metadata.sumw)) return 3;
      for (std::size_t selection = 0; selection < histograms.size(); ++selection)
      {
        if (!histograms[selection]->add(file.GetDirectory(kSelectionKeys[selection]))) return 4;
      }
      expected_begin = metadata.entry_end;
      ++partial_count;
    }
    if (expected_begin != sample_total_entries)
    {
      std::cerr << "Incomplete shard coverage for " << definition.name << ": " << expected_begin << "/" << sample_total_entries << std::endl;
      return 3;
    }
    sample_names.push_back(definition.name);
    map_counts.push_back(sample_map_count);
    sum_weights.push_back(sample_sumw);
  }

  TFile output((output_base + "/displaced_partner_mass_diagnostic.root").c_str(), "RECREATE");
  if (output.IsZombie()) return 5;
  for (std::size_t selection = 0; selection < histograms.size(); ++selection)
  {
    output.mkdir(kSelectionKeys[selection])->cd();
    histograms[selection]->write(kSelectionKeys[selection]);
  }
  output.cd();
  int schema_version = 1, source_schema = kMapSchema;
  std::string metadata_family = family, metadata_filter = sample_filter;
  bool metadata_require_complete = require_complete_partials;
  int metadata_shard_count_override = shard_count_override;
  TTree metadata("metadata", "Merged displaced-partner mass diagnostic metadata");
  metadata.Branch("schema_version", &schema_version);
  metadata.Branch("selection_settings", &reference.selection_settings);
  metadata.Branch("source_map_schema_version", &source_schema);
  metadata.Branch("family", &metadata_family);
  metadata.Branch("sample_filter", &metadata_filter);
  metadata.Branch("require_complete_partials", &metadata_require_complete);
  metadata.Branch("shard_count_override", &metadata_shard_count_override);
  metadata.Branch("topology_min_cluster_energy", &reference.topology_threshold);
  metadata.Branch("production_tagging_partner_min_energy", &reference.production_tag_threshold);
  metadata.Branch("diagnostic_min_cluster_energy", &reference.diagnostic_floor);
  metadata.Branch("emulated_tagging_partner_min_energy", &reference.emulated_tag_threshold);
  metadata.Branch("pi0_mass_min", &reference.mass_min);
  metadata.Branch("pi0_mass_max", &reference.mass_max);
  metadata.Branch("analysis_release", &reference.release);
  metadata.Branch("model_sha256", &reference.model);
  metadata.Branch("sample_names", &sample_names);
  metadata.Branch("sample_map_counts", &map_counts);
  metadata.Branch("sample_sum_generator_weights", &sum_weights);
  metadata.Branch("partial_count", &partial_count);
  metadata.Fill();
  metadata.Write();
  TTree summary("summary", "Merged per-selection displaced-pair summary");
  unsigned int index = 0;
  std::string key;
  unsigned long long selected_pairs = 0, window_pairs = 0;
  double selected_pb = 0, window_pb = 0, raw_fraction = 0, weighted_fraction = 0;
  summary.Branch("selection_index", &index);
  summary.Branch("selection_key", &key);
  summary.Branch("selected_pairs", &selected_pairs);
  summary.Branch("pi0_window_pairs", &window_pairs);
  summary.Branch("selected_cross_section_pb", &selected_pb);
  summary.Branch("pi0_window_cross_section_pb", &window_pb);
  summary.Branch("raw_pi0_window_fraction", &raw_fraction);
  summary.Branch("weighted_pi0_window_fraction", &weighted_fraction);
  for (index = 0; index < histograms.size(); ++index)
  {
    key = kSelectionKeys[index];
    auto* count = histograms[index]->mass_count.get();
    auto* weighted = histograms[index]->mass_pb.get();
    const int low = count->FindBin(std::nextafter(kPi0MassMin, kPi0MassMax));
    const int high = count->FindBin(std::nextafter(kPi0MassMax, kPi0MassMin));
    selected_pairs = std::llround(count->Integral(0, count->GetNbinsX() + 1));
    window_pairs = std::llround(count->Integral(low, high));
    selected_pb = weighted->Integral(0, weighted->GetNbinsX() + 1);
    window_pb = weighted->Integral(low, high);
    raw_fraction = selected_pairs != 0 ? static_cast<double>(window_pairs) / selected_pairs : 0;
    weighted_fraction = selected_pb != 0 ? window_pb / selected_pb : 0;
    summary.Fill();
  }
  summary.Write();
  output.Close();

  SetsPhenixStyle();
  gStyle->SetOptStat(0);
  for (std::size_t selection = 0; selection < histograms.size(); ++selection)
  {
    const std::string directory = output_base + "/" + kSelectionKeys[selection];
    gSystem->mkdir(directory.c_str(), true);
    std::ostringstream text;
    text << (family == "jet" ? "Pythia8 p+p Jet, " : "Pythia8 p+p PhotonJet, ") << kSelectionLabels[selection]
         << ", E_{partner} > " << reference.emulated_tag_threshold << " GeV";
    draw_mass(histograms[selection]->mass_pb.get(), directory + "/displaced_partner_mass.pdf", text.str());
    draw_2d(histograms[selection]->mass_vs_anchor_et_pb.get(), directory + "/displaced_partner_mass_vs_anchor_et.pdf", "Anchor cluster E_{T} [GeV]", text.str());
    draw_2d(histograms[selection]->mass_vs_truth_pt_pb.get(), directory + "/displaced_partner_mass_vs_truth_pi0_pt.pdf", "Truth #pi^{0} p_{T} [GeV]", text.str());
    draw_shapes(histograms[selection]->mass_vs_anchor_et_count.get(), kAnchorEtEdges, directory + "/displaced_partner_mass_by_anchor_et_shape.pdf", text.str());
    draw_shapes(histograms[selection]->mass_vs_truth_pt_count.get(), kTruthPtEdges, directory + "/displaced_partner_mass_by_truth_pi0_pt_shape.pdf", text.str());
  }
  std::cout << "MergeDisplacedPartnerMassDiagnostic - partials/output = " << partial_count << "/" << output_base << std::endl;
  return 0;
}
