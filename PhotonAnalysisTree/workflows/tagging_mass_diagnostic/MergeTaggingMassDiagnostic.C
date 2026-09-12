#include "TaggingMassDiagnostic.h"

int MergeTaggingMassDiagnostic(const std::string family, const std::string partial_root, const std::string output_base,
    const std::string sample_filter = "", const int shard_count_override = 0, const bool require_complete_partials = true, const bool make_plots = true)
{
  using namespace tagging_mass;
  auto definitions = samples(family);
  if (definitions.empty() || partial_root.empty() || output_base.empty() || shard_count_override < 0 || (shard_count_override > 0 && sample_filter.empty())) return 1;
  if (!sample_filter.empty()) definitions.erase(std::remove_if(definitions.begin(), definitions.end(), [&](const auto& s) { return s.name != sample_filter; }), definitions.end());
  if (definitions.empty()) return 1;
  const std::string output_path = output_base + "/tagging_mass_diagnostic.root";
  if (!gSystem->AccessPathName(output_path.c_str())) { std::cerr << "Output already exists: " << output_path << std::endl; return 2; }
  std::array<std::unique_ptr<Histograms>, 6> groups;
  for (std::size_t i = 0; i < groups.size(); ++i) groups[i] = std::make_unique<Histograms>(kSelectionKeys[i]);
  Flow flow;
  Metadata reference;
  std::vector<Metadata> sources;
  for (const auto& definition : definitions)
  {
    const int shards = shard_count_override > 0 ? shard_count_override : std::string(definition.name) == "jet12" ? 10 : 1;
    Metadata sample_reference;
    for (int shard = 0; shard < shards; ++shard)
    {
      const std::string path = partial_root + "/" + definition.name + "/shard_" + std::to_string(shard) + "/tagging_mass_diagnostic.root";
      TFile file(path.c_str(), "READ");
      Metadata value;
      if (file.IsZombie() || !value.read(file) || value.family != family || value.sample != definition.name || value.shard_index != shard || value.shard_count != shards ||
          value.max_events != 0 || (require_complete_partials && !value.require_complete) || value.total_entries <= 0 ||
          value.entry_begin != value.total_entries * shard / shards || value.entry_end != value.total_entries * (shard + 1) / shards || value.entry_end <= value.entry_begin)
      { std::cerr << "Invalid partial/coverage: " << path << std::endl; return 3; }
      if (sources.empty()) reference = value;
      if (!reference.compatible(value) || reference.require_complete != value.require_complete) return 3;
      if (shard == 0) sample_reference = value;
      if (value.total_entries != sample_reference.total_entries || value.map_count != sample_reference.map_count || !close(value.sumw, sample_reference.sumw)) return 3;
      for (std::size_t i = 0; i < groups.size(); ++i) if (!groups[i]->add(file.GetDirectory(kSelectionKeys[i]))) return 4;
      if (!flow.add(file)) return 4;
      sources.push_back(value);
    }
  }
  if (gSystem->mkdir(output_base.c_str(), true) != 0 && gSystem->AccessPathName(output_base.c_str())) return 5;
  const std::string temporary = output_path + ".tmp." + std::to_string(gSystem->GetPid());
  TFile output(temporary.c_str(), "CREATE");
  if (output.IsZombie()) return 5;
  for (std::size_t i = 0; i < groups.size(); ++i) { output.mkdir(kSelectionKeys[i])->cd(); groups[i]->write(); }
  output.cd(); flow.write();
  // Preserve the complete metadata of every contributing shard; sumw is per sample, never summed across shards.
  for (std::size_t i = 0; i < sources.size(); ++i)
  {
    output.mkdir(("sources/partial_" + std::to_string(i)).c_str(), "", true)->cd(); sources[i].write();
  }
  output.cd();
  TTree metadata("metadata", "Merged tagging mass diagnostic");
  int schema_version = 1, source_schema = 5, topology_version = 11;
  auto settings = reference.selection_settings, region_settings = reference.region_settings;
  std::string metadata_family = family, metadata_filter = sample_filter, release = reference.release, model = reference.model;
  unsigned long long partial_count = sources.size();
  bool complete_maps = reference.require_complete;
  metadata.Branch("schema_version", &schema_version); metadata.Branch("source_map_schema_version", &source_schema);
  metadata.Branch("pi0_topology_algorithm_version", &topology_version); metadata.Branch("selection_settings", &settings); metadata.Branch("region_settings", &region_settings);
  metadata.Branch("family", &metadata_family); metadata.Branch("sample_filter", &metadata_filter); metadata.Branch("partial_count", &partial_count);
  metadata.Branch("require_complete", &complete_maps); metadata.Branch("release", &release); metadata.Branch("model", &model);
  metadata.Fill(); metadata.Write(); output.Close();
  if (output.TestBit(TFile::kWriteError) || !gSystem->AccessPathName(output_path.c_str()) || gSystem->Rename(temporary.c_str(), output_path.c_str()) != 0) return 5;
  if (make_plots) plot_groups(groups, reference, output_base);
  std::cout << "Merged " << sources.size() << " partials: " << output_path << std::endl;
  return 0;
}
