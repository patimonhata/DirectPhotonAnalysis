#ifndef PHOTON_BDT_TRAINING_CLUSTERS_H
#define PHOTON_BDT_TRAINING_CLUSTERS_H

#include "../photon_candidate_selection/CandidateSelectionCommon.h"
#include "../../src/Pi0AnchorTopologyEvaluator.h"

#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <TUUID.h>

#include <cstdint>
#include <filesystem>
#include <stdexcept>

namespace bdt_training
{
using photon_tree::Pi0AnchorTopology;
using photon_tree::composition::Main;

struct ShapeFeature { const char* input; const char* output; };
constexpr std::array<ShapeFeature, 8> shapes = {{
    {"w_eta_cogx", "cluster_weta_cogx"}, {"w_phi_cogx", "cluster_wphi_cogx"}, {"e11_over_e33", "e11_over_e33"},
    {"et1", "cluster_et1"}, {"et2", "cluster_et2"}, {"et3", "cluster_et3"}, {"et4", "cluster_et4"}, {"e32_over_e35", "e32_over_e35"}}};

std::array<int, 3> labels(bool photonjet, Main main, int overview, int topology, bool isolated)
{
  if (photonjet) return main == Main::prompt && isolated ? std::array<int, 3>{1, 1, 1} : std::array<int, 3>{-1, -1, -1};
  if (main == Main::decay_photon) return {-1, 1, 1};
  if (overview == 3) return {0, 0, 0};
  if (main != Main::pi0 && main != Main::eta_gg) return {-1, -1, -1};
  switch (static_cast<Pi0AnchorTopology>(topology))
  {
    case Pi0AnchorTopology::missing: return {1, 1, 1};
    case Pi0AnchorTopology::merged:
    case Pi0AnchorTopology::single_contaminated: return {0, 0, 0};
    case Pi0AnchorTopology::separated: return {-1, -1, 1};
    default: return {-1, -1, -1};
  }
}

// Mix the sample name and the full segment/event UID reproducibly across ROOT sessions and shards.
int event_split(std::uint64_t sample_hash, std::uint64_t uid)
{
  std::uint64_t hash = (sample_hash ^ uid) + UINT64_C(0x9e3779b97f4a7c15);
  hash = (hash ^ (hash >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
  hash = (hash ^ (hash >> 27)) * UINT64_C(0x94d049bb133111eb);
  const auto bucket = (hash ^ (hash >> 31)) % 100;
  return bucket < 70 ? 0 : bucket < 85 ? 1 : 2;
}

struct ClusterRow
{
  std::string sample_name;
  ULong64_t event_uid = 0;
  UInt_t cluster_id = 0;
  int split = 0, common_main = 0, common_origin = 0, topology_applicability = 0, topology = -1;
  bool inputs_valid = false, event_eligible = false, truth_iso_valid = false, truth_pass_isolated = false;
  double weight_pb = 0, et = 0, eta = 0, z = 0;
  std::array<float, shapes.size()> shape{};
  std::array<int, 3> label{};

  void branch(TTree& tree)
  {
    tree.Branch("sample_name", &sample_name);
    tree.Branch("event_uid", &event_uid);
    tree.Branch("cluster_id", &cluster_id);
    tree.Branch("split", &split);
    tree.Branch("event_eligible", &event_eligible);
    tree.Branch("weight_pb", &weight_pb);
    tree.Branch("common_main", &common_main);
    tree.Branch("common_origin", &common_origin);
    tree.Branch("topology_applicability", &topology_applicability);
    tree.Branch("topology", &topology);
    tree.Branch("truth_iso_valid", &truth_iso_valid);
    tree.Branch("truth_pass_isolated", &truth_pass_isolated);
    tree.Branch("inputs_valid", &inputs_valid);
    tree.Branch("cluster_Et", &et);
    tree.Branch("cluster_Eta", &eta);
    tree.Branch("vertexz", &z);
    for (std::size_t i = 0; i < shapes.size(); ++i) tree.Branch(shapes[i].output, &shape[i]);
    tree.Branch("label_v1", &label[0]);
    tree.Branch("label_v2", &label[1]);
    tree.Branch("label_v3", &label[2]);
  }
};

void extract_map(const MapMetadata& map, bool photonjet, double sample_sumw, std::uint64_t sample_hash, ClusterRow& row, TTree& output)
{
  TFile file(map.path.c_str(), "READ");
  auto* tree = file.Get<TTree>("event_tree");
  if (file.IsZombie() || !tree || static_cast<ULong64_t>(tree->GetEntries()) != map.events_written) throw std::runtime_error("invalid event_tree");

  TTreeReader reader(tree);
  const auto required = [tree](const std::string& name) {
    if (!tree->GetBranch(name.c_str())) throw std::runtime_error("missing branch " + name);
    return name;
  };
  TTreeReaderValue<UInt_t> n(reader, required("split_ncluster").c_str());
  TTreeReaderValue<ULong64_t> uid(reader, required("event_uid").c_str());
  TTreeReaderValue<double> z(reader, required("mbd_vertex_z").c_str()), numerator(reader, required("weight_numerator_pb").c_str());
  TTreeReaderValue<UChar_t> weight_valid(reader, required("event_weight_valid").c_str()), stitch_valid(reader, required("sample_stitching_valid").c_str());
  TTreeReaderValue<UChar_t> stitch_pass(reader, required("sample_stitching_pass").c_str());
  const auto cluster = [&](const std::string& name) { return required("split_cluster_" + name); };
  TTreeReaderArray<UInt_t> id(reader, cluster("id").c_str());
  TTreeReaderArray<double> et(reader, cluster("et").c_str()), eta(reader, cluster("eta").c_str());
  TTreeReaderArray<int> main(reader, cluster("common_main").c_str()), pdg(reader, cluster("common_selected_pdg").c_str());
  TTreeReaderArray<int> applicability(reader, cluster("common_applicability").c_str());
  TTreeReaderArray<int> pi0_topology(reader, cluster("pi0_anchor_topology").c_str()), eta_topology(reader, cluster("eta_anchor_topology").c_str());
  TTreeReaderArray<UChar_t> shower_valid(reader, cluster("shower_valid").c_str());
  TTreeReaderArray<UChar_t> pi0_valid(reader, cluster("pi0_anchor_valid").c_str()), eta_valid(reader, cluster("eta_anchor_valid").c_str());
  using Flags = TTreeReaderArray<UChar_t>;
  std::unique_ptr<Flags> iso_valid, isolated;
  if (photonjet)
  {
    iso_valid = std::make_unique<Flags>(reader, cluster("truth_iso_valid").c_str());
    isolated = std::make_unique<Flags>(reader, cluster("truth_pass_isolated").c_str());
  }
  std::array<std::unique_ptr<TTreeReaderArray<float>>, shapes.size()> features;
  for (std::size_t i = 0; i < shapes.size(); ++i)
    features[i] = std::make_unique<TTreeReaderArray<float>>(reader, cluster(std::string("shower_") + shapes[i].input).c_str());

  while (reader.Next())
  {
    row.z = *z;
    if (!(std::abs(row.z) < 60.0)) continue;
    row.event_uid = *uid;
    row.split = event_split(sample_hash, row.event_uid);
    row.event_eligible = *weight_valid && *stitch_valid && *stitch_pass;
    row.weight_pb = row.event_eligible ? *numerator / sample_sumw : 0.0;
    if (!std::isfinite(row.weight_pb)) throw std::runtime_error("nonfinite event weight");

    const auto aligned = [&](const auto&... columns) { return ((columns.GetSize() == *n) && ...); };
    if (!aligned(id, et, eta, main, pdg, applicability, pi0_topology, eta_topology, pi0_valid, eta_valid, shower_valid) ||
        (photonjet && !aligned(*iso_valid, *isolated))) throw std::runtime_error("misaligned cluster branches");
    for (const auto& feature : features) if (!aligned(*feature)) throw std::runtime_error("misaligned shower inputs");

    for (UInt_t c = 0; c < *n; ++c)
    {
      row.et = et[c];
      row.eta = eta[c];
      if (!(row.et > 3.0 && row.et < 35.0 && std::abs(row.eta) < 0.7)) continue;
      row.cluster_id = id[c];
      row.common_main = main[c];
      if (main[c] < 0 || main[c] > static_cast<int>(Main::photon_unknown)) throw std::runtime_error("invalid common_main");
      const auto origin = static_cast<Main>(main[c]);
      row.common_origin = photon_tree::composition::overview(origin, pdg[c]);
      row.topology_applicability = applicability[c];
      row.topology = -1;
      if (applicability[c] == static_cast<int>(photon_tree::composition::Applicability::applicable))
      {
        if (origin == Main::pi0 && pi0_valid[c]) row.topology = pi0_topology[c];
        if (origin == Main::eta_gg && eta_valid[c]) row.topology = eta_topology[c];
      }
      row.truth_iso_valid = photonjet && (*iso_valid)[c];
      row.truth_pass_isolated = row.truth_iso_valid && (*isolated)[c];
      row.label = row.event_eligible ? labels(photonjet, origin, row.common_origin, row.topology, row.truth_pass_isolated) : std::array<int, 3>{-1, -1, -1};
      row.inputs_valid = shower_valid[c];
      for (std::size_t i = 0; i < shapes.size(); ++i)
      {
        row.shape[i] = (*features[i])[c];
        row.inputs_valid &= std::isfinite(row.shape[i]);
      }
      if (output.Fill() < 0) throw std::runtime_error("failed to write cluster");
    }
  }
  if (reader.GetEntryStatus() != TTreeReader::kEntryBeyondEnd) throw std::runtime_error("failed reading event_tree");
}

struct SampleInput
{
  std::string map_root, name;
  int shard_index = 0, shard_count = 1;
};

int write_clusters(const std::vector<SampleInput>& inputs, const std::string& output_path)
{
  namespace fs = std::filesystem;
  std::string location = output_path;
  fs::path temporary;
  try
  {
    if (inputs.empty() || output_path.empty()) throw std::runtime_error("missing inputs or output path");
    if (fs::exists(output_path)) throw std::runtime_error("output already exists");
    if (!make_output_directory(output_path)) throw std::runtime_error("cannot create output directory");
    temporary = output_path + ".tmp." + TUUID().AsString();
    TFile output(temporary.c_str(), "CREATE");
    if (output.IsZombie()) throw std::runtime_error("cannot create temporary output");
    ClusterRow row;
    TTree clusters("clusters", "One row per fiducial and kinematic split cluster");
    row.branch(clusters);

    int source_schema = kMbdMapSchema, shard_index = 0, shard_count = 0;
    Long64_t map_begin = 0, map_end = 0;
    double sample_sumw = 0;
    std::string input_root;
    TTree metadata("metadata", "One entry per sample and extraction shard");
    metadata.Branch("source_map_root", &input_root);
    metadata.Branch("source_schema", &source_schema);
    metadata.Branch("sample_name", &row.sample_name);
    metadata.Branch("sample_sum_generator_weights", &sample_sumw);
    metadata.Branch("shard_index", &shard_index);
    metadata.Branch("shard_count", &shard_count);
    metadata.Branch("map_begin", &map_begin);
    metadata.Branch("map_end", &map_end);

    for (const auto& input : inputs)
    {
      row.sample_name = input.name;
      shard_index = input.shard_index;
      shard_count = input.shard_count;
      location = input.map_root + "/" + input.name;
      const bool photonjet = input.name.rfind("photonjet", 0) == 0;
      const auto samples = sample_definitions(photonjet ? "photonjet" : "jet");
      const auto sample = std::find_if(samples.begin(), samples.end(), [&](const auto& s) { return s.name == input.name; });
      if (sample == samples.end() || input.map_root.empty() || shard_index < 0 || shard_count <= shard_index)
        throw std::runtime_error("invalid sample, map root or shard");
      input_root = fs::absolute(input.map_root).string();
      std::vector<MapMetadata> maps;
      if (!collect_maps(location + "/map_*.root", *sample, true, maps)) throw std::runtime_error("incomplete or incompatible input maps");
      if (static_cast<std::size_t>(shard_count) > maps.size()) throw std::runtime_error("more shards than input maps");
      sample_sumw = 0;
      for (const auto& map : maps) sample_sumw += map.sum_generator_weight_processed;
      if (!std::isfinite(sample_sumw) || sample_sumw <= 0) throw std::runtime_error("invalid sample normalization");
      map_begin = maps.size() * static_cast<std::size_t>(shard_index) / shard_count;
      map_end = maps.size() * static_cast<std::size_t>(shard_index + 1) / shard_count;
      std::uint64_t sample_hash = UINT64_C(14695981039346656037);
      for (unsigned char ch : input.name) sample_hash = (sample_hash ^ ch) * UINT64_C(1099511628211);
      std::cout << input.name << ": processing maps [" << map_begin << ", " << map_end << ")" << std::endl;
      for (Long64_t i = map_begin; i < map_end; ++i)
      {
        location = maps[i].path;
        extract_map(maps[i], photonjet, sample_sumw, sample_hash, row, clusters);
      }
      if (metadata.Fill() < 0) throw std::runtime_error("failed writing sample metadata");
    }
    location = output_path;
    output.cd();
    if (output.Write() <= 0 || output.TestBit(TFile::kWriteError)) throw std::runtime_error("failed writing ROOT output");
    const Long64_t rows = clusters.GetEntries();
    clusters.SetDirectory(nullptr);
    metadata.SetDirectory(nullptr);
    output.Close();
    if (output.TestBit(TFile::kWriteError)) throw std::runtime_error("failed closing ROOT output");
    fs::create_hard_link(temporary, output_path);
    fs::remove(temporary);
    std::cout << "Wrote " << rows << " rows to " << output_path << std::endl;
    return 0;
  }
  catch (const std::exception& e)
  {
    std::error_code ignored;
    if (!temporary.empty()) fs::remove(temporary, ignored);
    std::cerr << "BDT training extraction (" << location << "): " << e.what() << std::endl;
    return 1;
  }
}
}

#endif
