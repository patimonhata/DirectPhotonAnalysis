#ifndef PHOTON_CANDIDATE_CANDIDATESELECTIONCOMMON_H
#define PHOTON_CANDIDATE_CANDIDATESELECTIONCOMMON_H

#include "../../src/EtaOriginFractions.h"
#include "../common/plotting/sPhenixStyle.C"

#include "Pi0BdtHistograms.h"
#include <TCanvas.h>
#include <TDirectory.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <THStack.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TObjArray.h>
#include <TPad.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace
{
constexpr std::size_t kSpectrumCount = 13;
constexpr std::size_t kSelectionCount = 6;
constexpr std::array<const char*, kSelectionCount> kSelectionKeys = {
    "kinematic", "preselection", "preselection_tight", "preselection_isolation", "region_a", "region_a_tagging_veto"};
constexpr std::array<const char*, kSelectionCount> kSelectionLabels = {
    "Kinematic", "Pre-selection", "Pre-selection + TightBDT", "Pre-selection + Isolation", "Region A: isolated and tight", "Region A + Tagging veto"};
constexpr std::array<const char*, kSelectionCount> kSelectionDefinitions = {
    "pass_kinematics", "pass_kinematics_and_pass_preselection", "pass_kinematics_and_pass_preselection_and_pass_tight",
    "pass_kinematics_and_pass_preselection_and_pass_isolated", "pass_region_a", "pass_region_a_and_not_pi0_tag_and_not_eta_tag"};
constexpr std::array<const char*, kSpectrumCount> kKeys = {
    "prompt", "pi0_anchor", "separated", "merged", "single_contaminated", "missing",
    "missing_energy_band_taggable", "missing_energy_band_not_taggable", "missing_acceptance",
    "missing_low_energy", "missing_unclustered_or_no_cemc_deposit", "missing_other", "other"};
std::array<std::string, kSpectrumCount> kLabels = {
    "Prompt-#gamma cluster", "#pi^{0}-main anchor", "Separated", "Merged", "Single contaminated", "Missing (total)",
    "Missing: energy band, mass inside", "Missing: energy band, mass outside", "Missing: acceptance",
    "Missing: low energy", "Missing: unclustered or no CEMC deposit", "Missing: other", "Other"};
constexpr std::array<int, kSpectrumCount> kColors = {
    kRed + 1, kBlue + 1, kAzure + 7, kMagenta + 1, kCyan + 2, kGreen + 2,
    kOrange + 7, kPink + 7, kViolet + 1, kYellow + 2, kSpring + 5, kGreen + 3, kGray + 2};
constexpr std::array<std::size_t, 7> kSummarySpectrum = {0, 1, 2, 3, 4, 5, 12};
constexpr std::array<std::size_t, 10> kDetailedCategories = {2, 3, 4, 6, 7, 9, 10, 8, 11, 12};
constexpr std::array<std::size_t, 5> kSummaryCategories = {2, 3, 4, 5, 12};
void set_missing_labels(double low, double high)
{
  std::ostringstream band, lower;
  band << low << " < E #leq " << high << " GeV";
  lower << "Missing: E #leq " << low << " GeV";
  kLabels[6] = "Missing: " + band.str() + ", mass inside";
  kLabels[7] = "Missing: " + band.str() + ", mass outside";
  kLabels[9] = lower.str();
}
constexpr int kCanvasWidth = 1100;
constexpr int kCanvasHeight = 1000;
constexpr double kPlotEtMaximum = 40.0;
constexpr double kPlotPadHeight = 0.68;
constexpr double kCompositionPlotWidth = 0.70;
constexpr double kSpectrumYMinimum = 1e-2;
constexpr double kSpectrumYMaximum = 5e6;
constexpr Long64_t kEventTreeCacheSize = 64LL * 1024LL * 1024LL;

struct SampleDefinition
{
  const char* name;
  double window_min;
  double window_max;
  double cross_section_pb;
  bool upper_unbounded;
  long long expected_manifest_end;
};

std::vector<SampleDefinition> sample_definitions(const std::string& family)
{
  if (family == "jet")
  {
    return {
        {"jet3", 0.0, 5.0, 1.2147e9, false, 10000},
        {"jet5", 5.0, 9.0, 1.3878e8, false, 10000},
        {"jet8", 9.0, 14.0, 1.3013e7, false, 10000},
        {"jet12", 14.0, 21.0, 1.4903e6, false, 100000},
        {"jet20", 21.0, 32.0, 6.2623e4, false, 10000},
        {"jet30", 32.0, 42.0, 2.5298e3, false, 10000},
        {"jet40", 42.0, 0.0, 1.3553e2, true, 10000}};
  }
  if (family == "photonjet")
  {
    return {
        {"photonjet3", 0.0, 5.0, 1.0348e6, false, 10000},
        {"photonjet5", 5.0, 14.0, 1.4636e5, false, 10000},
        {"photonjet10", 14.0, 22.0, 6.9447e3, false, 10000},
        {"photonjet20", 22.0, 200.0, 1.3045e2, false, 10000}};
  }
  return {};
}

constexpr int kMbdMapSchema = 11;

struct MbdMapSettings
{
  double max_abs_mbd_vertex_z = -1.0;
  double shower_shape_min_cluster_et = -1.0;
  double candidate_abs_eta_max = -1.0;
  double topology_anchor_abs_eta_max = -1.0;
  double candidate_et_min = -1.0;
  double candidate_et_max = -1.0;
  std::string reconstruction_vertex_source;
  std::string topology_vertex_source;

  void create_branches(TTree& tree)
  {
    tree.Branch("max_abs_mbd_vertex_z", &max_abs_mbd_vertex_z);
    tree.Branch("shower_shape_min_cluster_et", &shower_shape_min_cluster_et);
    tree.Branch("candidate_abs_eta_max", &candidate_abs_eta_max);
    tree.Branch("topology_anchor_abs_eta_max", &topology_anchor_abs_eta_max);
    tree.Branch("candidate_et_min", &candidate_et_min);
    tree.Branch("candidate_et_max", &candidate_et_max);
    tree.Branch("reconstruction_vertex_source", &reconstruction_vertex_source);
    tree.Branch("topology_vertex_source", &topology_vertex_source);
  }
};

struct MapMetadata
{
  std::string path;
  int common_composition_algorithm_version = -1;
  int schema_version = -1;
  int decay_origin_algorithm_version = -1;
  double decay_origin_max_proper_length_cm = -1, decay_origin_max_flight_cm = -1;
  std::string input_manifest;
  std::string sample_name;
  std::string analysis_release;
  std::string model_sha256;
  long long manifest_begin = -1;
  long long manifest_end = -1;
  long long input_file_count = -1;
  unsigned int map_chunk_id = 0;
  double cross_section_pb = 0.0;
  double window_min = 0.0;
  double window_max = 0.0;
  bool upper_unbounded = false;
  double sum_generator_weight_processed = 0.0;
  unsigned long long events_processed = 0;
  double min_cluster_energy = -1.0;
  double partner_diagnostic_min_cluster_energy = -1.0;
  double pi0_partner_min_energy = -1.0;
  double eta_partner_min_energy = -1.0;
  double pi0_mass_min = -1.0;
  double pi0_mass_max = -1.0;
  double eta_mass_min = -1.0;
  double eta_mass_max = -1.0;
  double missing_energy_min = -1.0;
  double missing_energy_max = -1.0;
  double min_photon_energy_recovery = -1.0;
  int pi0_topology_algorithm_version = -1;
  int eta_topology_algorithm_version = -1;
  unsigned long long events_written = 0;
  unsigned long long events_vertex_rejected = 0;
  unsigned long long events_invalid = 0;
  double truth_isolation_radius = -1.0;
  MbdMapSettings mbd;
};

template <class T>
bool bind(TTree* tree, const char* name, T* address)
{
  return tree->GetBranch(name) && tree->SetBranchAddress(name, address) >= 0;
}

template <class T>
bool bind_active(TTree* tree, const char* name, T* address)
{
  if (!tree->GetBranch(name)) return false;
  tree->SetBranchStatus(name, true);
  tree->AddBranchToCache(name, true);
  return tree->SetBranchAddress(name, address) >= 0;
}

bool same_double(double left, double right)
{
  return std::abs(left - right) <= 1e-11 * std::max({1.0, std::abs(left), std::abs(right)});
}

bool read_mbd_map_settings(TTree* tree, MbdMapSettings& value)
{
  std::string* reconstruction_vertex_source = nullptr;
  std::string* topology_vertex_source = nullptr;
  const bool ok = bind(tree, "max_abs_mbd_vertex_z", &value.max_abs_mbd_vertex_z) &&
      bind(tree, "shower_shape_min_cluster_et", &value.shower_shape_min_cluster_et) &&
      bind(tree, "candidate_abs_eta_max", &value.candidate_abs_eta_max) &&
      bind(tree, "topology_anchor_abs_eta_max", &value.topology_anchor_abs_eta_max) &&
      bind(tree, "candidate_et_min", &value.candidate_et_min) && bind(tree, "candidate_et_max", &value.candidate_et_max) &&
      bind(tree, "reconstruction_vertex_source", &reconstruction_vertex_source) && bind(tree, "topology_vertex_source", &topology_vertex_source);
  const bool loaded = ok && tree->GetEntry(0) > 0 && reconstruction_vertex_source && topology_vertex_source;
  if (loaded)
  {
    value.reconstruction_vertex_source = *reconstruction_vertex_source;
    value.topology_vertex_source = *topology_vertex_source;
  }
  tree->ResetBranchAddresses();
  return loaded && std::isfinite(value.max_abs_mbd_vertex_z) && value.max_abs_mbd_vertex_z > 0.0 &&
      std::isfinite(value.shower_shape_min_cluster_et) && value.shower_shape_min_cluster_et >= 0.0 &&
      std::isfinite(value.candidate_abs_eta_max) && value.candidate_abs_eta_max > 0.0 &&
      std::isfinite(value.topology_anchor_abs_eta_max) && value.topology_anchor_abs_eta_max > 0.0 &&
      std::isfinite(value.candidate_et_min) && std::isfinite(value.candidate_et_max) && value.candidate_et_min >= 0.0 &&
      value.candidate_et_max > value.candidate_et_min && value.reconstruction_vertex_source == "MbdVertexMap" &&
      value.topology_vertex_source == "signal HepMC collision";
}

bool compatible(const MbdMapSettings& value, const MbdMapSettings& reference)
{
  return same_double(value.max_abs_mbd_vertex_z, reference.max_abs_mbd_vertex_z) &&
      same_double(value.shower_shape_min_cluster_et, reference.shower_shape_min_cluster_et) &&
      same_double(value.candidate_abs_eta_max, reference.candidate_abs_eta_max) &&
      same_double(value.topology_anchor_abs_eta_max, reference.topology_anchor_abs_eta_max) &&
      same_double(value.candidate_et_min, reference.candidate_et_min) && same_double(value.candidate_et_max, reference.candidate_et_max) &&
      value.reconstruction_vertex_source == reference.reconstruction_vertex_source && value.topology_vertex_source == reference.topology_vertex_source;
}

bool read_metadata(const std::string& path, MapMetadata& value)
{
  TFile file(path.c_str(), "READ");
  auto* tree = file.Get<TTree>("metadata");
  if (file.IsZombie() || !tree || tree->GetEntries() != 1) return false;

  std::string* input_manifest = nullptr;
  std::string* sample_name = nullptr;
  std::string* analysis_release = nullptr;
  std::string* model_sha256 = nullptr;
  bool ok = true;
  // Schema-11 maps with truth isolation but no explicit radius used R=0.4.
  value.truth_isolation_radius = tree->GetBranch("truth_isolation_definition") ? 0.4 : -1.0;
  if (tree->GetBranch("truth_isolation_radius")) ok &= bind(tree, "truth_isolation_radius", &value.truth_isolation_radius);
  ok &= bind(tree, "schema_version", &value.schema_version);
  ok &= bind(tree, "common_composition_algorithm_version", &value.common_composition_algorithm_version);
  ok &= bind(tree, "decay_origin_algorithm_version", &value.decay_origin_algorithm_version);
  ok &= bind(tree, "decay_origin_max_proper_length_cm", &value.decay_origin_max_proper_length_cm);
  ok &= bind(tree, "decay_origin_max_flight_cm", &value.decay_origin_max_flight_cm);
  ok &= bind(tree, "input_manifest", &input_manifest);
  ok &= bind(tree, "sample_name", &sample_name);
  ok &= bind(tree, "analysis_release", &analysis_release);
  ok &= bind(tree, "model_sha256", &model_sha256);
  ok &= bind(tree, "manifest_begin", &value.manifest_begin);
  ok &= bind(tree, "manifest_end", &value.manifest_end);
  ok &= bind(tree, "input_file_count", &value.input_file_count);
  ok &= bind(tree, "map_chunk_id", &value.map_chunk_id);
  ok &= bind(tree, "sample_cross_section_pb", &value.cross_section_pb);
  ok &= bind(tree, "sample_window_min", &value.window_min);
  ok &= bind(tree, "sample_window_max", &value.window_max);
  ok &= bind(tree, "min_cluster_energy", &value.min_cluster_energy);
  ok &= bind(tree, "partner_diagnostic_min_cluster_energy", &value.partner_diagnostic_min_cluster_energy);
  ok &= bind(tree, "pi0_partner_min_energy", &value.pi0_partner_min_energy);
  ok &= bind(tree, "eta_partner_min_energy", &value.eta_partner_min_energy);
  ok &= bind(tree, "pi0_mass_min", &value.pi0_mass_min);
  ok &= bind(tree, "pi0_mass_max", &value.pi0_mass_max);
  ok &= bind(tree, "eta_mass_min", &value.eta_mass_min);
  ok &= bind(tree, "eta_mass_max", &value.eta_mass_max);
  ok &= bind(tree, "missing_energy_min", &value.missing_energy_min);
  ok &= bind(tree, "missing_energy_max", &value.missing_energy_max);
  ok &= bind(tree, "min_photon_energy_recovery", &value.min_photon_energy_recovery);
  ok &= bind(tree, "pi0_topology_algorithm_version", &value.pi0_topology_algorithm_version);
  ok &= bind(tree, "eta_topology_algorithm_version", &value.eta_topology_algorithm_version);
  ok &= bind(tree, "sample_upper_unbounded", &value.upper_unbounded);
  ok &= bind(tree, "sum_generator_weight_processed", &value.sum_generator_weight_processed);
  ok &= bind(tree, "n_events_processed", &value.events_processed);
  ok &= bind(tree, "n_events_written", &value.events_written);
  ok &= bind(tree, "n_events_vertex_rejected", &value.events_vertex_rejected);
  ok &= bind(tree, "n_events_invalid", &value.events_invalid);
  if (!ok || tree->GetEntry(0) <= 0 || !input_manifest || !sample_name || !analysis_release || !model_sha256) return false;

  value.path = path;
  value.input_manifest = *input_manifest;
  value.sample_name = *sample_name;
  value.analysis_release = *analysis_release;
  value.model_sha256 = *model_sha256;
  if (value.schema_version != kMbdMapSchema)
  {
    std::cerr << "Unsupported map schema " << value.schema_version << " in " << path << ", expected " << kMbdMapSchema << std::endl;
    return false;
  }
  return read_mbd_map_settings(tree, value.mbd);
}

bool valid_metadata(const MapMetadata& value, const SampleDefinition& sample)
{
  return value.common_composition_algorithm_version == 1 && value.decay_origin_algorithm_version == 1 &&
      same_double(value.decay_origin_max_proper_length_cm, 1e-4) && same_double(value.decay_origin_max_flight_cm, 1e-2) && value.sample_name == sample.name && !value.input_manifest.empty() &&
      !value.analysis_release.empty() && !value.model_sha256.empty() && value.manifest_begin >= 0 &&
      std::isfinite(value.min_cluster_energy) && value.min_cluster_energy >= 0.0 &&
      std::isfinite(value.partner_diagnostic_min_cluster_energy) && same_double(value.partner_diagnostic_min_cluster_energy, 0.0) &&
      std::isfinite(value.pi0_partner_min_energy) && value.pi0_partner_min_energy >= 0.0 && value.pi0_topology_algorithm_version == 14 && value.eta_topology_algorithm_version == 4 &&
      std::isfinite(value.eta_partner_min_energy) && value.eta_partner_min_energy >= 0.0 &&
      std::isfinite(value.pi0_mass_min) && value.pi0_mass_min >= 0.0 &&
      std::isfinite(value.pi0_mass_max) && value.pi0_mass_max >= 0.0 &&
      std::isfinite(value.eta_mass_min) && value.eta_mass_min >= 0.0 &&
      std::isfinite(value.eta_mass_max) && value.eta_mass_max >= 0.0 &&
      std::isfinite(value.missing_energy_min) && value.missing_energy_min >= 0.0 &&
      std::isfinite(value.missing_energy_max) && value.missing_energy_max >= 0.0 &&
      std::isfinite(value.min_photon_energy_recovery) && value.min_photon_energy_recovery >= 0.0 &&
      value.pi0_mass_min < value.pi0_mass_max && value.eta_mass_min < value.eta_mass_max &&
      value.missing_energy_min < value.missing_energy_max && same_double(value.min_photon_energy_recovery, 0.0) &&
      value.manifest_end > value.manifest_begin && value.input_file_count == value.manifest_end - value.manifest_begin &&
      same_double(value.cross_section_pb, sample.cross_section_pb) && same_double(value.window_min, sample.window_min) &&
      same_double(value.window_max, sample.window_max) && value.upper_unbounded == sample.upper_unbounded &&
      std::isfinite(value.sum_generator_weight_processed) &&
      value.events_processed == value.events_written + value.events_vertex_rejected + value.events_invalid;
}

bool compatible(const MapMetadata& value, const MapMetadata& reference)
{
  return value.common_composition_algorithm_version == reference.common_composition_algorithm_version && value.schema_version == reference.schema_version && value.decay_origin_algorithm_version == reference.decay_origin_algorithm_version &&
      same_double(value.decay_origin_max_proper_length_cm, reference.decay_origin_max_proper_length_cm) &&
      same_double(value.decay_origin_max_flight_cm, reference.decay_origin_max_flight_cm) && value.input_manifest == reference.input_manifest &&
      value.sample_name == reference.sample_name && value.analysis_release == reference.analysis_release &&
      value.model_sha256 == reference.model_sha256 && same_double(value.cross_section_pb, reference.cross_section_pb) &&
      same_double(value.window_min, reference.window_min) && same_double(value.window_max, reference.window_max) &&
      value.upper_unbounded == reference.upper_unbounded &&
      same_double(value.min_cluster_energy, reference.min_cluster_energy) &&
      same_double(value.partner_diagnostic_min_cluster_energy, reference.partner_diagnostic_min_cluster_energy) &&
      same_double(value.pi0_partner_min_energy, reference.pi0_partner_min_energy) &&
      same_double(value.eta_partner_min_energy, reference.eta_partner_min_energy) &&
      same_double(value.pi0_mass_min, reference.pi0_mass_min) &&
      same_double(value.pi0_mass_max, reference.pi0_mass_max) &&
      same_double(value.eta_mass_min, reference.eta_mass_min) &&
      same_double(value.eta_mass_max, reference.eta_mass_max) &&
      same_double(value.missing_energy_min, reference.missing_energy_min) &&
      same_double(value.missing_energy_max, reference.missing_energy_max) &&
      same_double(value.min_photon_energy_recovery, reference.min_photon_energy_recovery) && value.pi0_topology_algorithm_version == reference.pi0_topology_algorithm_version &&
      value.eta_topology_algorithm_version == reference.eta_topology_algorithm_version &&
      compatible(value.mbd, reference.mbd);
}

bool collect_maps(const std::string& pattern, const SampleDefinition& sample, bool require_complete, std::vector<MapMetadata>& maps)
{
  TChain chain("metadata");
  const int matched = chain.Add(pattern.c_str());
  const TObjArray* files = chain.GetListOfFiles();
  if (matched <= 0 || !files || files->GetEntries() <= 0)
  {
    if (require_complete)
    {
      std::cerr << "No maps matched for " << sample.name << ": " << pattern << std::endl;
      return false;
    }
    return true;
  }

  std::set<std::string> unique_paths;
  for (int index = 0; index < files->GetEntries(); ++index)
  {
    const TObject* element = files->At(index);
    const std::string path = element ? element->GetTitle() : "";
    MapMetadata metadata;
    if (path.empty() || !unique_paths.insert(path).second || !read_metadata(path, metadata) || !valid_metadata(metadata, sample))
    {
      std::cerr << "Invalid map metadata: " << path << std::endl;
      return false;
    }
    maps.push_back(metadata);
  }
  std::sort(maps.begin(), maps.end(), [](const auto& left, const auto& right) { return left.manifest_begin < right.manifest_begin; });

  const MapMetadata& reference = maps.front();
  long long next_begin = 0;
  for (std::size_t index = 0; index < maps.size(); ++index)
  {
    const MapMetadata& map = maps[index];
    if (!compatible(map, reference) || map.manifest_begin != next_begin || map.map_chunk_id != index)
    {
      std::cerr << "Incompatible or noncontiguous map: " << map.path << ", expected begin/chunk " << next_begin << "/" << index << std::endl;
      return false;
    }
    next_begin = map.manifest_end;
  }
  if (next_begin > sample.expected_manifest_end || (require_complete && next_begin != sample.expected_manifest_end))
  {
    std::cerr << "Incomplete or over-complete sample " << sample.name << ": end=" << next_begin
              << ", expected=" << sample.expected_manifest_end << std::endl;
    return false;
  }
  return true;
}

struct Spectra
{
  std::array<std::unique_ptr<TH1D>, kSpectrumCount> counts;
  std::array<std::unique_ptr<TH1D>, kSpectrumCount> weighted_pb;

  Spectra(int n_bins, double et_max)
  {
    for (std::size_t index = 0; index < kSpectrumCount; ++index)
    {
      const std::string count_name = std::string("h_region_a_") + kKeys[index] + "_et_count";
      const std::string weighted_name = std::string("h_region_a_") + kKeys[index] + "_et_pb";
      counts[index] = std::make_unique<TH1D>(count_name.c_str(), "", n_bins, 0.0, et_max);
      weighted_pb[index] = std::make_unique<TH1D>(weighted_name.c_str(), "", n_bins, 0.0, et_max);
      counts[index]->SetDirectory(nullptr);
      weighted_pb[index]->SetDirectory(nullptr);
      counts[index]->Sumw2();
      weighted_pb[index]->Sumw2();
    }
  }

  void fill(std::size_t index, double et, double weight)
  {
    counts[index]->Fill(et);
    weighted_pb[index]->Fill(et, weight);
  }
};

bool pass_selection(std::size_t selection, unsigned char pass_kinematics, unsigned char pass_preselection,
                    unsigned char pass_tight, unsigned char pass_isolated, unsigned char pass_region_a,
                    unsigned char pi0_tag, unsigned char eta_tag)
{
  const bool kinematic = pass_kinematics != 0U;
  const bool preselection = kinematic && pass_preselection != 0U;
  switch (selection)
  {
    case 0: return kinematic;
    case 1: return preselection;
    case 2: return preselection && pass_tight != 0U;
    case 3: return preselection && pass_isolated != 0U;
    case 4: return pass_region_a != 0U;
    case 5: return pass_region_a != 0U && pi0_tag == 0U && eta_tag == 0U;
    default: return false;
  }
}

bool valid_partition(const std::array<std::unique_ptr<TH1D>, kSpectrumCount>& histograms)
{
  for (int bin = 0; bin <= histograms[0]->GetNbinsX() + 1; ++bin)
  {
    double missing_sum = 0.0;
    for (std::size_t index = 6; index <= 11; ++index) missing_sum += histograms[index]->GetBinContent(bin);
    const double category_sum = histograms[2]->GetBinContent(bin) + histograms[3]->GetBinContent(bin) +
        histograms[4]->GetBinContent(bin) + histograms[5]->GetBinContent(bin) + histograms[12]->GetBinContent(bin);
    const double scale = std::max({1.0, std::abs(histograms[1]->GetBinContent(bin)), std::abs(category_sum), std::abs(missing_sum)});
    if (std::abs(histograms[5]->GetBinContent(bin) - missing_sum) > 1e-9 * scale ||
        std::abs(histograms[1]->GetBinContent(bin) - category_sum) > 1e-9 * scale) return false;
    double missing_w2 = 0, category_w2 = 0;
    for (std::size_t i = 6; i <= 11; ++i) missing_w2 += std::pow(histograms[i]->GetBinError(bin),2);
    for (int i : {2,3,4,5,12}) category_w2 += std::pow(histograms[i]->GetBinError(bin),2);
    if (std::abs(missing_w2 - std::pow(histograms[5]->GetBinError(bin),2)) > 1e-9 * std::max(1.,missing_w2) ||
        std::abs(category_w2 - std::pow(histograms[1]->GetBinError(bin),2)) > 1e-9 * std::max(1.,category_w2)) return false;
  }
  return true;
}

bool make_output_directory(const std::string& output_base)
{
  const std::size_t slash = output_base.find_last_of('/');
  if (slash == std::string::npos) return true;
  const std::string directory = output_base.substr(0, slash);
  return directory.empty() || !gSystem->AccessPathName(directory.c_str()) || gSystem->mkdir(directory.c_str(), true) == 0;
}

}

#endif
