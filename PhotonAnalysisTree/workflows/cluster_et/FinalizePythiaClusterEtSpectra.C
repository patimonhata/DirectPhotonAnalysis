#include "../../macro/Utilities/sPhenixStyle.C"

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TObjArray.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace
{
constexpr std::size_t kHistogramCount = 8;
constexpr std::size_t kPlotHistogramCount = 5;

struct PartialMetadata
{
  std::string path;
  int schema_version = 0;
  std::string manifest_path;
  std::string cluster_collection;
  std::string prompt_selection;
  std::string pi0_selection;
  std::string topology_priority;
  std::string projection_scheme;
  std::string raw_truth_tower_node;
  std::string truth_cell_node;
  std::string truth_hit_node;
  std::string energy_topology_priority;
  std::string energy_matching_scheme;
  std::string energy_candidate_selection;
  long long manifest_begin = -1;
  long long manifest_end = -1;
  int signal_embedding_id = 0;
  int n_bins = 0;
  int pi0_truth_matching_algorithm_version = 0;
  double et_max = 0.0;
  double truth_eta_max = 0.0;
  double cluster_eta_max = 0.0;
  double min_cluster_energy = 0.0;
  double dominant_fraction_min = 0.0;
  double pi0_contributor_fraction_min = 0.0;
  double min_energy_contribution_fraction = 0.0;
  double separated_delta_r_cut = 0.0;
  double merged_delta_r_cut = 0.0;
  double response_min = 0.0;
  double response_max = 0.0;
  unsigned char bin_width_normalized = 1U;
  unsigned long long events_processed = 0;
  unsigned long long events_written = 0;
  unsigned long long events_invalid = 0;
  unsigned long long prompt_cluster_count = 0;
  unsigned long long pi0_cluster_count = 0;
  unsigned long long pi0_cluster_g4_decay_count = 0;
  unsigned long long pi0_cluster_generator_decay_count = 0;
  unsigned long long pi0_candidate_g4_decay_count = 0;
  unsigned long long pi0_candidate_generator_decay_count = 0;
  unsigned long long pi0_malformed_daughters_count = 0;
  unsigned long long pi0_projection_failure_count = 0;
  unsigned long long pi0_separated_count = 0;
  unsigned long long pi0_merged_count = 0;
  unsigned long long pi0_missing_count = 0;
  unsigned long long pi0_none_count = 0;
  unsigned long long pi0_ambiguous_count = 0;
  unsigned long long pi0_separated_cluster_fill_count = 0;
  unsigned long long pi0_merged_cluster_fill_count = 0;
  unsigned long long pi0_missing_cluster_fill_count = 0;
  unsigned long long pi0_energy_separated_count = 0;
  unsigned long long pi0_energy_merged_count = 0;
  unsigned long long pi0_energy_missing_count = 0;
  unsigned long long pi0_energy_none_count = 0;
  unsigned long long pi0_energy_match_invalid_count = 0;
  unsigned long long pi0_energy_separated_cluster_fill_count = 0;
  unsigned long long pi0_energy_merged_cluster_fill_count = 0;
  unsigned long long pi0_energy_missing_cluster_fill_count = 0;
};

template <class T>
bool bind(TTree* tree, const char* name, T* address)
{
  return tree->GetBranch(name) && tree->SetBranchAddress(name, address) >= 0;
}

bool read_metadata(const std::string& path, PartialMetadata& value)
{
  TFile input(path.c_str(), "READ");
  TTree* tree = nullptr;
  input.GetObject("metadata", tree);
  if (input.IsZombie() || !tree || tree->GetEntries() != 1)
  {
    return false;
  }

  std::string* manifest_path = nullptr;
  std::string* cluster_collection = nullptr;
  std::string* prompt_selection = nullptr;
  std::string* pi0_selection = nullptr;
  std::string* topology_priority = nullptr;
  std::string* projection_scheme = nullptr;
  std::string* raw_truth_tower_node = nullptr;
  std::string* truth_cell_node = nullptr;
  std::string* truth_hit_node = nullptr;
  std::string* energy_topology_priority = nullptr;
  std::string* energy_matching_scheme = nullptr;
  std::string* energy_candidate_selection = nullptr;

  bool ok = true;
  ok &= bind(tree, "schema_version", &value.schema_version);
  ok &= bind(tree, "manifest_path", &manifest_path);
  ok &= bind(tree, "manifest_begin", &value.manifest_begin);
  ok &= bind(tree, "manifest_end", &value.manifest_end);
  ok &= bind(tree, "cluster_collection", &cluster_collection);
  ok &= bind(tree, "prompt_selection", &prompt_selection);
  ok &= bind(tree, "pi0_selection", &pi0_selection);
  ok &= bind(tree, "topology_priority", &topology_priority);
  ok &= bind(tree, "projection_scheme", &projection_scheme);
  ok &= bind(tree, "raw_truth_tower_node", &raw_truth_tower_node);
  ok &= bind(tree, "truth_cell_node", &truth_cell_node);
  ok &= bind(tree, "truth_hit_node", &truth_hit_node);
  ok &= bind(tree, "energy_topology_priority", &energy_topology_priority);
  ok &= bind(tree, "energy_matching_scheme", &energy_matching_scheme);
  ok &= bind(tree, "energy_candidate_selection", &energy_candidate_selection);
  ok &= bind(tree, "signal_embedding_id", &value.signal_embedding_id);
  ok &= bind(tree, "n_bins", &value.n_bins);
  ok &= bind(tree, "pi0_truth_matching_algorithm_version",
             &value.pi0_truth_matching_algorithm_version);
  ok &= bind(tree, "et_max", &value.et_max);
  ok &= bind(tree, "truth_eta_max", &value.truth_eta_max);
  ok &= bind(tree, "cluster_eta_max", &value.cluster_eta_max);
  ok &= bind(tree, "min_cluster_energy", &value.min_cluster_energy);
  ok &= bind(tree, "dominant_fraction_min", &value.dominant_fraction_min);
  ok &= bind(tree, "pi0_contributor_fraction_min",
             &value.pi0_contributor_fraction_min);
  ok &= bind(tree, "min_energy_contribution_fraction",
             &value.min_energy_contribution_fraction);
  ok &= bind(tree, "separated_delta_r_cut", &value.separated_delta_r_cut);
  ok &= bind(tree, "merged_delta_r_cut", &value.merged_delta_r_cut);
  ok &= bind(tree, "response_min", &value.response_min);
  ok &= bind(tree, "response_max", &value.response_max);
  ok &= bind(tree, "bin_width_normalized", &value.bin_width_normalized);
  ok &= bind(tree, "events_processed", &value.events_processed);
  ok &= bind(tree, "events_written", &value.events_written);
  ok &= bind(tree, "events_invalid", &value.events_invalid);
  ok &= bind(tree, "prompt_cluster_count", &value.prompt_cluster_count);
  ok &= bind(tree, "pi0_cluster_count", &value.pi0_cluster_count);
  ok &= bind(tree, "pi0_cluster_g4_decay_count",
             &value.pi0_cluster_g4_decay_count);
  ok &= bind(tree, "pi0_cluster_generator_decay_count",
             &value.pi0_cluster_generator_decay_count);
  ok &= bind(tree, "pi0_candidate_g4_decay_count",
             &value.pi0_candidate_g4_decay_count);
  ok &= bind(tree, "pi0_candidate_generator_decay_count",
             &value.pi0_candidate_generator_decay_count);
  ok &= bind(tree, "pi0_malformed_daughters_count",
             &value.pi0_malformed_daughters_count);
  ok &= bind(tree, "pi0_projection_failure_count",
             &value.pi0_projection_failure_count);
  ok &= bind(tree, "pi0_separated_count", &value.pi0_separated_count);
  ok &= bind(tree, "pi0_merged_count", &value.pi0_merged_count);
  ok &= bind(tree, "pi0_missing_count", &value.pi0_missing_count);
  ok &= bind(tree, "pi0_none_count", &value.pi0_none_count);
  ok &= bind(tree, "pi0_ambiguous_count", &value.pi0_ambiguous_count);
  ok &= bind(tree, "pi0_separated_cluster_fill_count",
             &value.pi0_separated_cluster_fill_count);
  ok &= bind(tree, "pi0_merged_cluster_fill_count",
             &value.pi0_merged_cluster_fill_count);
  ok &= bind(tree, "pi0_missing_cluster_fill_count",
             &value.pi0_missing_cluster_fill_count);
  ok &= bind(tree, "pi0_energy_separated_count",
             &value.pi0_energy_separated_count);
  ok &= bind(tree, "pi0_energy_merged_count",
             &value.pi0_energy_merged_count);
  ok &= bind(tree, "pi0_energy_missing_count",
             &value.pi0_energy_missing_count);
  ok &= bind(tree, "pi0_energy_none_count",
             &value.pi0_energy_none_count);
  ok &= bind(tree, "pi0_energy_match_invalid_count",
             &value.pi0_energy_match_invalid_count);
  ok &= bind(tree, "pi0_energy_separated_cluster_fill_count",
             &value.pi0_energy_separated_cluster_fill_count);
  ok &= bind(tree, "pi0_energy_merged_cluster_fill_count",
             &value.pi0_energy_merged_cluster_fill_count);
  ok &= bind(tree, "pi0_energy_missing_cluster_fill_count",
             &value.pi0_energy_missing_cluster_fill_count);

  if (!ok || tree->GetEntry(0) <= 0 || !manifest_path || !cluster_collection ||
      !prompt_selection || !pi0_selection || !topology_priority ||
      !projection_scheme || !raw_truth_tower_node || !truth_cell_node ||
      !truth_hit_node || !energy_topology_priority || !energy_matching_scheme ||
      !energy_candidate_selection)
  {
    return false;
  }

  value.path = path;
  value.manifest_path = *manifest_path;
  value.cluster_collection = *cluster_collection;
  value.prompt_selection = *prompt_selection;
  value.pi0_selection = *pi0_selection;
  value.topology_priority = *topology_priority;
  value.projection_scheme = *projection_scheme;
  value.raw_truth_tower_node = *raw_truth_tower_node;
  value.truth_cell_node = *truth_cell_node;
  value.truth_hit_node = *truth_hit_node;
  value.energy_topology_priority = *energy_topology_priority;
  value.energy_matching_scheme = *energy_matching_scheme;
  value.energy_candidate_selection = *energy_candidate_selection;
  return true;
}

bool same_double(double first, double second)
{
  return std::abs(first - second) < 1e-12;
}

bool compatible(const PartialMetadata& value, const PartialMetadata& reference)
{
  const unsigned long long candidate_count =
      value.pi0_candidate_g4_decay_count +
      value.pi0_candidate_generator_decay_count;
  return value.schema_version == 2 &&
      value.manifest_path == reference.manifest_path &&
      value.cluster_collection == reference.cluster_collection &&
      value.prompt_selection == reference.prompt_selection &&
      value.pi0_selection == reference.pi0_selection &&
      value.topology_priority == reference.topology_priority &&
      value.projection_scheme == reference.projection_scheme &&
      value.raw_truth_tower_node == reference.raw_truth_tower_node &&
      value.truth_cell_node == reference.truth_cell_node &&
      value.truth_hit_node == reference.truth_hit_node &&
      value.energy_topology_priority == reference.energy_topology_priority &&
      value.energy_matching_scheme == reference.energy_matching_scheme &&
      value.energy_candidate_selection == reference.energy_candidate_selection &&
      value.signal_embedding_id == reference.signal_embedding_id &&
      value.pi0_truth_matching_algorithm_version ==
          reference.pi0_truth_matching_algorithm_version &&
      value.n_bins == reference.n_bins &&
      same_double(value.et_max, reference.et_max) &&
      same_double(value.truth_eta_max, reference.truth_eta_max) &&
      same_double(value.cluster_eta_max, reference.cluster_eta_max) &&
      same_double(value.min_cluster_energy, reference.min_cluster_energy) &&
      same_double(value.dominant_fraction_min, reference.dominant_fraction_min) &&
      same_double(value.pi0_contributor_fraction_min,
                  reference.pi0_contributor_fraction_min) &&
      same_double(value.min_energy_contribution_fraction,
                  reference.min_energy_contribution_fraction) &&
      same_double(value.separated_delta_r_cut,
                  reference.separated_delta_r_cut) &&
      same_double(value.merged_delta_r_cut, reference.merged_delta_r_cut) &&
      same_double(value.response_min, reference.response_min) &&
      same_double(value.response_max, reference.response_max) &&
      value.bin_width_normalized == 0U && value.events_processed > 0 &&
      value.events_written == value.events_processed && value.events_invalid == 0 &&
      value.pi0_cluster_count == value.pi0_cluster_g4_decay_count +
          value.pi0_cluster_generator_decay_count &&
      value.pi0_separated_cluster_fill_count ==
          2ULL * value.pi0_separated_count &&
      value.pi0_merged_cluster_fill_count == value.pi0_merged_count &&
      value.pi0_missing_cluster_fill_count == value.pi0_missing_count &&
      candidate_count == value.pi0_separated_count + value.pi0_merged_count +
          value.pi0_missing_count + value.pi0_none_count +
          value.pi0_ambiguous_count + value.pi0_projection_failure_count &&
      value.pi0_energy_separated_cluster_fill_count ==
          2ULL * value.pi0_energy_separated_count &&
      value.pi0_energy_merged_cluster_fill_count ==
          value.pi0_energy_merged_count &&
      value.pi0_energy_missing_cluster_fill_count ==
          value.pi0_energy_missing_count &&
      candidate_count == value.pi0_energy_separated_count +
          value.pi0_energy_merged_count + value.pi0_energy_missing_count +
          value.pi0_energy_none_count;
}

bool valid_histogram(const TH1D* histogram, const PartialMetadata& metadata,
                     unsigned long long entries)
{
  if (!histogram || histogram->GetNbinsX() != metadata.n_bins ||
      std::abs(histogram->GetXaxis()->GetXmin()) > 1e-12 ||
      !same_double(histogram->GetXaxis()->GetXmax(), metadata.et_max) ||
      histogram->GetSumw2N() == 0 ||
      std::abs(histogram->GetEntries() - static_cast<double>(entries)) > 0.5)
  {
    return false;
  }
  for (int bin = 0; bin <= histogram->GetNbinsX() + 1; ++bin)
  {
    if (!std::isfinite(histogram->GetBinContent(bin)) ||
        !std::isfinite(histogram->GetBinError(bin)))
    {
      return false;
    }
  }
  return true;
}

bool make_output_directory(const std::string& output_base)
{
  const std::size_t slash = output_base.find_last_of('/');
  if (slash == std::string::npos)
  {
    return true;
  }
  const std::string directory = output_base.substr(0, slash);
  return directory.empty() || !gSystem->AccessPathName(directory.c_str()) ||
      gSystem->mkdir(directory.c_str(), true) == 0;
}

double smallest_positive(
    const std::array<std::unique_ptr<TH1D>, kHistogramCount>& histograms,
    const std::array<std::size_t, kPlotHistogramCount>& indices)
{
  double result = std::numeric_limits<double>::infinity();
  for (const std::size_t index : indices)
  {
    const TH1D* histogram = histograms[index].get();
    for (int bin = 1; bin <= histogram->GetNbinsX(); ++bin)
    {
      const double value = histogram->GetBinContent(bin);
      if (value > 0.0 && value < result)
      {
        result = value;
      }
    }
  }
  return std::isfinite(result) ? result : 0.0;
}
}

int FinalizePythiaClusterEtSpectra(
    const std::string partial_pattern =
        "output/cluster_et_partial/prompt_pi0_eta07_energy_contribution/partial_*.root",
    const std::string output_base =
        "output/plots/newnewtempminbias_cluster_et_prompt_pi0_eta07",
    const long long expected_manifest_begin = 0,
    const long long expected_manifest_end = -1)
{
  if (partial_pattern.empty() || output_base.empty() ||
      expected_manifest_begin < 0 || expected_manifest_end < -1 ||
      (expected_manifest_end >= 0 &&
       expected_manifest_end <= expected_manifest_begin))
  {
    return 1;
  }

  TChain chain("metadata");
  const int matched = chain.Add(partial_pattern.c_str());
  const TObjArray* files = chain.GetListOfFiles();
  if (matched <= 0 || !files || files->GetEntries() <= 0)
  {
    std::cerr << "FinalizePythiaClusterEtSpectra - no partials matched"
              << std::endl;
    return 2;
  }

  std::vector<PartialMetadata> partials;
  std::set<std::string> unique_paths;
  for (int index = 0; index < files->GetEntries(); ++index)
  {
    const TObject* element = files->At(index);
    const std::string path = element ? element->GetTitle() : "";
    PartialMetadata metadata;
    if (path.empty() || !unique_paths.insert(path).second ||
        !read_metadata(path, metadata))
    {
      std::cerr << "FinalizePythiaClusterEtSpectra - invalid partial: "
                << path << std::endl;
      return 3;
    }
    partials.push_back(metadata);
  }

  std::sort(partials.begin(), partials.end(),
      [](const auto& left, const auto& right) {
        return left.manifest_begin < right.manifest_begin;
      });
  const PartialMetadata& reference = partials.front();
  long long next_begin = expected_manifest_begin;
  for (const PartialMetadata& partial : partials)
  {
    if (!compatible(partial, reference) ||
        partial.manifest_begin != next_begin)
    {
      std::cerr
          << "FinalizePythiaClusterEtSpectra - incompatible or noncontiguous partial: "
          << partial.path << ", expected begin " << next_begin << std::endl;
      return 4;
    }
    next_begin = partial.manifest_end;
  }
  if (expected_manifest_end >= 0 && next_begin != expected_manifest_end)
  {
    return 4;
  }

  const std::array<std::string, kHistogramCount> raw_names = {
      "h_prompt_cluster_et_raw",
      "h_pi0_cluster_et_raw",
      "h_pi0_separated_cluster_et_raw",
      "h_pi0_merged_cluster_et_raw",
      "h_pi0_missing_cluster_et_raw",
      "h_pi0_separated_energy_contribution_cluster_et_raw",
      "h_pi0_merged_energy_contribution_cluster_et_raw",
      "h_pi0_missing_energy_contribution_cluster_et_raw"};
  const std::array<std::string, kHistogramCount> density_names = {
      "h_prompt_cluster_et_density",
      "h_pi0_cluster_et_density",
      "h_pi0_separated_cluster_et_density",
      "h_pi0_merged_cluster_et_density",
      "h_pi0_missing_cluster_et_density",
      "h_pi0_separated_energy_contribution_cluster_et_density",
      "h_pi0_merged_energy_contribution_cluster_et_density",
      "h_pi0_missing_energy_contribution_cluster_et_density"};

  std::array<std::unique_ptr<TH1D>, kHistogramCount> raw;
  for (std::size_t index = 0; index < raw.size(); ++index)
  {
    raw[index] = std::make_unique<TH1D>(
        raw_names[index].c_str(), "", reference.n_bins, 0.0,
        reference.et_max);
    raw[index]->Sumw2();
  }

  PartialMetadata total = reference;
  total.events_processed = total.events_written = total.events_invalid = 0;
  total.prompt_cluster_count = total.pi0_cluster_count = 0;
  total.pi0_cluster_g4_decay_count =
      total.pi0_cluster_generator_decay_count = 0;
  total.pi0_candidate_g4_decay_count =
      total.pi0_candidate_generator_decay_count = 0;
  total.pi0_malformed_daughters_count =
      total.pi0_projection_failure_count = 0;
  total.pi0_separated_count = total.pi0_merged_count =
      total.pi0_missing_count = 0;
  total.pi0_none_count = total.pi0_ambiguous_count = 0;
  total.pi0_separated_cluster_fill_count =
      total.pi0_merged_cluster_fill_count =
      total.pi0_missing_cluster_fill_count = 0;
  total.pi0_energy_separated_count = total.pi0_energy_merged_count =
      total.pi0_energy_missing_count = total.pi0_energy_none_count = 0;
  total.pi0_energy_match_invalid_count = 0;
  total.pi0_energy_separated_cluster_fill_count =
      total.pi0_energy_merged_cluster_fill_count =
      total.pi0_energy_missing_cluster_fill_count = 0;

  for (const PartialMetadata& partial : partials)
  {
    TFile input(partial.path.c_str(), "READ");
    const std::array<unsigned long long, kHistogramCount> expected = {
        partial.prompt_cluster_count,
        partial.pi0_cluster_count,
        partial.pi0_separated_cluster_fill_count,
        partial.pi0_merged_cluster_fill_count,
        partial.pi0_missing_cluster_fill_count,
        partial.pi0_energy_separated_cluster_fill_count,
        partial.pi0_energy_merged_cluster_fill_count,
        partial.pi0_energy_missing_cluster_fill_count};
    for (std::size_t index = 0; index < raw.size(); ++index)
    {
      TH1D* histogram = nullptr;
      input.GetObject(raw_names[index].c_str(), histogram);
      if (input.IsZombie() ||
          !valid_histogram(histogram, partial, expected[index]) ||
          !raw[index]->Add(histogram))
      {
        std::cerr << "FinalizePythiaClusterEtSpectra - invalid histogram in "
                  << partial.path << std::endl;
        return 5;
      }
    }

    total.events_processed += partial.events_processed;
    total.events_written += partial.events_written;
    total.events_invalid += partial.events_invalid;
    total.prompt_cluster_count += partial.prompt_cluster_count;
    total.pi0_cluster_count += partial.pi0_cluster_count;
    total.pi0_cluster_g4_decay_count +=
        partial.pi0_cluster_g4_decay_count;
    total.pi0_cluster_generator_decay_count +=
        partial.pi0_cluster_generator_decay_count;
    total.pi0_candidate_g4_decay_count +=
        partial.pi0_candidate_g4_decay_count;
    total.pi0_candidate_generator_decay_count +=
        partial.pi0_candidate_generator_decay_count;
    total.pi0_malformed_daughters_count +=
        partial.pi0_malformed_daughters_count;
    total.pi0_projection_failure_count +=
        partial.pi0_projection_failure_count;
    total.pi0_separated_count += partial.pi0_separated_count;
    total.pi0_merged_count += partial.pi0_merged_count;
    total.pi0_missing_count += partial.pi0_missing_count;
    total.pi0_none_count += partial.pi0_none_count;
    total.pi0_ambiguous_count += partial.pi0_ambiguous_count;
    total.pi0_separated_cluster_fill_count +=
        partial.pi0_separated_cluster_fill_count;
    total.pi0_merged_cluster_fill_count +=
        partial.pi0_merged_cluster_fill_count;
    total.pi0_missing_cluster_fill_count +=
        partial.pi0_missing_cluster_fill_count;
    total.pi0_energy_separated_count +=
        partial.pi0_energy_separated_count;
    total.pi0_energy_merged_count += partial.pi0_energy_merged_count;
    total.pi0_energy_missing_count += partial.pi0_energy_missing_count;
    total.pi0_energy_none_count += partial.pi0_energy_none_count;
    total.pi0_energy_match_invalid_count +=
        partial.pi0_energy_match_invalid_count;
    total.pi0_energy_separated_cluster_fill_count +=
        partial.pi0_energy_separated_cluster_fill_count;
    total.pi0_energy_merged_cluster_fill_count +=
        partial.pi0_energy_merged_cluster_fill_count;
    total.pi0_energy_missing_cluster_fill_count +=
        partial.pi0_energy_missing_cluster_fill_count;
  }

  std::array<std::unique_ptr<TH1D>, kHistogramCount> density;
  const std::array<int, kHistogramCount> colors = {
      kRed + 1, kBlue + 1, kAzure + 7, kMagenta + 1, kGreen + 2,
      kAzure + 7, kMagenta + 1, kGreen + 2};
  for (std::size_t index = 0; index < raw.size(); ++index)
  {
    density[index].reset(static_cast<TH1D*>(
        raw[index]->Clone(density_names[index].c_str())));
    density[index]->Scale(1.0, "width");
    density[index]->SetStats(false);
    density[index]->SetFillStyle(0);
    density[index]->SetLineWidth(index < 2 ? 3 : 2);
    density[index]->SetLineColor(colors[index]);
    density[index]->GetXaxis()->SetTitle("Cluster E_{T} [GeV]");
    density[index]->GetYaxis()->SetTitle("Clusters / GeV");

    raw[index]->SetStats(false);
    raw[index]->SetFillStyle(0);
    raw[index]->SetLineWidth(index < 2 ? 3 : 2);
    raw[index]->SetLineColor(colors[index]);
    raw[index]->GetXaxis()->SetTitle("Cluster E_{T} [GeV]");
    raw[index]->GetYaxis()->SetTitle("Counts / bin");
  }

  if (!make_output_directory(output_base))
  {
    return 6;
  }
  SetsPhenixStyle();

  const auto draw_plot =
      [&](const std::array<std::size_t, kPlotHistogramCount>& indices,
          const std::string& canvas_name, const std::string& output_path,
          const char* matching_label)
      {
        double maximum = 0.0;
        for (const std::size_t index : indices)
        {
          maximum = std::max(maximum, raw[index]->GetMaximum());
        }
        const double minimum = smallest_positive(raw, indices);
        raw[indices[0]]->SetMinimum(minimum > 0.0 ? 0.5 * minimum : 0.5);
        raw[indices[0]]->SetMaximum(maximum > 0.0 ? 5.0 * maximum : 1.0);

        TCanvas canvas(canvas_name.c_str(), "Pythia cluster ET spectra",
                       1000, 800);
        canvas.SetLogy();
        raw[indices[0]]->Draw("HIST");
        for (std::size_t position = 1;
             position < indices.size(); ++position)
        {
          raw[indices[position]]->Draw("HIST SAME");
        }

        TLegend legend(0.47, 0.38, 0.89, 0.65);
        legend.AddEntry(raw[indices[0]].get(),
                        "Prompt-#gamma cluster", "l");
        legend.AddEntry(raw[indices[1]].get(),
                        "#pi^{0}-origin cluster", "l");
        legend.AddEntry(raw[indices[2]].get(),
                        "#pi^{0}: separated", "l");
        legend.AddEntry(raw[indices[3]].get(),
                        "#pi^{0}: merged", "l");
        legend.AddEntry(raw[indices[4]].get(),
                        "#pi^{0}: missing partner", "l");
        legend.Draw();

        TLatex label;
        label.SetNDC();
        label.SetTextAlign(13);
        label.DrawLatex(0.23, 0.92, "#it{#bf{sPHENIX}} Internal");
        label.DrawLatex(0.23, 0.84, "Pythia8 p+p minimum bias");
        std::ostringstream eta_label;
        eta_label << "|#eta^{truth}| < " << reference.truth_eta_max
                  << ", |#eta^{cluster}| < "
                  << reference.cluster_eta_max;
        label.DrawLatex(0.23, 0.76, eta_label.str().c_str());
        label.DrawLatex(0.23, 0.68, matching_label);
        canvas.RedrawAxis();
        canvas.SaveAs(output_path.c_str());
      };

  draw_plot({0U, 1U, 2U, 3U, 4U},
            "c_pythia_cluster_et_geometric", output_base + ".pdf",
            "Topology: geometrical matching");
  draw_plot({0U, 1U, 5U, 6U, 7U},
            "c_pythia_cluster_et_energy_contribution",
            output_base + "_energy_contribution.pdf",
            "Topology: energy-deposit matching");

  TFile output((output_base + ".root").c_str(), "RECREATE");
  if (output.IsZombie())
  {
    return 6;
  }
  for (std::size_t index = 0; index < raw.size(); ++index)
  {
    raw[index]->Write();
    density[index]->Write();
  }

  int output_schema_version = 2;
  long long manifest_begin = partials.front().manifest_begin;
  long long manifest_end = partials.back().manifest_end;
  long long partial_file_count =
      static_cast<long long>(partials.size());
  long long input_file_count = manifest_end - manifest_begin;
  unsigned char contains_raw_histograms = 1U;
  unsigned char contains_bin_width_normalized_histograms = 1U;
  TTree metadata("metadata", "Final Pythia cluster ET metadata");
  metadata.Branch("schema_version", &output_schema_version);
  metadata.Branch("manifest_path", &total.manifest_path);
  metadata.Branch("manifest_begin", &manifest_begin);
  metadata.Branch("manifest_end", &manifest_end);
  metadata.Branch("partial_file_count", &partial_file_count);
  metadata.Branch("input_file_count", &input_file_count);
  metadata.Branch("cluster_collection", &total.cluster_collection);
  metadata.Branch("prompt_selection", &total.prompt_selection);
  metadata.Branch("pi0_selection", &total.pi0_selection);
  metadata.Branch("topology_priority", &total.topology_priority);
  metadata.Branch("projection_scheme", &total.projection_scheme);
  metadata.Branch("raw_truth_tower_node", &total.raw_truth_tower_node);
  metadata.Branch("truth_cell_node", &total.truth_cell_node);
  metadata.Branch("truth_hit_node", &total.truth_hit_node);
  metadata.Branch("energy_topology_priority",
                  &total.energy_topology_priority);
  metadata.Branch("energy_matching_scheme",
                  &total.energy_matching_scheme);
  metadata.Branch("energy_candidate_selection",
                  &total.energy_candidate_selection);
  metadata.Branch("signal_embedding_id", &total.signal_embedding_id);
  metadata.Branch("n_bins", &total.n_bins);
  metadata.Branch("pi0_truth_matching_algorithm_version",
                  &total.pi0_truth_matching_algorithm_version);
  metadata.Branch("et_max", &total.et_max);
  metadata.Branch("truth_eta_max", &total.truth_eta_max);
  metadata.Branch("cluster_eta_max", &total.cluster_eta_max);
  metadata.Branch("min_cluster_energy", &total.min_cluster_energy);
  metadata.Branch("dominant_fraction_min",
                  &total.dominant_fraction_min);
  metadata.Branch("pi0_contributor_fraction_min",
                  &total.pi0_contributor_fraction_min);
  metadata.Branch("min_energy_contribution_fraction",
                  &total.min_energy_contribution_fraction);
  metadata.Branch("separated_delta_r_cut",
                  &total.separated_delta_r_cut);
  metadata.Branch("merged_delta_r_cut", &total.merged_delta_r_cut);
  metadata.Branch("response_min", &total.response_min);
  metadata.Branch("response_max", &total.response_max);
  metadata.Branch("contains_raw_histograms",
                  &contains_raw_histograms);
  metadata.Branch("contains_bin_width_normalized_histograms",
                  &contains_bin_width_normalized_histograms);
  metadata.Branch("events_processed", &total.events_processed);
  metadata.Branch("events_written", &total.events_written);
  metadata.Branch("events_invalid", &total.events_invalid);
  metadata.Branch("prompt_cluster_count",
                  &total.prompt_cluster_count);
  metadata.Branch("pi0_cluster_count", &total.pi0_cluster_count);
  metadata.Branch("pi0_cluster_g4_decay_count",
                  &total.pi0_cluster_g4_decay_count);
  metadata.Branch("pi0_cluster_generator_decay_count",
                  &total.pi0_cluster_generator_decay_count);
  metadata.Branch("pi0_candidate_g4_decay_count",
                  &total.pi0_candidate_g4_decay_count);
  metadata.Branch("pi0_candidate_generator_decay_count",
                  &total.pi0_candidate_generator_decay_count);
  metadata.Branch("pi0_malformed_daughters_count",
                  &total.pi0_malformed_daughters_count);
  metadata.Branch("pi0_projection_failure_count",
                  &total.pi0_projection_failure_count);
  metadata.Branch("pi0_separated_count", &total.pi0_separated_count);
  metadata.Branch("pi0_merged_count", &total.pi0_merged_count);
  metadata.Branch("pi0_missing_count", &total.pi0_missing_count);
  metadata.Branch("pi0_none_count", &total.pi0_none_count);
  metadata.Branch("pi0_ambiguous_count", &total.pi0_ambiguous_count);
  metadata.Branch("pi0_separated_cluster_fill_count",
                  &total.pi0_separated_cluster_fill_count);
  metadata.Branch("pi0_merged_cluster_fill_count",
                  &total.pi0_merged_cluster_fill_count);
  metadata.Branch("pi0_missing_cluster_fill_count",
                  &total.pi0_missing_cluster_fill_count);
  metadata.Branch("pi0_energy_separated_count",
                  &total.pi0_energy_separated_count);
  metadata.Branch("pi0_energy_merged_count",
                  &total.pi0_energy_merged_count);
  metadata.Branch("pi0_energy_missing_count",
                  &total.pi0_energy_missing_count);
  metadata.Branch("pi0_energy_none_count",
                  &total.pi0_energy_none_count);
  metadata.Branch("pi0_energy_match_invalid_count",
                  &total.pi0_energy_match_invalid_count);
  metadata.Branch("pi0_energy_separated_cluster_fill_count",
                  &total.pi0_energy_separated_cluster_fill_count);
  metadata.Branch("pi0_energy_merged_cluster_fill_count",
                  &total.pi0_energy_merged_cluster_fill_count);
  metadata.Branch("pi0_energy_missing_cluster_fill_count",
                  &total.pi0_energy_missing_cluster_fill_count);
  metadata.Fill();
  metadata.Write();
  output.Close();
  if (output.TestBit(TFile::kWriteError))
  {
    return 6;
  }

  std::cout << "FinalizePythiaClusterEtSpectra - partials/files/events/prompt/pi0"
            << "/geometric-sep/merged/missing/energy-sep/merged/missing = "
            << partial_file_count << "/" << input_file_count << "/"
            << total.events_processed << "/"
            << total.prompt_cluster_count << "/"
            << total.pi0_cluster_count << "/"
            << total.pi0_separated_cluster_fill_count << "/"
            << total.pi0_merged_cluster_fill_count << "/"
            << total.pi0_missing_cluster_fill_count << "/"
            << total.pi0_energy_separated_cluster_fill_count << "/"
            << total.pi0_energy_merged_cluster_fill_count << "/"
            << total.pi0_energy_missing_cluster_fill_count << std::endl;
  return 0;
}
