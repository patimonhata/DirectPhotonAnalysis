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
constexpr std::size_t kHistogramCount = 6;
constexpr std::size_t kCategoryCount = 4;

struct PartialMetadata
{
  std::string path;
  int schema_version = 0;
  std::string manifest_path;
  std::string cluster_collection;
  std::string classification_unit;
  std::string pi0_selection;
  std::string partner_selection;
  std::string topology_definition;
  std::string topology_priority;
  std::string response_policy;
  long long manifest_begin = -1;
  long long manifest_end = -1;
  int signal_embedding_id = 0;
  int n_bins = 0;
  int matcher_version = 0;
  double et_max = 0.0;
  double truth_eta_max = 0.0;
  double anchor_cluster_eta_max = 0.0;
  double partner_cluster_eta_max = 0.0;
  double min_cluster_energy = 0.0;
  double dominant_fraction_min = 0.0;
  double anchor_pi0_fraction_min = 0.0;
  double min_energy_contribution_fraction = 0.0;
  unsigned char bin_width_normalized = 1U;
  unsigned long long events_processed = 0;
  unsigned long long events_written = 0;
  unsigned long long events_invalid = 0;
  unsigned long long cluster_considered = 0;
  unsigned long long cluster_invalid_truth = 0;
  unsigned long long prompt_count = 0;
  unsigned long long candidate_g4 = 0;
  unsigned long long candidate_generator = 0;
  unsigned long long malformed_daughters = 0;
  unsigned long long anchor_count = 0;
  unsigned long long anchor_g4 = 0;
  unsigned long long anchor_generator = 0;
  unsigned long long ambiguous_main = 0;
  unsigned long long energy_match_invalid = 0;
  unsigned long long separated_count = 0;
  unsigned long long merged_count = 0;
  unsigned long long missing_count = 0;
  unsigned long long other_count = 0;
};

template <class T>
bool bind(TTree* tree, const char* name, T* address)
{
  return tree->GetBranch(name) &&
      tree->SetBranchAddress(name, address) >= 0;
}

bool same_double(double left, double right)
{
  return std::abs(left - right) <=
      1e-12 * std::max({1.0, std::abs(left), std::abs(right)});
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
  std::string* classification_unit = nullptr;
  std::string* pi0_selection = nullptr;
  std::string* partner_selection = nullptr;
  std::string* topology_definition = nullptr;
  std::string* topology_priority = nullptr;
  std::string* response_policy = nullptr;
  bool ok = true;
  ok &= bind(tree, "schema_version", &value.schema_version);
  ok &= bind(tree, "manifest_path", &manifest_path);
  ok &= bind(tree, "manifest_begin", &value.manifest_begin);
  ok &= bind(tree, "manifest_end", &value.manifest_end);
  ok &= bind(tree, "cluster_collection", &cluster_collection);
  ok &= bind(tree, "classification_unit", &classification_unit);
  ok &= bind(tree, "pi0_selection", &pi0_selection);
  ok &= bind(tree, "partner_selection", &partner_selection);
  ok &= bind(tree, "topology_definition", &topology_definition);
  ok &= bind(tree, "topology_priority", &topology_priority);
  ok &= bind(tree, "response_policy", &response_policy);
  ok &= bind(tree, "signal_embedding_id", &value.signal_embedding_id);
  ok &= bind(tree, "n_bins", &value.n_bins);
  ok &= bind(tree, "et_max", &value.et_max);
  ok &= bind(tree, "truth_eta_max", &value.truth_eta_max);
  ok &= bind(tree, "anchor_cluster_eta_max",
             &value.anchor_cluster_eta_max);
  ok &= bind(tree, "partner_cluster_eta_max",
             &value.partner_cluster_eta_max);
  ok &= bind(tree, "min_cluster_energy", &value.min_cluster_energy);
  ok &= bind(tree, "dominant_fraction_min",
             &value.dominant_fraction_min);
  ok &= bind(tree, "anchor_pi0_fraction_min",
             &value.anchor_pi0_fraction_min);
  ok &= bind(tree, "min_energy_contribution_fraction",
             &value.min_energy_contribution_fraction);
  ok &= bind(tree, "pi0_truth_matching_algorithm_version",
             &value.matcher_version);
  ok &= bind(tree, "bin_width_normalized",
             &value.bin_width_normalized);
  ok &= bind(tree, "events_processed", &value.events_processed);
  ok &= bind(tree, "events_written", &value.events_written);
  ok &= bind(tree, "events_invalid", &value.events_invalid);
  ok &= bind(tree, "cluster_considered_count",
             &value.cluster_considered);
  ok &= bind(tree, "cluster_invalid_truth_count",
             &value.cluster_invalid_truth);
  ok &= bind(tree, "prompt_cluster_count", &value.prompt_count);
  ok &= bind(tree, "pi0_candidate_g4_decay_count",
             &value.candidate_g4);
  ok &= bind(tree, "pi0_candidate_generator_decay_count",
             &value.candidate_generator);
  ok &= bind(tree, "pi0_malformed_daughters_count",
             &value.malformed_daughters);
  ok &= bind(tree, "anchor_cluster_count", &value.anchor_count);
  ok &= bind(tree, "anchor_g4_decay_count", &value.anchor_g4);
  ok &= bind(tree, "anchor_generator_decay_count",
             &value.anchor_generator);
  ok &= bind(tree, "anchor_ambiguous_main_count",
             &value.ambiguous_main);
  ok &= bind(tree, "energy_match_invalid_count",
             &value.energy_match_invalid);
  ok &= bind(tree, "separated_count", &value.separated_count);
  ok &= bind(tree, "merged_count", &value.merged_count);
  ok &= bind(tree, "missing_count", &value.missing_count);
  ok &= bind(tree, "other_count", &value.other_count);
  if (!ok || tree->GetEntry(0) <= 0 || !manifest_path ||
      !cluster_collection || !classification_unit || !pi0_selection ||
      !partner_selection || !topology_definition || !topology_priority ||
      !response_policy)
  {
    return false;
  }

  value.path = path;
  value.manifest_path = *manifest_path;
  value.cluster_collection = *cluster_collection;
  value.classification_unit = *classification_unit;
  value.pi0_selection = *pi0_selection;
  value.partner_selection = *partner_selection;
  value.topology_definition = *topology_definition;
  value.topology_priority = *topology_priority;
  value.response_policy = *response_policy;
  return true;
}

bool valid_metadata(const PartialMetadata& value)
{
  return value.schema_version == 1 &&
      !value.manifest_path.empty() &&
      value.manifest_begin >= 0 &&
      value.manifest_end > value.manifest_begin &&
      value.cluster_collection == "split" &&
      value.classification_unit ==
          "every_cluster_with_selected_pi0_as_grouped_main_contributor" &&
      value.pi0_selection ==
          "signal_g4_primary_pi0_or_generator_pi0_with_exactly_two_g4_photons" &&
      value.partner_selection ==
          "same_energy_cut_as_anchor_partner_eta_cut_configurable" &&
      value.topology_definition ==
          "anchor_membership_in_direct_daughter_maximum_deposit_clusters" &&
      value.topology_priority ==
          "ambiguous_main_to_other_then_merged_then_separated_then_missing_then_other" &&
      value.response_policy == "not_used_for_classification" &&
      value.signal_embedding_id > 0 && value.n_bins > 0 &&
      value.et_max > 0.0 && value.truth_eta_max > 0.0 &&
      value.anchor_cluster_eta_max > 0.0 &&
      std::isfinite(value.partner_cluster_eta_max) &&
      value.min_cluster_energy >= 0.0 &&
      value.dominant_fraction_min >= 0.0 &&
      value.dominant_fraction_min <= 1.0 &&
      value.anchor_pi0_fraction_min >= 0.0 &&
      value.anchor_pi0_fraction_min <= 1.0 &&
      value.min_energy_contribution_fraction >= 0.0 &&
      value.min_energy_contribution_fraction < 1.0 &&
      value.matcher_version > 0 &&
      value.bin_width_normalized == 0U &&
      value.events_processed > 0 &&
      value.events_written + value.events_invalid ==
          value.events_processed &&
      value.cluster_invalid_truth <= value.cluster_considered &&
      value.anchor_count == value.anchor_g4 + value.anchor_generator &&
      value.anchor_count ==
          value.separated_count + value.merged_count +
          value.missing_count + value.other_count &&
      value.ambiguous_main <= value.other_count;
}

bool compatible(const PartialMetadata& value,
                const PartialMetadata& reference)
{
  return valid_metadata(value) &&
      value.schema_version == reference.schema_version &&
      value.manifest_path == reference.manifest_path &&
      value.cluster_collection == reference.cluster_collection &&
      value.classification_unit == reference.classification_unit &&
      value.pi0_selection == reference.pi0_selection &&
      value.partner_selection == reference.partner_selection &&
      value.topology_definition == reference.topology_definition &&
      value.topology_priority == reference.topology_priority &&
      value.response_policy == reference.response_policy &&
      value.signal_embedding_id == reference.signal_embedding_id &&
      value.n_bins == reference.n_bins &&
      value.matcher_version == reference.matcher_version &&
      same_double(value.et_max, reference.et_max) &&
      same_double(value.truth_eta_max, reference.truth_eta_max) &&
      same_double(value.anchor_cluster_eta_max,
                  reference.anchor_cluster_eta_max) &&
      same_double(value.partner_cluster_eta_max,
                  reference.partner_cluster_eta_max) &&
      same_double(value.min_cluster_energy,
                  reference.min_cluster_energy) &&
      same_double(value.dominant_fraction_min,
                  reference.dominant_fraction_min) &&
      same_double(value.anchor_pi0_fraction_min,
                  reference.anchor_pi0_fraction_min) &&
      same_double(value.min_energy_contribution_fraction,
                  reference.min_energy_contribution_fraction);
}

bool valid_histogram(const TH1D* histogram,
                     const PartialMetadata& metadata,
                     unsigned long long expected_entries)
{
  if (!histogram || histogram->GetNbinsX() != metadata.n_bins ||
      std::abs(histogram->GetXaxis()->GetXmin()) > 1e-12 ||
      !same_double(histogram->GetXaxis()->GetXmax(), metadata.et_max) ||
      histogram->GetSumw2N() == 0 ||
      std::abs(histogram->GetEntries() -
               static_cast<double>(expected_entries)) > 0.5)
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
  return directory.empty() ||
      !gSystem->AccessPathName(directory.c_str()) ||
      gSystem->mkdir(directory.c_str(), true) == 0;
}

double smallest_positive(
    const std::array<std::unique_ptr<TH1D>, kHistogramCount>& histograms)
{
  double result = std::numeric_limits<double>::infinity();
  for (const auto& histogram : histograms)
  {
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

int FinalizePythiaPi0AnchorClusterSpectra(
    const std::string partial_pattern =
        "output/pi0_anchor_topology_partial/eta07_full_partner_fgamma0p3_jet5/partial_*.root",
    const std::string output_base =
        "output/plots/pi0_anchor_topology/jet5/eta07_full_partner_fgamma0p3",
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
    std::cerr
        << "FinalizePythiaPi0AnchorClusterSpectra - no partials matched"
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
      std::cerr
          << "FinalizePythiaPi0AnchorClusterSpectra - invalid partial: "
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
          << "FinalizePythiaPi0AnchorClusterSpectra - incompatible or noncontiguous partial: "
          << partial.path << ", expected begin " << next_begin
          << std::endl;
      return 4;
    }
    next_begin = partial.manifest_end;
  }
  if (expected_manifest_end >= 0 && next_begin != expected_manifest_end)
  {
    return 4;
  }

  const std::array<std::string, kHistogramCount> raw_names = {
      "h_prompt_cluster_et_raw", "h_pi0_anchor_cluster_et_raw",
      "h_pi0_anchor_separated_cluster_et_raw",
      "h_pi0_anchor_merged_cluster_et_raw",
      "h_pi0_anchor_missing_cluster_et_raw",
      "h_pi0_anchor_other_cluster_et_raw"};
  const std::array<std::string, kHistogramCount> density_names = {
      "h_prompt_cluster_et_density",
      "h_pi0_anchor_cluster_et_density",
      "h_pi0_anchor_separated_cluster_et_density",
      "h_pi0_anchor_merged_cluster_et_density",
      "h_pi0_anchor_missing_cluster_et_density",
      "h_pi0_anchor_other_cluster_et_density"};

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
  total.cluster_considered = total.cluster_invalid_truth = 0;
  total.prompt_count = total.candidate_g4 = total.candidate_generator = 0;
  total.malformed_daughters = total.anchor_count = 0;
  total.anchor_g4 = total.anchor_generator = total.ambiguous_main = 0;
  total.energy_match_invalid = 0;
  total.separated_count = total.merged_count = 0;
  total.missing_count = total.other_count = 0;

  for (const PartialMetadata& partial : partials)
  {
    TFile input(partial.path.c_str(), "READ");
    const std::array<unsigned long long, kHistogramCount> counts = {
        partial.prompt_count, partial.anchor_count,
        partial.separated_count, partial.merged_count,
        partial.missing_count, partial.other_count};
    std::array<TH1D*, kHistogramCount> partial_histograms{};
    for (std::size_t index = 0; index < raw.size(); ++index)
    {
      input.GetObject(raw_names[index].c_str(),
                      partial_histograms[index]);
      if (input.IsZombie() ||
          !valid_histogram(
              partial_histograms[index], partial, counts[index]) ||
          !raw[index]->Add(partial_histograms[index]))
      {
        std::cerr
            << "FinalizePythiaPi0AnchorClusterSpectra - invalid histogram in "
            << partial.path << std::endl;
        return 5;
      }
    }
    for (int bin = 0; bin <= partial.n_bins + 1; ++bin)
    {
      const double categories =
          partial_histograms[2]->GetBinContent(bin) +
          partial_histograms[3]->GetBinContent(bin) +
          partial_histograms[4]->GetBinContent(bin) +
          partial_histograms[5]->GetBinContent(bin);
      if (std::abs(
              partial_histograms[1]->GetBinContent(bin) - categories) >
          1e-9)
      {
        return 5;
      }
    }

    total.events_processed += partial.events_processed;
    total.events_written += partial.events_written;
    total.events_invalid += partial.events_invalid;
    total.cluster_considered += partial.cluster_considered;
    total.cluster_invalid_truth += partial.cluster_invalid_truth;
    total.prompt_count += partial.prompt_count;
    total.candidate_g4 += partial.candidate_g4;
    total.candidate_generator += partial.candidate_generator;
    total.malformed_daughters += partial.malformed_daughters;
    total.anchor_count += partial.anchor_count;
    total.anchor_g4 += partial.anchor_g4;
    total.anchor_generator += partial.anchor_generator;
    total.ambiguous_main += partial.ambiguous_main;
    total.energy_match_invalid += partial.energy_match_invalid;
    total.separated_count += partial.separated_count;
    total.merged_count += partial.merged_count;
    total.missing_count += partial.missing_count;
    total.other_count += partial.other_count;
  }

  for (int bin = 0; bin <= reference.n_bins + 1; ++bin)
  {
    const double categories = raw[2]->GetBinContent(bin) +
        raw[3]->GetBinContent(bin) + raw[4]->GetBinContent(bin) +
        raw[5]->GetBinContent(bin);
    if (std::abs(raw[1]->GetBinContent(bin) - categories) > 1e-9)
    {
      return 5;
    }
  }

  std::array<std::unique_ptr<TH1D>, kHistogramCount> density;
  const std::array<int, kHistogramCount> colors = {
      kRed + 1, kBlue + 1, kAzure + 7, kMagenta + 1,
      kGreen + 2, kGray + 2};
  for (std::size_t index = 0; index < raw.size(); ++index)
  {
    density[index].reset(static_cast<TH1D*>(
        raw[index]->Clone(density_names[index].c_str())));
    density[index]->SetDirectory(nullptr);
    density[index]->Scale(1.0, "width");
    density[index]->SetStats(false);
    density[index]->SetLineColor(colors[index]);
    density[index]->SetLineWidth(index < 2 ? 3 : 2);
    density[index]->GetXaxis()->SetTitle("Cluster E_{T} [GeV]");
    density[index]->GetYaxis()->SetTitle("Clusters / GeV");

    raw[index]->SetStats(false);
    raw[index]->SetLineColor(colors[index]);
    raw[index]->SetLineWidth(index < 2 ? 3 : 2);
    raw[index]->GetXaxis()->SetTitle("Cluster E_{T} [GeV]");
    raw[index]->GetYaxis()->SetTitle("Counts / bin");
  }

  const std::array<std::string, kCategoryCount> fraction_names = {
      "h_pi0_anchor_separated_fraction",
      "h_pi0_anchor_merged_fraction",
      "h_pi0_anchor_missing_fraction",
      "h_pi0_anchor_other_fraction"};
  std::array<std::unique_ptr<TH1D>, kCategoryCount> fractions;
  for (std::size_t index = 0; index < fractions.size(); ++index)
  {
    fractions[index].reset(static_cast<TH1D*>(
        raw[index + 2]->Clone(fraction_names[index].c_str())));
    fractions[index]->SetDirectory(nullptr);
    fractions[index]->Divide(
        raw[index + 2].get(), raw[1].get(), 1.0, 1.0, "B");
    fractions[index]->SetStats(false);
    fractions[index]->SetLineColor(colors[index + 2]);
    fractions[index]->SetMarkerColor(colors[index + 2]);
    fractions[index]->SetMarkerStyle(20 + static_cast<int>(index));
    fractions[index]->SetLineWidth(2);
    fractions[index]->GetXaxis()->SetTitle("Anchor cluster E_{T} [GeV]");
    fractions[index]->GetYaxis()->SetTitle("Category / all anchors");
  }

  if (!make_output_directory(output_base))
  {
    return 6;
  }
  SetsPhenixStyle();

  TCanvas spectrum_canvas(
      "c_pythia_pi0_anchor_cluster_et",
      "Pythia pi0 anchor cluster spectra", 1000, 800);
  spectrum_canvas.SetLogy();
  double maximum = 0.0;
  for (const auto& histogram : raw)
  {
    maximum = std::max(maximum, histogram->GetMaximum());
  }
  const double minimum = smallest_positive(raw);
  raw[0]->SetMinimum(minimum > 0.0 ? 0.5 * minimum : 0.5);
  raw[0]->SetMaximum(maximum > 0.0 ? 5.0 * maximum : 1.0);
  raw[0]->Draw("HIST");
  for (std::size_t index = 1; index < raw.size(); ++index)
  {
    raw[index]->Draw("HIST SAME");
  }
  TLegend spectrum_legend(0.49, 0.36, 0.89, 0.68);
  spectrum_legend.SetBorderSize(0);
  spectrum_legend.AddEntry(raw[0].get(), "Prompt-#gamma cluster", "l");
  spectrum_legend.AddEntry(raw[1].get(), "#pi^{0}-main anchor", "l");
  spectrum_legend.AddEntry(raw[2].get(), "Separated", "l");
  spectrum_legend.AddEntry(raw[3].get(), "Merged", "l");
  spectrum_legend.AddEntry(raw[4].get(), "Missing partner", "l");
  spectrum_legend.AddEntry(raw[5].get(), "Other", "l");
  spectrum_legend.Draw();

  TLatex spectrum_label;
  spectrum_label.SetNDC();
  spectrum_label.SetTextAlign(13);
  spectrum_label.DrawLatex(0.20, 0.92, "#it{#bf{sPHENIX}} Internal");
  spectrum_label.DrawLatex(0.20, 0.84, "Pythia8 p+p minimum bias");
  std::ostringstream anchor_label;
  anchor_label << "Anchor |#eta| < "
               << reference.anchor_cluster_eta_max;
  spectrum_label.DrawLatex(
      0.20, 0.76, anchor_label.str().c_str());
  const std::string partner_label =
      reference.partner_cluster_eta_max <= 0.0
          ? "Partner: full CEMC cluster container"
          : "Partner: software #eta cut applied";
  spectrum_label.DrawLatex(0.20, 0.68, partner_label.c_str());
  spectrum_label.DrawLatex(
      0.20, 0.60, "Topology: direct-daughter energy deposit");
  spectrum_canvas.RedrawAxis();
  spectrum_canvas.SaveAs((output_base + ".pdf").c_str());

  TCanvas fraction_canvas(
      "c_pythia_pi0_anchor_category_fractions",
      "Pythia pi0 anchor category fractions", 1000, 800);
  fractions[0]->SetMinimum(0.0);
  fractions[0]->SetMaximum(1.05);
  fractions[0]->Draw("E1");
  for (std::size_t index = 1; index < fractions.size(); ++index)
  {
    fractions[index]->Draw("E1 SAME");
  }
  TLegend fraction_legend(0.58, 0.62, 0.89, 0.86);
  fraction_legend.SetBorderSize(0);
  fraction_legend.AddEntry(fractions[0].get(), "Separated", "lep");
  fraction_legend.AddEntry(fractions[1].get(), "Merged", "lep");
  fraction_legend.AddEntry(fractions[2].get(), "Missing partner", "lep");
  fraction_legend.AddEntry(fractions[3].get(), "Other", "lep");
  fraction_legend.Draw();
  TLatex fraction_label;
  fraction_label.SetNDC();
  fraction_label.SetTextAlign(13);
  fraction_label.DrawLatex(0.20, 0.92, "#it{#bf{sPHENIX}} Internal");
  fraction_label.DrawLatex(0.20, 0.84, "Pythia8 p+p minimum bias");
  fraction_label.DrawLatex(
      0.20, 0.76, "Denominator: all #pi^{0}-main anchors");
  fraction_canvas.RedrawAxis();
  fraction_canvas.SaveAs(
      (output_base + "_category_fractions.pdf").c_str());

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
  for (auto& histogram : fractions)
  {
    histogram->Write();
  }

  int output_schema_version = 1;
  long long manifest_begin = partials.front().manifest_begin;
  long long manifest_end = partials.back().manifest_end;
  long long partial_file_count =
      static_cast<long long>(partials.size());
  long long input_file_count = manifest_end - manifest_begin;
  unsigned char contains_raw_histograms = 1U;
  unsigned char contains_bin_width_normalized_histograms = 1U;
  unsigned char contains_category_fractions = 1U;
  TTree metadata(
      "metadata", "Final Pythia pi0 anchor-cluster metadata");
  metadata.Branch("schema_version", &output_schema_version);
  metadata.Branch("manifest_path", &total.manifest_path);
  metadata.Branch("manifest_begin", &manifest_begin);
  metadata.Branch("manifest_end", &manifest_end);
  metadata.Branch("partial_file_count", &partial_file_count);
  metadata.Branch("input_file_count", &input_file_count);
  metadata.Branch("cluster_collection", &total.cluster_collection);
  metadata.Branch("classification_unit", &total.classification_unit);
  metadata.Branch("pi0_selection", &total.pi0_selection);
  metadata.Branch("partner_selection", &total.partner_selection);
  metadata.Branch("topology_definition", &total.topology_definition);
  metadata.Branch("topology_priority", &total.topology_priority);
  metadata.Branch("response_policy", &total.response_policy);
  metadata.Branch("signal_embedding_id", &total.signal_embedding_id);
  metadata.Branch("n_bins", &total.n_bins);
  metadata.Branch("et_max", &total.et_max);
  metadata.Branch("truth_eta_max", &total.truth_eta_max);
  metadata.Branch(
      "anchor_cluster_eta_max", &total.anchor_cluster_eta_max);
  metadata.Branch(
      "partner_cluster_eta_max", &total.partner_cluster_eta_max);
  metadata.Branch("min_cluster_energy", &total.min_cluster_energy);
  metadata.Branch(
      "dominant_fraction_min", &total.dominant_fraction_min);
  metadata.Branch(
      "anchor_pi0_fraction_min", &total.anchor_pi0_fraction_min);
  metadata.Branch(
      "min_energy_contribution_fraction",
      &total.min_energy_contribution_fraction);
  metadata.Branch(
      "pi0_truth_matching_algorithm_version", &total.matcher_version);
  metadata.Branch(
      "contains_raw_histograms", &contains_raw_histograms);
  metadata.Branch(
      "contains_bin_width_normalized_histograms",
      &contains_bin_width_normalized_histograms);
  metadata.Branch(
      "contains_category_fractions", &contains_category_fractions);
  metadata.Branch("events_processed", &total.events_processed);
  metadata.Branch("events_written", &total.events_written);
  metadata.Branch("events_invalid", &total.events_invalid);
  metadata.Branch(
      "cluster_considered_count", &total.cluster_considered);
  metadata.Branch(
      "cluster_invalid_truth_count", &total.cluster_invalid_truth);
  metadata.Branch("prompt_cluster_count", &total.prompt_count);
  metadata.Branch(
      "pi0_candidate_g4_decay_count", &total.candidate_g4);
  metadata.Branch(
      "pi0_candidate_generator_decay_count",
      &total.candidate_generator);
  metadata.Branch(
      "pi0_malformed_daughters_count", &total.malformed_daughters);
  metadata.Branch("anchor_cluster_count", &total.anchor_count);
  metadata.Branch("anchor_g4_decay_count", &total.anchor_g4);
  metadata.Branch(
      "anchor_generator_decay_count", &total.anchor_generator);
  metadata.Branch(
      "anchor_ambiguous_main_count", &total.ambiguous_main);
  metadata.Branch(
      "energy_match_invalid_count", &total.energy_match_invalid);
  metadata.Branch("separated_count", &total.separated_count);
  metadata.Branch("merged_count", &total.merged_count);
  metadata.Branch("missing_count", &total.missing_count);
  metadata.Branch("other_count", &total.other_count);
  metadata.Fill();
  metadata.Write();
  output.Close();
  if (output.TestBit(TFile::kWriteError))
  {
    return 6;
  }

  std::cout
      << "FinalizePythiaPi0AnchorClusterSpectra - partials/files/events"
      << "/anchor/separated/merged/missing/other = "
      << partial_file_count << "/" << input_file_count << "/"
      << total.events_processed << "/" << total.anchor_count << "/"
      << total.separated_count << "/" << total.merged_count << "/"
      << total.missing_count << "/" << total.other_count << std::endl;
  return 0;
}
