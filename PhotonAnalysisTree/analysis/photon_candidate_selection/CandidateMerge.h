#ifndef PHOTON_CANDIDATE_CANDIDATEMERGE_H
#define PHOTON_CANDIDATE_CANDIDATEMERGE_H

#include "CandidateSelection.h"


#include <TGraphErrors.h>

namespace candidate_composition_merge
{
// Kinematic, preselection, tight, isolation, region A, region A + tagging veto.
constexpr std::array<std::size_t, kSelectionCount> kPreviousCutSelections = {0, 0, 1, 1, 1, 4};
constexpr std::array<const char*, kSelectionCount> kPreviousCutLabels = {
    "Kinematic", "Kinematic", "Pre-selection", "Pre-selection", "Pre-selection", "Region A"};
struct Metadata
{
  int common_composition_algorithm_version = -1;
  int schema_version = -1;
  int source_map_schema_version = -1;
  std::string family;
  std::string map_root;
  bool require_complete = false;
  int n_bins = -1;
  double et_max = -1.0;
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
  int other_origin_algorithm_version = -1;
  int decay_origin_algorithm_version = -1;
  double decay_origin_max_proper_length_cm = -1, decay_origin_max_flight_cm = -1;
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

constexpr std::size_t kCandidateOriginCount = 4;
constexpr std::array<const char*, kCandidateOriginCount> kCandidateOriginKeys = {"prompt", "pi0", "eta", "other"};
constexpr std::array<const char*, kCandidateOriginCount> kCandidateOriginLabels = {
    "Prompt-#gamma clusters", "#pi^{0}-origin clusters", "#eta-origin clusters", "Other"};
constexpr std::array<int, kCandidateOriginCount> kCandidateOriginColors = {
    candidate_composition::kColors[candidate_composition::prompt], candidate_composition::kColors[candidate_composition::pi0_separated],
    candidate_composition::kColors[candidate_composition::eta], candidate_composition::kColors[candidate_composition::other]};

std::array<std::unique_ptr<TH1D>, kCandidateOriginCount> make_candidate_origin_density(
    const candidate_composition::Histograms& histograms, const std::string& prefix)
{
  using namespace candidate_composition;
  std::array<std::unique_ptr<TH1D>, kCandidateOriginCount> density;
  const std::array<std::size_t, kCandidateOriginCount> source_indices = {prompt, pi0_separated, eta, other};
  for (std::size_t index = 0; index < kCandidateOriginCount; ++index)
  {
    density[index].reset(static_cast<TH1D*>(histograms.weighted[source_indices[index]]->Clone((prefix + kCandidateOriginKeys[index] + "_et_pb_per_gev").c_str())));
    density[index]->SetDirectory(nullptr);
  }
  density[1]->Reset("ICES");
  for (std::size_t index = pi0_separated; index <= pi0_unknown; ++index) density[1]->Add(histograms.weighted[index].get());
  for (std::size_t index = 0; index < kCandidateOriginCount; ++index)
  {
    density[index]->Scale(1.0, "width");
    density[index]->SetLineColor(kCandidateOriginColors[index]);
    style_axes(density[index].get(), "Weighted Counts [a.u.]");
  }
  return density;
}

void draw_candidate_origin_spectrum(const std::array<std::unique_ptr<TH1D>, kCandidateOriginCount>& density,
                                    const std::string& output_path, const std::string& selection_key)
{
  SetsPhenixStyle();
  TCanvas canvas(("c_" + selection_key + "_candidate_origin_spectrum").c_str(), "", kCanvasWidth, kCanvasHeight);
  auto plot_pad = make_plot_pad(selection_key + "_candidate_origin_spectrum_pad", true, kPlotPadHeight);
  auto frame = make_plot_frame(*density.front(), "h_" + selection_key + "_candidate_origin_frame", "Weighted Counts [a.u.]",
                                           kSpectrumYMinimum, kSpectrumYMaximum);
  scale_axes_to_canvas(frame.get(), kPlotPadHeight);
  frame->GetYaxis()->SetTitleOffset(1.0);
  frame->Draw("AXIS");
  for (const auto& histogram : density) histogram->Draw("HIST SAME");
  TLegend legend(0.42, 0.58, 0.94, 0.94);
  scale_legend_to_canvas(legend, kPlotPadHeight);
  legend.SetFillStyle(0);
  for (std::size_t index = 0; index < kCandidateOriginCount; ++index) legend.AddEntry(density[index].get(), kCandidateOriginLabels[index], "l");
  legend.Draw();
  plot_pad->RedrawAxis();
  canvas.cd();
  draw_annotations();
  canvas.SaveAs(output_path.c_str());
}

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
  ok &= bind(tree, "common_composition_algorithm_version", &value.common_composition_algorithm_version);
  ok &= bind(tree, "source_map_schema_version", &value.source_map_schema_version);
  ok &= bind(tree, "family", &family);
  ok &= bind(tree, "map_root", &map_root);
  ok &= bind(tree, "require_complete", &value.require_complete);
  ok &= bind(tree, "n_bins", &value.n_bins);
  ok &= bind(tree, "et_max", &value.et_max);
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
  ok &= bind(tree, "other_origin_algorithm_version", &value.other_origin_algorithm_version);
  ok &= bind(tree, "decay_origin_algorithm_version", &value.decay_origin_algorithm_version);
  ok &= bind(tree, "decay_origin_max_proper_length_cm", &value.decay_origin_max_proper_length_cm);
  ok &= bind(tree, "decay_origin_max_flight_cm", &value.decay_origin_max_flight_cm);
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
  if (!ok || tree->GetEntry(0) <= 0 || value.common_composition_algorithm_version != 1 || value.decay_origin_algorithm_version != decay_origin::kAlgorithmVersion ||
      !same_double(value.decay_origin_max_proper_length_cm, decay_origin::kMaxProperLengthCm) ||
      !same_double(value.decay_origin_max_flight_cm, decay_origin::kMaxFlightCm) || value.pi0_topology_algorithm_version != 14 || value.eta_topology_algorithm_version != 4 || value.other_origin_algorithm_version != other_origin::kAlgorithmVersion || !family || !map_root || !majority_comparison || !eta_definition ||
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
  return value.common_composition_algorithm_version == reference.common_composition_algorithm_version && value.schema_version == reference.schema_version && value.source_map_schema_version == reference.source_map_schema_version &&
      value.family == reference.family && value.map_root == reference.map_root &&
      value.require_complete == reference.require_complete && value.n_bins == reference.n_bins && same_double(value.et_max, reference.et_max) &&
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
      same_double(value.min_photon_energy_recovery, reference.min_photon_energy_recovery) &&
      value.pi0_topology_algorithm_version == reference.pi0_topology_algorithm_version &&
      value.eta_topology_algorithm_version == reference.eta_topology_algorithm_version &&
      value.other_origin_algorithm_version == reference.other_origin_algorithm_version && value.decay_origin_algorithm_version == reference.decay_origin_algorithm_version &&
      same_double(value.decay_origin_max_proper_length_cm, reference.decay_origin_max_proper_length_cm) &&
      same_double(value.decay_origin_max_flight_cm, reference.decay_origin_max_flight_cm) &&
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

struct SurvivalCurve
{
  std::unique_ptr<TH1D> histogram;
  std::unique_ptr<TGraphErrors> graph;
  std::string label;
};

bool valid_subset(const TH1D& numerator_counts, const TH1D& denominator_counts,
                  const TH1D& numerator_weighted, const TH1D& denominator_weighted)
{
  for (int bin = 0; bin <= numerator_counts.GetNbinsX() + 1; ++bin)
  {
    const double numerator_count = numerator_counts.GetBinContent(bin);
    const double denominator_count = denominator_counts.GetBinContent(bin);
    const double count_scale = std::max({1.0, std::abs(numerator_count), std::abs(denominator_count)});
    const double numerator_w2 = std::pow(numerator_weighted.GetBinError(bin), 2);
    const double denominator_w2 = std::pow(denominator_weighted.GetBinError(bin), 2);
    const double w2_scale = std::max({1.0, numerator_w2, denominator_w2});
    if (numerator_count < -1e-9 * count_scale || numerator_count > denominator_count + 1e-9 * count_scale ||
        numerator_w2 > denominator_w2 + 1e-9 * w2_scale) return false;
  }
  return true;
}

std::unique_ptr<TH1D> sum_histograms(const std::vector<const TH1D*>& inputs, const std::string& name)
{
  if (inputs.empty()) return nullptr;
  auto result = std::unique_ptr<TH1D>(static_cast<TH1D*>(inputs.front()->Clone(name.c_str())));
  result->SetDirectory(nullptr);
  result->Reset("ICES");
  for (const TH1D* input : inputs)
  {
    if (!input || !result->Add(input)) return nullptr;
  }
  return result;
}

std::unique_ptr<SurvivalCurve> make_survival_curve(const TH1D& numerator, const TH1D& denominator,
                                                   const std::string& name, const std::string& label,
                                                   int color, const std::string& denominator_label = "Kinematic")
{
  auto curve = std::make_unique<SurvivalCurve>();
  curve->histogram = candidate_composition::fraction_histogram(numerator, denominator, name);
  curve->histogram->SetLineColor(color);
  curve->histogram->SetMarkerColor(color);
  style_axes(curve->histogram.get(), ("#splitline{Survival fraction}{relative to " + denominator_label + "}").c_str());
  curve->graph = std::make_unique<TGraphErrors>();
  curve->graph->SetName((name + "_valid_bins").c_str());
  curve->graph->SetLineColor(color);
  curve->graph->SetMarkerColor(color);
  curve->label = label;
  for (int bin = 1; bin <= denominator.GetNbinsX(); ++bin)
  {
    if (denominator.GetBinContent(bin) == 0.0) continue;
    const double value = curve->histogram->GetBinContent(bin);
    const double error = curve->histogram->GetBinError(bin);
    if (!std::isfinite(value) || !std::isfinite(error)) continue;
    const int point = curve->graph->GetN();
    curve->graph->SetPoint(point, denominator.GetXaxis()->GetBinCenter(bin), value);
    curve->graph->SetPointError(point, 0.0, error);
  }
  return curve;
}

void draw_survival_curves(const std::vector<const SurvivalCurve*>& curves, const TH1D& axis_source,
                          const std::string& output_path, const std::string& canvas_name, bool detailed, const std::string& denominator_label = "Kinematic")
{
  SetsPhenixStyle();
  const bool side_legend = detailed || curves.size() > 4;
  TCanvas canvas(canvas_name.c_str(), "", side_legend ? kDetailedCanvasWidth : kCanvasWidth, kCanvasHeight);
  constexpr double plot_height = kPlotPadHeight;
  auto plot_pad = make_plot_pad(canvas_name + "_pad", false, plot_height, side_legend ? kDetailedPlotWidth : 1.0);
  auto frame = make_plot_frame(axis_source, canvas_name + "_frame", ("#splitline{Survival fraction}{relative to " + denominator_label + "}").c_str(), 0.0, 1.05);
  scale_axes_to_canvas(frame.get(), plot_height);
  if (side_legend) frame->GetXaxis()->SetTitleOffset(1.0);
  frame->GetYaxis()->SetTitleOffset(1.0);
  frame->Draw("AXIS");
  for (const SurvivalCurve* curve : curves) curve->graph->Draw("PE1 SAME");
  std::vector<std::unique_ptr<TLegend>> legends;
  if (side_legend) canvas.cd();
  legends.push_back(side_legend
      ? std::make_unique<TLegend>(kDetailedLegendLeft, detailed ? 0.05 : 0.25, 0.99, detailed ? 0.95 : 0.75)
      : std::make_unique<TLegend>(0.60, 0.68, 0.955, 0.94));
  for (const SurvivalCurve* curve : curves)
  {
    if (side_legend) add_readable_legend_entry(*legends[0], curve->graph.get(), curve->label, "lep");
    else legends[0]->AddEntry(curve->graph.get(), curve->label.c_str(), "lep");
  }
  for (auto& legend : legends)
  {
    if (!side_legend) scale_legend_to_canvas(*legend, plot_height);
    legend->SetFillStyle(0);
    legend->Draw();
  }
  plot_pad->RedrawAxis();
  canvas.cd();
  draw_annotations();
  const bool embed_fonts = output_path.size() >= 4 && output_path.compare(output_path.size() - 4, 4, ".pdf") == 0 &&
      std::any_of(curves.begin(), curves.end(), [](const SurvivalCurve* curve) { return curve->label.find("#mu") != std::string::npos; });
  if (embed_fonts) canvas.Print(output_path.c_str(), "EmbedFonts");
  else canvas.SaveAs(output_path.c_str());
}
}

#endif
