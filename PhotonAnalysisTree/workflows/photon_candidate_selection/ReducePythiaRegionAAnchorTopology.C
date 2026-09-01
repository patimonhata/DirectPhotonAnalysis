#include "../../macro/Utilities/sPhenixStyle.C"

#include <TCanvas.h>
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
#include <limits>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace
{
constexpr std::size_t kSpectrumCount = 14;
constexpr std::array<const char*, kSpectrumCount> kKeys = {
    "prompt", "pi0_anchor", "separated", "merged", "single_contaminated", "missing",
    "missing_energy_threshold", "missing_displaced_partner", "missing_acceptance",
    "missing_no_cemc_deposit", "missing_unclustered_deposit", "missing_match_incomplete",
    "missing_other", "other"};
constexpr std::array<const char*, kSpectrumCount> kLabels = {
    "Prompt-#gamma cluster", "#pi^{0}-main anchor", "Separated", "Merged", "Single contaminated", "Missing (total)",
    "Missing: energy threshold", "Missing: displaced partner cluster", "Missing: acceptance",
    "Missing: no CEMC deposit", "Missing: unclustered deposit", "Missing: match incomplete",
    "Missing: other", "Other"};
constexpr std::array<int, kSpectrumCount> kColors = {
    kRed + 1, kBlue + 1, kAzure + 7, kMagenta + 1, kCyan + 2, kGreen + 2,
    kOrange + 7, kPink + 7, kViolet + 1, kYellow + 2, kSpring + 5, kBlue + 3, kGreen + 3, kGray + 2};
constexpr std::array<std::size_t, 7> kSummarySpectrum = {0, 1, 2, 3, 4, 5, 13};
constexpr std::array<std::size_t, 11> kDetailedCategories = {2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13};
constexpr std::array<std::size_t, 5> kSummaryCategories = {2, 3, 4, 5, 13};
constexpr int kCanvasWidth = 1100;
constexpr int kCanvasHeight = 900;

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
        {"jet3", 0.0, 5.0, 1.2147e9, false, 10001},
        {"jet5", 5.0, 9.0, 1.3878e8, false, 10001},
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

struct MapMetadata
{
  std::string path;
  int schema_version = -1;
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
  unsigned long long events_written = 0;
  unsigned long long events_vertex_rejected = 0;
  unsigned long long events_invalid = 0;
};

template <class T>
bool bind(TTree* tree, const char* name, T* address)
{
  return tree->GetBranch(name) && tree->SetBranchAddress(name, address) >= 0;
}

bool same_double(double left, double right)
{
  return std::abs(left - right) <= 1e-11 * std::max({1.0, std::abs(left), std::abs(right)});
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
  ok &= bind(tree, "schema_version", &value.schema_version);
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
  return true;
}

bool valid_metadata(const MapMetadata& value, const SampleDefinition& sample)
{
  return value.schema_version == 3 && value.sample_name == sample.name && !value.input_manifest.empty() &&
      !value.analysis_release.empty() && !value.model_sha256.empty() && value.manifest_begin >= 0 &&
      value.manifest_end > value.manifest_begin && value.input_file_count == value.manifest_end - value.manifest_begin &&
      same_double(value.cross_section_pb, sample.cross_section_pb) && same_double(value.window_min, sample.window_min) &&
      same_double(value.window_max, sample.window_max) && value.upper_unbounded == sample.upper_unbounded &&
      std::isfinite(value.sum_generator_weight_processed) && value.sum_generator_weight_processed > 0.0 &&
      value.events_processed == value.events_written + value.events_vertex_rejected + value.events_invalid;
}

bool compatible(const MapMetadata& value, const MapMetadata& reference)
{
  return value.schema_version == reference.schema_version && value.input_manifest == reference.input_manifest &&
      value.sample_name == reference.sample_name && value.analysis_release == reference.analysis_release &&
      value.model_sha256 == reference.model_sha256 && same_double(value.cross_section_pb, reference.cross_section_pb) &&
      same_double(value.window_min, reference.window_min) && same_double(value.window_max, reference.window_max) &&
      value.upper_unbounded == reference.upper_unbounded;
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

bool valid_partition(const std::array<std::unique_ptr<TH1D>, kSpectrumCount>& histograms)
{
  for (int bin = 0; bin <= histograms[0]->GetNbinsX() + 1; ++bin)
  {
    double missing_sum = 0.0;
    for (std::size_t index = 6; index <= 12; ++index) missing_sum += histograms[index]->GetBinContent(bin);
    const double category_sum = histograms[2]->GetBinContent(bin) + histograms[3]->GetBinContent(bin) +
        histograms[4]->GetBinContent(bin) + histograms[5]->GetBinContent(bin) + histograms[13]->GetBinContent(bin);
    const double scale = std::max({1.0, std::abs(histograms[1]->GetBinContent(bin)), std::abs(category_sum), std::abs(missing_sum)});
    if (std::abs(histograms[5]->GetBinContent(bin) - missing_sum) > 1e-9 * scale ||
        std::abs(histograms[1]->GetBinContent(bin) - category_sum) > 1e-9 * scale) return false;
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

std::unique_ptr<TPad> make_plot_pad(const std::string& name, bool log_y = false)
{
  auto pad = std::make_unique<TPad>(name.c_str(), "", 0.0, 0.0, 1.0, 0.62);
  pad->SetLeftMargin(0.13);
  pad->SetRightMargin(0.04);
  pad->SetBottomMargin(0.16);
  pad->SetTopMargin(0.04);
  pad->SetTicks(1, 1);
  pad->SetLogy(log_y);
  pad->Draw();
  pad->cd();
  return pad;
}

void style_axes(TH1* histogram, const char* y_title)
{
  histogram->SetStats(false);
  histogram->GetXaxis()->SetTitle("Cluster E_{T} [GeV]");
  histogram->GetYaxis()->SetTitle(y_title);
  for (TAxis* axis : {histogram->GetXaxis(), histogram->GetYaxis()})
  {
    axis->SetLabelSize(0.05);
    axis->SetTitleSize(0.05);
    axis->CenterTitle();
  }
  histogram->GetXaxis()->SetTitleOffset(1.4);
  histogram->GetYaxis()->SetTitleOffset(1.15);
}

void draw_annotations(const std::string& family_label)
{
  TLatex label;
  label.SetNDC();
  label.SetTextAlign(13);
  label.SetTextSize(0.026);
  label.DrawLatex(0.06, 0.96, "#it{#bf{sPHENIX}} Internal");
  label.DrawLatex(0.06, 0.91, family_label.c_str());
  label.DrawLatex(0.06, 0.86, "Region A: isolated and tight");
  label.DrawLatex(0.06, 0.81, "5 < E_{T} < 35 GeV, |#eta| < 0.7");
  label.DrawLatex(0.06, 0.76, "Stored topology: E_{cluster} > 0.1 GeV");
  label.DrawLatex(0.06, 0.71, "|z_{vtx}^{truth}| < 60 cm");
}

double smallest_positive(const std::array<std::unique_ptr<TH1D>, kSpectrumCount>& histograms, const std::vector<std::size_t>& indices)
{
  double result = std::numeric_limits<double>::infinity();
  for (std::size_t index : indices)
  {
    for (int bin = 1; bin <= histograms[index]->GetNbinsX(); ++bin)
    {
      const double value = histograms[index]->GetBinContent(bin);
      if (value > 0.0) result = std::min(result, value);
    }
  }
  return std::isfinite(result) ? result : 0.0;
}

void draw_spectrum(const std::array<std::unique_ptr<TH1D>, kSpectrumCount>& density, const std::vector<std::size_t>& indices,
                   const std::string& output_path, const std::string& family_label, bool detailed)
{
  TCanvas canvas(("c_" + std::string(detailed ? "detailed" : "summary") + "_region_a_topology_spectrum").c_str(), "", kCanvasWidth, kCanvasHeight);
  auto plot_pad = make_plot_pad(std::string(detailed ? "detailed" : "summary") + "_region_a_topology_spectrum_pad", true);
  double maximum = 0.0;
  for (std::size_t index : indices) maximum = std::max(maximum, density[index]->GetMaximum());
  const double minimum = smallest_positive(density, indices);
  density[indices.front()]->SetMinimum(minimum > 0.0 ? 0.5 * minimum : 1e-12);
  density[indices.front()]->SetMaximum(maximum > 0.0 ? 2.0 * maximum : 1.0);
  density[indices.front()]->Draw("HIST");
  for (std::size_t position = 1; position < indices.size(); ++position) density[indices[position]]->Draw("HIST SAME");
  canvas.cd();
  TLegend legend(detailed ? 0.48 : 0.56, detailed ? 0.62 : 0.68, 0.95, 0.97);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetTextSize(detailed ? 0.016 : 0.024);
  for (std::size_t index : indices) legend.AddEntry(density[index].get(), kLabels[index], "l");
  legend.Draw();
  draw_annotations(family_label);
  plot_pad->cd();
  plot_pad->RedrawAxis();
  canvas.cd();
  canvas.SaveAs(output_path.c_str());
}

std::vector<std::unique_ptr<TH1D>> make_fractions(const std::array<std::unique_ptr<TH1D>, kSpectrumCount>& weighted_pb,
                                                  const std::vector<std::size_t>& indices, const std::string& prefix)
{
  std::vector<std::unique_ptr<TH1D>> result;
  for (std::size_t index : indices)
  {
    const std::string name = prefix + kKeys[index] + "_fraction";
    auto histogram = std::unique_ptr<TH1D>(static_cast<TH1D*>(weighted_pb[index]->Clone(name.c_str())));
    histogram->SetDirectory(nullptr);
    histogram->Divide(weighted_pb[index].get(), weighted_pb[1].get(), 1.0, 1.0, "B");
    histogram->SetLineColor(kColors[index]);
    histogram->SetMarkerColor(kColors[index]);
    histogram->SetMarkerStyle(20 + static_cast<int>(result.size()));
    histogram->SetMarkerSize(0.9);
    histogram->SetLineWidth(2);
    style_axes(histogram.get(), "Fraction");
    result.push_back(std::move(histogram));
  }
  return result;
}

void draw_fraction_lines(const std::vector<std::unique_ptr<TH1D>>& fractions, const std::vector<std::size_t>& indices,
                         const std::string& output_path, const std::string& family_label, bool detailed)
{
  TCanvas canvas(("c_" + std::string(detailed ? "detailed" : "summary") + "_region_a_topology_fraction").c_str(), "", kCanvasWidth, kCanvasHeight);
  auto plot_pad = make_plot_pad(std::string(detailed ? "detailed" : "summary") + "_region_a_topology_fraction_pad");
  fractions.front()->SetMinimum(0.0);
  fractions.front()->SetMaximum(1.05);
  fractions.front()->Draw("E1");
  for (std::size_t index = 1; index < fractions.size(); ++index) fractions[index]->Draw("E1 SAME");
  canvas.cd();
  TLegend legend(detailed ? 0.48 : 0.56, detailed ? 0.62 : 0.73, 0.95, 0.97);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetTextSize(detailed ? 0.018 : 0.024);
  for (std::size_t index = 0; index < fractions.size(); ++index) legend.AddEntry(fractions[index].get(), kLabels[indices[index]], "lep");
  legend.Draw();
  draw_annotations(family_label);
  plot_pad->cd();
  plot_pad->RedrawAxis();
  canvas.cd();
  canvas.SaveAs(output_path.c_str());
}

void draw_fraction_stack(std::vector<std::unique_ptr<TH1D>>& fractions, const std::vector<std::size_t>& indices,
                         const std::string& output_path, const std::string& family_label, bool detailed)
{
  TCanvas canvas(("c_" + std::string(detailed ? "detailed" : "summary") + "_region_a_topology_fraction_stack").c_str(), "", kCanvasWidth, kCanvasHeight);
  auto plot_pad = make_plot_pad(std::string(detailed ? "detailed" : "summary") + "_region_a_topology_fraction_stack_pad");
  THStack stack(("stack_" + std::string(detailed ? "detailed" : "summary") + "_region_a_topology_fraction").c_str(), "");
  for (std::size_t index = 0; index < fractions.size(); ++index)
  {
    fractions[index]->SetFillColor(kColors[indices[index]]);
    fractions[index]->SetLineColor(kBlack);
    fractions[index]->SetMarkerStyle(0);
    stack.Add(fractions[index].get());
  }
  stack.SetMinimum(0.0);
  stack.SetMaximum(1.05);
  stack.Draw("HIST");
  stack.GetXaxis()->SetTitle("Cluster E_{T} [GeV]");
  stack.GetYaxis()->SetTitle("Fraction");
  stack.GetXaxis()->SetLabelSize(0.05);
  stack.GetYaxis()->SetLabelSize(0.05);
  stack.GetXaxis()->SetTitleSize(0.05);
  stack.GetYaxis()->SetTitleSize(0.05);
  stack.GetXaxis()->CenterTitle();
  stack.GetYaxis()->CenterTitle();
  canvas.cd();
  TLegend legend(detailed ? 0.48 : 0.56, detailed ? 0.62 : 0.73, 0.95, 0.97);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetTextSize(detailed ? 0.018 : 0.024);
  for (std::size_t index = 0; index < fractions.size(); ++index) legend.AddEntry(fractions[index].get(), kLabels[indices[index]], "f");
  legend.Draw();
  draw_annotations(family_label);
  plot_pad->cd();
  plot_pad->RedrawAxis();
  canvas.cd();
  canvas.SaveAs(output_path.c_str());
}
}

int ReducePythiaRegionAAnchorTopology(
    const std::string family = "jet",
    const std::string map_root = "/sphenix/user/ryotaro/DirectPhotonAnalysis/PhotonAnalysisTree/output/intermediate_files/photon_candidate_selection",
    const std::string output_base = "",
    const bool require_complete = true,
    const int n_bins = 200,
    const double et_max = 40.0)
{
  const std::vector<SampleDefinition> samples = sample_definitions(family);
  if (samples.empty() || map_root.empty() || n_bins <= 0 || !std::isfinite(et_max) || et_max <= 0.0)
  {
    std::cerr << "Family must be jet or photonjet, with valid paths and binning" << std::endl;
    return 1;
  }
  const std::string resolved_output_base = output_base.empty()
      ? "/sphenix/user/ryotaro/DirectPhotonAnalysis/PhotonAnalysisTree/output/plots/photon_candidate_selection/region_a_pi0_anchor_topology/" + family + "/region_a_pi0_anchor_topology"
      : output_base;
  if (!make_output_directory(resolved_output_base)) return 2;

  Spectra spectra(n_bins, et_max);
  std::vector<std::string> sample_names;
  std::vector<unsigned long long> sample_map_counts;
  std::vector<unsigned long long> sample_events_written;
  std::vector<unsigned long long> sample_events_stitch_pass;
  std::vector<unsigned long long> sample_region_a_clusters;
  std::vector<unsigned long long> sample_region_a_prompt_clusters;
  std::vector<unsigned long long> sample_region_a_anchor_clusters;
  std::vector<double> sample_cross_sections_pb;
  std::vector<double> sample_sum_generator_weights;
  std::string family_analysis_release;
  std::string family_model_sha256;

  for (const SampleDefinition& sample : samples)
  {
    const std::string pattern = map_root + "/" + sample.name + "/map_*.root";
    std::vector<MapMetadata> maps;
    if (!collect_maps(pattern, sample, require_complete, maps)) return 3;
    if (maps.empty()) continue;
    if (family_analysis_release.empty())
    {
      family_analysis_release = maps.front().analysis_release;
      family_model_sha256 = maps.front().model_sha256;
    }
    else if (maps.front().analysis_release != family_analysis_release || maps.front().model_sha256 != family_model_sha256)
    {
      std::cerr << "Analysis release or BDT model differs across samples at " << sample.name << std::endl;
      return 3;
    }
    double sample_sumw = 0.0;
    unsigned long long expected_written = 0;
    for (const MapMetadata& map : maps)
    {
      sample_sumw += map.sum_generator_weight_processed;
      expected_written += map.events_written;
    }
    if (!std::isfinite(sample_sumw) || sample_sumw <= 0.0) return 4;

    unsigned long long events_written = 0;
    unsigned long long events_stitch_pass = 0;
    unsigned long long region_a_clusters = 0;
    unsigned long long region_a_prompt_clusters = 0;
    unsigned long long region_a_anchor_clusters = 0;
    for (const MapMetadata& map : maps)
    {
      TFile file(map.path.c_str(), "READ");
      auto* tree = file.Get<TTree>("event_tree");
      if (file.IsZombie() || !tree) return 5;
      unsigned int ncluster = 0;
      unsigned char event_weight_valid = 0;
      unsigned char stitching_valid = 0;
      unsigned char stitching_pass = 0;
      double weight_numerator_pb = 0.0;
      std::vector<double>* cluster_et = nullptr;
      std::vector<unsigned char>* region_a = nullptr;
      std::vector<unsigned char>* prompt = nullptr;
      std::vector<unsigned char>* anchor_valid = nullptr;
      std::vector<int>* topology = nullptr;
      std::vector<int>* missing_category = nullptr;
      bool ok = true;
      ok &= bind(tree, "split_ncluster", &ncluster);
      ok &= bind(tree, "event_weight_valid", &event_weight_valid);
      ok &= bind(tree, "sample_stitching_valid", &stitching_valid);
      ok &= bind(tree, "sample_stitching_pass", &stitching_pass);
      ok &= bind(tree, "weight_numerator_pb", &weight_numerator_pb);
      ok &= bind(tree, "split_cluster_et", &cluster_et);
      ok &= bind(tree, "split_cluster_pass_region_a", &region_a);
      ok &= bind(tree, "split_cluster_truth_prompt_cluster", &prompt);
      ok &= bind(tree, "split_cluster_pi0_anchor_valid", &anchor_valid);
      ok &= bind(tree, "split_cluster_pi0_anchor_topology", &topology);
      ok &= bind(tree, "split_cluster_pi0_anchor_missing_category", &missing_category);
      if (!ok) return 5;

      events_written += static_cast<unsigned long long>(tree->GetEntries());
      for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry)
      {
        tree->GetEntry(entry);
        if (!cluster_et || !region_a || !prompt || !anchor_valid || !topology || !missing_category ||
            cluster_et->size() != ncluster || region_a->size() != ncluster || prompt->size() != ncluster ||
            anchor_valid->size() != ncluster || topology->size() != ncluster || missing_category->size() != ncluster)
        {
          std::cerr << "Malformed cluster vectors in " << map.path << " entry " << entry << std::endl;
          return 6;
        }
        if (!event_weight_valid || !stitching_valid || !stitching_pass) continue;
        if (!std::isfinite(weight_numerator_pb)) return 6;
        ++events_stitch_pass;
        const double weight_pb = weight_numerator_pb / sample_sumw;
        for (std::size_t cluster = 0; cluster < ncluster; ++cluster)
        {
          const double et = (*cluster_et)[cluster];
          if (!std::isfinite(et)) return 6;
          if (!(*region_a)[cluster]) continue;
          ++region_a_clusters;
          if ((*prompt)[cluster])
          {
            spectra.fill(0, et, weight_pb);
            ++region_a_prompt_clusters;
          }
          if (!(*anchor_valid)[cluster]) continue;
          spectra.fill(1, et, weight_pb);
          ++region_a_anchor_clusters;
          const int topology_value = (*topology)[cluster];
          if (topology_value == 1) spectra.fill(2, et, weight_pb);
          else if (topology_value == 2) spectra.fill(3, et, weight_pb);
          else if (topology_value == 4) spectra.fill(4, et, weight_pb);
          else if (topology_value == 0) spectra.fill(13, et, weight_pb);
          else if (topology_value == 3)
          {
            spectra.fill(5, et, weight_pb);
            const int missing_value = (*missing_category)[cluster];
            if (missing_value == 1) spectra.fill(6, et, weight_pb);
            else if (missing_value == 4) spectra.fill(7, et, weight_pb);
            else if (missing_value == 2) spectra.fill(8, et, weight_pb);
            else if (missing_value == 5) spectra.fill(9, et, weight_pb);
            else if (missing_value == 6) spectra.fill(10, et, weight_pb);
            else if (missing_value == 7) spectra.fill(11, et, weight_pb);
            else if (missing_value == 3) spectra.fill(12, et, weight_pb);
            else
            {
              std::cerr << "Invalid stored missing category in " << map.path << " entry " << entry << std::endl;
              return 6;
            }
          }
          else
          {
            std::cerr << "Invalid stored topology in " << map.path << " entry " << entry << std::endl;
            return 6;
          }
        }
      }
    }
    if (events_written != expected_written) return 6;
    sample_names.emplace_back(sample.name);
    sample_map_counts.push_back(maps.size());
    sample_events_written.push_back(events_written);
    sample_events_stitch_pass.push_back(events_stitch_pass);
    sample_region_a_clusters.push_back(region_a_clusters);
    sample_region_a_prompt_clusters.push_back(region_a_prompt_clusters);
    sample_region_a_anchor_clusters.push_back(region_a_anchor_clusters);
    sample_cross_sections_pb.push_back(sample.cross_section_pb);
    sample_sum_generator_weights.push_back(sample_sumw);
    std::cout << "Reduce sample/maps/events/stitch/regionA/prompt/anchor/sumw = " << sample.name << "/" << maps.size() << "/"
              << events_written << "/" << events_stitch_pass << "/" << region_a_clusters << "/" << region_a_prompt_clusters
              << "/" << region_a_anchor_clusters << "/" << sample_sumw << std::endl;
  }

  if (sample_names.empty())
  {
    std::cerr << "No sample maps were reduced" << std::endl;
    return 3;
  }
  if (!valid_partition(spectra.counts) || !valid_partition(spectra.weighted_pb))
  {
    std::cerr << "Topology partition failed" << std::endl;
    return 7;
  }

  std::array<std::unique_ptr<TH1D>, kSpectrumCount> density;
  for (std::size_t index = 0; index < kSpectrumCount; ++index)
  {
    density[index].reset(static_cast<TH1D*>(spectra.weighted_pb[index]->Clone((std::string("h_region_a_") + kKeys[index] + "_et_pb_per_gev").c_str())));
    density[index]->SetDirectory(nullptr);
    density[index]->Scale(1.0, "width");
    density[index]->SetLineColor(kColors[index]);
    density[index]->SetLineWidth(index < 2 ? 3 : 2);
    style_axes(density[index].get(), "d#sigma/dE_{T} [pb/GeV]");
    spectra.counts[index]->SetLineColor(kColors[index]);
    spectra.weighted_pb[index]->SetLineColor(kColors[index]);
  }
  const std::vector<std::size_t> detailed_spectrum_indices = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
  const std::vector<std::size_t> summary_spectrum_indices(kSummarySpectrum.begin(), kSummarySpectrum.end());
  const std::vector<std::size_t> detailed_category_indices(kDetailedCategories.begin(), kDetailedCategories.end());
  const std::vector<std::size_t> summary_category_indices(kSummaryCategories.begin(), kSummaryCategories.end());
  auto detailed_fractions = make_fractions(spectra.weighted_pb, detailed_category_indices, "h_region_a_detailed_");
  auto summary_fractions = make_fractions(spectra.weighted_pb, summary_category_indices, "h_region_a_summary_");

  SetsPhenixStyle();
  const std::string family_label = family == "jet" ? "Pythia8 p+p Jet samples" : "Pythia8 p+p PhotonJet samples";
  draw_spectrum(density, detailed_spectrum_indices, resolved_output_base + "_detailed.pdf", family_label, true);
  draw_spectrum(density, summary_spectrum_indices, resolved_output_base + ".pdf", family_label, false);
  draw_fraction_lines(detailed_fractions, detailed_category_indices, resolved_output_base + "_category_fractions_detailed.pdf", family_label, true);
  draw_fraction_lines(summary_fractions, summary_category_indices, resolved_output_base + "_category_fractions.pdf", family_label, false);
  draw_fraction_stack(detailed_fractions, detailed_category_indices, resolved_output_base + "_category_fraction_stack_detailed.pdf", family_label, true);
  draw_fraction_stack(summary_fractions, summary_category_indices, resolved_output_base + "_category_fraction_stack.pdf", family_label, false);

  TFile output((resolved_output_base + ".root").c_str(), "RECREATE");
  if (output.IsZombie()) return 8;
  for (std::size_t index = 0; index < kSpectrumCount; ++index)
  {
    spectra.counts[index]->Write();
    spectra.weighted_pb[index]->Write();
    density[index]->Write();
  }
  for (auto& histogram : detailed_fractions) histogram->Write();
  for (auto& histogram : summary_fractions) histogram->Write();

  int output_schema_version = 1;
  int source_map_schema_version = 3;
  std::string selection_definition = "sample_stitching_valid_and_pass_and_region_a_and_pi0_anchor_valid";
  std::string prompt_definition = "sample_stitching_valid_and_pass_and_region_a_and_truth_prompt_cluster";
  std::string topology_source = "stored_map_topology_evaluated_with_strict_cluster_energy_gt_0p1_GeV";
  std::string weight_definition = "sample_cross_section_pb_times_generator_weight_divided_by_full_sample_sum_generator_weight_processed";
  std::string metadata_family = family;
  std::string metadata_map_root = map_root;
  bool metadata_require_complete = require_complete;
  int metadata_n_bins = n_bins;
  double metadata_et_max = et_max;
  unsigned int sample_count = sample_names.size();
  TTree metadata("metadata", "Region-A pi0-anchor topology reduce metadata");
  metadata.Branch("schema_version", &output_schema_version);
  metadata.Branch("source_map_schema_version", &source_map_schema_version);
  metadata.Branch("family", &metadata_family);
  metadata.Branch("map_root", &metadata_map_root);
  metadata.Branch("require_complete", &metadata_require_complete);
  metadata.Branch("n_bins", &metadata_n_bins);
  metadata.Branch("et_max", &metadata_et_max);
  metadata.Branch("selection_definition", &selection_definition);
  metadata.Branch("prompt_definition", &prompt_definition);
  metadata.Branch("topology_source", &topology_source);
  metadata.Branch("analysis_release", &family_analysis_release);
  metadata.Branch("model_sha256", &family_model_sha256);
  metadata.Branch("weight_definition", &weight_definition);
  metadata.Branch("sample_count", &sample_count);
  metadata.Branch("sample_names", &sample_names);
  metadata.Branch("sample_map_counts", &sample_map_counts);
  metadata.Branch("sample_events_written", &sample_events_written);
  metadata.Branch("sample_events_stitch_pass", &sample_events_stitch_pass);
  metadata.Branch("sample_region_a_clusters", &sample_region_a_clusters);
  metadata.Branch("sample_region_a_prompt_clusters", &sample_region_a_prompt_clusters);
  metadata.Branch("sample_region_a_anchor_clusters", &sample_region_a_anchor_clusters);
  metadata.Branch("sample_cross_sections_pb", &sample_cross_sections_pb);
  metadata.Branch("sample_sum_generator_weights", &sample_sum_generator_weights);
  metadata.Fill();
  metadata.Write();
  output.Close();

  std::cout << "ReducePythiaRegionAAnchorTopology - family/samples/anchor/output = " << family << "/" << sample_count << "/"
            << spectra.counts[1]->GetEntries() << "/" << resolved_output_base << std::endl;
  return 0;
}
