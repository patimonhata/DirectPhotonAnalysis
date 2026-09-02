#include "../../macro/Utilities/sPhenixStyle.C"

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
#include <limits>
#include <numeric>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace
{
constexpr std::size_t kSpectrumCount = 14;
constexpr std::size_t kSelectionCount = 5;
constexpr std::array<const char*, kSelectionCount> kSelectionKeys = {
    "kinematic", "preselection", "preselection_tight", "preselection_isolation", "region_a_tagging_veto"};
constexpr std::array<const char*, kSelectionCount> kSelectionLabels = {
    "Kinematic", "Pre-selection", "Pre-selection + TightBDT", "Pre-selection + Isolation", "Region A + Tagging veto"};
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
  double min_cluster_energy = -1.0;
  double partner_diagnostic_min_cluster_energy = -1.0;
  double meson_partner_min_energy = -1.0;
  int pi0_topology_algorithm_version = -1;
  unsigned long long events_written = 0;
  unsigned long long events_vertex_rejected = 0;
  unsigned long long events_invalid = 0;
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
  ok &= bind(tree, "min_cluster_energy", &value.min_cluster_energy);
  ok &= bind(tree, "partner_diagnostic_min_cluster_energy", &value.partner_diagnostic_min_cluster_energy);
  ok &= bind(tree, "meson_partner_min_energy", &value.meson_partner_min_energy);
  ok &= bind(tree, "pi0_topology_algorithm_version", &value.pi0_topology_algorithm_version);
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
  return value.schema_version == 4 && value.sample_name == sample.name && !value.input_manifest.empty() &&
      !value.analysis_release.empty() && !value.model_sha256.empty() && value.manifest_begin >= 0 &&
      std::isfinite(value.min_cluster_energy) && value.min_cluster_energy >= 0.0 &&
      std::isfinite(value.partner_diagnostic_min_cluster_energy) && same_double(value.partner_diagnostic_min_cluster_energy, 0.1) &&
      std::isfinite(value.meson_partner_min_energy) && value.meson_partner_min_energy >= 0.0 && value.pi0_topology_algorithm_version == 8 &&
      value.manifest_end > value.manifest_begin && value.input_file_count == value.manifest_end - value.manifest_begin &&
      same_double(value.cross_section_pb, sample.cross_section_pb) && same_double(value.window_min, sample.window_min) &&
      same_double(value.window_max, sample.window_max) && value.upper_unbounded == sample.upper_unbounded &&
      std::isfinite(value.sum_generator_weight_processed) &&
      value.events_processed == value.events_written + value.events_vertex_rejected + value.events_invalid;
}

bool compatible(const MapMetadata& value, const MapMetadata& reference)
{
  return value.schema_version == reference.schema_version && value.input_manifest == reference.input_manifest &&
      value.sample_name == reference.sample_name && value.analysis_release == reference.analysis_release &&
      value.model_sha256 == reference.model_sha256 && same_double(value.cross_section_pb, reference.cross_section_pb) &&
      same_double(value.window_min, reference.window_min) && same_double(value.window_max, reference.window_max) &&
      value.upper_unbounded == reference.upper_unbounded &&
      same_double(value.min_cluster_energy, reference.min_cluster_energy) &&
      same_double(value.partner_diagnostic_min_cluster_energy, reference.partner_diagnostic_min_cluster_energy) &&
      same_double(value.meson_partner_min_energy, reference.meson_partner_min_energy) && value.pi0_topology_algorithm_version == reference.pi0_topology_algorithm_version;
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
    case 4: return pass_region_a != 0U && pi0_tag == 0U && eta_tag == 0U;
    default: return false;
  }
}

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
  label.DrawLatex(0.06, 0.86, "5 < E_{T} < 35 GeV, |#eta| < 0.7");
  label.DrawLatex(0.06, 0.81, "|z_{vtx}^{truth}| < 60 cm");
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
    histogram->Reset("ICES");
    for (int bin = 0; bin <= histogram->GetNbinsX() + 1; ++bin)
    {
      const double denominator = weighted_pb[1]->GetBinContent(bin);
      if (denominator == 0.0) continue;
      const double numerator = weighted_pb[index]->GetBinContent(bin);
      const double value = numerator / denominator;
      const double numerator_w2 = std::pow(weighted_pb[index]->GetBinError(bin), 2);
      const double denominator_w2 = std::pow(weighted_pb[1]->GetBinError(bin), 2);
      const double complement_w2 = std::max(0.0, denominator_w2 - numerator_w2);
      const double variance = (std::pow(1.0 - value, 2) * numerator_w2 + std::pow(value, 2) * complement_w2) / std::pow(denominator, 2);
      histogram->SetBinContent(bin, value);
      histogram->SetBinError(bin, std::sqrt(std::max(0.0, variance)));
    }
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



namespace candidate_composition
{
enum Category : std::size_t
{
  denominator,
  prompt,
  pi0_separated,
  pi0_merged,
  pi0_single_contaminated,
  pi0_missing,
  pi0_other,
  eta,
  other,
  category_count
};

constexpr double kMajorityThreshold = 0.5;
constexpr int kSignalEmbeddingId = 1;
constexpr std::array<const char*, category_count> kKeys = {
    "denominator", "prompt", "pi0_separated", "pi0_merged", "pi0_single_contaminated", "pi0_missing", "pi0_other", "eta", "other"};
constexpr std::array<const char*, category_count> kLabels = {
    "Selected candidates", "Prompt #gamma", "#pi^{0}: separated", "#pi^{0}: merged", "#pi^{0}: single contaminated",
    "#pi^{0}: missing", "#pi^{0}: other", "#eta", "Other"};
constexpr std::array<int, category_count> kColors = {
    kBlack, kRed + 1, kAzure + 7, kMagenta + 1, kCyan + 2, kGreen + 2, kGray + 1, kOrange + 7, kGray + 2};

struct SelectionDefinition
{
  const char* key;
  const char* label;
};

constexpr std::array<SelectionDefinition, 7> kSelectionDefinitions = {{
    {"kinematic", "Kinematic"},
    {"preselection", "Pre-selection"},
    {"preselection_tight", "Pre-selection + TightBDT"},
    {"preselection_isolation", "Pre-selection + Isolation"},
    {"region_a", "Region A: isolated and tight"},
    {"region_a_tagging_veto", "Region A after #pi^{0}/#eta tag veto"},
    {"final_photon", "Region A after #pi^{0}/#eta tag veto"},
}};

int required_shard_count(const std::string& sample_name)
{
  return sample_name == "jet12" ? 10 : 1;
}

const SelectionDefinition* find_selection(const std::string& key)
{
  const auto found = std::find_if(kSelectionDefinitions.begin(), kSelectionDefinitions.end(), [&](const SelectionDefinition& definition) { return key == definition.key; });
  return found == kSelectionDefinitions.end() ? nullptr : &*found;
}

bool passes_selection(const std::string& selection, unsigned char kinematic, unsigned char preselection, unsigned char tight,
                      unsigned char isolated, unsigned char region_a, unsigned char final_photon)
{
  if (selection == "kinematic") return kinematic;
  if (selection == "preselection") return kinematic && preselection;
  if (selection == "preselection_tight") return kinematic && preselection && tight;
  if (selection == "preselection_isolation") return kinematic && preselection && isolated;
  if (selection == "region_a") return region_a;
  return (selection == "region_a_tagging_veto" || selection == "final_photon") && final_photon;
}

struct Histograms
{
  std::array<std::unique_ptr<TH1D>, category_count> counts;
  std::array<std::unique_ptr<TH1D>, category_count> weighted;

  Histograms(int n_bins, double et_max)
  {
    for (std::size_t index = 0; index < category_count; ++index)
    {
      counts[index] = std::make_unique<TH1D>((std::string("h_candidate_") + kKeys[index] + "_et_count").c_str(), "", n_bins, 0.0, et_max);
      weighted[index] = std::make_unique<TH1D>((std::string("h_candidate_") + kKeys[index] + "_et_pb").c_str(), "", n_bins, 0.0, et_max);
      counts[index]->SetDirectory(nullptr);
      weighted[index]->SetDirectory(nullptr);
      counts[index]->Sumw2();
      weighted[index]->Sumw2();
    }
  }

  void fill(std::size_t index, double et, double weight)
  {
    counts[index]->Fill(et);
    weighted[index]->Fill(et, weight);
  }
};

bool valid_partition(const std::array<std::unique_ptr<TH1D>, category_count>& histograms)
{
  for (int bin = 0; bin <= histograms.front()->GetNbinsX() + 1; ++bin)
  {
    double sum = 0.0;
    for (std::size_t index = 1; index < category_count; ++index) sum += histograms[index]->GetBinContent(bin);
    const double scale = std::max({1.0, std::abs(histograms[denominator]->GetBinContent(bin)), std::abs(sum)});
    if (std::abs(histograms[denominator]->GetBinContent(bin) - sum) > 1e-9 * scale) return false;
  }
  return true;
}

std::unique_ptr<TH1D> fraction_histogram(const TH1D& numerator, const TH1D& total, const std::string& name)
{
  auto result = std::unique_ptr<TH1D>(static_cast<TH1D*>(numerator.Clone(name.c_str())));
  result->Reset("ICES");
  result->SetDirectory(nullptr);
  for (int bin = 0; bin <= result->GetNbinsX() + 1; ++bin)
  {
    const double den = total.GetBinContent(bin);
    if (den == 0.0) continue;
    const double num = numerator.GetBinContent(bin);
    const double value = num / den;
    const double num_w2 = std::pow(numerator.GetBinError(bin), 2);
    const double den_w2 = std::pow(total.GetBinError(bin), 2);
    const double complement_w2 = std::max(0.0, den_w2 - num_w2);
    const double variance = (std::pow(1.0 - value, 2) * num_w2 + std::pow(value, 2) * complement_w2) / std::pow(den, 2);
    result->SetBinContent(bin, value);
    result->SetBinError(bin, std::sqrt(std::max(0.0, variance)));
  }
  return result;
}

bool valid_contributors(unsigned int ncluster, const std::vector<unsigned int>& offset, const std::vector<int>& g4_pdg,
                        const std::vector<int>& embedding, const std::vector<float>& fraction,
                        const std::vector<int>& source, const std::vector<int>& parent)
{
  const std::size_t size = fraction.size();
  return offset.size() == static_cast<std::size_t>(ncluster) + 1U && !offset.empty() && offset.front() == 0U && offset.back() == size &&
      g4_pdg.size() == size && embedding.size() == size && source.size() == size && parent.size() == size &&
      std::adjacent_find(offset.begin(), offset.end(), std::greater<unsigned int>()) == offset.end();
}

double eta_fraction(std::size_t cluster, const std::vector<unsigned int>& offset, const std::vector<int>& g4_pdg,
                    const std::vector<int>& embedding, const std::vector<float>& fraction,
                    const std::vector<int>& source, const std::vector<int>& parent)
{
  double result = 0.0;
  for (std::size_t contributor = offset[cluster]; contributor < offset[cluster + 1U]; ++contributor)
  {
    if (embedding[contributor] != kSignalEmbeddingId) continue;
    const bool g4_eta = std::abs(g4_pdg[contributor]) == 221;
    const bool generator_eta_photon = g4_pdg[contributor] == 22 && (parent[contributor] == 221 || source[contributor] == 3);
    if (g4_eta || generator_eta_photon) result += fraction[contributor];
  }
  return result;
}

std::size_t pi0_category(int topology)
{
  if (topology == 1) return pi0_separated;
  if (topology == 2) return pi0_merged;
  if (topology == 4) return pi0_single_contaminated;
  if (topology == 3) return pi0_missing;
  return topology == 0 ? pi0_other : category_count;
}

void draw_stack(std::array<std::unique_ptr<TH1D>, category_count>& fractions, const std::string& output,
                const std::string& family, const std::string& selection, double min_cluster_energy, bool detailed)
{
  SetsPhenixStyle();
  const std::string detail = detailed ? "detailed" : "summary";
  TCanvas canvas(("c_" + detail + "_photon_candidate_composition").c_str(), "", kCanvasWidth, kCanvasHeight);
  auto plot_pad = make_plot_pad(detail + "_photon_candidate_composition_pad");
  THStack stack(("stack_" + detail + "_photon_candidate_composition").c_str(), "");
  std::unique_ptr<TH1D> pi0_fraction;
  if (detailed)
  {
    for (std::size_t index = 1; index < category_count; ++index)
    {
      fractions[index]->SetFillColor(kColors[index]);
      fractions[index]->SetLineColor(kBlack);
      stack.Add(fractions[index].get());
    }
  }
  else
  {
    pi0_fraction.reset(static_cast<TH1D*>(fractions[pi0_separated]->Clone("h_candidate_pi0_fraction")));
    pi0_fraction->Reset("ICES");
    pi0_fraction->SetDirectory(nullptr);
    for (std::size_t index = pi0_separated; index <= pi0_other; ++index) pi0_fraction->Add(fractions[index].get());
    pi0_fraction->SetFillColor(kColors[pi0_separated]);
    pi0_fraction->SetLineColor(kBlack);
    for (std::size_t index : {prompt, eta, other})
    {
      fractions[index]->SetFillColor(kColors[index]);
      fractions[index]->SetLineColor(kBlack);
    }
    stack.Add(fractions[prompt].get());
    stack.Add(pi0_fraction.get());
    stack.Add(fractions[eta].get());
    stack.Add(fractions[other].get());
  }
  stack.SetMinimum(0.0);
  stack.SetMaximum(1.05);
  stack.Draw("HIST");
  stack.GetXaxis()->SetTitle("Cluster E_{T} [GeV]");
  stack.GetYaxis()->SetTitle("Fraction of selected candidates");
  for (TAxis* axis : {stack.GetXaxis(), stack.GetYaxis()})
  {
    axis->SetLabelSize(0.05);
    axis->SetTitleSize(0.05);
    axis->CenterTitle();
  }
  stack.GetXaxis()->SetTitleOffset(1.4);
  stack.GetYaxis()->SetTitleOffset(1.15);
  canvas.cd();
  TLegend legend(detailed ? 0.48 : 0.56, detailed ? 0.60 : 0.73, 0.95, 0.97);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetTextSize(detailed ? 0.018 : 0.024);
  if (detailed)
  {
    for (std::size_t index = 1; index < category_count; ++index) legend.AddEntry(fractions[index].get(), kLabels[index], "f");
  }
  else
  {
    legend.AddEntry(fractions[prompt].get(), kLabels[prompt], "f");
    legend.AddEntry(pi0_fraction.get(), "#pi^{0}", "f");
    legend.AddEntry(fractions[eta].get(), kLabels[eta], "f");
    legend.AddEntry(fractions[other].get(), kLabels[other], "f");
  }
  legend.Draw();
  TLatex label;
  label.SetNDC();
  label.SetTextAlign(13);
  label.SetTextSize(0.026);
  label.DrawLatex(0.06, 0.96, "#it{#bf{sPHENIX}} Internal");
  label.DrawLatex(0.06, 0.91, (family == "jet" ? "Pythia8 p+p Jet samples" : "Pythia8 p+p PhotonJet samples"));
  const SelectionDefinition* selected_definition = find_selection(selection);
  label.DrawLatex(0.06, 0.86, selected_definition ? selected_definition->label : selection.c_str());
  std::ostringstream threshold;
  threshold << "E_{cluster} > " << min_cluster_energy << " GeV; truth contribution > 50%";
  label.DrawLatex(0.06, 0.81, threshold.str().c_str());
  plot_pad->cd();
  plot_pad->RedrawAxis();
  canvas.cd();
  canvas.SaveAs(output.c_str());
}
}

int ReducePythiaPhotonCandidateSelection(
    const std::string family,
    const std::string map_root,
    const std::string output_base,
    const std::string selection,
    const bool require_complete,
    const int n_bins,
    const double et_max,
    const std::string sample_name,
    const int shard_index)
{
  using namespace candidate_composition;
  const SelectionDefinition* selected_definition = find_selection(selection);
  const std::vector<SampleDefinition> family_samples = sample_definitions(family);
  const auto selected_sample = std::find_if(family_samples.begin(), family_samples.end(), [&](const SampleDefinition& sample) { return sample.name == sample_name; });
  const int shard_count = required_shard_count(sample_name);
  if (selected_sample == family_samples.end() || shard_index < 0 || shard_index >= shard_count || map_root.empty() || !selected_definition ||
      n_bins <= 0 || !std::isfinite(et_max) || et_max <= 0.0) return 1;
  std::string normalized_map_root = map_root;
  while (normalized_map_root.size() > 1U && normalized_map_root.back() == '/') normalized_map_root.pop_back();
  const std::string configuration = normalized_map_root.substr(normalized_map_root.find_last_of('/') + 1U);
  const std::string default_output_suffix = "/partial/" + sample_name + "/shard_" + std::to_string(shard_index) + "/photon_candidate_selection";
  const std::string resolved_output_base = output_base.empty()
      ? "/sphenix/user/ryotaro/DirectPhotonAnalysis/PhotonAnalysisTree/output/plots/photon_candidate_selection/reduce/" +
            configuration + "/" + selection + "/" + family + default_output_suffix
      : output_base;
  if (!make_output_directory(resolved_output_base)) return 2;

  Histograms histograms(n_bins, et_max);
  std::array<std::unique_ptr<Spectra>, kSelectionCount> topology_histograms;
  for (std::size_t topology_selection = 0; topology_selection < kSelectionCount; ++topology_selection)
    topology_histograms[topology_selection] = std::make_unique<Spectra>(n_bins, et_max);
  std::vector<std::string> sample_names;
  std::vector<unsigned long long> sample_map_counts;
  std::vector<double> sample_sum_generator_weights;
  std::string analysis_release, model_sha256;
  double min_cluster_energy = -1.0;
  double partner_diagnostic_min_cluster_energy = -1.0;
  double meson_partner_min_energy = -1.0;
  int pi0_topology_algorithm_version = -1;
  unsigned long long selected_count = 0, prompt_count = 0, pi0_count = 0, eta_count = 0, other_count = 0;
  unsigned long long overlap_count = 0, half_boundary_count = 0, invalid_truth_count = 0;
  unsigned long long events_written = 0, expected_events_written = 0, events_stitch_pass = 0;
  unsigned long long region_a_clusters = 0, region_a_prompt_clusters = 0, region_a_anchor_clusters = 0;
  std::array<unsigned long long, kSelectionCount> topology_selected_clusters = {};
  std::array<unsigned long long, kSelectionCount> topology_selected_anchor_clusters = {};

  const SampleDefinition& sample = *selected_sample;
  unsigned long long shard_map_begin = 0, shard_map_end = 0, total_map_count = 0;
  std::vector<MapMetadata> maps;
  if (!collect_maps(map_root + "/" + sample.name + "/map_*.root", sample, require_complete, maps)) return 3;
  if (maps.empty()) return 3;
  analysis_release = maps.front().analysis_release;
  model_sha256 = maps.front().model_sha256;
  min_cluster_energy = maps.front().min_cluster_energy;
  partner_diagnostic_min_cluster_energy = maps.front().partner_diagnostic_min_cluster_energy;
  meson_partner_min_energy = maps.front().meson_partner_min_energy;
  pi0_topology_algorithm_version = maps.front().pi0_topology_algorithm_version;
  double sample_sumw = 0.0;
  for (const auto& map : maps) sample_sumw += map.sum_generator_weight_processed;
  if (!std::isfinite(sample_sumw) || sample_sumw <= 0.0) return 4;
  total_map_count = maps.size();
  shard_map_begin = total_map_count * static_cast<unsigned long long>(shard_index) / static_cast<unsigned long long>(shard_count);
  shard_map_end = total_map_count * static_cast<unsigned long long>(shard_index + 1) / static_cast<unsigned long long>(shard_count);
  if (shard_map_begin >= shard_map_end) return 4;
  for (unsigned long long map_index = shard_map_begin; map_index < shard_map_end; ++map_index)
    expected_events_written += maps[map_index].events_written;
  sample_names.push_back(sample.name);
  sample_map_counts.push_back(maps.size());
  sample_sum_generator_weights.push_back(sample_sumw);

  for (unsigned long long map_index = shard_map_begin; map_index < shard_map_end; ++map_index)
  {
    const MapMetadata& map = maps[map_index];
    TFile file(map.path.c_str(), "READ");
    auto* tree = file.Get<TTree>("event_tree");
    if (file.IsZombie() || !tree) return 5;
    unsigned int ncluster = 0;
    unsigned char weight_valid = 0, stitch_valid = 0, stitch_pass = 0;
    double weight_numerator = 0.0;
    std::vector<double>* et = nullptr;
    std::vector<unsigned char>* pass_kinematics = nullptr;
    std::vector<unsigned char>* pass_preselection = nullptr;
    std::vector<unsigned char>* pass_tight = nullptr;
    std::vector<unsigned char>* pass_isolated = nullptr;
    std::vector<unsigned char>* pass_region_a = nullptr;
    std::vector<unsigned char>* pi0_tag = nullptr;
    std::vector<unsigned char>* eta_tag = nullptr;
    std::vector<unsigned char>* pass_final_photon = nullptr;
    std::vector<unsigned char>* truth_valid = nullptr;
    std::vector<unsigned char>* prompt_flag = nullptr;
    std::vector<float>* dominant_fraction = nullptr;
    std::vector<unsigned char>* pi0_valid = nullptr;
    std::vector<float>* pi0_fraction = nullptr;
    std::vector<int>* topology = nullptr;
    std::vector<int>* missing_category = nullptr;
    std::vector<unsigned int>* offset = nullptr;
    std::vector<int>* g4_pdg = nullptr;
    std::vector<int>* embedding = nullptr;
    std::vector<float>* fraction = nullptr;
    std::vector<int>* source = nullptr;
    std::vector<int>* parent = nullptr;
    tree->SetBranchStatus("*", false);
    tree->SetCacheSize(kEventTreeCacheSize);
    bool ok = bind_active(tree, "split_ncluster", &ncluster) &&
        bind_active(tree, "event_weight_valid", &weight_valid) &&
        bind_active(tree, "sample_stitching_valid", &stitch_valid) &&
        bind_active(tree, "sample_stitching_pass", &stitch_pass) &&
        bind_active(tree, "weight_numerator_pb", &weight_numerator) &&
        bind_active(tree, "split_cluster_et", &et) &&
        bind_active(tree, "split_cluster_pass_kinematics", &pass_kinematics) &&
        bind_active(tree, "split_cluster_pass_preselection", &pass_preselection) &&
        bind_active(tree, "split_cluster_pass_tight", &pass_tight) &&
        bind_active(tree, "split_cluster_pass_isolated", &pass_isolated) &&
        bind_active(tree, "split_cluster_pass_region_a", &pass_region_a) &&
        bind_active(tree, "split_cluster_pi0_tag", &pi0_tag) &&
        bind_active(tree, "split_cluster_eta_tag", &eta_tag) &&
        bind_active(tree, "split_cluster_pass_final_photon", &pass_final_photon) &&
        bind_active(tree, "split_cluster_truth_valid", &truth_valid) &&
        bind_active(tree, "split_cluster_truth_prompt_cluster", &prompt_flag) &&
        bind_active(tree, "split_cluster_truth_dominant_fraction", &dominant_fraction) &&
        bind_active(tree, "split_cluster_pi0_anchor_valid", &pi0_valid) &&
        bind_active(tree, "split_cluster_pi0_anchor_main_fraction", &pi0_fraction) &&
        bind_active(tree, "split_cluster_pi0_anchor_topology", &topology) &&
        bind_active(tree, "split_cluster_pi0_anchor_missing_category", &missing_category) &&
        bind_active(tree, "split_cluster_truth_contributor_offset", &offset) &&
        bind_active(tree, "split_cluster_truth_contributor_g4_pdg_id", &g4_pdg) &&
        bind_active(tree, "split_cluster_truth_contributor_embedding_id", &embedding) &&
        bind_active(tree, "split_cluster_truth_contributor_fraction", &fraction) &&
        bind_active(tree, "split_cluster_truth_contributor_photon_source", &source) &&
        bind_active(tree, "split_cluster_truth_contributor_classification_parent_pdg", &parent);
    if (!ok) return 5;
    tree->StopCacheLearningPhase();
    events_written += static_cast<unsigned long long>(tree->GetEntries());

    for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry)
    {
      tree->GetEntry(entry);
      if (!et || !pass_kinematics || !pass_preselection || !pass_tight || !pass_isolated || !pass_region_a || !pi0_tag || !eta_tag ||
          !pass_final_photon || !truth_valid || !prompt_flag || !dominant_fraction || !pi0_valid || !pi0_fraction || !topology || !missing_category ||
          !offset || !g4_pdg || !embedding || !fraction || !source || !parent || et->size() != ncluster || pass_kinematics->size() != ncluster ||
          pass_preselection->size() != ncluster || pass_tight->size() != ncluster || pass_isolated->size() != ncluster ||
          pass_region_a->size() != ncluster || pi0_tag->size() != ncluster || eta_tag->size() != ncluster || pass_final_photon->size() != ncluster ||
          truth_valid->size() != ncluster || prompt_flag->size() != ncluster || dominant_fraction->size() != ncluster ||
          pi0_valid->size() != ncluster || pi0_fraction->size() != ncluster || topology->size() != ncluster || missing_category->size() != ncluster ||
          !valid_contributors(ncluster, *offset, *g4_pdg, *embedding, *fraction, *source, *parent)) return 6;
      if (!weight_valid || !stitch_valid || !stitch_pass) continue;
      ++events_stitch_pass;
      const double weight = weight_numerator / sample_sumw;
      if (!std::isfinite(weight)) return 6;
      for (std::size_t cluster = 0; cluster < ncluster; ++cluster)
      {
        const double cluster_et = (*et)[cluster];
        if (!std::isfinite(cluster_et)) return 6;
        const bool computed_final_photon = (*pass_region_a)[cluster] && !(*pi0_tag)[cluster] && !(*eta_tag)[cluster];
        if (static_cast<bool>((*pass_final_photon)[cluster]) != computed_final_photon) return 6;
        if ((*pass_region_a)[cluster])
        {
          ++region_a_clusters;
          if ((*prompt_flag)[cluster]) ++region_a_prompt_clusters;
          if ((*pi0_valid)[cluster]) ++region_a_anchor_clusters;
        }
        for (std::size_t topology_selection = 0; topology_selection < kSelectionCount; ++topology_selection)
        {
          if (!pass_selection(topology_selection, (*pass_kinematics)[cluster], (*pass_preselection)[cluster], (*pass_tight)[cluster],
                              (*pass_isolated)[cluster], (*pass_region_a)[cluster], (*pi0_tag)[cluster], (*eta_tag)[cluster])) continue;
          ++topology_selected_clusters[topology_selection];
          Spectra& spectra = *topology_histograms[topology_selection];
          if ((*prompt_flag)[cluster]) spectra.fill(0, cluster_et, weight);
          if (!(*pi0_valid)[cluster]) continue;
          ++topology_selected_anchor_clusters[topology_selection];
          spectra.fill(1, cluster_et, weight);
          const int topology_value = (*topology)[cluster];
          if (topology_value == 1) spectra.fill(2, cluster_et, weight);
          else if (topology_value == 2) spectra.fill(3, cluster_et, weight);
          else if (topology_value == 4) spectra.fill(4, cluster_et, weight);
          else if (topology_value == 0) spectra.fill(13, cluster_et, weight);
          else if (topology_value == 3)
          {
            spectra.fill(5, cluster_et, weight);
            const int missing_value = (*missing_category)[cluster];
            if (missing_value == 1) spectra.fill(6, cluster_et, weight);
            else if (missing_value == 4) spectra.fill(7, cluster_et, weight);
            else if (missing_value == 2) spectra.fill(8, cluster_et, weight);
            else if (missing_value == 5) spectra.fill(9, cluster_et, weight);
            else if (missing_value == 6) spectra.fill(10, cluster_et, weight);
            else if (missing_value == 7) spectra.fill(11, cluster_et, weight);
            else if (missing_value == 3) spectra.fill(12, cluster_et, weight);
            else return 6;
          }
          else return 6;
        }
        const bool passes_composition = candidate_composition::passes_selection(
            selection, (*pass_kinematics)[cluster], (*pass_preselection)[cluster], (*pass_tight)[cluster],
            (*pass_isolated)[cluster], (*pass_region_a)[cluster], (*pass_final_photon)[cluster]);
        if (!passes_composition) continue;
        const double eta_contribution = eta_fraction(cluster, *offset, *g4_pdg, *embedding, *fraction, *source, *parent);
        if (!std::isfinite(eta_contribution) || eta_contribution < -1e-6 || eta_contribution > 1.0 + 1e-5) return 6;
        ++selected_count;
        histograms.fill(denominator, cluster_et, weight);
        if (!(*truth_valid)[cluster]) ++invalid_truth_count;
        const bool is_prompt = (*prompt_flag)[cluster] && (*dominant_fraction)[cluster] > kMajorityThreshold;
        const bool is_pi0 = (*pi0_valid)[cluster] && (*pi0_fraction)[cluster] > kMajorityThreshold;
        const bool is_eta = eta_contribution > kMajorityThreshold;
        const int majority_count = static_cast<int>(is_prompt) + static_cast<int>(is_pi0) + static_cast<int>(is_eta);
        if (((*prompt_flag)[cluster] && (*dominant_fraction)[cluster] == kMajorityThreshold) ||
            ((*pi0_valid)[cluster] && (*pi0_fraction)[cluster] == kMajorityThreshold) || eta_contribution == kMajorityThreshold) ++half_boundary_count;
        std::size_t category = other;
        if (majority_count > 1) ++overlap_count;
        else if (is_prompt) category = prompt;
        else if (is_pi0) category = pi0_category((*topology)[cluster]);
        else if (is_eta) category = eta;
        if (category == category_count) return 6;
        histograms.fill(category, cluster_et, weight);
        if (category == prompt) ++prompt_count;
        else if (category >= pi0_separated && category <= pi0_other) ++pi0_count;
        else if (category == eta) ++eta_count;
        else ++other_count;
      }
    }
  }

  if (events_written != expected_events_written || sample_names.empty() || selected_count != prompt_count + pi0_count + eta_count + other_count ||
      !valid_partition(histograms.counts) || !valid_partition(histograms.weighted)) return 7;
  for (std::size_t topology_selection = 0; topology_selection < kSelectionCount; ++topology_selection)
    if (!::valid_partition(topology_histograms[topology_selection]->counts) || !::valid_partition(topology_histograms[topology_selection]->weighted_pb)) return 7;
  TFile output((resolved_output_base + ".root").c_str(), "RECREATE");
  if (output.IsZombie()) return 8;
  TDirectory* composition_directory = output.mkdir("composition");
  if (!composition_directory) return 8;
  composition_directory->cd();
  for (std::size_t index = 0; index < category_count; ++index)
  {
    histograms.counts[index]->Write();
    histograms.weighted[index]->Write();
  }
  output.cd();
  TDirectory* topology_root_directory = output.mkdir("anchor_topology");
  if (!topology_root_directory) return 8;
  for (std::size_t topology_selection = 0; topology_selection < kSelectionCount; ++topology_selection)
  {
    TDirectory* topology_directory = topology_root_directory->mkdir(kSelectionKeys[topology_selection]);
    if (!topology_directory) return 8;
    topology_directory->cd();
    for (std::size_t index = 0; index < kSpectrumCount; ++index)
    {
      topology_histograms[topology_selection]->counts[index]->Write();
      topology_histograms[topology_selection]->weighted_pb[index]->Write();
    }
  }
  output.cd();
  int schema_version = 3, source_schema_version = 4, signal_embedding_id = kSignalEmbeddingId;
  double majority_threshold = kMajorityThreshold;
  std::string majority_comparison = "strictly_greater_than";
  std::string eta_definition = "sum_signal_embedding_g4_eta_or_generator_eta_decay_photon_contributor_fraction";
  std::string other_definition = "no_unique_prompt_pi0_eta_strict_majority_including_exact_0p5";
  std::string weight_definition = "cross_section_pb_times_generator_weight_divided_by_full_sample_sum_generator_weight_processed";
  std::string metadata_family = family;
  std::string metadata_selection = selection;
  bool metadata_require_complete = require_complete;
  int metadata_n_bins = n_bins;
  double metadata_et_max = et_max;
  double sample_cross_section_pb = sample.cross_section_pb;
  TTree metadata("metadata", "Photon-candidate selection partial metadata");
  metadata.Branch("schema_version", &schema_version);
  metadata.Branch("source_map_schema_version", &source_schema_version);
  metadata.Branch("family", &metadata_family);
  metadata.Branch("selection", &metadata_selection);
  metadata.Branch("map_root", &normalized_map_root);
  metadata.Branch("require_complete", &metadata_require_complete);
  metadata.Branch("n_bins", &metadata_n_bins);
  metadata.Branch("et_max", &metadata_et_max);
  metadata.Branch("min_cluster_energy", &min_cluster_energy);
  metadata.Branch("partner_diagnostic_min_cluster_energy", &partner_diagnostic_min_cluster_energy);
  metadata.Branch("meson_partner_min_energy", &meson_partner_min_energy);
  metadata.Branch("pi0_topology_algorithm_version", &pi0_topology_algorithm_version);
  metadata.Branch("signal_embedding_id", &signal_embedding_id);
  metadata.Branch("majority_threshold", &majority_threshold);
  metadata.Branch("majority_comparison", &majority_comparison);
  metadata.Branch("eta_definition", &eta_definition);
  metadata.Branch("other_definition", &other_definition);
  metadata.Branch("weight_definition", &weight_definition);
  metadata.Branch("analysis_release", &analysis_release);
  metadata.Branch("model_sha256", &model_sha256);
  metadata.Branch("sample_names", &sample_names);
  metadata.Branch("sample_map_counts", &sample_map_counts);
  metadata.Branch("sample_sum_generator_weights", &sample_sum_generator_weights);
  metadata.Branch("sample_cross_section_pb", &sample_cross_section_pb);
  metadata.Branch("events_written", &events_written);
  metadata.Branch("events_stitch_pass", &events_stitch_pass);
  metadata.Branch("region_a_clusters", &region_a_clusters);
  metadata.Branch("region_a_prompt_clusters", &region_a_prompt_clusters);
  metadata.Branch("region_a_anchor_clusters", &region_a_anchor_clusters);
  metadata.Branch("selected_cluster_count", &selected_count);
  metadata.Branch("prompt_cluster_count", &prompt_count);
  metadata.Branch("pi0_cluster_count", &pi0_count);
  metadata.Branch("eta_cluster_count", &eta_count);
  metadata.Branch("other_cluster_count", &other_count);
  metadata.Branch("overlap_cluster_count", &overlap_count);
  metadata.Branch("half_boundary_cluster_count", &half_boundary_count);
  metadata.Branch("invalid_truth_cluster_count", &invalid_truth_count);
  metadata.Fill();
  metadata.Write();
  TTree topology_summary("topology_summary", "Per-selection anchor-topology partial summary");
  unsigned int topology_selection_index = 0;
  std::string topology_selection_key;
  unsigned long long topology_selection_selected = 0, topology_selection_anchor = 0;
  topology_summary.Branch("selection_index", &topology_selection_index);
  topology_summary.Branch("selection_key", &topology_selection_key);
  topology_summary.Branch("selected_clusters", &topology_selection_selected);
  topology_summary.Branch("selected_anchor_clusters", &topology_selection_anchor);
  for (topology_selection_index = 0; topology_selection_index < kSelectionCount; ++topology_selection_index)
  {
    topology_selection_key = kSelectionKeys[topology_selection_index];
    topology_selection_selected = topology_selected_clusters[topology_selection_index];
    topology_selection_anchor = topology_selected_anchor_clusters[topology_selection_index];
    topology_summary.Fill();
  }
  topology_summary.Write();
  std::string shard_sample_name = sample_name;
  int metadata_shard_index = shard_index, metadata_shard_count = shard_count;
  TTree shard_metadata("shard_metadata", "Photon-candidate selection partial shard metadata");
  shard_metadata.Branch("sample_name", &shard_sample_name);
  shard_metadata.Branch("shard_index", &metadata_shard_index);
  shard_metadata.Branch("shard_count", &metadata_shard_count);
  shard_metadata.Branch("map_index_begin", &shard_map_begin);
  shard_metadata.Branch("map_index_end", &shard_map_end);
  shard_metadata.Branch("total_map_count", &total_map_count);
  shard_metadata.Fill();
  shard_metadata.Write();
  output.Close();
  std::cout << "ReducePythiaPhotonCandidateSelection - family/selection/selected/prompt/pi0/eta/other/overlap/half/output = "
            << family << "/" << selection << "/" << selected_count << "/" << prompt_count << "/" << pi0_count << "/" << eta_count << "/"
            << other_count << "/" << overlap_count << "/" << half_boundary_count << "/" << resolved_output_base;
  std::cout << " (partial " << sample_name << " shard " << shard_index << "/" << shard_count << ", maps " << shard_map_begin << ":" << shard_map_end << ")";
  std::cout << std::endl;
  return overlap_count == 0ULL ? 0 : 9;
}
