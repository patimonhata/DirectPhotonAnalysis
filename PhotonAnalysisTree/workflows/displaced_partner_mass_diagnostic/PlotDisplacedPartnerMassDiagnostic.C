#include "../../macro/Utilities/sPhenixStyle.C"

#include <TCanvas.h>
#include <TChain.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

namespace
{
constexpr int kMapSchema = 4;
constexpr int kTopologyVersion = 8;
constexpr double kTopologyThreshold = 0.5;
constexpr double kDiagnosticFloor = 0.1;
constexpr double kPi0MassMin = 0.10;
constexpr double kPi0MassMax = 0.20;
constexpr int kMassBins = 240;
constexpr double kMassMax = 1.2;
constexpr std::array<const char*, 6> kSelectionKeys = {
    "inclusive", "kinematic", "preselection", "preselection_tight", "preselection_isolation", "region_a"};
constexpr std::array<const char*, 6> kSelectionLabels = {
    "All topology anchors", "Kinematic", "Pre-selection", "Pre-selection + TightBDT", "Pre-selection + Isolation", "Region A"};
constexpr std::array<double, 10> kAnchorEtEdges = {0, 5, 6, 8, 10, 15, 20, 35, 50, 100};
constexpr std::array<double, 11> kTruthPtEdges = {0, 3, 5, 6, 8, 10, 15, 20, 35, 50, 100};
constexpr std::array<int, 10> kColors = {kBlack, kBlue + 1, kRed + 1, kGreen + 2, kMagenta + 1,
                                         kOrange + 7, kCyan + 2, kViolet + 1, kAzure + 7, kPink + 7};

struct Sample
{
  const char* name;
  long long expected_end;
};

std::vector<Sample> samples(const std::string& family)
{
  if (family == "jet") return {{"jet3", 10000}, {"jet5", 10000}, {"jet8", 10000}, {"jet12", 100000},
                                {"jet20", 10000}, {"jet30", 10000}, {"jet40", 10000}};
  if (family == "photonjet") return {{"photonjet3", 10000}, {"photonjet5", 10000}, {"photonjet10", 10000}, {"photonjet20", 10000}};
  return {};
}

template <class T>
bool bind_branch(TTree& tree, const char* name, T* address)
{
  return tree.GetBranch(name) && tree.SetBranchAddress(name, address) >= 0;
}


template <class T>
bool bind_active(TTree& tree, const char* name, T* address)
{
  if (!tree.GetBranch(name)) return false;
  tree.SetBranchStatus(name, true);
  tree.AddBranchToCache(name, true);
  return tree.SetBranchAddress(name, address) >= 0;
}
bool close(double left, double right)
{
  return std::abs(left - right) <= 1e-11 * std::max({1.0, std::abs(left), std::abs(right)});
}

struct Input
{
  Sample sample;
  std::string pattern;
  unsigned long long maps = 0;
  double sumw = 0;
};

bool inspect(const Sample& sample, const std::string& pattern, bool require_complete, Input& input,
             std::string& common_release, std::string& common_model, double& common_tag_threshold)
{
  TChain tree("metadata");
  if (tree.Add(pattern.c_str()) <= 0) return false;
  int schema = -1, topology_version = -1;
  unsigned int chunk = 0;
  long long begin = -1, end = -1;
  double topology_threshold = -1, diagnostic_floor = -1, tag_threshold = -1, mass_min = -1, mass_max = -1, sumw = 0;
  std::string* sample_name = nullptr;
  std::string* release = nullptr;
  std::string* model = nullptr;
  const bool ok = bind_branch(tree, "schema_version", &schema) && bind_branch(tree, "sample_name", &sample_name) && bind_branch(tree, "analysis_release", &release) &&
      bind_branch(tree, "model_sha256", &model) && bind_branch(tree, "manifest_begin", &begin) && bind_branch(tree, "manifest_end", &end) && bind_branch(tree, "map_chunk_id", &chunk) &&
      bind_branch(tree, "sum_generator_weight_processed", &sumw) && bind_branch(tree, "min_cluster_energy", &topology_threshold) &&
      bind_branch(tree, "partner_diagnostic_min_cluster_energy", &diagnostic_floor) && bind_branch(tree, "meson_partner_min_energy", &tag_threshold) &&
      bind_branch(tree, "pi0_mass_min", &mass_min) && bind_branch(tree, "pi0_mass_max", &mass_max) && bind_branch(tree, "pi0_topology_algorithm_version", &topology_version);
  if (!ok) return false;
  long long expected_begin = 0;
  for (Long64_t entry = 0; entry < tree.GetEntries(); ++entry)
  {
    if (tree.GetEntry(entry) <= 0 || !sample_name || !release || !model || schema != kMapSchema || *sample_name != sample.name || begin != expected_begin ||
        end <= begin || chunk != static_cast<unsigned int>(entry) || !close(topology_threshold, kTopologyThreshold) || !close(diagnostic_floor, kDiagnosticFloor) ||
        !close(mass_min, kPi0MassMin) || !close(mass_max, kPi0MassMax) || topology_version != kTopologyVersion || !std::isfinite(sumw)) return false;
    if (common_release.empty())
    {
      common_release = *release;
      common_model = *model;
      common_tag_threshold = tag_threshold;
    }
    else if (common_release != *release || common_model != *model || !close(common_tag_threshold, tag_threshold)) return false;
    input.sumw += sumw;
    expected_begin = end;
  }
  if (expected_begin > sample.expected_end || (require_complete && expected_begin != sample.expected_end) || !(input.sumw > 0)) return false;
  input.sample = sample;
  input.pattern = pattern;
  input.maps = tree.GetEntries();
  return true;
}

struct Histograms
{
  std::unique_ptr<TH1D> mass_count;
  std::unique_ptr<TH1D> mass_pb;
  std::unique_ptr<TH1D> partner_e_count;
  std::unique_ptr<TH1D> partner_e_pb;
  std::unique_ptr<TH2D> mass_vs_anchor_et_count;
  std::unique_ptr<TH2D> mass_vs_anchor_et_pb;
  std::unique_ptr<TH2D> mass_vs_truth_pt_count;
  std::unique_ptr<TH2D> mass_vs_truth_pt_pb;

  explicit Histograms(const std::string& key)
    : mass_count(std::make_unique<TH1D>(("h_" + key + "_mass_count").c_str(), "", kMassBins, 0, kMassMax)),
      mass_pb(std::make_unique<TH1D>(("h_" + key + "_mass_pb").c_str(), "", kMassBins, 0, kMassMax)),
      partner_e_count(std::make_unique<TH1D>(("h_" + key + "_partner_e_count").c_str(), "", 120, 0, 0.6)),
      partner_e_pb(std::make_unique<TH1D>(("h_" + key + "_partner_e_pb").c_str(), "", 120, 0, 0.6)),
      mass_vs_anchor_et_count(std::make_unique<TH2D>(("h_" + key + "_mass_vs_anchor_et_count").c_str(), "", kAnchorEtEdges.size() - 1,
                                                     kAnchorEtEdges.data(), kMassBins, 0, kMassMax)),
      mass_vs_anchor_et_pb(std::make_unique<TH2D>(("h_" + key + "_mass_vs_anchor_et_pb").c_str(), "", kAnchorEtEdges.size() - 1,
                                                  kAnchorEtEdges.data(), kMassBins, 0, kMassMax)),
      mass_vs_truth_pt_count(std::make_unique<TH2D>(("h_" + key + "_mass_vs_truth_pi0_pt_count").c_str(), "", kTruthPtEdges.size() - 1,
                                                    kTruthPtEdges.data(), kMassBins, 0, kMassMax)),
      mass_vs_truth_pt_pb(std::make_unique<TH2D>(("h_" + key + "_mass_vs_truth_pi0_pt_pb").c_str(), "", kTruthPtEdges.size() - 1,
                                                 kTruthPtEdges.data(), kMassBins, 0, kMassMax))
  {
    for (TH1* histogram : {static_cast<TH1*>(mass_count.get()), static_cast<TH1*>(mass_pb.get()), static_cast<TH1*>(partner_e_count.get()), static_cast<TH1*>(partner_e_pb.get()), static_cast<TH1*>(mass_vs_anchor_et_count.get()),
                           static_cast<TH1*>(mass_vs_anchor_et_pb.get()), static_cast<TH1*>(mass_vs_truth_pt_count.get()), static_cast<TH1*>(mass_vs_truth_pt_pb.get())}) histogram->Sumw2();
  }

  void fill(double mass, double partner_e, double anchor_et, double truth_pt, double weight)
  {
    mass_count->Fill(mass);
    mass_pb->Fill(mass, weight);
    partner_e_count->Fill(partner_e);
    partner_e_pb->Fill(partner_e, weight);
    mass_vs_anchor_et_count->Fill(anchor_et, mass);
    mass_vs_anchor_et_pb->Fill(anchor_et, mass, weight);
    mass_vs_truth_pt_count->Fill(truth_pt, mass);
    mass_vs_truth_pt_pb->Fill(truth_pt, mass, weight);
  }

  bool add(TDirectory* directory)
  {
    if (!directory) return false;
    const std::array<TH1*, 8> targets = {mass_count.get(), mass_pb.get(), partner_e_count.get(), partner_e_pb.get(), mass_vs_anchor_et_count.get(),
                                         mass_vs_anchor_et_pb.get(), mass_vs_truth_pt_count.get(), mass_vs_truth_pt_pb.get()};
    for (TH1* target : targets)
    {
      auto* source = directory->Get<TH1>(target->GetName());
      if (!source || !target->Add(source)) return false;
    }
    return true;
  }

  void write(const std::string& key)
  {
    for (TH1* histogram : {static_cast<TH1*>(mass_count.get()), static_cast<TH1*>(mass_pb.get()), static_cast<TH1*>(partner_e_count.get()), static_cast<TH1*>(partner_e_pb.get()), static_cast<TH1*>(mass_vs_anchor_et_count.get()),
                           static_cast<TH1*>(mass_vs_anchor_et_pb.get()), static_cast<TH1*>(mass_vs_truth_pt_count.get()), static_cast<TH1*>(mass_vs_truth_pt_pb.get())}) histogram->Write();
    for (int bin = 1; bin <= mass_vs_anchor_et_count->GetNbinsX(); ++bin)
    {
      std::unique_ptr<TH1D> count(mass_vs_anchor_et_count->ProjectionY(("h_" + key + "_anchor_et_bin" + std::to_string(bin) + "_mass_count").c_str(), bin, bin));
      std::unique_ptr<TH1D> weighted(mass_vs_anchor_et_pb->ProjectionY(("h_" + key + "_anchor_et_bin" + std::to_string(bin) + "_mass_pb").c_str(), bin, bin));
      count->Write();
      weighted->Write();
    }
    for (int bin = 1; bin <= mass_vs_truth_pt_count->GetNbinsX(); ++bin)
    {
      std::unique_ptr<TH1D> count(mass_vs_truth_pt_count->ProjectionY(("h_" + key + "_truth_pi0_pt_bin" + std::to_string(bin) + "_mass_count").c_str(), bin, bin));
      std::unique_ptr<TH1D> weighted(mass_vs_truth_pt_pb->ProjectionY(("h_" + key + "_truth_pi0_pt_bin" + std::to_string(bin) + "_mass_pb").c_str(), bin, bin));
      count->Write();
      weighted->Write();
    }
  }
};

void lines(double maximum)
{
  for (double x : {kPi0MassMin, kPi0MassMax})
  {
    TLine line(x, 0, x, maximum);
    line.SetLineColor(kRed + 1);
    line.SetLineStyle(2);
    line.SetLineWidth(2);
    line.DrawClone();
  }
}

void label(const std::string& text_value)
{
  TLatex text;
  text.SetNDC();
  text.SetTextSize(0.031);
  text.DrawLatex(0.15, 0.92, "#it{#bf{sPHENIX}} Internal");
  text.DrawLatex(0.15, 0.87, text_value.c_str());
}

void draw_mass(TH1D* histogram, const std::string& path, const std::string& text)
{
  TCanvas canvas("c_mass", "", 1000, 800);
  histogram->SetStats(false);
  histogram->SetLineColor(kBlue + 1);
  histogram->SetLineWidth(2);
  histogram->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
  histogram->GetYaxis()->SetTitle("Cross section [pb / 5 MeV]");
  histogram->Draw("HIST");
  lines(std::max(1.0, 1.05 * histogram->GetMaximum()));
  label(text);
  canvas.SaveAs(path.c_str());
}

void draw_2d(TH2D* histogram, const std::string& path, const char* x_title, const std::string& text)
{
  TCanvas canvas(("c_" + std::string(histogram->GetName())).c_str(), "", 1000, 800);
  canvas.SetRightMargin(0.15);
  histogram->SetStats(false);
  histogram->GetXaxis()->SetTitle(x_title);
  histogram->GetYaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
  histogram->GetZaxis()->SetTitle("Cross section [pb]");
  histogram->Draw("COLZ");
  for (double y : {kPi0MassMin, kPi0MassMax})
  {
    TLine line(histogram->GetXaxis()->GetXmin(), y, histogram->GetXaxis()->GetXmax(), y);
    line.SetLineColor(kRed + 1);
    line.SetLineStyle(2);
    line.SetLineWidth(2);
    line.DrawClone();
  }
  label(text);
  canvas.SaveAs(path.c_str());
}

template <std::size_t N>
void draw_shapes(TH2D* source, const std::array<double, N>& edges, const std::string& path, const std::string& text)
{
  TCanvas canvas(("c_shape_" + std::string(source->GetName())).c_str(), "", 1100, 850);
  TLegend legend(0.54, 0.44, 0.94, 0.90);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetTextSize(0.024);
  std::vector<std::unique_ptr<TH1D>> shapes;
  double maximum = 0;
  for (int bin = 1; bin <= source->GetNbinsX(); ++bin)
  {
    auto shape = std::unique_ptr<TH1D>(source->ProjectionY(("shape_" + std::to_string(bin)).c_str(), bin, bin));
    const double entries = shape->Integral(0, shape->GetNbinsX() + 1);
    if (entries > 0) shape->Scale(1.0 / entries);
    shape->SetLineColor(kColors[(bin - 1) % kColors.size()]);
    shape->SetLineWidth(2);
    maximum = std::max(maximum, shape->GetMaximum());
    std::ostringstream item;
    item << std::fixed << std::setprecision(0) << edges[bin - 1] << "-" << edges[bin] << " GeV (N=" << entries << ")";
    legend.AddEntry(shape.get(), item.str().c_str(), "l");
    shapes.push_back(std::move(shape));
  }
  shapes.front()->SetStats(false);
  shapes.front()->SetMaximum(std::max(0.02, 1.15 * maximum));
  shapes.front()->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
  shapes.front()->GetYaxis()->SetTitle("Unit-normalized pairs");
  shapes.front()->Draw("HIST");
  for (std::size_t index = 1; index < shapes.size(); ++index) shapes[index]->Draw("HIST SAME");
  lines(shapes.front()->GetMaximum());
  legend.Draw();
  label(text);
  canvas.SaveAs(path.c_str());
}
}

int PlotDisplacedPartnerMassDiagnostic(
    const std::string family = "jet",
    const std::string map_root = "/sphenix/tg/tg01/coldqcd/ryotaro/DirectPhotonAnalysis/photon_candidate_selection/cluster_e_gt_0p5",
    const std::string output_base = "",
    const bool require_complete = true,
    const double tagging_partner_min_energy = 0.2,
    const std::string sample_filter = "",
    const Long64_t max_events_per_sample = 0,
    const int shard_index = 0,
    const int shard_count = 1,
    const bool make_plots = true)
{
  auto definitions = samples(family);
  if (definitions.empty() || map_root.empty() || !std::isfinite(tagging_partner_min_energy) || !(tagging_partner_min_energy > kDiagnosticFloor) ||
      !(tagging_partner_min_energy < kTopologyThreshold) || max_events_per_sample < 0 || shard_index < 0 || shard_count <= 0 || shard_index >= shard_count ||
      (shard_count != 1 && (sample_filter.empty() || max_events_per_sample != 0))) return 1;
  if (!sample_filter.empty())
  {
    definitions.erase(std::remove_if(definitions.begin(), definitions.end(), [&](const auto& sample) { return sample_filter != sample.name; }), definitions.end());
    if (definitions.empty()) return 1;
  }
  std::string root = map_root;
  while (root.size() > 1 && root.back() == '/') root.pop_back();
  const std::string configuration = root.substr(root.find_last_of('/') + 1);
  const std::string output_directory = output_base.empty()
      ? "/sphenix/user/ryotaro/DirectPhotonAnalysis/PhotonAnalysisTree/output/plots/displaced_partner_mass_diagnostic/" + configuration + "/" + family
      : output_base;
  if (gSystem->mkdir(output_directory.c_str(), true) != 0 && gSystem->AccessPathName(output_directory.c_str())) return 2;

  std::vector<Input> inputs;
  std::string release, model;
  double production_tag_threshold = -1;
  for (const auto& definition : definitions)
  {
    Input input;
    if (!inspect(definition, root + "/" + definition.name + "/map_*.root", require_complete, input, release, model, production_tag_threshold)) return 3;
    inputs.push_back(input);
  }

  gROOT->cd();
  std::array<std::unique_ptr<Histograms>, kSelectionKeys.size()> histograms;
  for (std::size_t selection = 0; selection < histograms.size(); ++selection) histograms[selection] = std::make_unique<Histograms>(kSelectionKeys[selection]);
  std::array<unsigned long long, kSelectionKeys.size()> counts = {}, window_counts = {};
  std::array<double, kSelectionKeys.size()> yields = {}, window_yields = {};
  unsigned long long events_considered = 0;
  for (const auto& input : inputs)
  {
    TChain tree("event_tree");
    if (tree.Add(input.pattern.c_str()) <= 0) return 4;
    tree.SetBranchStatus("*", false);
    tree.SetCacheSize(64LL * 1024LL * 1024LL);
    unsigned char weight_valid = 0, stitch_valid = 0, stitch_pass = 0;
    double weight_numerator = 0;
    unsigned int ncluster = 0;
    std::vector<double>* et = nullptr;
    std::vector<unsigned char>* kinematic = nullptr;
    std::vector<unsigned char>* preselection = nullptr;
    std::vector<unsigned char>* tight = nullptr;
    std::vector<unsigned char>* isolated = nullptr;
    std::vector<unsigned char>* region_a = nullptr;
    std::vector<unsigned char>* anchor_valid = nullptr;
    std::vector<int>* candidate_index = nullptr;
    std::vector<int>* topology = nullptr;
    std::vector<int>* missing_category = nullptr;
    std::vector<int>* alignment = nullptr;
    std::vector<int>* tag_status = nullptr;
    std::vector<float>* diagnostic_mass = nullptr;
    std::vector<float>* truth_mass = nullptr;
    std::vector<float>* partner_e = nullptr;
    std::vector<float>* candidate_pt = nullptr;
    const bool ok = bind_active(tree, "event_weight_valid", &weight_valid) && bind_active(tree, "sample_stitching_valid", &stitch_valid) &&
        bind_active(tree, "sample_stitching_pass", &stitch_pass) && bind_active(tree, "weight_numerator_pb", &weight_numerator) &&
        bind_active(tree, "split_ncluster", &ncluster) && bind_active(tree, "split_cluster_et", &et) &&
        bind_active(tree, "split_cluster_pass_kinematics", &kinematic) && bind_active(tree, "split_cluster_pass_preselection", &preselection) &&
        bind_active(tree, "split_cluster_pass_tight", &tight) && bind_active(tree, "split_cluster_pass_isolated", &isolated) &&
        bind_active(tree, "split_cluster_pass_region_a", &region_a) && bind_active(tree, "split_cluster_pi0_anchor_valid", &anchor_valid) &&
        bind_active(tree, "split_cluster_pi0_anchor_candidate_index", &candidate_index) && bind_active(tree, "split_cluster_pi0_anchor_topology", &topology) &&
        bind_active(tree, "split_cluster_pi0_anchor_missing_category", &missing_category) &&
        bind_active(tree, "split_cluster_pi0_anchor_partner_alignment", &alignment) &&
        bind_active(tree, "split_cluster_pi0_anchor_truth_partner_tag_status", &tag_status) &&
        bind_active(tree, "split_cluster_pi0_anchor_partner_diagnostic_mass", &diagnostic_mass) &&
        bind_active(tree, "split_cluster_pi0_anchor_truth_partner_mass", &truth_mass) &&
        bind_active(tree, "split_cluster_pi0_anchor_truth_partner_cluster_e", &partner_e) && bind_active(tree, "pi0_topology_candidate_pt", &candidate_pt);
    if (!ok) return 4;
    const Long64_t total_entries = tree.GetEntries();
    const Long64_t entry_begin = total_entries * shard_index / shard_count;
    Long64_t entry_end = total_entries * (shard_index + 1) / shard_count;
    if (max_events_per_sample > 0) entry_end = std::min(entry_end, entry_begin + max_events_per_sample);
    events_considered += entry_end - entry_begin;
    for (Long64_t entry = entry_begin; entry < entry_end; ++entry)
    {
      if (tree.GetEntry(entry) <= 0) return 5;
      const std::array<std::size_t, 15> sizes = {et ? et->size() : 0, kinematic ? kinematic->size() : 0, preselection ? preselection->size() : 0,
          tight ? tight->size() : 0, isolated ? isolated->size() : 0, region_a ? region_a->size() : 0, anchor_valid ? anchor_valid->size() : 0,
          candidate_index ? candidate_index->size() : 0, topology ? topology->size() : 0, missing_category ? missing_category->size() : 0,
          alignment ? alignment->size() : 0, tag_status ? tag_status->size() : 0, diagnostic_mass ? diagnostic_mass->size() : 0,
          truth_mass ? truth_mass->size() : 0, partner_e ? partner_e->size() : 0};
      if (!candidate_pt || std::any_of(sizes.begin(), sizes.end(), [&](std::size_t size) { return size != ncluster; })) return 5;
      if (!weight_valid || !stitch_valid || !stitch_pass) continue;
      const double weight = weight_numerator / input.sumw;
      if (!std::isfinite(weight)) return 5;
      for (std::size_t cluster = 0; cluster < ncluster; ++cluster)
      {
        if (!(*anchor_valid)[cluster] || (*topology)[cluster] != 3 || (*missing_category)[cluster] != 4 || (*alignment)[cluster] != 2 ||
            (*tag_status)[cluster] != 4 || !((*partner_e)[cluster] > tagging_partner_min_energy)) continue;
        const int candidate = (*candidate_index)[cluster];
        const double mass = (*diagnostic_mass)[cluster];
        if (candidate < 0 || static_cast<std::size_t>(candidate) >= candidate_pt->size() || !std::isfinite(mass) || mass < 0 ||
            !std::isfinite((*truth_mass)[cluster]) || std::abs(mass - (*truth_mass)[cluster]) > 1e-5) return 5;
        const std::array<bool, 6> pass = {true, (*kinematic)[cluster] != 0, (*kinematic)[cluster] && (*preselection)[cluster],
            (*kinematic)[cluster] && (*preselection)[cluster] && (*tight)[cluster],
            (*kinematic)[cluster] && (*preselection)[cluster] && (*isolated)[cluster], (*region_a)[cluster] != 0};
        const bool in_window = mass > kPi0MassMin && mass < kPi0MassMax;
        for (std::size_t selection = 0; selection < histograms.size(); ++selection)
        {
          if (!pass[selection]) continue;
          histograms[selection]->fill(mass, (*partner_e)[cluster], (*et)[cluster], (*candidate_pt)[candidate], weight);
          ++counts[selection];
          if (in_window) ++window_counts[selection];
        }
      }
    }
  }
  for (std::size_t selection = 0; selection < histograms.size(); ++selection)
  {
    yields[selection] = histograms[selection]->mass_pb->Integral(0, kMassBins + 1);
    const int low = histograms[selection]->mass_pb->FindBin(std::nextafter(kPi0MassMin, kPi0MassMax));
    const int high = histograms[selection]->mass_pb->FindBin(std::nextafter(kPi0MassMax, kPi0MassMin));
    window_yields[selection] = histograms[selection]->mass_pb->Integral(low, high);
  }

  TFile output((output_directory + "/displaced_partner_mass_diagnostic.root").c_str(), "RECREATE");
  if (output.IsZombie()) return 6;
  for (std::size_t selection = 0; selection < histograms.size(); ++selection)
  {
    output.mkdir(kSelectionKeys[selection])->cd();
    histograms[selection]->write(kSelectionKeys[selection]);
  }
  output.cd();
  int schema_version = 1, source_schema = kMapSchema;
  double topology_threshold = kTopologyThreshold, diagnostic_floor = kDiagnosticFloor, mass_min = kPi0MassMin, mass_max = kPi0MassMax;
  bool complete = require_complete;
  int metadata_shard_index = shard_index, metadata_shard_count = shard_count;
  Long64_t metadata_total_entries = 0, metadata_entry_begin = 0, metadata_entry_end = 0;
  if (!inputs.empty())
  {
    TChain coverage("event_tree");
    coverage.Add(inputs.front().pattern.c_str());
    metadata_total_entries = coverage.GetEntries();
    metadata_entry_begin = metadata_total_entries * shard_index / shard_count;
    metadata_entry_end = metadata_total_entries * (shard_index + 1) / shard_count;
    if (max_events_per_sample > 0) metadata_entry_end = std::min(metadata_entry_end, metadata_entry_begin + max_events_per_sample);
  }
  Long64_t metadata_max_events = max_events_per_sample;
  double emulated_tag_threshold = tagging_partner_min_energy;
  std::string metadata_family = family, metadata_filter = sample_filter;
  std::vector<std::string> sample_names;
  std::vector<unsigned long long> map_counts;
  std::vector<double> sum_weights;
  for (const auto& input : inputs)
  {
    sample_names.push_back(input.sample.name);
    map_counts.push_back(input.maps);
    sum_weights.push_back(input.sumw);
  }
  TTree metadata("metadata", "Displaced-partner mass diagnostic metadata");
  metadata.Branch("schema_version", &schema_version);
  metadata.Branch("source_map_schema_version", &source_schema);
  metadata.Branch("family", &metadata_family);
  metadata.Branch("sample_filter", &metadata_filter);
  metadata.Branch("require_complete", &complete);
  metadata.Branch("max_events_per_sample", &metadata_max_events);
  metadata.Branch("shard_index", &metadata_shard_index);
  metadata.Branch("total_entries", &metadata_total_entries);
  metadata.Branch("shard_count", &metadata_shard_count);
  metadata.Branch("entry_begin", &metadata_entry_begin);
  metadata.Branch("entry_end", &metadata_entry_end);
  metadata.Branch("topology_min_cluster_energy", &topology_threshold);
  metadata.Branch("production_tagging_partner_min_energy", &production_tag_threshold);
  metadata.Branch("diagnostic_min_cluster_energy", &diagnostic_floor);
  metadata.Branch("emulated_tagging_partner_min_energy", &emulated_tag_threshold);
  metadata.Branch("pi0_mass_min", &mass_min);
  metadata.Branch("pi0_mass_max", &mass_max);
  metadata.Branch("analysis_release", &release);
  metadata.Branch("model_sha256", &model);
  metadata.Branch("sample_names", &sample_names);
  metadata.Branch("sample_map_counts", &map_counts);
  metadata.Branch("sample_sum_generator_weights", &sum_weights);
  metadata.Branch("events_considered", &events_considered);
  metadata.Fill();
  metadata.Write();
  TTree summary("summary", "Per-selection displaced-pair summary");
  unsigned int index = 0;
  std::string key;
  unsigned long long selected_pairs = 0, pi0_window_pairs = 0;
  double selected_pb = 0, pi0_window_pb = 0, raw_fraction = 0, weighted_fraction = 0;
  summary.Branch("selection_index", &index);
  summary.Branch("selection_key", &key);
  summary.Branch("selected_pairs", &selected_pairs);
  summary.Branch("pi0_window_pairs", &pi0_window_pairs);
  summary.Branch("selected_cross_section_pb", &selected_pb);
  summary.Branch("pi0_window_cross_section_pb", &pi0_window_pb);
  summary.Branch("raw_pi0_window_fraction", &raw_fraction);
  summary.Branch("weighted_pi0_window_fraction", &weighted_fraction);
  for (index = 0; index < histograms.size(); ++index)
  {
    key = kSelectionKeys[index];
    selected_pairs = counts[index];
    pi0_window_pairs = window_counts[index];
    selected_pb = yields[index];
    pi0_window_pb = window_yields[index];
    raw_fraction = selected_pairs ? static_cast<double>(pi0_window_pairs) / selected_pairs : 0;
    weighted_fraction = selected_pb != 0 ? pi0_window_pb / selected_pb : 0;
    summary.Fill();
  }
  summary.Write();
  output.Close();

  if (!make_plots)
  {
    std::cout << "ReduceDisplacedPartnerMassDiagnostic - sample/shard/pairs/output = " << sample_filter << "/" << shard_index << "/" << shard_count
              << "/" << counts[0] << "/" << output_directory << std::endl;
    return 0;
  }
  SetsPhenixStyle();
  gStyle->SetOptStat(0);
  for (std::size_t selection = 0; selection < histograms.size(); ++selection)
  {
    const std::string directory = output_directory + "/" + kSelectionKeys[selection];
    gSystem->mkdir(directory.c_str(), true);
    std::ostringstream text;
    text << (family == "jet" ? "Pythia8 p+p Jet, " : "Pythia8 p+p PhotonJet, ") << kSelectionLabels[selection]
         << ", E_{partner} > " << tagging_partner_min_energy << " GeV";
    draw_mass(histograms[selection]->mass_pb.get(), directory + "/displaced_partner_mass.pdf", text.str());
    draw_2d(histograms[selection]->mass_vs_anchor_et_pb.get(), directory + "/displaced_partner_mass_vs_anchor_et.pdf", "Anchor cluster E_{T} [GeV]", text.str());
    draw_2d(histograms[selection]->mass_vs_truth_pt_pb.get(), directory + "/displaced_partner_mass_vs_truth_pi0_pt.pdf", "Truth #pi^{0} p_{T} [GeV]", text.str());
    draw_shapes(histograms[selection]->mass_vs_anchor_et_count.get(), kAnchorEtEdges, directory + "/displaced_partner_mass_by_anchor_et_shape.pdf", text.str());
    draw_shapes(histograms[selection]->mass_vs_truth_pt_count.get(), kTruthPtEdges, directory + "/displaced_partner_mass_by_truth_pi0_pt_shape.pdf", text.str());
  }
  std::cout << "PlotDisplacedPartnerMassDiagnostic - samples/maps/pairs/output = " << inputs.size() << "/"
            << std::accumulate(map_counts.begin(), map_counts.end(), 0ULL) << "/" << counts[0] << "/" << output_directory << std::endl;
  return 0;
}
