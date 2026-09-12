#ifndef TAGGING_MASS_DIAGNOSTIC_H
#define TAGGING_MASS_DIAGNOSTIC_H
#include "../photon_candidate_selection/PhotonCandidateSettings.h"
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
#include <TLeaf.h>
#include <set>

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

namespace tagging_mass
{
constexpr int kMapSchema = 5;
constexpr int kTopologyVersion = 10;

constexpr double kDiagnosticFloor = 0.0;
double kPi0MassMin = 0.10;
double kPi0MassMax = 0.20;
constexpr int kMassBins = 240;
constexpr double kMassMax = 1.2;
constexpr std::array<const char*, 6> kSelectionKeys = {
    "separated_truth_pair", "separated_truth_pair_below_threshold", "prompt_pi0_all_pairs", "prompt_pi0_window_pairs", "prompt_eta_all_pairs", "prompt_eta_window_pairs"};
constexpr std::array<const char*, 6> kSelectionLabels = {
    "Separated truth pair", "Truth partner below tagging threshold", "Prompt: all #pi^{0} partners", "Prompt: #pi^{0} window pairs", "Prompt: all #eta partners", "Prompt: #eta window pairs"};
constexpr std::array<double, 10> kAnchorEtEdges = {0, 5, 6, 8, 10, 15, 20, 35, 50, 100};
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
             std::string& common_release, std::string& common_model, double& common_tag_threshold, std::vector<double>& common_settings, std::vector<double>& region_settings)
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
      bind_branch(tree, "partner_diagnostic_min_cluster_energy", &diagnostic_floor) && bind_branch(tree, "pi0_partner_min_energy", &tag_threshold) &&
      bind_branch(tree, "pi0_mass_min", &mass_min) && bind_branch(tree, "pi0_mass_max", &mass_max) && bind_branch(tree, "pi0_topology_algorithm_version", &topology_version);
  if (!ok) return false;
  const char* region_names[] = {"shower_shape_min_tower_energy", "candidate_et_min", "candidate_et_max", "candidate_abs_eta_max",
      "max_abs_vertex_z", "isolation_radius", "isolation_scale", "isolation_offset", "nonisolation_gap"};
  std::vector<double> region_values(9);
  for (std::size_t i = 0; i < region_values.size(); ++i)
    if (!bind_branch(tree, region_names[i], &region_values[i])) return false;
  std::vector<double> sample_settings;
  long long expected_begin = 0;
  for (Long64_t entry = 0; entry < tree.GetEntries(); ++entry)
  {
    if (tree.GetEntry(entry) <= 0 || !sample_name || !release || !model || schema != kMapSchema || *sample_name != sample.name || begin != expected_begin ||
        end <= begin || chunk != static_cast<unsigned int>(entry) || !close(diagnostic_floor, kDiagnosticFloor) ||
        topology_version != kTopologyVersion || !std::isfinite(sumw)) return false;
    if (region_settings.empty()) region_settings = region_values;
    if (region_settings.size() != region_values.size()) return false;
    for (std::size_t i = 0; i < region_values.size(); ++i)
      if (!std::isfinite(region_values[i]) || !close(region_values[i], region_settings[i])) return false;
    // Check cross section and stitching bounds within each sample as well as full-sample sumw.
    std::vector<double> current_sample;
    for (const char* name : {"sample_cross_section_pb", "sample_window_min", "sample_window_max", "sample_upper_unbounded"})
    {
      auto* leaf = tree.GetLeaf(name);
      if (!leaf || !std::isfinite(leaf->GetValue())) return false;
      current_sample.push_back(leaf->GetValue());
    }
    if (sample_settings.empty()) sample_settings = current_sample;
    for (std::size_t i = 0; i < current_sample.size(); ++i)
      if (!close(current_sample[i], sample_settings[i])) return false;
    std::vector<double> settings;
    if (!photon_candidate_settings::read(tree.GetFile()->GetName(), settings)) return false;
    if (common_settings.empty()) common_settings = settings;
    else if (!photon_candidate_settings::same(common_settings, settings)) return false;
    if (settings[1] < settings[0] || settings[2] < settings[0])
    {
      std::cerr << "Map does not store every eligible tagging partner: " << pattern << std::endl;
      return false;
    }
    kPi0MassMin = settings[3];
    kPi0MassMax = settings[4];
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

// Group histograms are pair-weighted; flow histograms below are candidate-weighted.
struct Histograms
{
  std::unique_ptr<TH1D> mass_count, mass_pb;
  std::unique_ptr<TH2D> mass_vs_anchor_et_count, mass_vs_anchor_et_pb, window_count, window_pb;
  explicit Histograms(const std::string& key)
  {
    mass_count = std::make_unique<TH1D>(("h_" + key + "_mass_count").c_str(), "", kMassBins, 0, kMassMax);
    mass_pb = std::make_unique<TH1D>(("h_" + key + "_mass_pb").c_str(), "", kMassBins, 0, kMassMax);
    mass_vs_anchor_et_count = std::make_unique<TH2D>(("h_" + key + "_mass_vs_anchor_et_count").c_str(), "", 9, kAnchorEtEdges.data(), kMassBins, 0, kMassMax);
    mass_vs_anchor_et_pb = std::make_unique<TH2D>(("h_" + key + "_mass_vs_anchor_et_pb").c_str(), "", 9, kAnchorEtEdges.data(), kMassBins, 0, kMassMax);
    window_count = std::make_unique<TH2D>(("h_" + key + "_window_count").c_str(), "", 9, kAnchorEtEdges.data(), 2, -0.5, 1.5);
    window_pb = std::make_unique<TH2D>(("h_" + key + "_window_pb").c_str(), "", 9, kAnchorEtEdges.data(), 2, -0.5, 1.5);
    for (auto* h : objects()) { h->SetDirectory(nullptr); h->Sumw2(); }
  }
  std::array<TH1*, 6> objects() const { return {mass_count.get(), mass_pb.get(), mass_vs_anchor_et_count.get(), mass_vs_anchor_et_pb.get(), window_count.get(), window_pb.get()}; }
  void fill(double mass, double et, double weight, double low, double high)
  {
    mass_count->Fill(mass); mass_pb->Fill(mass, weight);
    mass_vs_anchor_et_count->Fill(et, mass); mass_vs_anchor_et_pb->Fill(et, mass, weight);
    const bool inside = mass > low && mass < high;
    window_count->Fill(et, inside); window_pb->Fill(et, inside, weight);
  }
  static bool add_checked(TH1* target, TH1* source)
  {
    if (!source || source->IsA() != target->IsA()) return false;
    for (auto axes : {std::make_pair(target->GetXaxis(), source->GetXaxis()), std::make_pair(target->GetYaxis(), source->GetYaxis())})
    {
      if (axes.first->GetNbins() != axes.second->GetNbins()) return false;
      for (int i = 1; i <= axes.first->GetNbins() + 1; ++i)
        if (!close(axes.first->GetBinLowEdge(i), axes.second->GetBinLowEdge(i))) return false;
    }
    return target->Add(source);
  }
  bool add(TDirectory* directory)
  {
    if (!directory) return false;
    for (auto* h : objects()) if (!add_checked(h, directory->Get<TH1>(h->GetName()))) return false;
    return true;
  }
  void write()
  {
    for (auto* h : objects()) h->Write();
    for (int bin = 1; bin <= 9; ++bin)
    {
      for (auto* h : {mass_vs_anchor_et_count.get(), mass_vs_anchor_et_pb.get()})
      {
        std::unique_ptr<TH1D> projection(h->ProjectionY((std::string(h->GetName()) + "_et_bin" + std::to_string(bin)).c_str(), bin, bin));
        projection->SetDirectory(nullptr); projection->Write();
      }
    }
    // Exact strict-window classification, independent of mass bin edges and overflow.
    for (auto* h : {window_count.get(), window_pb.get()})
    {
      TH1D fraction((std::string(h->GetName()) + "_fraction").c_str(), "", 9, kAnchorEtEdges.data());
      fraction.SetDirectory(nullptr);
      for (int bin = 0; bin <= 10; ++bin)
      {
        const double no = h->GetBinContent(bin, 1), yes = h->GetBinContent(bin, 2), total = no + yes;
        if (total == 0) continue;
        const double value = yes / total;
        fraction.SetBinContent(bin, value);
        fraction.SetBinError(bin, std::sqrt(std::pow((1 - value) * h->GetBinError(bin, 2), 2) + std::pow(value * h->GetBinError(bin, 1), 2)) / std::abs(total));
      }
      fraction.Write();
    }
  }
};

struct Flow
{
  TH2D count{"candidate_flow_count", "", 9, kAnchorEtEdges.data(), 8, -0.5, 7.5};
  TH2D weighted{"candidate_flow_pb", "", 9, kAnchorEtEdges.data(), 8, -0.5, 7.5};
  Flow()
  {
    const char* labels[] = {"separated", "truth_pair_valid", "truth_partner_below_threshold", "truth_pair_invalid", "prompt", "prompt_pi0_veto", "prompt_eta_veto", "prompt_any_veto"};
    for (auto* h : {&count, &weighted})
    {
      h->SetDirectory(nullptr); h->Sumw2();
      for (int i = 0; i < 8; ++i) h->GetYaxis()->SetBinLabel(i + 1, labels[i]);
    }
  }
  void fill(int code, double et, double weight) { count.Fill(et, code); weighted.Fill(et, code, weight); }
  void write() { count.Write(); weighted.Write(); }
  bool add(TFile& file) { return Histograms::add_checked(&count, file.Get<TH1>(count.GetName())) && Histograms::add_checked(&weighted, file.Get<TH1>(weighted.GetName())); }
};

struct Metadata
{
  int schema_version = 1, source_map_schema_version = 5, pi0_topology_algorithm_version = 10;
  int shard_index = 0, shard_count = 1;
  Long64_t total_entries = 0, entry_begin = 0, entry_end = 0, max_events = 0;
  bool require_complete = true;
  unsigned long long map_count = 0;
  double sumw = 0;
  std::string family, sample, release, model;
  std::vector<double> selection_settings, region_settings;
  void write()
  {
    TTree tree("metadata", "Tagging mass diagnostic metadata");
#define B(x) tree.Branch(#x, &x)
    B(schema_version); B(source_map_schema_version); B(pi0_topology_algorithm_version); B(shard_index); B(shard_count);
    B(total_entries); B(entry_begin); B(entry_end); B(max_events); B(require_complete); B(map_count); B(sumw);
    B(family); B(sample); B(release); B(model); B(selection_settings); B(region_settings);
#undef B
    tree.Fill(); tree.Write();
  }
  bool read(TFile& file)
  {
    auto* t = file.Get<TTree>("metadata");
    if (!t || t->GetEntries() != 1) return false;
#define B(x) if (!bind_branch(*t, #x, &x)) return false
    B(schema_version); B(source_map_schema_version); B(pi0_topology_algorithm_version); B(shard_index); B(shard_count);
    B(total_entries); B(entry_begin); B(entry_end); B(max_events); B(require_complete); B(map_count); B(sumw);
#undef B
    std::string *f = nullptr, *s = nullptr, *r = nullptr, *m = nullptr;
    std::vector<double> *settings = nullptr, *region = nullptr;
    if (!bind_branch(*t, "family", &f) || !bind_branch(*t, "sample", &s) || !bind_branch(*t, "release", &r) || !bind_branch(*t, "model", &m) ||
        !bind_branch(*t, "selection_settings", &settings) || !bind_branch(*t, "region_settings", &region) || t->GetEntry(0) <= 0 || !f || !s || !r || !m || !settings || !region) return false;
    family = *f; sample = *s; release = *r; model = *m; selection_settings = *settings; region_settings = *region;
    t->ResetBranchAddresses();
    return schema_version == 1 && source_map_schema_version == 5 && pi0_topology_algorithm_version == 10 && selection_settings.size() == 11 && region_settings.size() == 9 &&
        photon_candidate_settings::same(selection_settings, selection_settings) && selection_settings[1] >= selection_settings[0] && selection_settings[2] >= selection_settings[0] &&
        selection_settings[3] < selection_settings[4] && selection_settings[5] < selection_settings[6] && selection_settings[9] == 0 && selection_settings[10] == 0 && std::isfinite(sumw) && sumw > 0 && map_count > 0;
  }
  bool compatible(const Metadata& other) const
  {
    if (family != other.family || release != other.release || model != other.model || !photon_candidate_settings::same(selection_settings, other.selection_settings) || region_settings.size() != other.region_settings.size()) return false;
    for (std::size_t i = 0; i < region_settings.size(); ++i)
      if (!std::isfinite(region_settings[i]) || !close(region_settings[i], other.region_settings[i])) return false;
    return true;
  }
};

inline double pair_mass(double e, double eta, double phi, double pe, double peta, double pphi)
{
  const double pt = e / std::cosh(eta), partner_pt = pe / std::cosh(peta);
  const double delta_phi = std::atan2(std::sin(phi - pphi), std::cos(phi - pphi));
  const double mass2 = 2 * pt * partner_pt * (std::cosh(eta - peta) - std::cos(delta_phi));
  return std::isfinite(mass2) && mass2 >= 0 ? std::sqrt(mass2) : -1.0;
}

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
  text.DrawLatex(0.15, 0.96, "#it{#bf{sPHENIX}} Internal");
  std::istringstream stream(text_value);
  std::string line;
  double y = 0.91;
  while (std::getline(stream, line)) { text.DrawLatex(0.15, y, line.c_str()); y -= 0.04; }
}

void draw_mass(TH1D* histogram, const std::string& path, const std::string& text)
{
  TCanvas canvas("c_mass", "", 1000, 800);
  canvas.SetTopMargin(0.23);
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
  canvas.SetTopMargin(0.23);
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
  canvas.SetTopMargin(0.23);
  canvas.SetRightMargin(0.30);
  TLegend legend(0.72, 0.30, 0.99, 0.76);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetTextSize(0.024);
  std::vector<std::unique_ptr<TH1D>> shapes;
  double maximum = 0;
  for (int bin = 1; bin <= source->GetNbinsX(); ++bin)
  {
    auto shape = std::unique_ptr<TH1D>(source->ProjectionY(("shape_" + std::to_string(bin)).c_str(), bin, bin));
    const double entries = shape->Integral(0, shape->GetNbinsX() + 1);
    if (entries != 0) shape->Scale(1.0 / entries);
    else shape->Reset();
    shape->SetLineColor(kColors[(bin - 1) % kColors.size()]);
    shape->SetLineWidth(2);
    maximum = std::max(maximum, shape->GetMaximum());
    std::ostringstream item;
    item << std::fixed << std::setprecision(0) << edges[bin - 1] << "-" << edges[bin] << " GeV";
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
void plot_groups(const std::array<std::unique_ptr<Histograms>, 6>& groups, const Metadata& metadata, const std::string& output)
{
  SetsPhenixStyle(); gStyle->SetOptStat(0);
  for (std::size_t i = 0; i < groups.size(); ++i)
  {
    const bool eta = i >= 4;
    kPi0MassMin = metadata.selection_settings[eta ? 5 : 3];
    kPi0MassMax = metadata.selection_settings[eta ? 6 : 4];
    const std::string directory = output + "/" + kSelectionKeys[i];
    gSystem->mkdir(directory.c_str(), true);
    std::ostringstream caption;
    caption << (metadata.family == "jet" ? "Pythia8 p+p Jet" : "Pythia8 p+p PhotonJet") << ", Region A before veto\n" << kSelectionLabels[i] << "\n";
    if (i >= 2) caption << "E_{partner} > " << metadata.selection_settings[eta ? 2 : 1] << " GeV";
    else if (i == 1) caption << "E_{representative} #leq " << metadata.selection_settings[1] << " GeV (subset)";
    else caption << "No representative partner cut; topology partner E > " << metadata.selection_settings[1] << " GeV";
    caption << "; window " << kPi0MassMin << "-" << kPi0MassMax << " GeV";
    draw_mass(groups[i]->mass_pb.get(), directory + "/mass.pdf", caption.str());
    draw_2d(groups[i]->mass_vs_anchor_et_pb.get(), directory + "/mass_vs_candidate_et.pdf", "Candidate cluster E_{T} [GeV]", caption.str());
    draw_shapes(groups[i]->mass_vs_anchor_et_count.get(), kAnchorEtEdges, directory + "/mass_by_candidate_et_raw_shape.pdf", caption.str());
    draw_shapes(groups[i]->mass_vs_anchor_et_pb.get(), kAnchorEtEdges, directory + "/mass_by_candidate_et_weighted_shape.pdf", caption.str());
  }
}
}
#endif
