#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TObjArray.h>
#include <TObject.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace
{
struct PartialMetadata
{
  std::string path;
  std::string manifest_path;
  std::string first_suffix;
  std::string last_suffix;
  std::string cluster_collection;
  std::string pi0_selection;
  std::string cluster_selection;
  std::string fraction_definition;
  std::string zero_threshold_definition;
  long long manifest_begin = -1;
  long long manifest_end = -1;
  int schema_version = 0;
  int signal_embedding_id = 0;
  int truth_matcher_algorithm_version = 0;
  int pt_bins = 0;
  int multiplicity_max = 0;
  int cluster_energy_bins = 0;
  double pt_max = 0.0;
  double cluster_energy_max = 0.0;
  double truth_eta_max = 0.0;
  double cluster_eta_max = 0.0;
  unsigned char cluster_energy_cut_applied = 1U;
  std::array<double, 4> thresholds{};
  unsigned long long events_processed = 0;
  unsigned long long events_written = 0;
  unsigned long long events_invalid = 0;
  unsigned long long cluster_considered = 0;
  unsigned long long cluster_invalid_truth = 0;
  unsigned long long candidate_count = 0;
  unsigned long long candidate_g4_primary = 0;
  unsigned long long candidate_generator = 0;
  unsigned long long malformed_daughters = 0;
  unsigned long long pair_evaluated = 0;
  unsigned long long pair_positive = 0;
};

template <class T>
bool bind(TTree* tree, const char* name, T* address)
{
  return tree->GetBranch(name) &&
      tree->SetBranchAddress(name, address) >= 0;
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
  std::string* first_suffix = nullptr;
  std::string* last_suffix = nullptr;
  std::string* cluster_collection = nullptr;
  std::string* pi0_selection = nullptr;
  std::string* cluster_selection = nullptr;
  std::string* fraction_definition = nullptr;
  std::string* zero_threshold_definition = nullptr;
  bool ok = true;
  ok &= bind(tree, "schema_version", &value.schema_version);
  ok &= bind(tree, "manifest_path", &manifest_path);
  ok &= bind(tree, "manifest_begin", &value.manifest_begin);
  ok &= bind(tree, "manifest_end", &value.manifest_end);
  ok &= bind(tree, "first_suffix", &first_suffix);
  ok &= bind(tree, "last_suffix", &last_suffix);
  ok &= bind(tree, "cluster_collection", &cluster_collection);
  ok &= bind(tree, "pi0_selection", &pi0_selection);
  ok &= bind(tree, "cluster_selection", &cluster_selection);
  ok &= bind(tree, "fraction_definition", &fraction_definition);
  ok &= bind(tree, "zero_threshold_definition", &zero_threshold_definition);
  ok &= bind(tree, "signal_embedding_id", &value.signal_embedding_id);
  ok &= bind(tree, "truth_matcher_algorithm_version",
             &value.truth_matcher_algorithm_version);
  ok &= bind(tree, "pt_bins", &value.pt_bins);
  ok &= bind(tree, "pt_max", &value.pt_max);
  ok &= bind(tree, "multiplicity_max", &value.multiplicity_max);
  ok &= bind(tree, "cluster_energy_bins", &value.cluster_energy_bins);
  ok &= bind(tree, "cluster_energy_max", &value.cluster_energy_max);
  ok &= bind(tree, "truth_eta_max", &value.truth_eta_max);
  ok &= bind(tree, "cluster_eta_max", &value.cluster_eta_max);
  ok &= bind(tree, "cluster_energy_cut_applied",
             &value.cluster_energy_cut_applied);
  ok &= bind(tree, "fraction_thresholds", value.thresholds.data());
  ok &= bind(tree, "events_processed", &value.events_processed);
  ok &= bind(tree, "events_written", &value.events_written);
  ok &= bind(tree, "events_invalid", &value.events_invalid);
  ok &= bind(tree, "cluster_considered_count", &value.cluster_considered);
  ok &= bind(tree, "cluster_invalid_truth_count",
             &value.cluster_invalid_truth);
  ok &= bind(tree, "pi0_candidate_count", &value.candidate_count);
  ok &= bind(tree, "pi0_candidate_g4_primary_count",
             &value.candidate_g4_primary);
  ok &= bind(tree, "pi0_candidate_generator_count",
             &value.candidate_generator);
  ok &= bind(tree, "pi0_malformed_daughters_count",
             &value.malformed_daughters);
  ok &= bind(tree, "pi0_cluster_pair_evaluated_count",
             &value.pair_evaluated);
  ok &= bind(tree, "pi0_cluster_pair_positive_count",
             &value.pair_positive);
  if (!ok || tree->GetEntry(0) <= 0 || value.schema_version != 2 ||
      !manifest_path || !first_suffix ||
      !last_suffix || !cluster_collection || !pi0_selection ||
      !cluster_selection || !fraction_definition || !zero_threshold_definition)
  {
    return false;
  }
  value.path = path;
  value.manifest_path = *manifest_path;
  value.first_suffix = *first_suffix;
  value.last_suffix = *last_suffix;
  value.cluster_collection = *cluster_collection;
  value.pi0_selection = *pi0_selection;
  value.cluster_selection = *cluster_selection;
  value.fraction_definition = *fraction_definition;
  value.zero_threshold_definition = *zero_threshold_definition;
  return true;
}

bool same_double(double left, double right)
{
  return std::abs(left - right) <=
      1e-12 * std::max({1.0, std::abs(left), std::abs(right)});
}

bool compatible(const PartialMetadata& value, const PartialMetadata& reference)
{
  if (value.schema_version != reference.schema_version ||
      value.manifest_path != reference.manifest_path ||
      value.cluster_collection != reference.cluster_collection ||
      value.pi0_selection != reference.pi0_selection ||
      value.cluster_selection != reference.cluster_selection ||
      value.fraction_definition != reference.fraction_definition ||
      value.zero_threshold_definition != reference.zero_threshold_definition ||
      value.signal_embedding_id != reference.signal_embedding_id ||
      value.truth_matcher_algorithm_version !=
          reference.truth_matcher_algorithm_version ||
      value.pt_bins != reference.pt_bins ||
      value.multiplicity_max != reference.multiplicity_max ||
      value.cluster_energy_bins != reference.cluster_energy_bins ||
      value.cluster_energy_cut_applied !=
          reference.cluster_energy_cut_applied ||
      !same_double(value.pt_max, reference.pt_max) ||
      !same_double(value.cluster_energy_max, reference.cluster_energy_max) ||
      !same_double(value.truth_eta_max, reference.truth_eta_max) ||
      !same_double(value.cluster_eta_max, reference.cluster_eta_max))
  {
    return false;
  }
  for (std::size_t i = 0; i < value.thresholds.size(); ++i)
  {
    if (!same_double(value.thresholds[i], reference.thresholds[i]))
    {
      return false;
    }
  }
  return true;
}

std::vector<std::string> histogram_names()
{
  std::vector<std::string> result = {
      "h_pi0_maximum_compatible_fraction_raw",
      "h_pi0_second_compatible_fraction_raw",
      "h_pi0_compatible_fraction_vs_cluster_energy_raw"};
  const std::array<std::string, 4> thresholds = {
      "0p0", "0p1", "0p3", "0p5"};
  const std::array<std::string, 2> pathways = {
      "g4_primary", "generator"};
  for (const std::string& threshold : thresholds)
  {
    result.push_back(
        "h_pi0_cluster_multiplicity_fmin_" + threshold + "_raw");
    result.push_back(
        "h_pi0_cluster_multiplicity_vs_truth_pt_fmin_" + threshold + "_raw");
    for (const std::string& pathway : pathways)
    {
      result.push_back("h_pi0_cluster_multiplicity_fmin_" + threshold + "_" +
                       pathway + "_raw");
    }
  }
  return result;
}

unsigned long long expected_entries(const std::string& name,
                                    const PartialMetadata& metadata)
{
  if (name == "h_pi0_compatible_fraction_vs_cluster_energy_raw")
  {
    return metadata.pair_positive;
  }
  if (name.find("_g4_primary_raw") != std::string::npos)
  {
    return metadata.candidate_g4_primary;
  }
  if (name.find("_generator_raw") != std::string::npos)
  {
    return metadata.candidate_generator;
  }
  return metadata.candidate_count;
}

bool make_output_directory(const std::string& output_base)
{
  const std::size_t slash = output_base.find_last_of('/');
  if (slash == std::string::npos)
  {
    return true;
  }
  const std::string directory = output_base.substr(0, slash);
  return directory.empty() || gSystem->mkdir(directory.c_str(), true) == 0 ||
      !gSystem->AccessPathName(directory.c_str());
}
}

int FinalizePythiaPi0ClusterMultiplicity(
    const std::string partial_pattern =
        "output/pi0_cluster_multiplicity_partial/pilot_primary_generator_eta07_no_cluster_energy_cut/partial_*.root",
    const std::string output_base =
        "output/plots/pi0_cluster_multiplicity_primary_generator_eta07_no_cluster_energy_cut",
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
    std::cerr << "FinalizePythiaPi0ClusterMultiplicity - no partials matched"
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
      std::cerr << "FinalizePythiaPi0ClusterMultiplicity - invalid partial: "
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
      std::cerr << "FinalizePythiaPi0ClusterMultiplicity - incompatible or noncontiguous partial: "
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

  const std::vector<std::string> names = histogram_names();
  std::map<std::string, std::unique_ptr<TH1>> histograms;
  PartialMetadata total = reference;
  total.events_processed = total.events_written = total.events_invalid = 0;
  total.cluster_considered = total.cluster_invalid_truth = 0;
  total.candidate_count = total.candidate_g4_primary =
      total.candidate_generator = 0;
  total.malformed_daughters = total.pair_evaluated = total.pair_positive = 0;

  for (const PartialMetadata& partial : partials)
  {
    TFile input(partial.path.c_str(), "READ");
    if (input.IsZombie())
    {
      return 5;
    }
    for (const std::string& name : names)
    {
      TH1* source = nullptr;
      input.GetObject(name.c_str(), source);
      if (!source ||
          std::abs(source->GetEntries() -
                   static_cast<double>(expected_entries(name, partial))) > 0.5)
      {
        std::cerr << "FinalizePythiaPi0ClusterMultiplicity - invalid histogram "
                  << name << " in " << partial.path << std::endl;
        return 5;
      }
      auto found = histograms.find(name);
      if (found == histograms.end())
      {
        std::unique_ptr<TH1> clone(
            static_cast<TH1*>(source->Clone(name.c_str())));
        clone->SetDirectory(nullptr);
        histograms.emplace(name, std::move(clone));
      }
      else if (!found->second->Add(source))
      {
        return 5;
      }
    }
    total.events_processed += partial.events_processed;
    total.events_written += partial.events_written;
    total.events_invalid += partial.events_invalid;
    total.cluster_considered += partial.cluster_considered;
    total.cluster_invalid_truth += partial.cluster_invalid_truth;
    total.candidate_count += partial.candidate_count;
    total.candidate_g4_primary += partial.candidate_g4_primary;
    total.candidate_generator += partial.candidate_generator;
    total.malformed_daughters += partial.malformed_daughters;
    total.pair_evaluated += partial.pair_evaluated;
    total.pair_positive += partial.pair_positive;
  }

  if (!make_output_directory(output_base))
  {
    return 6;
  }

  const std::array<std::string, 4> tags = {"0p0", "0p1", "0p3", "0p5"};
  const std::array<int, 4> colors = {kBlack, kBlue + 1, kGreen + 2, kRed + 1};
  std::array<std::unique_ptr<TH1D>, 4> probability;
  double maximum_probability = 0.0;
  for (std::size_t index = 0; index < tags.size(); ++index)
  {
    const std::string raw_name =
        "h_pi0_cluster_multiplicity_fmin_" + tags[index] + "_raw";
    TH1D* raw = dynamic_cast<TH1D*>(histograms.at(raw_name).get());
    probability[index].reset(static_cast<TH1D*>(raw->Clone(
        ("h_pi0_cluster_multiplicity_fmin_" + tags[index] +
         "_probability").c_str())));
    probability[index]->SetDirectory(nullptr);
    const double integral = probability[index]->Integral(
        0, probability[index]->GetNbinsX() + 1);
    if (integral > 0.0)
    {
      probability[index]->Scale(1.0 / integral);
    }
    probability[index]->SetStats(false);
    probability[index]->SetLineColor(colors[index]);
    probability[index]->SetMarkerColor(colors[index]);
    probability[index]->SetLineWidth(2);
    probability[index]->GetXaxis()->SetTitle(
        "N_{cluster} compatible with selected #pi^{0}");
    probability[index]->GetYaxis()->SetTitle("Probability");
    maximum_probability =
        std::max(maximum_probability, probability[index]->GetMaximum());
  }

  gStyle->SetOptStat(0);
  TCanvas multiplicity_canvas(
      "c_pi0_cluster_multiplicity", "Pi0 cluster multiplicity", 1000, 800);
  multiplicity_canvas.SetLogy();
  probability[0]->SetMinimum(1e-7);
  probability[0]->SetMaximum(
      maximum_probability > 0.0 ? 5.0 * maximum_probability : 1.0);
  probability[0]->Draw("HIST");
  for (std::size_t index = 1; index < probability.size(); ++index)
  {
    probability[index]->Draw("HIST SAME");
  }
  TLegend legend(0.58, 0.62, 0.88, 0.86);
  legend.SetBorderSize(0);
  legend.AddEntry(probability[0].get(), "f_{#pi^{0}} > 0", "l");
  legend.AddEntry(probability[1].get(), "f_{#pi^{0}} #geq 0.1", "l");
  legend.AddEntry(probability[2].get(), "f_{#pi^{0}} #geq 0.3", "l");
  legend.AddEntry(probability[3].get(), "f_{#pi^{0}} #geq 0.5", "l");
  legend.Draw();
  TLatex label;
  label.SetNDC();
  label.SetTextAlign(13);
  label.DrawLatex(0.18, 0.92, "#it{#bf{sPHENIX}} Internal");
  label.DrawLatex(0.18, 0.84, "Pythia8 p+p minimum bias");
  std::ostringstream selection;
  selection << "|#eta_{#pi^{0}}^{truth}| < " << reference.truth_eta_max
            << ", |#eta_{cluster}| < " << reference.cluster_eta_max;
  label.DrawLatex(0.18, 0.76, selection.str().c_str());
  label.DrawLatex(0.18, 0.68, "No cluster-energy threshold");
  multiplicity_canvas.RedrawAxis();
  multiplicity_canvas.SaveAs((output_base + ".pdf").c_str());

  TCanvas pt_canvas(
      "c_pi0_cluster_multiplicity_vs_pt", "Multiplicity versus pi0 pT",
      1200, 1000);
  pt_canvas.Divide(2, 2);
  for (std::size_t index = 0; index < tags.size(); ++index)
  {
    pt_canvas.cd(static_cast<int>(index + 1));
    gPad->SetLogz();
    const std::string name =
        "h_pi0_cluster_multiplicity_vs_truth_pt_fmin_" + tags[index] + "_raw";
    TH2D* histogram = dynamic_cast<TH2D*>(histograms.at(name).get());
    histogram->SetTitle(index == 0 ? "f_{#pi^{0}} > 0" :
        (std::string("f_{#pi^{0}} #geq ") +
         std::to_string(reference.thresholds[index])).c_str());
    histogram->GetXaxis()->SetTitle("p_{T}^{truth,#pi^{0}} [GeV]");
    histogram->GetYaxis()->SetTitle("N_{cluster}");
    histogram->Draw("COLZ");
  }
  pt_canvas.SaveAs((output_base + "_vs_truth_pt.pdf").c_str());

  TCanvas fraction_canvas(
      "c_pi0_fraction_vs_cluster_energy", "Fraction versus cluster energy",
      1000, 800);
  fraction_canvas.SetLogz();
  TH2D* fraction_vs_energy = dynamic_cast<TH2D*>(
      histograms.at("h_pi0_compatible_fraction_vs_cluster_energy_raw").get());
  fraction_vs_energy->GetXaxis()->SetTitle("Cluster energy [GeV]");
  fraction_vs_energy->GetYaxis()->SetTitle("f_{#pi^{0}}^{cluster}");
  fraction_vs_energy->Draw("COLZ");
  fraction_canvas.SaveAs(
      (output_base + "_fraction_vs_cluster_energy.pdf").c_str());

  TFile output((output_base + ".root").c_str(), "RECREATE");
  if (output.IsZombie())
  {
    return 6;
  }
  for (auto& item : histograms)
  {
    item.second->Write();
  }
  for (auto& histogram : probability)
  {
    histogram->Write();
  }

  double summary_threshold = 0.0;
  double summary_mean = 0.0;
  double summary_probability_ge2 = 0.0;
  double summary_probability_ge3 = 0.0;
  double summary_overflow_probability = 0.0;
  TTree summary("threshold_summary", "Pi0 cluster multiplicity summary");
  summary.Branch("fraction_threshold", &summary_threshold);
  summary.Branch("mean_multiplicity", &summary_mean);
  summary.Branch("probability_ge2", &summary_probability_ge2);
  summary.Branch("probability_ge3", &summary_probability_ge3);
  summary.Branch("overflow_probability", &summary_overflow_probability);
  for (std::size_t index = 0; index < probability.size(); ++index)
  {
    TH1D* histogram = probability[index].get();
    summary_threshold = reference.thresholds[index];
    summary_mean = histogram->GetMean();
    summary_probability_ge2 = histogram->Integral(
        histogram->FindFixBin(2.0), histogram->GetNbinsX() + 1);
    summary_probability_ge3 = histogram->Integral(
        histogram->FindFixBin(3.0), histogram->GetNbinsX() + 1);
    summary_overflow_probability =
        histogram->GetBinContent(histogram->GetNbinsX() + 1);
    summary.Fill();
    std::cout << "fmin=" << summary_threshold
              << " mean=" << summary_mean
              << " P(N>=2)=" << summary_probability_ge2
              << " P(N>=3)=" << summary_probability_ge3
              << " overflow=" << summary_overflow_probability << std::endl;
  }
  summary.Write();

  int output_schema_version = 2;
  long long manifest_begin = partials.front().manifest_begin;
  long long manifest_end = partials.back().manifest_end;
  long long partial_file_count = static_cast<long long>(partials.size());
  long long input_file_count = manifest_end - manifest_begin;
  TTree metadata("metadata", "Final Pythia pi0 cluster multiplicity metadata");
  metadata.Branch("schema_version", &output_schema_version);
  metadata.Branch("manifest_path", &total.manifest_path);
  metadata.Branch("manifest_begin", &manifest_begin);
  metadata.Branch("manifest_end", &manifest_end);
  metadata.Branch("partial_file_count", &partial_file_count);
  metadata.Branch("input_file_count", &input_file_count);
  metadata.Branch("cluster_collection", &total.cluster_collection);
  metadata.Branch("pi0_selection", &total.pi0_selection);
  metadata.Branch("cluster_selection", &total.cluster_selection);
  metadata.Branch("fraction_definition", &total.fraction_definition);
  metadata.Branch("zero_threshold_definition",
                  &total.zero_threshold_definition);
  metadata.Branch("signal_embedding_id", &total.signal_embedding_id);
  metadata.Branch("truth_matcher_algorithm_version",
                  &total.truth_matcher_algorithm_version);
  metadata.Branch("pt_bins", &total.pt_bins);
  metadata.Branch("pt_max", &total.pt_max);
  metadata.Branch("multiplicity_max", &total.multiplicity_max);
  metadata.Branch("cluster_energy_bins", &total.cluster_energy_bins);
  metadata.Branch("cluster_energy_max", &total.cluster_energy_max);
  metadata.Branch("truth_eta_max", &total.truth_eta_max);
  metadata.Branch("cluster_eta_max", &total.cluster_eta_max);
  metadata.Branch("cluster_energy_cut_applied",
                  &total.cluster_energy_cut_applied);
  metadata.Branch("fraction_thresholds", total.thresholds.data(),
                  "fraction_thresholds[4]/D");
  metadata.Branch("events_processed", &total.events_processed);
  metadata.Branch("events_written", &total.events_written);
  metadata.Branch("events_invalid", &total.events_invalid);
  metadata.Branch("cluster_considered_count", &total.cluster_considered);
  metadata.Branch("cluster_invalid_truth_count",
                  &total.cluster_invalid_truth);
  metadata.Branch("pi0_candidate_count", &total.candidate_count);
  metadata.Branch("pi0_candidate_g4_primary_count",
                  &total.candidate_g4_primary);
  metadata.Branch("pi0_candidate_generator_count",
                  &total.candidate_generator);
  metadata.Branch("pi0_malformed_daughters_count",
                  &total.malformed_daughters);
  metadata.Branch("pi0_cluster_pair_evaluated_count",
                  &total.pair_evaluated);
  metadata.Branch("pi0_cluster_pair_positive_count",
                  &total.pair_positive);
  metadata.Fill();
  metadata.Write();
  output.Close();
  if (output.TestBit(TFile::kWriteError))
  {
    return 6;
  }

  std::cout << "FinalizePythiaPi0ClusterMultiplicity - partials/files/events/pi0/positive pairs = "
            << partial_file_count << "/" << input_file_count << "/"
            << total.events_written << "/" << total.candidate_count << "/"
            << total.pair_positive << std::endl;
  return 0;
}
