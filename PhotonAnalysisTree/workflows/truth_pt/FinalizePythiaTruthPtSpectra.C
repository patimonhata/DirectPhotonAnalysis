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
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <set>
#include <string>
#include <vector>

namespace
{
struct PartialMetadata
{
  std::string path;
  Int_t schema_version = 0;
  std::string photon_selection;
  std::string pi0_decay_photon_selection;
  std::string manifest_path;
  std::string input_file_prefix;
  Long64_t manifest_begin = -1;
  Long64_t manifest_end = -1;
  Long64_t files_added = 0;
  Long64_t events_processed = 0;
  Int_t n_bins = 0;
  double pt_max = 0.0;
  double max_abs_eta = 0.0;
  UChar_t use_event_weight = 0U;
  UChar_t bin_width_normalized = 1U;
  ULong64_t prompt_photon_count = 0ULL;
  ULong64_t pi0_count = 0ULL;
  ULong64_t pi0_decay_photon_count = 0ULL;
  ULong64_t hepmc_pi0_decay_photon_count = 0ULL;
  ULong64_t g4_pi0_decay_photon_count = 0ULL;
  ULong64_t malformed_event_count = 0ULL;
  ULong64_t invalid_weight_event_count = 0ULL;
};

template <class T>
bool bind_branch(TTree* tree, const char* name, T* address)
{
  return tree->GetBranch(name) && tree->SetBranchAddress(name, address) >= 0;
}

bool read_metadata(const std::string& path, PartialMetadata& result)
{
  TFile input(path.c_str(), "READ");
  TTree* tree = nullptr;
  input.GetObject("metadata", tree);
  if (input.IsZombie() || !tree || tree->GetEntries() != 1)
  {
    return false;
  }
  std::string* photon_selection = nullptr;
  std::string* pi0_decay_photon_selection = nullptr;
  std::string* manifest_path = nullptr;
  std::string* input_file_prefix = nullptr;
  bool ok = true;
  ok &= bind_branch(tree, "schema_version", &result.schema_version);
  ok &= bind_branch(tree, "photon_selection", &photon_selection);
  ok &= bind_branch(tree, "pi0_decay_photon_selection", &pi0_decay_photon_selection);
  ok &= bind_branch(tree, "manifest_path", &manifest_path);
  ok &= bind_branch(tree, "input_file_prefix", &input_file_prefix);
  ok &= bind_branch(tree, "manifest_begin", &result.manifest_begin);
  ok &= bind_branch(tree, "manifest_end", &result.manifest_end);
  ok &= bind_branch(tree, "files_added", &result.files_added);
  ok &= bind_branch(tree, "events_processed", &result.events_processed);
  ok &= bind_branch(tree, "n_bins", &result.n_bins);
  ok &= bind_branch(tree, "pt_max", &result.pt_max);
  ok &= bind_branch(tree, "max_abs_eta", &result.max_abs_eta);
  ok &= bind_branch(tree, "use_event_weight", &result.use_event_weight);
  ok &= bind_branch(tree, "bin_width_normalized", &result.bin_width_normalized);
  ok &= bind_branch(tree, "prompt_photon_count", &result.prompt_photon_count);
  ok &= bind_branch(tree, "pi0_count", &result.pi0_count);
  ok &= bind_branch(tree, "pi0_decay_photon_count", &result.pi0_decay_photon_count);
  ok &= bind_branch(tree, "hepmc_pi0_decay_photon_count", &result.hepmc_pi0_decay_photon_count);
  ok &= bind_branch(tree, "g4_pi0_decay_photon_count", &result.g4_pi0_decay_photon_count);
  ok &= bind_branch(tree, "malformed_event_count", &result.malformed_event_count);
  ok &= bind_branch(tree, "invalid_weight_event_count", &result.invalid_weight_event_count);
  if (!ok || tree->GetEntry(0) <= 0 || !photon_selection ||
      !pi0_decay_photon_selection || !manifest_path || !input_file_prefix)
  {
    return false;
  }
  result.path = path;
  result.photon_selection = *photon_selection;
  result.pi0_decay_photon_selection = *pi0_decay_photon_selection;
  result.manifest_path = *manifest_path;
  result.input_file_prefix = *input_file_prefix;
  return true;
}

bool compatible(const PartialMetadata& value, const PartialMetadata& reference)
{
  return value.schema_version == 3 &&
      value.photon_selection == "prompt_category_1_or_2" &&
      value.photon_selection == reference.photon_selection &&
      value.pi0_decay_photon_selection ==
          "hepmc_final_photon_with_valid_single_pi0_origin_plus_"
          "g4_immediate_photon_daughter_of_signal_primary_pi0" &&
      value.pi0_decay_photon_selection == reference.pi0_decay_photon_selection &&
      value.manifest_path == reference.manifest_path &&
      value.input_file_prefix == "G4Hits_" && value.input_file_prefix == reference.input_file_prefix &&
      value.n_bins == reference.n_bins &&
      std::abs(value.pt_max - reference.pt_max) < 1e-12 &&
      std::abs(value.max_abs_eta - reference.max_abs_eta) < 1e-12 &&
      value.use_event_weight == reference.use_event_weight &&
      value.bin_width_normalized == 0U &&
      value.pi0_decay_photon_count ==
          value.hepmc_pi0_decay_photon_count + value.g4_pi0_decay_photon_count &&
      value.files_added == value.manifest_end - value.manifest_begin &&
      value.files_added > 0 && value.events_processed > 0 &&
      value.malformed_event_count == 0ULL && value.invalid_weight_event_count == 0ULL;
}

bool compatible_histogram(const TH1D* histogram, const PartialMetadata& metadata,
    const ULong64_t expected_entries)
{
  const bool valid = histogram && histogram->GetNbinsX() == metadata.n_bins &&
      std::abs(histogram->GetXaxis()->GetXmin()) < 1e-12 &&
      std::abs(histogram->GetXaxis()->GetXmax() - metadata.pt_max) < 1e-12 &&
      histogram->GetSumw2N() > 0 &&
      std::abs(histogram->GetEntries() - static_cast<double>(expected_entries)) < 0.5;
  if (!valid)
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

bool histogram_sum_matches(const TH1D* total, const TH1D* hepmc, const TH1D* g4)
{
  if (!total || !hepmc || !g4 ||
      std::abs(total->GetEntries() - hepmc->GetEntries() - g4->GetEntries()) > 0.5)
  {
    return false;
  }
  for (int bin = 0; bin <= total->GetNbinsX() + 1; ++bin)
  {
    const double expected_content =
        hepmc->GetBinContent(bin) + g4->GetBinContent(bin);
    const double expected_error2 =
        std::pow(hepmc->GetBinError(bin), 2) + std::pow(g4->GetBinError(bin), 2);
    const double content_scale = std::max(1.0, std::abs(expected_content));
    const double error2_scale = std::max(1.0, std::abs(expected_error2));
    if (std::abs(total->GetBinContent(bin) - expected_content) >
            1e-10 * content_scale ||
        std::abs(std::pow(total->GetBinError(bin), 2) - expected_error2) >
            1e-10 * error2_scale)
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

double smallest_positive_bin(const std::vector<TH1D*>& histograms)
{
  double result = std::numeric_limits<double>::infinity();
  for (const TH1D* histogram : histograms)
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

int FinalizePythiaTruthPtSpectra(
    const std::string partial_pattern = "output/truth_pt_partial/minimum_bias/prompt_eta07_unweighted_inclusive_pi0_decay/partial_*.root",
    const std::string output_base = "output/plots/truth_pT/minbias/prompt_eta07_inclusive_pi0_decay",
    const Long64_t expected_manifest_begin = 0,
    const Long64_t expected_manifest_end = -1)
{
  if (partial_pattern.empty() || output_base.empty() ||
      expected_manifest_begin < 0 || expected_manifest_end < -1 ||
      (expected_manifest_end >= 0 && expected_manifest_end <= expected_manifest_begin)) {
    std::cerr << "FinalizePythiaTruthPtSpectra - invalid argument" << std::endl;
    return 1;
  }

  TChain partial_chain("metadata");
  const int matched_files = partial_chain.Add(partial_pattern.c_str());
  const TObjArray* file_elements = partial_chain.GetListOfFiles();
  if (matched_files <= 0 || !file_elements || file_elements->GetEntries() <= 0) {
    std::cerr << "FinalizePythiaTruthPtSpectra - no partial files matched " << partial_pattern << std::endl;
    return 2;
  }

  std::vector<PartialMetadata> partials;
  std::set<std::string> unique_paths;
  for (int index = 0; index < file_elements->GetEntries(); ++index) {
    const TObject* element = file_elements->At(index);
    const std::string path = element ? element->GetTitle() : "";
    PartialMetadata metadata;
    if (path.empty() || !unique_paths.insert(path).second || !read_metadata(path, metadata)) {
      std::cerr << "FinalizePythiaTruthPtSpectra - invalid or duplicate partial: " << path << std::endl;
      return 3;
    }
    partials.push_back(metadata);
  }
  std::sort(partials.begin(), partials.end(), [](const auto& left, const auto& right) {
    return left.manifest_begin < right.manifest_begin;
  });

  const PartialMetadata& reference = partials.front();
  Long64_t next_begin = expected_manifest_begin;
  for (const PartialMetadata& partial : partials) {
    if (!compatible(partial, reference) || partial.manifest_begin != next_begin) {
      std::cerr << "FinalizePythiaTruthPtSpectra - incompatible, missing, or overlapping range at "
                << partial.path << "; expected begin " << next_begin << ", observed ["
                << partial.manifest_begin << ":" << partial.manifest_end << "]" << std::endl;
      return 4;
    }
    next_begin = partial.manifest_end;
  }
  if (expected_manifest_end >= 0 && next_begin != expected_manifest_end) {
    std::cerr << "FinalizePythiaTruthPtSpectra - final manifest end mismatch: " << next_begin << " != " << expected_manifest_end << std::endl;
    return 4;
  }

  TH1::AddDirectory(false);
  TH1D prompt_raw("h_prompt_photon_truth_pt_raw", "", reference.n_bins, 0.0, reference.pt_max);
  TH1D pi0_raw("h_pi0_truth_pt_raw", "", reference.n_bins, 0.0, reference.pt_max);
  TH1D pi0_decay_raw("h_pi0_decay_photon_truth_pt_raw", "", reference.n_bins, 0.0, reference.pt_max);
  TH1D hepmc_pi0_decay_raw("h_hepmc_pi0_decay_photon_truth_pt_raw", "", reference.n_bins, 0.0, reference.pt_max);
  TH1D g4_pi0_decay_raw("h_g4_pi0_decay_photon_truth_pt_raw", "", reference.n_bins, 0.0, reference.pt_max);
  for (TH1D* histogram : {&prompt_raw, &pi0_raw, &pi0_decay_raw, &hepmc_pi0_decay_raw, &g4_pi0_decay_raw}) {
    histogram->Sumw2();
  }

  Long64_t total_input_files = 0;
  Long64_t total_events = 0;
  ULong64_t total_prompt_photons = 0ULL;
  ULong64_t total_pi0 = 0ULL;
  ULong64_t total_pi0_decay_photons = 0ULL;
  ULong64_t total_hepmc_pi0_decay_photons = 0ULL;
  ULong64_t total_g4_pi0_decay_photons = 0ULL;
  for (const PartialMetadata& partial : partials) {
    TFile input(partial.path.c_str(), "READ");
    TH1D* prompt = nullptr;
    TH1D* pi0 = nullptr;
    TH1D* pi0_decay = nullptr;
    TH1D* hepmc_pi0_decay = nullptr;
    TH1D* g4_pi0_decay = nullptr;
    input.GetObject("h_prompt_photon_truth_pt_raw", prompt);
    input.GetObject("h_pi0_truth_pt_raw", pi0);
    input.GetObject("h_pi0_decay_photon_truth_pt_raw", pi0_decay);
    input.GetObject("h_hepmc_pi0_decay_photon_truth_pt_raw", hepmc_pi0_decay);
    input.GetObject("h_g4_pi0_decay_photon_truth_pt_raw", g4_pi0_decay);
    if (input.IsZombie() ||
        !compatible_histogram(prompt, partial, partial.prompt_photon_count) ||
        !compatible_histogram(pi0, partial, partial.pi0_count) ||
        !compatible_histogram(pi0_decay, partial, partial.pi0_decay_photon_count) ||
        !compatible_histogram(hepmc_pi0_decay, partial,
            partial.hepmc_pi0_decay_photon_count) ||
        !compatible_histogram(g4_pi0_decay, partial,
            partial.g4_pi0_decay_photon_count) ||
        !histogram_sum_matches(pi0_decay, hepmc_pi0_decay, g4_pi0_decay) ||
        prompt_raw.Add(prompt) == false || pi0_raw.Add(pi0) == false ||
        pi0_decay_raw.Add(pi0_decay) == false ||
        hepmc_pi0_decay_raw.Add(hepmc_pi0_decay) == false ||
        g4_pi0_decay_raw.Add(g4_pi0_decay) == false)
    {
      std::cerr << "FinalizePythiaTruthPtSpectra - invalid histogram in " << partial.path << std::endl;
      return 5;
    }
    total_input_files += partial.files_added;
    total_events += partial.events_processed;
    total_prompt_photons += partial.prompt_photon_count;
    total_pi0 += partial.pi0_count;
    total_pi0_decay_photons += partial.pi0_decay_photon_count;
    total_hepmc_pi0_decay_photons += partial.hepmc_pi0_decay_photon_count;
    total_g4_pi0_decay_photons += partial.g4_pi0_decay_photon_count;
  }
  if (total_pi0_decay_photons != total_hepmc_pi0_decay_photons + total_g4_pi0_decay_photons) {
    std::cerr << "FinalizePythiaTruthPtSpectra - inconsistent pi0 decay totals" << std::endl;
    return 5;
  }

  std::unique_ptr<TH1D> prompt_density(static_cast<TH1D*>(prompt_raw.Clone("h_prompt_photon_truth_pt_density")));
  std::unique_ptr<TH1D> pi0_density(static_cast<TH1D*>(pi0_raw.Clone("h_pi0_truth_pt_density")));
  std::unique_ptr<TH1D> pi0_decay_density(static_cast<TH1D*>(pi0_decay_raw.Clone("h_pi0_decay_photon_truth_pt_density")));
  std::unique_ptr<TH1D> hepmc_pi0_decay_density(static_cast<TH1D*>(hepmc_pi0_decay_raw.Clone("h_hepmc_pi0_decay_photon_truth_pt_density")));
  std::unique_ptr<TH1D> g4_pi0_decay_density(static_cast<TH1D*>(g4_pi0_decay_raw.Clone("h_g4_pi0_decay_photon_truth_pt_density")));
  for (TH1D* histogram : {prompt_density.get(), pi0_density.get(), pi0_decay_density.get(), hepmc_pi0_decay_density.get(), g4_pi0_decay_density.get()}) {
    histogram->Scale(1.0, "width");
    histogram->SetStats(false);
    histogram->SetFillStyle(0);
    histogram->SetLineWidth(3);
    histogram->GetXaxis()->SetTitle("Truth p_{T} [GeV/#it{c}]");
    histogram->GetYaxis()->SetTitle(reference.use_event_weight ? "Weighted particles / (GeV/#it{c})" : "Particles / (GeV/#it{c})");
  }
  prompt_density->SetLineColor(kRed + 1);
  pi0_density->SetLineColor(kBlue + 1);
  pi0_decay_density->SetLineColor(kGreen + 2);

  for (TH1D* histogram : {&prompt_raw, &pi0_raw, &pi0_decay_raw, &hepmc_pi0_decay_raw, &g4_pi0_decay_raw}) {
    histogram->SetStats(false);
    histogram->SetFillStyle(0);
    histogram->SetLineWidth(3);
    histogram->GetXaxis()->SetTitle("Truth p_{T} [GeV/#it{c}]");
    histogram->GetYaxis()->SetTitle(reference.use_event_weight ? "Weighted counts / bin" : "Counts / bin");
  }
  prompt_raw.SetLineColor(kRed + 1);
  pi0_raw.SetLineColor(kBlue + 1);
  pi0_decay_raw.SetLineColor(kGreen + 2);

  const double maximum = std::max({prompt_raw.GetMaximum(), pi0_raw.GetMaximum(), pi0_decay_raw.GetMaximum()});
  const double smallest = smallest_positive_bin({&prompt_raw, &pi0_raw, &pi0_decay_raw});
  prompt_raw.SetMinimum(smallest > 0.0 ? 0.5 * smallest : 0.5);
  prompt_raw.SetMaximum(maximum > 0.0 ? 5.0 * maximum : 1.0);

  if (!make_output_directory(output_base)) {
    std::cerr << "FinalizePythiaTruthPtSpectra - cannot create output directory" << std::endl;
    return 6;
  }
  SetsPhenixStyle();
  TCanvas canvas("c_pythia_truth_pt_spectra", "Pythia truth pT spectra", 1000, 800);
  canvas.SetLogy();
  prompt_raw.Draw("HIST");
  pi0_raw.Draw("HIST SAME");
  pi0_decay_raw.Draw("HIST SAME");
  TLegend legend(0.40, 0.45, 0.79, 0.65);
  legend.AddEntry(&prompt_raw, "Prompt #gamma (direct + frag.)", "l");
  legend.AddEntry(&pi0_raw, "Inclusive #pi^{0}", "l");
  legend.AddEntry(&pi0_decay_raw, "#gamma from #pi^{0}", "l");
  legend.Draw();
  TLatex label;
  label.SetNDC();
  label.SetTextAlign(13);
  label.DrawLatex(0.23, 0.92, "#it{#bf{sPHENIX}} Internal");
  label.DrawLatex(0.23, 0.84, "Pythia8 p+p MB");
  const std::string selection = reference.max_abs_eta < 0.0
      ? "No #eta selection" : "|#eta^{truth}| < " + std::to_string(reference.max_abs_eta);
  label.DrawLatex(0.23, 0.76, selection.c_str());
  canvas.RedrawAxis();
  canvas.SaveAs((output_base + ".pdf").c_str());

  TFile output((output_base + ".root").c_str(), "RECREATE");
  if (output.IsZombie()) {
    return 6;
  }
  prompt_raw.Write();
  pi0_raw.Write();
  pi0_decay_raw.Write();
  hepmc_pi0_decay_raw.Write();
  g4_pi0_decay_raw.Write();
  prompt_density->Write();
  pi0_density->Write();
  pi0_decay_density->Write();
  hepmc_pi0_decay_density->Write();
  g4_pi0_decay_density->Write();

  Int_t output_schema_version = 2;
  std::string photon_selection = reference.photon_selection;
  std::string pi0_decay_photon_selection = reference.pi0_decay_photon_selection;
  std::string manifest_path = reference.manifest_path;
  std::string input_file_prefix = reference.input_file_prefix;
  Long64_t manifest_begin = partials.front().manifest_begin;
  Long64_t manifest_end = partials.back().manifest_end;
  Long64_t partial_file_count = static_cast<Long64_t>(partials.size());
  Int_t n_bins = reference.n_bins;
  double pt_max = reference.pt_max;
  double max_abs_eta = reference.max_abs_eta;
  UChar_t use_event_weight = reference.use_event_weight;
  UChar_t contains_raw_histograms = 1U;
  UChar_t contains_bin_width_normalized_histograms = 1U;
  TTree metadata("metadata", "Final Pythia truth pT metadata");
  metadata.Branch("schema_version", &output_schema_version);
  metadata.Branch("photon_selection", &photon_selection);
  metadata.Branch("pi0_decay_photon_selection", &pi0_decay_photon_selection);
  metadata.Branch("manifest_path", &manifest_path);
  metadata.Branch("input_file_prefix", &input_file_prefix);
  metadata.Branch("manifest_begin", &manifest_begin);
  metadata.Branch("manifest_end", &manifest_end);
  metadata.Branch("partial_file_count", &partial_file_count);
  metadata.Branch("input_file_count", &total_input_files);
  metadata.Branch("events_processed", &total_events);
  metadata.Branch("n_bins", &n_bins);
  metadata.Branch("pt_max", &pt_max);
  metadata.Branch("max_abs_eta", &max_abs_eta);
  metadata.Branch("use_event_weight", &use_event_weight);
  metadata.Branch("contains_raw_histograms", &contains_raw_histograms);
  metadata.Branch("contains_bin_width_normalized_histograms", &contains_bin_width_normalized_histograms);
  metadata.Branch("prompt_photon_count", &total_prompt_photons);
  metadata.Branch("pi0_count", &total_pi0);
  metadata.Branch("pi0_decay_photon_count", &total_pi0_decay_photons);
  metadata.Branch("hepmc_pi0_decay_photon_count", &total_hepmc_pi0_decay_photons);
  metadata.Branch("g4_pi0_decay_photon_count", &total_g4_pi0_decay_photons);
  metadata.Fill();
  metadata.Write();
  output.Close();
  if (output.TestBit(TFile::kWriteError)) {
    return 6;
  }

  std::cout << "FinalizePythiaTruthPtSpectra - partials/files/events/prompt/pi0/"
               "pi0 decay (total/HepMC/G4) = "
            << partial_file_count << "/" << total_input_files << "/" << total_events << "/"
            << total_prompt_photons << "/" << total_pi0 << "/"
            << total_pi0_decay_photons << "/" << total_hepmc_pi0_decay_photons << "/"
            << total_g4_pi0_decay_photons << std::endl;
  std::cout << "Wrote " << output_base << ".pdf and " << output_base << ".root" << std::endl;
  return 0;
}
