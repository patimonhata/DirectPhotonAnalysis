#include "Utilities/sPhenixStyle.C"

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace
{
template <class T>
bool bind_branch(TTree* tree, const char* name, T* address)
{
  if (!tree->GetBranch(name))
  {
    std::cerr << "PlotPythiaTruthPtSpectra - missing branch: " << name << std::endl;
    return false;
  }
  tree->SetBranchStatus(name, true);
  return tree->SetBranchAddress(name, address) >= 0;
}

template <class T>
bool aligned(const std::vector<T>* values, const std::size_t expected)
{
  return values && values->size() == expected;
}

bool in_acceptance(const float eta, const double max_abs_eta)
{
  return max_abs_eta < 0.0 || std::abs(static_cast<double>(eta)) < max_abs_eta;
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
}

int PlotPythiaTruthPtSpectra(
    const std::string input_pattern = "output/truth_root/pythia_truth_spectrum_tree_*.root",
    const std::string output_base = "output/plots/pythia_truth_pt_spectra",
    const int n_bins = 100,
    const double pt_max = 20.0,
    const double max_abs_eta = -1.0,
    const bool use_event_weight = false)
{
  if (input_pattern.empty() || output_base.empty() || n_bins <= 0 || pt_max <= 0.0 ||
      !std::isfinite(pt_max) || !std::isfinite(max_abs_eta))
  {
    std::cerr << "PlotPythiaTruthPtSpectra - invalid argument" << std::endl;
    return 1;
  }
  SetsPhenixStyle();
  TH1::AddDirectory(false);
  if (!make_output_directory(output_base))
  {
    std::cerr << "PlotPythiaTruthPtSpectra - failed to create output directory" << std::endl;
    return 1;
  }

  TChain tree("event_tree");
  const int n_files = tree.Add(input_pattern.c_str());
  if (n_files <= 0 || tree.GetEntries() <= 0)
  {
    std::cerr << "PlotPythiaTruthPtSpectra - no event_tree entries matched "
              << input_pattern << std::endl;
    return 2;
  }
  tree.SetBranchStatus("*", false);

  UChar_t event_weight_valid = 0U;
  double event_weight = 1.0;
  UInt_t truth_photon_n = 0U;
  std::vector<unsigned char>* photon_kinematics_valid = nullptr;
  std::vector<float>* photon_pt = nullptr;
  std::vector<float>* photon_eta = nullptr;
  std::vector<unsigned char>* photon_classification_valid = nullptr;
  std::vector<int>* photon_category = nullptr;
  UInt_t truth_pi0_n = 0U;
  std::vector<unsigned char>* pi0_kinematics_valid = nullptr;
  std::vector<float>* pi0_pt = nullptr;
  std::vector<float>* pi0_eta = nullptr;
  UInt_t truth_pi0_decay_photon_n = 0U;
  std::vector<unsigned char>* pi0_decay_photon_kinematics_valid = nullptr;
  std::vector<float>* pi0_decay_photon_pt = nullptr;
  std::vector<float>* pi0_decay_photon_eta = nullptr;

  bool ok = true;
  ok &= bind_branch(&tree, "event_weight_valid", &event_weight_valid);
  ok &= bind_branch(&tree, "event_weight", &event_weight);
  ok &= bind_branch(&tree, "truth_photon_n", &truth_photon_n);
  ok &= bind_branch(&tree, "truth_photon_kinematics_valid", &photon_kinematics_valid);
  ok &= bind_branch(&tree, "truth_photon_pt", &photon_pt);
  ok &= bind_branch(&tree, "truth_photon_eta", &photon_eta);
  ok &= bind_branch(&tree, "truth_photon_classification_valid", &photon_classification_valid);
  ok &= bind_branch(&tree, "truth_photon_category", &photon_category);
  ok &= bind_branch(&tree, "truth_pi0_n", &truth_pi0_n);
  ok &= bind_branch(&tree, "truth_pi0_kinematics_valid", &pi0_kinematics_valid);
  ok &= bind_branch(&tree, "truth_pi0_pt", &pi0_pt);
  ok &= bind_branch(&tree, "truth_pi0_eta", &pi0_eta);
  ok &= bind_branch(&tree, "truth_pi0_decay_photon_n", &truth_pi0_decay_photon_n);
  ok &= bind_branch(&tree, "truth_pi0_decay_photon_kinematics_valid", &pi0_decay_photon_kinematics_valid);
  ok &= bind_branch(&tree, "truth_pi0_decay_photon_pt", &pi0_decay_photon_pt);
  ok &= bind_branch(&tree, "truth_pi0_decay_photon_eta", &pi0_decay_photon_eta);
  if (!ok)
  {
    return 3;
  }

  TH1D prompt_photon("h_prompt_photon_truth_pt", "", n_bins, 0.0, pt_max);
  TH1D pi0("h_pi0_truth_pt", "", n_bins, 0.0, pt_max);
  TH1D pi0_decay_photon("h_pi0_decay_photon_truth_pt", "", n_bins, 0.0, pt_max);
  for (TH1D* histogram : {&prompt_photon, &pi0, &pi0_decay_photon})
  {
    histogram->Sumw2();
    histogram->SetStats(false);
    histogram->SetFillStyle(0);
    histogram->GetXaxis()->SetTitle("Truth p_{T} [GeV/#it{c}]");
    histogram->GetYaxis()->SetTitle(
        use_event_weight ? "Weighted counts / bin" : "Counts / bin");
  }
  prompt_photon.SetLineColor(kRed + 1);
  prompt_photon.SetLineWidth(3);
  pi0.SetLineColor(kBlue + 1);
  pi0.SetLineWidth(3);
  pi0_decay_photon.SetLineColor(kGreen + 2);
  pi0_decay_photon.SetLineWidth(3);

  ULong64_t n_prompt_photon = 0ULL;
  ULong64_t n_pi0 = 0ULL;
  ULong64_t n_pi0_decay_photon = 0ULL;
  ULong64_t malformed_events = 0ULL;
  ULong64_t invalid_weight_events = 0ULL;
  for (Long64_t entry = 0; entry < tree.GetEntries(); ++entry)
  {
    if (tree.GetEntry(entry) <= 0)
    {
      ++malformed_events;
      continue;
    }
    const std::size_t n_photon = truth_photon_n;
    const std::size_t n_pi0_event = truth_pi0_n;
    const std::size_t n_pi0_decay_photon_event = truth_pi0_decay_photon_n;
    const bool photon_aligned =
        aligned(photon_kinematics_valid, n_photon) && aligned(photon_pt, n_photon) &&
        aligned(photon_eta, n_photon) && aligned(photon_classification_valid, n_photon) &&
        aligned(photon_category, n_photon);
    const bool pi0_aligned = aligned(pi0_kinematics_valid, n_pi0_event) &&
        aligned(pi0_pt, n_pi0_event) && aligned(pi0_eta, n_pi0_event);
    const bool pi0_decay_photon_aligned =
        aligned(pi0_decay_photon_kinematics_valid, n_pi0_decay_photon_event) &&
        aligned(pi0_decay_photon_pt, n_pi0_decay_photon_event) &&
        aligned(pi0_decay_photon_eta, n_pi0_decay_photon_event);
    if (!photon_aligned || !pi0_aligned || !pi0_decay_photon_aligned)
    {
      ++malformed_events;
      continue;
    }
    if (use_event_weight && (!event_weight_valid || !std::isfinite(event_weight)))
    {
      ++invalid_weight_events;
      continue;
    }
    const double weight = use_event_weight ? event_weight : 1.0;

    for (std::size_t particle = 0; particle < n_photon; ++particle)
    {
      if (!(*photon_kinematics_valid)[particle] ||
          !in_acceptance((*photon_eta)[particle], max_abs_eta))
      {
        continue;
      }
      if ((*photon_classification_valid)[particle] &&
          ((*photon_category)[particle] == 1 || (*photon_category)[particle] == 2))
      {
        prompt_photon.Fill((*photon_pt)[particle], weight);
        ++n_prompt_photon;
      }
    }
    for (std::size_t particle = 0; particle < n_pi0_event; ++particle)
    {
      if ((*pi0_kinematics_valid)[particle] &&
          in_acceptance((*pi0_eta)[particle], max_abs_eta))
      {
        pi0.Fill((*pi0_pt)[particle], weight);
        ++n_pi0;
      }
    }
    for (std::size_t particle = 0; particle < n_pi0_decay_photon_event; ++particle)
    {
      if ((*pi0_decay_photon_kinematics_valid)[particle] &&
          in_acceptance((*pi0_decay_photon_eta)[particle], max_abs_eta))
      {
        pi0_decay_photon.Fill((*pi0_decay_photon_pt)[particle], weight);
        ++n_pi0_decay_photon;
      }
    }
  }
  if (malformed_events != 0ULL || (use_event_weight && invalid_weight_events != 0ULL))
  {
    std::cerr << "PlotPythiaTruthPtSpectra - malformed/invalid-weight events = "
              << malformed_events << "/" << invalid_weight_events << std::endl;
    return 4;
  }

  const double maximum = std::max(
      {prompt_photon.GetMaximum(), pi0.GetMaximum(), pi0_decay_photon.GetMaximum()});
  prompt_photon.SetMinimum(0.5);
  prompt_photon.SetMaximum(maximum > 0.0 ? 5.0 * maximum : 1.0);

  TCanvas canvas("c_pythia_truth_pt_spectra", "Pythia truth pT spectra", 1000, 800);
  canvas.SetLogy();
  prompt_photon.Draw("HIST");
  pi0.Draw("HIST SAME");
  pi0_decay_photon.Draw("HIST SAME");
  TLegend legend(0.55, 0.64, 0.89, 0.84);
  legend.AddEntry(&prompt_photon, "Prompt #gamma (direct + fragmentation)", "l");
  legend.AddEntry(&pi0, "Last-copy #pi^{0}", "l");
  legend.AddEntry(&pi0_decay_photon, "#gamma from #pi^{0}", "l");
  legend.Draw();
  TLatex label;
  label.SetNDC();
  label.SetTextAlign(13);
  label.DrawLatex(0.18, 0.92, "#it{#bf{sPHENIX}} Internal");
  label.DrawLatex(0.18, 0.84, "Pythia8 p+p minimum bias");
  const std::string selection = max_abs_eta < 0.0
      ? "No #eta selection" : "|#eta^{truth}| < " + std::to_string(max_abs_eta);
  label.DrawLatex(0.18, 0.76, selection.c_str());
  canvas.RedrawAxis();
  canvas.SaveAs((output_base + ".pdf").c_str());

  TFile histogram_file((output_base + ".root").c_str(), "RECREATE");
  prompt_photon.Write();
  pi0.Write();
  pi0_decay_photon.Write();
  histogram_file.Close();

  std::cout << "PlotPythiaTruthPtSpectra - files/events/prompt photons/pi0/pi0 decay photons = "
            << n_files << "/" << tree.GetEntries() << "/" << n_prompt_photon << "/"
            << n_pi0 << "/" << n_pi0_decay_photon << std::endl;
  std::cout << "Wrote " << output_base << ".pdf and " << output_base << ".root" << std::endl;
  return 0;
}
