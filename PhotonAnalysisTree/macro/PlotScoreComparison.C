#include "Utilities/sPhenixStyle.C"

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TPad.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace
{
constexpr int score_nbins = 50;
constexpr double score_min = 0.0;
constexpr double score_max = 1.0;

struct HistogramPair
{
  std::unique_ptr<TH1D> bdt;
  std::unique_ptr<TH1D> mlp;
};

HistogramPair make_histogram_pair(const std::string& suffix)
{
  HistogramPair histograms;
  histograms.bdt = std::make_unique<TH1D>(
      ("h_bdt_" + suffix).c_str(), "", score_nbins, score_min, score_max);
  histograms.mlp = std::make_unique<TH1D>(
      ("h_mlp_" + suffix).c_str(), "", score_nbins, score_min, score_max);

  for (TH1D* histogram : {histograms.bdt.get(), histograms.mlp.get()})
  {
    histogram->SetDirectory(nullptr);
    histogram->Sumw2();
    histogram->SetStats(false);
    histogram->SetFillStyle(0);
  }

  histograms.bdt->SetLineColor(kBlue + 1);
  histograms.bdt->SetLineStyle(1);
  histograms.mlp->SetLineColor(kRed + 1);
  histograms.mlp->SetLineStyle(2);

  return histograms;
}

void normalize(HistogramPair& histograms)
{
  for (TH1D* histogram : {histograms.bdt.get(), histograms.mlp.get()})
  {
    if (histogram->GetEntries() > 0.0)
    {
      histogram->Scale(1.0 / histogram->GetEntries());
    }
  }
}

void configure_axes(HistogramPair& histograms)
{
  TH1D* frame = histograms.bdt.get();
  frame->GetXaxis()->SetTitle("Classifier score");
  frame->GetYaxis()->SetTitle("Fraction of clusters / 0.02");
  frame->GetXaxis()->SetRangeUser(score_min, score_max);
  frame->SetMinimum(0.0);

  const double maximum =
      std::max(histograms.bdt->GetMaximum(), histograms.mlp->GetMaximum());
  frame->SetMaximum(maximum > 0.0 ? 1.25 * maximum : 1.0);
}

void draw_histograms(HistogramPair& histograms)
{
  configure_axes(histograms);
  histograms.bdt->Draw("HIST");
  histograms.mlp->Draw("HIST SAME");
}

void draw_legend(HistogramPair& histograms,
                 const std::string& collection_label,
                 const double x1 = 0.20,
                 const double y1 = 0.67,
                 const double x2 = 0.88,
                 const double y2 = 0.84)
{
  TLegend legend(x1, y1, x2, y2);
  legend.AddEntry(histograms.bdt.get(), (collection_label + " BDT").c_str(), "l");
  legend.AddEntry(histograms.mlp.get(), (collection_label + " MLP").c_str(), "l");
  legend.DrawClone();
}

void draw_sphenix_label(const char* sample_text)
{
  TLatex label;
  label.SetNDC();
  label.SetTextAlign(13);
  label.DrawLatex(0.20, 0.91, "#it{#bf{sPHENIX}} Internal");
  label.DrawLatex(0.20, 0.83, sample_text);
}

void draw_panel_label(const std::string& text)
{
  TLatex label;
  label.SetNDC();
  label.SetTextAlign(13);
  label.DrawLatex(0.22, 0.92, text.c_str());
}

void draw_information_panel(HistogramPair& histograms,
                            const std::string& collection_label,
                            const char* binning_text,
                            const bool conditions)
{
  if (conditions)
  {
    TLatex label;
    label.SetNDC();
    label.SetTextAlign(13);
    label.DrawLatex(0.20, 0.91, "#it{#bf{sPHENIX}} Internal");
    label.DrawLatex(0.20, 0.82, "Single #pi^{0} gun");
    label.DrawLatex(
        0.20, 0.73, (collection_label + " common valid clusters").c_str());
    label.DrawLatex(0.20, 0.64, binning_text);
    return;
  }

  draw_legend(histograms, collection_label, 0.15, 0.35, 0.90, 0.65);
}

std::size_t find_bin(const double value, const std::vector<double>& edges)
{
  if (!std::isfinite(value) || value < edges.front() || value > edges.back())
  {
    return edges.size();
  }

  if (value == edges.back())
  {
    return edges.size() - 2U;
  }

  const auto upper = std::upper_bound(edges.begin(), edges.end(), value);
  return static_cast<std::size_t>(std::distance(edges.begin(), upper) - 1);
}

bool make_output_directory(const std::string& output_base)
{
  const std::size_t slash = output_base.find_last_of('/');
  if (slash == std::string::npos)
  {
    return true;
  }

  const std::string directory = output_base.substr(0, slash);
  if (directory.empty() || !gSystem->AccessPathName(directory.c_str()))
  {
    return true;
  }

  if (gSystem->mkdir(directory.c_str(), true) != 0)
  {
    std::cerr << "Failed to create output directory " << directory << std::endl;
    return false;
  }
  return true;
}

std::string collection_output_base(const std::string& output_base,
                                   const std::string& collection)
{
  const std::size_t slash = output_base.find_last_of("/");
  const std::string directory =
      slash == std::string::npos ? "" : output_base.substr(0, slash + 1U);
  const std::string stem =
      slash == std::string::npos ? output_base : output_base.substr(slash + 1U);
  return directory + collection + "/" + stem + "_" + collection;
}

int plot_collection(TTree* tree,
                    const std::string& collection,
                    const std::string& output_base)
{
  const std::string collection_label =
      collection == "nosplit" ? "NO_SPLIT" : "SPLIT";

  double truth_pt = 0.0;
  std::vector<double>* cluster_et = nullptr;
  std::vector<float>* bdt_score = nullptr;
  std::vector<unsigned char>* bdt_valid = nullptr;
  std::vector<float>* mlp_score = nullptr;
  std::vector<unsigned char>* mlp_valid = nullptr;

  tree->ResetBranchAddresses();
  tree->SetBranchStatus("*", false);

  bool branches_ok = true;
  const auto bind = [&](const std::string& branch_name, auto* address)
  {
    if (!tree->GetBranch(branch_name.c_str()))
    {
      std::cerr << "Missing branch " << branch_name << std::endl;
      branches_ok = false;
      return;
    }
    tree->SetBranchStatus(branch_name.c_str(), true);
    tree->SetBranchAddress(branch_name.c_str(), address);
  };

  bind("truth_pt", &truth_pt);
  bind(collection + "_cluster_et", &cluster_et);
  bind(collection + "_cluster_bdt_base_v3E_score", &bdt_score);
  bind(collection + "_cluster_bdt_base_v3E_valid", &bdt_valid);
  bind(collection + "_cluster_p_gamma", &mlp_score);
  bind(collection + "_cluster_p_gamma_valid", &mlp_valid);
  if (!branches_ok)
  {
    return 1;
  }

  HistogramPair inclusive = make_histogram_pair(collection + "_inclusive");

  const std::vector<double> truth_pt_edges = {
      3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0};
  const std::array<std::string, 6> truth_pt_labels = {
      "3 #leq p_{T}^{truth} < 5 GeV",
      "5 #leq p_{T}^{truth} < 7 GeV",
      "7 #leq p_{T}^{truth} < 9 GeV",
      "9 #leq p_{T}^{truth} < 11 GeV",
      "11 #leq p_{T}^{truth} < 13 GeV",
      "13 #leq p_{T}^{truth} #leq 15 GeV"};
  std::vector<HistogramPair> truth_pt_histograms;
  truth_pt_histograms.reserve(truth_pt_labels.size());
  for (std::size_t bin = 0; bin < truth_pt_labels.size(); ++bin)
  {
    truth_pt_histograms.push_back(
        make_histogram_pair(collection + "_truth_pt_" + std::to_string(bin)));
  }

  const std::vector<double> cluster_et_edges = {
      3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0};
  const std::array<std::string, 6> cluster_et_labels = {
      "3 #leq E_{T}^{cluster} < 5 GeV",
      "5 #leq E_{T}^{cluster} < 7 GeV",
      "7 #leq E_{T}^{cluster} < 9 GeV",
      "9 #leq E_{T}^{cluster} < 11 GeV",
      "11 #leq E_{T}^{cluster} < 13 GeV",
      "13 #leq E_{T}^{cluster} #leq 15 GeV"};
  std::vector<HistogramPair> cluster_et_histograms;
  cluster_et_histograms.reserve(cluster_et_labels.size());
  for (std::size_t bin = 0; bin < cluster_et_labels.size(); ++bin)
  {
    cluster_et_histograms.push_back(
        make_histogram_pair(collection + "_cluster_et_" + std::to_string(bin)));
  }

  Long64_t total_clusters = 0;
  Long64_t common_valid_clusters = 0;
  Long64_t malformed_events = 0;
  Long64_t truth_pt_outside_bins = 0;

  const Long64_t entries = tree->GetEntries();
  for (Long64_t entry = 0; entry < entries; ++entry)
  {
    if (entry % 1000 == 0)
    {
      std::cout << collection_label << " entry: " << entry << std::endl;
    }
    tree->GetEntry(entry);

    if (!cluster_et || !bdt_score || !bdt_valid || !mlp_score || !mlp_valid)
    {
      ++malformed_events;
      continue;
    }

    const std::size_t ncluster = cluster_et->size();
    total_clusters += static_cast<Long64_t>(ncluster);
    if (bdt_score->size() != ncluster || bdt_valid->size() != ncluster ||
        mlp_score->size() != ncluster || mlp_valid->size() != ncluster)
    {
      ++malformed_events;
      continue;
    }

    const std::size_t truth_pt_bin = find_bin(truth_pt, truth_pt_edges);
    for (std::size_t cluster = 0; cluster < ncluster; ++cluster)
    {
      if (!(*bdt_valid)[cluster] || !(*mlp_valid)[cluster])
      {
        continue;
      }

      const double bdt = (*bdt_score)[cluster];
      const double mlp = (*mlp_score)[cluster];
      const double et = (*cluster_et)[cluster];
      if (!std::isfinite(bdt) || !std::isfinite(mlp) || !std::isfinite(et))
      {
        continue;
      }

      ++common_valid_clusters;
      inclusive.bdt->Fill(bdt);
      inclusive.mlp->Fill(mlp);

      if (truth_pt_bin < truth_pt_histograms.size())
      {
        truth_pt_histograms[truth_pt_bin].bdt->Fill(bdt);
        truth_pt_histograms[truth_pt_bin].mlp->Fill(mlp);
      }
      else
      {
        ++truth_pt_outside_bins;
      }

      const std::size_t cluster_et_bin = find_bin(et, cluster_et_edges);
      if (cluster_et_bin < cluster_et_histograms.size())
      {
        cluster_et_histograms[cluster_et_bin].bdt->Fill(bdt);
        cluster_et_histograms[cluster_et_bin].mlp->Fill(mlp);
      }
    }
  }

  normalize(inclusive);
  for (HistogramPair& histograms : truth_pt_histograms)
  {
    normalize(histograms);
  }
  for (HistogramPair& histograms : cluster_et_histograms)
  {
    normalize(histograms);
  }

  if (!make_output_directory(output_base))
  {
    return 1;
  }

  TCanvas inclusive_canvas(
      ("c_score_" + collection + "_inclusive").c_str(),
      (collection_label + " inclusive score comparison").c_str(), 1000, 800);
  draw_histograms(inclusive);
  draw_legend(inclusive, collection_label, 0.55, 0.55, 0.90, 0.72);
  draw_sphenix_label(
      ("Single #pi^{0} gun, " + collection_label + " common valid clusters").c_str());
  inclusive_canvas.RedrawAxis();
  inclusive_canvas.SaveAs((output_base + "_inclusive.pdf").c_str());

  TCanvas truth_pt_canvas(
      ("c_score_" + collection + "_truth_pt").c_str(),
      (collection_label + " score comparison by truth pi0 pT").c_str(), 1350,
      1350);
  truth_pt_canvas.Divide(3, 3);
  for (std::size_t bin = 0; bin < truth_pt_histograms.size(); ++bin)
  {
    truth_pt_canvas.cd(static_cast<int>(bin + 1U));
    draw_histograms(truth_pt_histograms[bin]);
    draw_panel_label(truth_pt_labels[bin]);
    gPad->RedrawAxis();
  }
  truth_pt_canvas.cd(7);
  draw_information_panel(
      truth_pt_histograms.front(), collection_label,
      "Binned in truth #pi^{0} p_{T}", true);
  truth_pt_canvas.cd(8);
  draw_information_panel(
      truth_pt_histograms.front(), collection_label,
      "Binned in truth #pi^{0} p_{T}", false);
  truth_pt_canvas.SaveAs((output_base + "_truth_pt.pdf").c_str());

  TCanvas cluster_et_canvas(
      ("c_score_" + collection + "_cluster_et").c_str(),
      (collection_label + " score comparison by cluster ET").c_str(), 1350,
      1350);
  cluster_et_canvas.Divide(3, 3);
  for (std::size_t bin = 0; bin < cluster_et_histograms.size(); ++bin)
  {
    cluster_et_canvas.cd(static_cast<int>(bin + 1U));
    draw_histograms(cluster_et_histograms[bin]);
    draw_panel_label(cluster_et_labels[bin]);
    gPad->RedrawAxis();
  }
  cluster_et_canvas.cd(7);
  const std::string cluster_binning_text =
      "Binned in " + collection_label + " cluster E_{T}";
  draw_information_panel(
      cluster_et_histograms.front(), collection_label,
      cluster_binning_text.c_str(), true);
  cluster_et_canvas.cd(8);
  draw_information_panel(
      cluster_et_histograms.front(), collection_label,
      cluster_binning_text.c_str(), false);
  cluster_et_canvas.SaveAs((output_base + "_cluster_et.pdf").c_str());

  std::cout << "PlotScoreComparison (" << collection_label
            << ") - events/total clusters/common valid clusters/"
               "malformed events/truth-pT out-of-range clusters = "
            << entries << "/" << total_clusters << "/" << common_valid_clusters << "/"
            << malformed_events << "/" << truth_pt_outside_bins << std::endl;
  std::cout << "Wrote " << output_base << "_inclusive.pdf, "
            << output_base << "_truth_pt.pdf, and "
            << output_base << "_cluster_et.pdf" << std::endl;

  return malformed_events == 0 ? 0 : 1;
}
}  // namespace

int PlotScoreComparison(
    const std::string input_path = "PhotonAnalysisTree/output/merged/100kevents_pi0_3to15GeV_etapm1_vertexpm60_newBDTprediction.root",
    // const std::string input_path = "PhotonAnalysisTree/output/merged/100kevents_pi0_5to15GeV_etapm1.root",
    // const std::string input_path = "PhotonAnalysisTree/output/merged/100segments_Jet5.root",
    const std::string output_base =
        "PhotonAnalysisTree/output/plots/score_comparison/score_comparison")
{
  SetsPhenixStyle();
  TH1::AddDirectory(false);

  std::unique_ptr<TFile> input(TFile::Open(input_path.c_str(), "READ"));
  if (!input || input->IsZombie())
  {
    std::cerr << "Failed to open " << input_path << std::endl;
    return 1;
  }

  TTree* tree = dynamic_cast<TTree*>(input->Get("event_tree"));
  if (!tree)
  {
    std::cerr << "Missing TTree event_tree in " << input_path << std::endl;
    return 1;
  }

  const int nosplit_status =
      plot_collection(tree, "nosplit",
                      collection_output_base(output_base, "nosplit"));
  const int split_status =
      plot_collection(tree, "split",
                      collection_output_base(output_base, "split"));

  return nosplit_status == 0 && split_status == 0 ? 0 : 1;
}
