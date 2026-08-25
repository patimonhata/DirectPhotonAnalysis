#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TNamed.h>
#include <TStyle.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace {
struct EventData {
  unsigned long long event_uid = 0;
  std::vector<double> *cluster_e = nullptr;
};

struct ClusterCountPair {
  unsigned int no_split = 0;
  unsigned int split = 0;
};

bool SetBranches(TTree *tree, EventData &data, const std::string &label) {
  if (!tree || !tree->GetBranch("event_uid") || !tree->GetBranch("cluster_e")) {
    std::cerr << label << " tree must contain event_uid and cluster_e branches" << std::endl;
    return false;
  }

  tree->SetBranchAddress("event_uid", &data.event_uid);
  tree->SetBranchAddress("cluster_e", &data.cluster_e);
  return true;
}

unsigned int CountClusters(const std::vector<double> *cluster_energies, bool use_energy_threshold, double min_cluster_energy)
{
  if (!cluster_energies) {
    return 0;
  }

  if (!use_energy_threshold) {
    return static_cast<unsigned int>(cluster_energies->size());
  }

  return static_cast<unsigned int>(std::count_if(
      cluster_energies->begin(),
      cluster_energies->end(),
      [min_cluster_energy](double energy)
      {
        return std::isfinite(energy) && energy >= min_cluster_energy;
      }));
}

std::string PdfFileName(const std::string &root_file_name) {
  const std::string suffix = ".root";
  if (root_file_name.size() >= suffix.size() &&
      root_file_name.compare(
          root_file_name.size() - suffix.size(),
          suffix.size(),
          suffix) == 0)
  {
    return root_file_name.substr(0, root_file_name.size() - suffix.size()) + ".pdf";
  }
  return root_file_name + ".pdf";
}
}  // namespace

int CompareSplitClusterCounts(
    const std::string no_split_file =
        "/sphenix/user/ryotaro/DirectPhotonAnalysis/Pi0Reconstruction/output/100kevents_5GeV_pi0_eta0_towerinfo_NO_SPLIT_CLUSTERS_UID_tree.root",
    const std::string split_file =
        "/sphenix/user/ryotaro/DirectPhotonAnalysis/Pi0Reconstruction/output/100kevents_5GeV_pi0_eta0_towerinfo_SPLIT_CLUSTERS_UID_tree.root",
    const std::string output_file = "/sphenix/user/ryotaro/DirectPhotonAnalysis/Pi0Reconstruction/output/compare_split_cluster_counts_5GeV.root",
    const double min_cluster_energy = std::numeric_limits<double>::quiet_NaN())
{
  const bool use_energy_threshold = !std::isnan(min_cluster_energy);

  TFile *no_split_input = TFile::Open(no_split_file.c_str(), "READ");
  TFile *split_input = TFile::Open(split_file.c_str(), "READ");
  if (!no_split_input || no_split_input->IsZombie()) {
    std::cerr << "Failed to open " << no_split_file << std::endl;
    if (split_input) {
      split_input->Close();
    }
    return 1;
  }
  if (!split_input || split_input->IsZombie()) {
    std::cerr << "Failed to open " << split_file << std::endl;
    no_split_input->Close();
    return 1;
  }

  TTree *no_split_tree = dynamic_cast<TTree *>(no_split_input->Get("event_tree"));
  TTree *split_tree = dynamic_cast<TTree *>(split_input->Get("event_tree"));
  if (!no_split_tree || !split_tree) {
    std::cerr << "Missing event_tree in one of the input files" << std::endl;
    no_split_input->Close();
    split_input->Close();
    return 1;
  }

  EventData no_split;
  EventData split;
  if (!SetBranches(no_split_tree, no_split, "NO_SPLIT") || !SetBranches(split_tree, split, "SPLIT")) {
    no_split_input->Close();
    split_input->Close();
    return 1;
  }

  std::map<unsigned long long, unsigned int> split_count_by_uid;
  const Long64_t split_entries = split_tree->GetEntries();
  for (Long64_t entry = 0; entry < split_entries; ++entry) {
    split_tree->GetEntry(entry);
    if (!split.cluster_e) {
      std::cerr << "Null cluster_e in SPLIT entry " << entry << std::endl;
      no_split_input->Close();
      split_input->Close();
      return 1;
    }

    const unsigned int count = CountClusters( split.cluster_e, use_energy_threshold, min_cluster_energy);
    const auto result = split_count_by_uid.emplace(split.event_uid, count);
    if (!result.second) {
      std::cerr << "Duplicate event_uid in SPLIT tree: " << split.event_uid << std::endl;
      no_split_input->Close();
      split_input->Close();
      return 1;
    }
  }

  std::set<unsigned long long> no_split_uids;
  std::vector<ClusterCountPair> matched_pairs;
  unsigned long long unmatched_no_split = 0;
  const Long64_t no_split_entries = no_split_tree->GetEntries();
  matched_pairs.reserve(static_cast<std::size_t>(std::min(no_split_entries, split_entries)));

  for (Long64_t entry = 0; entry < no_split_entries; ++entry) {
    no_split_tree->GetEntry(entry);
    if (!no_split.cluster_e) {
      std::cerr << "Null cluster_e in NO_SPLIT entry " << entry << std::endl;
      no_split_input->Close();
      split_input->Close();
      return 1;
    }
    if (!no_split_uids.insert(no_split.event_uid).second) {
      std::cerr << "Duplicate event_uid in NO_SPLIT tree: " << no_split.event_uid << std::endl;
      no_split_input->Close();
      split_input->Close();
      return 1;
    }

    const auto split_iter = split_count_by_uid.find(no_split.event_uid);
    if (split_iter == split_count_by_uid.end()) {
      ++unmatched_no_split;
      continue;
    }

    const unsigned int no_split_count = CountClusters(no_split.cluster_e, use_energy_threshold, min_cluster_energy);
    matched_pairs.push_back({no_split_count, split_iter->second});
    split_count_by_uid.erase(split_iter);
  }

  const unsigned long long unmatched_split = split_count_by_uid.size();
  no_split_input->Close();
  split_input->Close();

  if (matched_pairs.empty()) {
    std::cerr << "No matching event_uid values were found" << std::endl;
    return 1;
  }

  unsigned int max_cluster_count = 0;
  unsigned int max_abs_delta = 0;
  unsigned long long unchanged_events = 0;
  unsigned long long increased_events = 0;
  unsigned long long decreased_events = 0;
  double sum_no_split = 0.0;
  double sum_split = 0.0;

  for (const ClusterCountPair &pair : matched_pairs) {
    max_cluster_count = std::max(max_cluster_count, std::max(pair.no_split, pair.split));
    const int delta = static_cast<int>(pair.split) - static_cast<int>(pair.no_split);
    max_abs_delta = std::max( max_abs_delta, static_cast<unsigned int>(std::abs(delta)));

    sum_no_split += pair.no_split;
    sum_split += pair.split;
    if (delta > 0) {
      ++increased_events;
    } else if (delta < 0) {
      ++decreased_events;
    } else {
      ++unchanged_events;
    }
  }

  TH1::AddDirectory(false);
  const int ncluster_bins = static_cast<int>(max_cluster_count) + 1;
  TH1D *h_ncluster_no_split = new TH1D(
      "h_ncluster_no_split",
      "Cluster count;N_{cluster};Matched events",
      ncluster_bins,
      -0.5,
      static_cast<double>(max_cluster_count) + 0.5);
  TH1D *h_ncluster_split = new TH1D(
      "h_ncluster_split",
      "Cluster count;N_{cluster};Matched events",
      ncluster_bins,
      -0.5,
      static_cast<double>(max_cluster_count) + 0.5);
  TH2D *h_ncluster_migration = new TH2D(
      "h_ncluster_migration",
      "Cluster-count migration;NO_SPLIT N_{cluster};SPLIT N_{cluster}",
      ncluster_bins,
      -0.5,
      static_cast<double>(max_cluster_count) + 0.5,
      ncluster_bins,
      -0.5,
      static_cast<double>(max_cluster_count) + 0.5);

  const int delta_bins = 2 * static_cast<int>(max_abs_delta) + 1;
  TH1D *h_delta_ncluster = new TH1D(
      "h_delta_ncluster",
      "Change in cluster count;N_{cluster}^{SPLIT} - N_{cluster}^{NO_SPLIT};Matched events",
      delta_bins,
      -static_cast<double>(max_abs_delta) - 0.5,
      static_cast<double>(max_abs_delta) + 0.5);

  for (const ClusterCountPair &pair : matched_pairs) {
    h_ncluster_no_split->Fill(pair.no_split);
    h_ncluster_split->Fill(pair.split);
    h_ncluster_migration->Fill(pair.no_split, pair.split);
    h_delta_ncluster->Fill( static_cast<int>(pair.split) - static_cast<int>(pair.no_split));
  }

  TH2D *h_ncluster_migration_row_percent = static_cast<TH2D *>(h_ncluster_migration->Clone("h_ncluster_migration_row_percent"));
  h_ncluster_migration_row_percent->Reset("ICES");
  h_ncluster_migration_row_percent->SetTitle("Migration probability for each NO_SPLIT count;NO_SPLIT N_{cluster};SPLIT N_{cluster}");

  for (int xbin = 1; xbin <= h_ncluster_migration->GetNbinsX(); ++xbin) {
    double row_sum = 0.0;
    for (int ybin = 1; ybin <= h_ncluster_migration->GetNbinsY(); ++ybin) {
      row_sum += h_ncluster_migration->GetBinContent(xbin, ybin);
    }
    if (row_sum <= 0.0) {
      continue;
    }
    for (int ybin = 1; ybin <= h_ncluster_migration->GetNbinsY(); ++ybin){
      const double percent = 100.0 * h_ncluster_migration->GetBinContent(xbin, ybin) / row_sum;
      h_ncluster_migration_row_percent->SetBinContent(xbin, ybin, percent);
    }
  }

  h_ncluster_no_split->SetLineColor(kBlue + 1);
  h_ncluster_no_split->SetLineWidth(2);
  h_ncluster_split->SetLineColor(kRed + 1);
  h_ncluster_split->SetLineWidth(2);
  h_ncluster_split->SetLineStyle(2);
  h_delta_ncluster->SetLineColor(kBlack);
  h_delta_ncluster->SetLineWidth(2);
  h_ncluster_migration->SetMarkerSize(max_cluster_count > 12 ? 0.55 : 0.8);
  h_ncluster_migration_row_percent->SetMarkerSize( max_cluster_count > 12 ? 0.55 : 0.8);

  std::ostringstream selection;
  if (use_energy_threshold) {
    selection << "E_{cluster} #geq " << std::fixed << std::setprecision(3) << min_cluster_energy << " GeV";
  } else {
    selection << "all saved clusters (negative energy included)";
  }

  TCanvas *canvas = new TCanvas("c_split_cluster_counts", "split cluster counts", 1400, 1000);
  canvas->Divide(2, 2);

  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat(".1f");

  canvas->cd(1);
  gPad->SetRightMargin(0.14);
  h_ncluster_migration->SetTitle(("Cluster-count migration, " + selection.str() + ";NO_SPLIT N_{cluster};SPLIT N_{cluster}").c_str());
  h_ncluster_migration->Draw("COLZ TEXT");

  canvas->cd(2);
  const double overlay_max = 1.15 * std::max(h_ncluster_no_split->GetMaximum(), h_ncluster_split->GetMaximum());
  h_ncluster_no_split->SetMaximum(overlay_max);
  h_ncluster_no_split->SetTitle(("Cluster-count distributions, " + selection.str() + ";N_{cluster};Matched events").c_str());
  h_ncluster_no_split->Draw("HIST");
  h_ncluster_split->Draw("HIST SAME");
  TLegend *legend = new TLegend(0.62, 0.72, 0.88, 0.88);
  legend->AddEntry(h_ncluster_no_split, "NO_SPLIT", "l");
  legend->AddEntry(h_ncluster_split, "SPLIT", "l");
  legend->Draw();

  canvas->cd(3);
  h_delta_ncluster->SetTitle(("Cluster-count change, " + selection.str() + ";N_{cluster}^{SPLIT} - N_{cluster}^{NO_SPLIT};Matched events").c_str());
  h_delta_ncluster->Draw("HIST");

  canvas->cd(4);
  gPad->SetRightMargin(0.14);
  h_ncluster_migration_row_percent->SetTitle(("Row-normalized migration [%], " + selection.str() + ";NO_SPLIT N_{cluster};SPLIT N_{cluster}").c_str());
  h_ncluster_migration_row_percent->Draw("COLZ TEXT");

  TFile *output = TFile::Open(output_file.c_str(), "RECREATE");
  if (!output || output->IsZombie()) {
    std::cerr << "Failed to open output " << output_file << std::endl;
    return 1;
  }

  output->cd();
  h_ncluster_no_split->Write();
  h_ncluster_split->Write();
  h_ncluster_migration->Write();
  h_ncluster_migration_row_percent->Write();
  h_delta_ncluster->Write();
  canvas->Write();

  const double matched_count = static_cast<double>(matched_pairs.size());
  std::ostringstream summary;
  summary << std::fixed << std::setprecision(6);
  summary << "no_split_file: " << no_split_file << "\n";
  summary << "split_file: " << split_file << "\n";
  summary << "output_file: " << output_file << "\n";
  if (use_energy_threshold) {
    summary << "min_cluster_energy: " << min_cluster_energy << "\n";
  } else {
    summary << "min_cluster_energy: all saved clusters"
            << " (negative energy included)\n";
  }
  summary << "no_split_entries: " << no_split_entries << "\n";
  summary << "split_entries: " << split_entries << "\n";
  summary << "matched_events: " << matched_pairs.size() << "\n";
  summary << "unmatched_no_split_events: " << unmatched_no_split << "\n";
  summary << "unmatched_split_events: " << unmatched_split << "\n";
  summary << "unchanged_cluster_count_events: " << unchanged_events << "\n";
  summary << "increased_cluster_count_events: " << increased_events << "\n";
  summary << "decreased_cluster_count_events: " << decreased_events << "\n";
  summary << "mean_no_split_cluster_count: " << sum_no_split / matched_count << "\n";
  summary << "mean_split_cluster_count: " << sum_split / matched_count << "\n";
  summary << "mean_delta_cluster_count: " << (sum_split - sum_no_split) / matched_count << "\n";

  TNamed summary_object("summary", summary.str().c_str());
  summary_object.Write();

  const std::string pdf_file = PdfFileName(output_file);
  canvas->SaveAs(pdf_file.c_str());
  output->Close();

  std::cout << summary.str();
  std::cout << "pdf_file: " << pdf_file << std::endl;
  return 0;
}
