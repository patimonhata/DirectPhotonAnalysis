#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TMath.h>
#include <TNamed.h>
#include <TStyle.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <limits>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace
{
struct EventData
{
  unsigned long long event_uid = 0;
  unsigned int ncluster = 0;
  unsigned int ncluster_all = 0;
  std::vector<double> *cluster_e = nullptr;
  std::vector<double> *cluster_px = nullptr;
  std::vector<double> *cluster_py = nullptr;
  std::vector<double> *cluster_pz = nullptr;
};

struct PairInfo
{
  double mass = std::numeric_limits<double>::quiet_NaN();
  double total_energy = std::numeric_limits<double>::quiet_NaN();
  double opening_angle = std::numeric_limits<double>::quiet_NaN();
  double energy_asymmetry = std::numeric_limits<double>::quiet_NaN();
};

bool SetBranches(TTree *tree, EventData &data)
{
  if (!tree || !tree->GetBranch("event_uid"))
  {
    std::cerr << "Missing required event_uid branch" << std::endl;
    return false;
  }

  tree->SetBranchAddress("event_uid", &data.event_uid);
  tree->SetBranchAddress("ncluster", &data.ncluster);
  tree->SetBranchAddress("ncluster_all", &data.ncluster_all);
  tree->SetBranchAddress("cluster_e", &data.cluster_e);
  tree->SetBranchAddress("cluster_px", &data.cluster_px);
  tree->SetBranchAddress("cluster_py", &data.cluster_py);
  tree->SetBranchAddress("cluster_pz", &data.cluster_pz);
  return true;
}

PairInfo BuildOnlyPair(const EventData &data)
{
  PairInfo info;
  if (!data.cluster_e || !data.cluster_px || !data.cluster_py || !data.cluster_pz)
  {
    return info;
  }
  if (data.cluster_e->size() != 2 ||
      data.cluster_px->size() != 2 ||
      data.cluster_py->size() != 2 ||
      data.cluster_pz->size() != 2)
  {
    return info;
  }

  const double e0 = data.cluster_e->at(0);
  const double e1 = data.cluster_e->at(1);
  const double px0 = data.cluster_px->at(0);
  const double py0 = data.cluster_py->at(0);
  const double pz0 = data.cluster_pz->at(0);
  const double px1 = data.cluster_px->at(1);
  const double py1 = data.cluster_py->at(1);
  const double pz1 = data.cluster_pz->at(1);

  info.total_energy = e0 + e1;
  const double px = px0 + px1;
  const double py = py0 + py1;
  const double pz = pz0 + pz1;
  const double mass2 = info.total_energy * info.total_energy - px * px - py * py - pz * pz;
  info.mass = std::sqrt(std::max(0.0, mass2));
  info.energy_asymmetry = info.total_energy > 0.0 ? std::abs(e0 - e1) / info.total_energy : std::numeric_limits<double>::quiet_NaN();

  const double p0 = std::sqrt(px0 * px0 + py0 * py0 + pz0 * pz0);
  const double p1 = std::sqrt(px1 * px1 + py1 * py1 + pz1 * pz1);
  if (p0 > 0.0 && p1 > 0.0)
  {
    const double cos_angle = std::clamp((px0 * px1 + py0 * py1 + pz0 * pz1) / (p0 * p1), -1.0, 1.0);
    info.opening_angle = std::acos(cos_angle);
  }

  return info;
}

void ScaleToUnitArea(TH1D *hist)
{
  if (!hist)
  {
    return;
  }
  const double integral = hist->Integral();
  if (integral > 0.0)
  {
    hist->Scale(1.0 / integral);
  }
}
}  // namespace

int CompareSplitClusterMass(
    const std::string no_split_file = "/sphenix/user/ryotaro/DirectPhotonAnalysis/Pi0Reconstruction/output/100kevents_5GeV_pi0_reconstruction_NO_SPLIT_CLUSTERS_tree.root",
    const std::string split_file = "/sphenix/user/ryotaro/DirectPhotonAnalysis/Pi0Reconstruction/output/100kevents_5GeV_pi0_reconstruction_SPLIT_CLUSTERS_tree.root",
    const std::string output_file = "/sphenix/user/ryotaro/DirectPhotonAnalysis/Pi0Reconstruction/output/compare_split_cluster_mass_5GeV.root")
{
  TFile *no_split_input = TFile::Open(no_split_file.c_str(), "READ");
  TFile *split_input = TFile::Open(split_file.c_str(), "READ");
  if (!no_split_input || no_split_input->IsZombie())
  {
    std::cerr << "Failed to open " << no_split_file << std::endl;
    return 1;
  }
  if (!split_input || split_input->IsZombie())
  {
    std::cerr << "Failed to open " << split_file << std::endl;
    return 1;
  }

  TTree *no_split_tree = dynamic_cast<TTree *>(no_split_input->Get("event_tree"));
  TTree *split_tree = dynamic_cast<TTree *>(split_input->Get("event_tree"));
  if (!no_split_tree || !split_tree)
  {
    std::cerr << "Missing event_tree in one of the input files" << std::endl;
    return 1;
  }

  EventData no_split;
  EventData split;
  if (!SetBranches(no_split_tree, no_split) || !SetBranches(split_tree, split))
  {
    std::cerr << "Both input trees must contain unique event_uid values" << std::endl;
    no_split_input->Close();
    split_input->Close();
    return 1;
  }

  std::map<unsigned long long, Long64_t> split_entry_by_uid;
  const Long64_t split_entries = split_tree->GetEntries();
  for (Long64_t entry = 0; entry < split_entries; ++entry)
  {
    split_tree->GetEntry(entry);
    const auto result = split_entry_by_uid.emplace(split.event_uid, entry);
    if (!result.second)
    {
      std::cerr << "Duplicate event_uid in SPLIT tree: " << split.event_uid << std::endl;
      no_split_input->Close();
      split_input->Close();
      return 1;
    }
  }

  TH1D *h_m_no_split_n2 = new TH1D("h_m_no_split_n2", "NO_SPLIT, N_{cluster}=2;M_{#gamma#gamma} [GeV];Events", 120, 0.0, 0.30);
  TH1D *h_m_split_n2_given_no_split_n2 = new TH1D("h_m_split_n2_given_no_split_n2", "SPLIT, N_{cluster}=2, same NO_SPLIT N_{cluster}=2 events;M_{#gamma#gamma} [GeV];Events", 120, 0.0, 0.30);
  TH1D *h_delta_m_split_minus_no_split = new TH1D("h_delta_m_split_minus_no_split", "SPLIT - NO_SPLIT for events with both N_{cluster}=2;#Delta M_{#gamma#gamma} [GeV];Events", 120, -0.12, 0.12);
  TH1D *h_ncluster_split_given_no_split_n2 = new TH1D("h_ncluster_split_given_no_split_n2", "SPLIT N_{cluster} for NO_SPLIT N_{cluster}=2 events;SPLIT N_{cluster};Events", 20, -0.5, 19.5);
  TH1D *h_energy_sum_ratio = new TH1D("h_energy_sum_ratio", "E_{#gamma#gamma}^{SPLIT}/E_{#gamma#gamma}^{NO_SPLIT} for both N_{cluster}=2;Energy sum ratio;Events", 120, 0.5, 1.5);
  TH1D *h_opening_angle_ratio = new TH1D("h_opening_angle_ratio", "#theta_{#gamma#gamma}^{SPLIT}/#theta_{#gamma#gamma}^{NO_SPLIT} for both N_{cluster}=2;Opening-angle ratio;Events", 120, 0.5, 1.5);
  TH2D *h_m_split_vs_no_split = new TH2D("h_m_split_vs_no_split", "Both N_{cluster}=2;NO_SPLIT M_{#gamma#gamma} [GeV];SPLIT M_{#gamma#gamma} [GeV]", 120, 0.0, 0.30, 120, 0.0, 0.30);
  TH2D *h_delta_m_vs_no_split_m = new TH2D("h_delta_m_vs_no_split_m", "Both N_{cluster}=2;NO_SPLIT M_{#gamma#gamma} [GeV];SPLIT - NO_SPLIT [GeV]", 120, 0.0, 0.30, 120, -0.12, 0.12);

  unsigned long long no_split_n2 = 0;
  unsigned long long matched_events = 0;
  unsigned long long split_n2_given_no_split_n2 = 0;
  unsigned long long split_not_n2_given_no_split_n2 = 0;
  double sum_no_split_mass = 0.0;
  double sum_split_mass = 0.0;
  double sum_delta_mass = 0.0;
  double sum_energy_ratio = 0.0;
  double sum_opening_ratio = 0.0;
  std::set<unsigned long long> no_split_uids;

  const Long64_t no_split_entries = no_split_tree->GetEntries();
  for (Long64_t entry = 0; entry < no_split_entries; ++entry)
  {
    no_split_tree->GetEntry(entry);
    if (!no_split_uids.insert(no_split.event_uid).second)
    {
      std::cerr << "Duplicate event_uid in NO_SPLIT tree: " << no_split.event_uid << std::endl;
      no_split_input->Close();
      split_input->Close();
      return 1;
    }
    if (no_split.ncluster != 2)
    {
      continue;
    }

    ++no_split_n2;
    const PairInfo no_split_pair = BuildOnlyPair(no_split);
    if (!std::isfinite(no_split_pair.mass))
    {
      continue;
    }
    h_m_no_split_n2->Fill(no_split_pair.mass);
    sum_no_split_mass += no_split_pair.mass;

    const auto split_entry_iter = split_entry_by_uid.find(no_split.event_uid);
    if (split_entry_iter == split_entry_by_uid.end())
    {
      continue;
    }

    ++matched_events;
    split_tree->GetEntry(split_entry_iter->second);
    h_ncluster_split_given_no_split_n2->Fill(static_cast<double>(split.ncluster));
    if (split.ncluster != 2)
    {
      ++split_not_n2_given_no_split_n2;
      continue;
    }

    const PairInfo split_pair = BuildOnlyPair(split);
    if (!std::isfinite(split_pair.mass))
    {
      continue;
    }

    ++split_n2_given_no_split_n2;
    const double delta_mass = split_pair.mass - no_split_pair.mass;
    h_m_split_n2_given_no_split_n2->Fill(split_pair.mass);
    h_delta_m_split_minus_no_split->Fill(delta_mass);
    h_m_split_vs_no_split->Fill(no_split_pair.mass, split_pair.mass);
    h_delta_m_vs_no_split_m->Fill(no_split_pair.mass, delta_mass);

    sum_split_mass += split_pair.mass;
    sum_delta_mass += delta_mass;

    if (no_split_pair.total_energy > 0.0)
    {
      const double ratio = split_pair.total_energy / no_split_pair.total_energy;
      h_energy_sum_ratio->Fill(ratio);
      sum_energy_ratio += ratio;
    }
    if (no_split_pair.opening_angle > 0.0)
    {
      const double ratio = split_pair.opening_angle / no_split_pair.opening_angle;
      h_opening_angle_ratio->Fill(ratio);
      sum_opening_ratio += ratio;
    }
  }

  TFile *output = TFile::Open(output_file.c_str(), "RECREATE");
  if (!output || output->IsZombie())
  {
    std::cerr << "Failed to open output " << output_file << std::endl;
    return 1;
  }

  output->cd();
  h_m_no_split_n2->Write();
  h_m_split_n2_given_no_split_n2->Write();
  h_delta_m_split_minus_no_split->Write();
  h_ncluster_split_given_no_split_n2->Write();
  h_energy_sum_ratio->Write();
  h_opening_angle_ratio->Write();
  h_m_split_vs_no_split->Write();
  h_delta_m_vs_no_split_m->Write();

  TH1D *h_m_no_split_n2_norm = static_cast<TH1D *>(h_m_no_split_n2->Clone("h_m_no_split_n2_norm"));
  TH1D *h_m_split_n2_given_no_split_n2_norm = static_cast<TH1D *>(h_m_split_n2_given_no_split_n2->Clone("h_m_split_n2_given_no_split_n2_norm"));
  ScaleToUnitArea(h_m_no_split_n2_norm);
  ScaleToUnitArea(h_m_split_n2_given_no_split_n2_norm);
  h_m_no_split_n2_norm->SetLineColor(kBlue + 1);
  h_m_no_split_n2_norm->SetLineWidth(2);
  h_m_split_n2_given_no_split_n2_norm->SetLineColor(kRed + 1);
  h_m_split_n2_given_no_split_n2_norm->SetLineWidth(2);
  h_m_split_n2_given_no_split_n2_norm->SetLineStyle(2);

  TCanvas *canvas = new TCanvas("c_mass_overlay", "mass overlay", 900, 700);
  h_m_no_split_n2_norm->Draw("hist");
  h_m_split_n2_given_no_split_n2_norm->Draw("hist same");
  TLegend *legend = new TLegend(0.53, 0.72, 0.88, 0.88);
  legend->AddEntry(h_m_no_split_n2_norm, "NO_SPLIT N=2", "l");
  legend->AddEntry(h_m_split_n2_given_no_split_n2_norm, "SPLIT N=2 in same events", "l");
  legend->Draw();
  canvas->Write();
  const std::string pdf_file = output_file.substr(0, output_file.find_last_of('.')) + "_mass_overlay.pdf";
  canvas->SaveAs(pdf_file.c_str());

  std::ostringstream summary;
  summary << std::fixed << std::setprecision(6);
  summary << "no_split_file: " << no_split_file << "\n";
  summary << "split_file: " << split_file << "\n";
  summary << "output_file: " << output_file << "\n";
  summary << "no_split_entries: " << no_split_entries << "\n";
  summary << "split_entries: " << split_entries << "\n";
  summary << "no_split_ncluster_eq_2: " << no_split_n2 << "\n";
  summary << "matched_no_split_n2_events: " << matched_events << "\n";
  summary << "split_ncluster_eq_2_given_no_split_n2: " << split_n2_given_no_split_n2 << "\n";
  summary << "split_ncluster_not_2_given_no_split_n2: " << split_not_n2_given_no_split_n2 << "\n";
  summary << "mean_no_split_mass_for_no_split_n2: " << (no_split_n2 > 0 ? sum_no_split_mass / no_split_n2 : 0.0) << "\n";
  summary << "mean_split_mass_for_both_n2: " << (split_n2_given_no_split_n2 > 0 ? sum_split_mass / split_n2_given_no_split_n2 : 0.0) << "\n";
  summary << "mean_delta_mass_split_minus_no_split_for_both_n2: " << (split_n2_given_no_split_n2 > 0 ? sum_delta_mass / split_n2_given_no_split_n2 : 0.0) << "\n";
  summary << "mean_energy_sum_ratio_for_both_n2: " << (split_n2_given_no_split_n2 > 0 ? sum_energy_ratio / split_n2_given_no_split_n2 : 0.0) << "\n";
  summary << "mean_opening_angle_ratio_for_both_n2: " << (split_n2_given_no_split_n2 > 0 ? sum_opening_ratio / split_n2_given_no_split_n2 : 0.0) << "\n";

  TNamed summary_object("summary", summary.str().c_str());
  summary_object.Write();
  std::cout << summary.str();

  output->Close();
  no_split_input->Close();
  split_input->Close();

  return 0;
}
