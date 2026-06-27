// ROOT macro: interactive CEMC cluster viewer from CemcClusterDumper output
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace
{
  bool isFinitePosition(float eta, float phi)
  {
    return std::isfinite(eta) && std::isfinite(phi) && std::fabs(eta) < 100.0f && std::fabs(phi) < 100.0f;
  }

  int findHighestEnergyCluster(const std::vector<float>& energy)
  {
    int selected = -1;
    float best_energy = -std::numeric_limits<float>::infinity();
    for (size_t i = 0; i < energy.size(); ++i) {
      if (std::isfinite(energy[i]) && energy[i] > best_energy) {
        best_energy = energy[i];
        selected = static_cast<int>(i);
      }
    }
    return selected;
  }

  float getClusterZ(const std::string& zmode,
                    size_t index,
                    const std::vector<float>& energy,
                    const std::vector<float>& et)
  {
    if (zmode == "energy") {
      return energy[index];
    }
    return et[index];
  }

  float getMaxFinite(const std::vector<float>& values)
  {
    float max_value = 0.0f;
    for (const float value : values) {
      if (std::isfinite(value) && value > max_value) {
        max_value = value;
      }
    }
    return max_value;
  }
}

void DrawCemcClusterLego(const char* file_name = "",
                         const char* tree_name = "cemc_clusters",
                         const char* mode = "vtx",
                         const char* zmode = "et",
                         int cluster_index = -1)
{
  if (!file_name || std::string(file_name).empty()) {
    std::cout << "Usage: DrawCemcClusterLego(\"clusters.root\", \"cemc_clusters\", \"vtx\", \"et\", -1)" << std::endl;
    return;
  }

  TFile* file = TFile::Open(file_name, "READ");
  if (!file || file->IsZombie()) {
    std::cout << "Failed to open file: " << file_name << std::endl;
    return;
  }

  TTree* tree = dynamic_cast<TTree*>(file->Get(tree_name));
  if (!tree) {
    std::cout << "Failed to find TTree: " << tree_name << std::endl;
    return;
  }

  int event = 0;
  int has_vertex = 0;
  float vtx_x = 0.0f;
  float vtx_y = 0.0f;
  float vtx_z = 0.0f;

  std::vector<unsigned int>* cluster_id = nullptr;
  std::vector<float>* cluster_energy = nullptr;
  std::vector<float>* cluster_ecore = nullptr;
  std::vector<float>* cluster_eta0 = nullptr;
  std::vector<float>* cluster_phi0 = nullptr;
  std::vector<float>* cluster_et0 = nullptr;
  std::vector<float>* cluster_etavtx = nullptr;
  std::vector<float>* cluster_phivtx = nullptr;
  std::vector<float>* cluster_etvtx = nullptr;
  std::vector<int>* cluster_n_towers = nullptr;
  std::vector<int>* cluster_lead_tower_ieta = nullptr;
  std::vector<int>* cluster_lead_tower_iphi = nullptr;

  std::vector<int>* member_cluster_index = nullptr;
  std::vector<int>* member_tower_ieta = nullptr;
  std::vector<int>* member_tower_iphi = nullptr;
  std::vector<float>* member_energy = nullptr;
  std::vector<float>* member_energy_fraction = nullptr;

  tree->SetBranchStatus("*", 0);
  tree->SetBranchStatus("event", 1);
  tree->SetBranchStatus("has_vertex", 1);
  tree->SetBranchStatus("vtx_x", 1);
  tree->SetBranchStatus("vtx_y", 1);
  tree->SetBranchStatus("vtx_z", 1);
  tree->SetBranchStatus("cluster_id", 1);
  tree->SetBranchStatus("cluster_energy", 1);
  tree->SetBranchStatus("cluster_ecore", 1);
  tree->SetBranchStatus("cluster_eta0", 1);
  tree->SetBranchStatus("cluster_phi0", 1);
  tree->SetBranchStatus("cluster_et0", 1);
  tree->SetBranchStatus("cluster_etavtx", 1);
  tree->SetBranchStatus("cluster_phivtx", 1);
  tree->SetBranchStatus("cluster_etvtx", 1);
  tree->SetBranchStatus("cluster_n_towers", 1);
  tree->SetBranchStatus("cluster_lead_tower_ieta", 1);
  tree->SetBranchStatus("cluster_lead_tower_iphi", 1);
  tree->SetBranchStatus("member_cluster_index", 1);
  tree->SetBranchStatus("member_tower_ieta", 1);
  tree->SetBranchStatus("member_tower_iphi", 1);
  tree->SetBranchStatus("member_energy", 1);
  tree->SetBranchStatus("member_energy_fraction", 1);

  tree->SetBranchAddress("event", &event);
  tree->SetBranchAddress("has_vertex", &has_vertex);
  tree->SetBranchAddress("vtx_x", &vtx_x);
  tree->SetBranchAddress("vtx_y", &vtx_y);
  tree->SetBranchAddress("vtx_z", &vtx_z);
  tree->SetBranchAddress("cluster_id", &cluster_id);
  tree->SetBranchAddress("cluster_energy", &cluster_energy);
  tree->SetBranchAddress("cluster_ecore", &cluster_ecore);
  tree->SetBranchAddress("cluster_eta0", &cluster_eta0);
  tree->SetBranchAddress("cluster_phi0", &cluster_phi0);
  tree->SetBranchAddress("cluster_et0", &cluster_et0);
  tree->SetBranchAddress("cluster_etavtx", &cluster_etavtx);
  tree->SetBranchAddress("cluster_phivtx", &cluster_phivtx);
  tree->SetBranchAddress("cluster_etvtx", &cluster_etvtx);
  tree->SetBranchAddress("cluster_n_towers", &cluster_n_towers);
  tree->SetBranchAddress("cluster_lead_tower_ieta", &cluster_lead_tower_ieta);
  tree->SetBranchAddress("cluster_lead_tower_iphi", &cluster_lead_tower_iphi);
  tree->SetBranchAddress("member_cluster_index", &member_cluster_index);
  tree->SetBranchAddress("member_tower_ieta", &member_tower_ieta);
  tree->SetBranchAddress("member_tower_iphi", &member_tower_iphi);
  tree->SetBranchAddress("member_energy", &member_energy);
  tree->SetBranchAddress("member_energy_fraction", &member_energy_fraction);

  const std::string requested_mode(mode);
  const std::string requested_zmode(zmode);
  const bool use_energy_z = (requested_zmode == "energy");
  if (requested_zmode != "et" && requested_zmode != "energy") {
    std::cout << "Unknown zmode: " << zmode << ". Use zmode=\"et\" or \"energy\"." << std::endl;
    return;
  }
  if (requested_mode != "vtx" && requested_mode != "raw") {
    std::cout << "Unknown mode: " << mode << ". Use mode=\"vtx\" or \"raw\"." << std::endl;
    return;
  }

  gStyle->SetOptStat(0);

  TCanvas* canvas = new TCanvas("cemc_cluster_lego", "CEMC cluster viewer", 1500, 1000);
  canvas->Divide(2, 2);

  TH2F* h_cluster_map = new TH2F("h_cluster_map",
                                 ";#eta;#phi;cluster E_{T} [GeV]",
                                 96, -1.2, 1.2,
                                 256, -TMath::Pi(), TMath::Pi());
  TH1F* h_cluster_energy = new TH1F("h_cluster_energy",
                                    ";cluster energy [GeV];clusters",
                                    100, 0.0, 50.0);
  TH2F* h_member_energy = new TH2F("h_member_energy",
                                   ";tower i#eta;tower i#phi;member energy [GeV]",
                                   96, -0.5, 95.5,
                                   256, -0.5, 255.5);
  TH2F* h_member_fraction = new TH2F("h_member_fraction",
                                     ";tower i#eta;tower i#phi;member energy fraction",
                                     96, -0.5, 95.5,
                                     256, -0.5, 255.5);

  const Long64_t n_entries = tree->GetEntries();
  for (Long64_t ievt = 0; ievt < n_entries; ++ievt) {
    tree->GetEntry(ievt);

    h_cluster_map->Reset("ICES");
    h_cluster_energy->Reset("ICES");
    h_member_energy->Reset("ICES");
    h_member_fraction->Reset("ICES");

    const bool use_vtx = (requested_mode == "vtx" && has_vertex);
    const std::vector<float>& cluster_eta = use_vtx ? *cluster_etavtx : *cluster_eta0;
    const std::vector<float>& cluster_phi = use_vtx ? *cluster_phivtx : *cluster_phi0;
    const std::vector<float>& cluster_et = use_vtx ? *cluster_etvtx : *cluster_et0;
    const std::string used_mode = use_vtx ? "vtx" : "raw";
    const std::string z_label = use_energy_z ? "E" : "E_{T}";

    const int n_clusters = static_cast<int>(cluster_energy->size());
    int selected_cluster = cluster_index;
    if (selected_cluster < 0 || selected_cluster >= n_clusters) {
      selected_cluster = findHighestEnergyCluster(*cluster_energy);
    }

    h_cluster_map->GetZaxis()->SetTitle(use_energy_z ? "cluster energy [GeV]" : "cluster E_{T} [GeV]");
    const float max_z = use_energy_z ? getMaxFinite(*cluster_energy) : getMaxFinite(cluster_et);
    if (max_z > 0.0f) {
      h_cluster_map->SetMaximum(1.05f * max_z);
    }

    const float max_energy = getMaxFinite(*cluster_energy);
    h_cluster_energy->SetBins(100, 0.0, std::max(1.0f, 1.1f * max_energy));

    for (int iclus = 0; iclus < n_clusters; ++iclus) {
      if (!isFinitePosition(cluster_eta[iclus], cluster_phi[iclus])) {
        continue;
      }
      const float z = getClusterZ(requested_zmode, iclus, *cluster_energy, cluster_et);
      if (std::isfinite(z)) {
        h_cluster_map->Fill(cluster_eta[iclus], cluster_phi[iclus], z);
      }
      if (std::isfinite((*cluster_energy)[iclus])) {
        h_cluster_energy->Fill((*cluster_energy)[iclus]);
      }
    }

    float selected_member_sum = 0.0f;
    if (selected_cluster >= 0) {
      for (size_t imember = 0; imember < member_cluster_index->size(); ++imember) {
        if ((*member_cluster_index)[imember] != selected_cluster) {
          continue;
        }
        const float energy = (*member_energy)[imember];
        if (std::isfinite(energy)) {
          selected_member_sum += energy;
        }
      }

      for (size_t imember = 0; imember < member_cluster_index->size(); ++imember) {
        if ((*member_cluster_index)[imember] != selected_cluster) {
          continue;
        }
        const int ieta = (*member_tower_ieta)[imember];
        const int iphi = (*member_tower_iphi)[imember];
        const float energy = (*member_energy)[imember];
        if (!std::isfinite(energy)) {
          continue;
        }
        const float fraction = (selected_member_sum > 0.0f) ? energy / selected_member_sum : (*member_energy_fraction)[imember];
        h_member_energy->Fill(ieta, iphi, energy);
        h_member_fraction->Fill(ieta, iphi, fraction);
      }
    }

    int selected_id = -1;
    float selected_energy = 0.0f;
    float selected_et = 0.0f;
    int selected_ntowers = 0;
    int lead_ieta = -1;
    int lead_iphi = -1;
    if (selected_cluster >= 0) {
      selected_id = static_cast<int>((*cluster_id)[selected_cluster]);
      selected_energy = (*cluster_energy)[selected_cluster];
      selected_et = cluster_et[selected_cluster];
      selected_ntowers = (*cluster_n_towers)[selected_cluster];
      lead_ieta = (*cluster_lead_tower_ieta)[selected_cluster];
      lead_iphi = (*cluster_lead_tower_iphi)[selected_cluster];
    }

    h_cluster_map->SetTitle(Form("Event %d  clusters=%d  mode=%s  z=%s  selected index=%d",
                                 event, n_clusters, used_mode.c_str(), requested_zmode.c_str(), selected_cluster));
    h_cluster_energy->SetTitle(Form("Event %d cluster energy spectrum", event));
    h_member_energy->SetTitle(Form("Selected cluster id=%d  E=%.3f GeV  E_{T}=%.3f GeV  towers=%d  lead=(%d,%d)",
                                   selected_id, selected_energy, selected_et, selected_ntowers, lead_ieta, lead_iphi));
    h_member_fraction->SetTitle(Form("Selected cluster tower energy fraction  sum(member E)=%.3f GeV",
                                     selected_member_sum));

    canvas->cd(1);
    h_cluster_map->Draw("lego2");

    canvas->cd(2);
    h_cluster_energy->Draw("hist");

    canvas->cd(3);
    h_member_energy->Draw("lego2");

    canvas->cd(4);
    h_member_fraction->SetMinimum(0.0);
    h_member_fraction->SetMaximum(1.0);
    h_member_fraction->Draw("lego2");

    canvas->Update();
    canvas->SaveAs("cluster_event.png");

    std::cout << "Event " << event
              << " displayed. selected_cluster=" << selected_cluster
              << " id=" << selected_id
              << " " << z_label << "=" << (use_energy_z ? selected_energy : selected_et)
              << " GeV. Press Enter to advance (q + Enter to quit): " << std::endl;
    std::string line;
    std::getline(std::cin, line);
    if (line == "q" || line == "Q") {
      break;
    }

    gSystem->ProcessEvents();
  }
}
