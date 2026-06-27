#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

int MakePi0HistogramsFromTree(
    const double min_cluster_energy = 0.07,
    const double mass_window_min = 0.10,
    const double mass_window_max = 0.18,
    const int processID=0)
{
  std::ostringstream pid;
  pid << std::setw(6) << std::setfill('0') << processID;
  std::string pid_str = pid.str();

  // const std::string input_file = Form("/sphenix/u/ryotaro/DirectPhotonAnalysis/Pi0Reconstruction/output/pi0_reconstruction_%s.root", pid_str.c_str());
  const std::string input_file = "/sphenix/user/ryotaro/DirectPhotonAnalysis/Pi0Reconstruction/output/100kevents_5GeV_pi0_reconstruction_SPLIT_CLUSTERS_tree.root";
  const std::string output_file = Form("/sphenix/user/ryotaro/DirectPhotonAnalysis/Pi0Reconstruction/output/5GeV_pi0_over70MeV_80-160MeV_rehist_%s.root", pid_str.c_str());  

  TFile *input = TFile::Open(input_file.c_str(), "READ");
  if (!input || input->IsZombie())
  {
    std::cout << "Failed to open input file: " << input_file << std::endl;
    return 1;
  }

  TTree *tree = dynamic_cast<TTree *>(input->Get("event_tree"));
  if (!tree)
  {
    std::cout << "Failed to find event_tree in: " << input_file << std::endl;
    input->Close();
    return 1;
  }

  unsigned int event = 0;
  unsigned int ncluster = 0;
  unsigned int ncluster_all = 0;
  std::vector<double> *cluster_e = nullptr;
  std::vector<double> *cluster_px = nullptr;
  std::vector<double> *cluster_py = nullptr;
  std::vector<double> *cluster_pz = nullptr;

  tree->SetBranchAddress("event", &event);
  tree->SetBranchAddress("ncluster", &ncluster);
  tree->SetBranchAddress("ncluster_all", &ncluster_all);
  tree->SetBranchAddress("cluster_e", &cluster_e);
  tree->SetBranchAddress("cluster_px", &cluster_px);
  tree->SetBranchAddress("cluster_py", &cluster_py);
  tree->SetBranchAddress("cluster_pz", &cluster_pz);

  TFile *output = TFile::Open(output_file.c_str(), "RECREATE");
  if (!output || output->IsZombie())
  {
    std::cout << "Failed to open output file: " << output_file << std::endl;
    input->Close();
    return 1;
  }

  TH1D *h_ncluster = new TH1D("h_ncluster", "CEMC clusters per event;N_{cluster};Events", 100, 0.0, 100.0);
  TH1D *h_cluster_e = new TH1D("h_cluster_e", "CEMC cluster energy;E_{cluster} [GeV];Clusters", 200, 0.0, 20.0);
  TH1D *h_m_gg = new TH1D("h_m_gg", "CEMC cluster pair invariant mass;M_{#gamma#gamma} [GeV];Pairs", 100, 0.0, 1.0);
  TH1D *h_pair_e_asym = new TH1D("h_pair_e_asym", "CEMC cluster pair energy asymmetry after mass window;(|E_{1}-E_{2}|)/(E_{1}+E_{2});Pairs", 100, -1.0, 1.0);

  const Long64_t entries = tree->GetEntries();
  for (Long64_t entry = 0; entry < entries; ++entry)
  {
    if (entry%1000 == 0) {
      std:cout << "event: " << entry << std::endl;
    }

    tree->GetEntry(entry);
    if (!cluster_e || !cluster_px || !cluster_py || !cluster_pz)
    {
      continue;
    }

    std::vector<unsigned int> selected_clusters;
    selected_clusters.reserve(cluster_e->size());
    for (std::size_t i = 0; i < cluster_e->size(); ++i)
    {
      const double energy = cluster_e->at(i);
      if (!std::isfinite(energy) || energy < min_cluster_energy)
      {
        continue;
      }

      h_cluster_e->Fill(energy);
      selected_clusters.push_back(static_cast<unsigned int>(i));
    }

    h_ncluster->Fill(static_cast<double>(selected_clusters.size()));

    for (std::size_t i = 0; i < selected_clusters.size(); ++i)
    {
      for (std::size_t j = i + 1; j < selected_clusters.size(); ++j)
      {
        const unsigned int first_index = selected_clusters[i];
        const unsigned int second_index = selected_clusters[j];
        const double first_energy = cluster_e->at(first_index);
        const double second_energy = cluster_e->at(second_index);
        const double total_energy = first_energy + second_energy;
        const double px = cluster_px->at(first_index) + cluster_px->at(second_index);
        const double py = cluster_py->at(first_index) + cluster_py->at(second_index);
        const double pz = cluster_pz->at(first_index) + cluster_pz->at(second_index);
        const double mass2 = total_energy * total_energy - px * px - py * py - pz * pz;
        const double mass = std::sqrt(std::max(0.0, mass2));

        h_m_gg->Fill(mass);

        if (mass < mass_window_min || mass > mass_window_max || total_energy <= 0.0)
        {
          continue;
        }

        const double energy_asymmetry = std::abs(first_energy - second_energy) / total_energy;
        h_pair_e_asym->Fill(energy_asymmetry);
      }
    }
  }

  output->cd();
  h_ncluster->Write();
  h_cluster_e->Write();
  h_m_gg->Write();
  h_pair_e_asym->Write();
  output->Close();
  input->Close();

  std::cout << "Wrote histograms to " << output_file << std::endl;
  std::cout << "min_cluster_energy = " << min_cluster_energy << std::endl;
  std::cout << "mass window = [" << mass_window_min << ", " << mass_window_max << "] GeV" << std::endl;

  return 0;
}
