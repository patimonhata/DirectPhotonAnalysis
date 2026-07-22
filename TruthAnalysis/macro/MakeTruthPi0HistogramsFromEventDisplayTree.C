#include <TFile.h>
#include <TH1D.h>
#include <TParameter.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace TruthPi0Hist
{
  static const double kInvalidDouble = -999.0;

  struct Gamma
  {
    int track_id = -999;
    int parent_id = -999;
    double e = kInvalidDouble;
    double px = kInvalidDouble;
    double py = kInvalidDouble;
    double pz = kInvalidDouble;
    double eta1 = kInvalidDouble;
    double phi1 = kInvalidDouble;
    bool has_projection = false;
  };

  struct Pi0
  {
    int event = -999;
    int track_id = -999;
    std::vector<Gamma> gammas;
  };

  struct Projection
  {
    double eta1 = kInvalidDouble;
    double phi1 = kInvalidDouble;
  };

  bool is_valid(double value)
  {
    return std::isfinite(value) && value != kInvalidDouble;
  }

  bool in_acceptance(const Gamma& gamma, double eta_max)
  {
    return gamma.has_projection && is_valid(gamma.eta1) && std::abs(gamma.eta1) < eta_max;
  }

  void create_output_directory(const std::string& filename)
  {
    const std::string::size_type slash_position = filename.find_last_of('/');
    if (slash_position == std::string::npos)
    {
      return;
    }

    const std::string directory = filename.substr(0, slash_position);
    if (!directory.empty())
    {
      gSystem->mkdir(directory.c_str(), true);
    }
  }

  double invariant_mass(const Gamma& first, const Gamma& second)
  {
    const double total_e = first.e + second.e;
    const double px = first.px + second.px;
    const double py = first.py + second.py;
    const double pz = first.pz + second.pz;
    const double mass2 = total_e * total_e - px * px - py * py - pz * pz;
    return std::sqrt(std::max(0.0, mass2));
  }

  double energy_asymmetry(const Gamma& first, const Gamma& second)
  {
    const double total_e = first.e + second.e;
    return total_e > 0.0 ? std::abs(first.e - second.e) / total_e : kInvalidDouble;
  }
}

int MakeTruthPi0HistogramsFromEventDisplayTree(
    const int processID=0,
    const double acceptance_eta_max = 1.1,
    const Long64_t progress_interval = 100000)
{
  std::ostringstream prid;
  prid << std::setw(6) << std::setfill('0') << processID;
  std::string pid_str = prid.str();

  const std::string input_file = Form("/sphenix/user/ryotaro/DirectPhotonAnalysis/EventDisplay/output/root/event_display_%s.root", pid_str.c_str());
  const std::string output_file = Form("/sphenix/user/ryotaro/DirectPhotonAnalysis/TruthAnalysis/output/root/100events_5GeVpi0_eta0_truth_histograms_%s_NOSPLIT_CLUSTERS.root", pid_str.c_str());


  using Pi0Key = std::pair<int, int>;
  using TruthPi0Hist::Gamma;
  using TruthPi0Hist::Pi0;
  using TruthPi0Hist::Projection;

  TFile* input = TFile::Open(input_file.c_str(), "READ");
  if (!input || input->IsZombie())
  {
    std::cout << "MakeTruthPi0HistogramsFromEventDisplayTree - failed to open input file: " << input_file << std::endl;
    return 1;
  }

  TTree* truth_particles = dynamic_cast<TTree*>(input->Get("truth_particles"));
  TTree* truth_segments = dynamic_cast<TTree*>(input->Get("truth_segments"));
  if (!truth_particles)
  {
    std::cout << "MakeTruthPi0HistogramsFromEventDisplayTree - failed to find truth_particles in: " << input_file << std::endl;
    input->Close();
    return 1;
  }
  if (!truth_segments)
  {
    std::cout << "MakeTruthPi0HistogramsFromEventDisplayTree - failed to find truth_segments in: " << input_file << std::endl;
    input->Close();
    return 1;
  }

  int event = 0;
  int track_id = 0;
  int pid = 0;
  int parent_id = 0;
  int is_primary = 0;
  double e = 0.0;
  double px = 0.0;
  double py = 0.0;
  double pz = 0.0;

  truth_particles->SetBranchStatus("*", 0);
  truth_particles->SetBranchStatus("event", 1);
  truth_particles->SetBranchStatus("track_id", 1);
  truth_particles->SetBranchStatus("pid", 1);
  truth_particles->SetBranchStatus("parent_id", 1);
  truth_particles->SetBranchStatus("is_primary", 1);
  truth_particles->SetBranchStatus("e", 1);
  truth_particles->SetBranchStatus("px", 1);
  truth_particles->SetBranchStatus("py", 1);
  truth_particles->SetBranchStatus("pz", 1);
  truth_particles->SetBranchAddress("event", &event);
  truth_particles->SetBranchAddress("track_id", &track_id);
  truth_particles->SetBranchAddress("pid", &pid);
  truth_particles->SetBranchAddress("parent_id", &parent_id);
  truth_particles->SetBranchAddress("is_primary", &is_primary);
  truth_particles->SetBranchAddress("e", &e);
  truth_particles->SetBranchAddress("px", &px);
  truth_particles->SetBranchAddress("py", &py);
  truth_particles->SetBranchAddress("pz", &pz);

  std::map<Pi0Key, Pi0> pi0_by_key;
  std::map<Pi0Key, std::vector<Gamma> > pending_gamma_by_parent;

  const Long64_t particle_entries = truth_particles->GetEntries();
  for (Long64_t entry = 0; entry < particle_entries; ++entry)
  {
    if (progress_interval > 0 && entry > 0 && entry % progress_interval == 0)
    {
      std::cout << "MakeTruthPi0HistogramsFromEventDisplayTree - read truth_particles " << entry << " / " << particle_entries << std::endl;
    }

    truth_particles->GetEntry(entry);
    if (pid == 111 && is_primary == 1)
    {
      const Pi0Key key(event, track_id);
      Pi0& pi0 = pi0_by_key[key];
      pi0.event = event;
      pi0.track_id = track_id;

      std::map<Pi0Key, std::vector<Gamma> >::iterator pending_iter = pending_gamma_by_parent.find(key);
      if (pending_iter != pending_gamma_by_parent.end())
      {
        pi0.gammas.insert(pi0.gammas.end(), pending_iter->second.begin(), pending_iter->second.end());
        pending_gamma_by_parent.erase(pending_iter);
      }
      continue;
    }

    if (pid == 22)
    {
      Gamma gamma;
      gamma.track_id = track_id;
      gamma.parent_id = parent_id;
      gamma.e = e;
      gamma.px = px;
      gamma.py = py;
      gamma.pz = pz;

      const Pi0Key parent_key(event, parent_id);
      std::map<Pi0Key, Pi0>::iterator pi0_iter = pi0_by_key.find(parent_key);
      if (pi0_iter != pi0_by_key.end())
      {
        pi0_iter->second.gammas.push_back(gamma);
      }
      else
      {
        pending_gamma_by_parent[parent_key].push_back(gamma);
      }
    }
  }

  int segment_event = 0;
  int segment_track_id = 0;
  int segment_pid = 0;
  int segment_type = 0;
  double eta1 = 0.0;
  double phi1 = 0.0;

  truth_segments->SetBranchStatus("*", 0);
  truth_segments->SetBranchStatus("event", 1);
  truth_segments->SetBranchStatus("track_id", 1);
  truth_segments->SetBranchStatus("pid", 1);
  truth_segments->SetBranchStatus("eta1", 1);
  truth_segments->SetBranchStatus("phi1", 1);
  truth_segments->SetBranchStatus("segment_type", 1);
  truth_segments->SetBranchAddress("event", &segment_event);
  truth_segments->SetBranchAddress("track_id", &segment_track_id);
  truth_segments->SetBranchAddress("pid", &segment_pid);
  truth_segments->SetBranchAddress("eta1", &eta1);
  truth_segments->SetBranchAddress("phi1", &phi1);
  truth_segments->SetBranchAddress("segment_type", &segment_type);

  std::map<Pi0Key, Projection> projection_by_gamma_key;
  const Long64_t segment_entries = truth_segments->GetEntries();
  for (Long64_t entry = 0; entry < segment_entries; ++entry)
  {
    if (progress_interval > 0 && entry > 0 && entry % progress_interval == 0)
    {
      std::cout << "MakeTruthPi0HistogramsFromEventDisplayTree - read truth_segments " << entry << " / " << segment_entries << std::endl;
    }

    truth_segments->GetEntry(entry);
    if (segment_pid != 22 || segment_type != 2)
    {
      continue;
    }

    Projection projection;
    projection.eta1 = eta1;
    projection.phi1 = phi1;
    projection_by_gamma_key[Pi0Key(segment_event, segment_track_id)] = projection;
  }

  for (std::map<Pi0Key, Pi0>::iterator pi0_iter = pi0_by_key.begin(); pi0_iter != pi0_by_key.end(); ++pi0_iter)
  {
    std::vector<Gamma>& gammas = pi0_iter->second.gammas;
    std::sort(gammas.begin(), gammas.end(), [](const Gamma& first, const Gamma& second) {
      return first.track_id < second.track_id;
    });

    for (Gamma& gamma : gammas)
    {
      const Pi0Key gamma_key(pi0_iter->second.event, gamma.track_id);
      std::map<Pi0Key, Projection>::const_iterator projection_iter = projection_by_gamma_key.find(gamma_key);
      if (projection_iter == projection_by_gamma_key.end())
      {
        continue;
      }

      gamma.eta1 = projection_iter->second.eta1;
      gamma.phi1 = projection_iter->second.phi1;
      gamma.has_projection = TruthPi0Hist::is_valid(gamma.eta1);
    }
  }

  TruthPi0Hist::create_output_directory(output_file);
  TFile* output = TFile::Open(output_file.c_str(), "RECREATE");
  if (!output || output->IsZombie())
  {
    std::cout << "MakeTruthPi0HistogramsFromEventDisplayTree - failed to open output file: " << output_file << std::endl;
    input->Close();
    return 1;
  }

  TH1D* h_truth_decay_ngamma = new TH1D("h_truth_decay_ngamma", "Direct truth #gamma daughters per primary #pi^{0};N_{#gamma};Primary #pi^{0}", 10, -0.5, 9.5);
  TH1D* h_truth_gamma_e = new TH1D("h_truth_gamma_e", "Truth #gamma energy from primary #pi^{0} #rightarrow 2#gamma;E_{#gamma}^{truth} [GeV];Photons", 200, 0.0, 20.0);
  TH1D* h_truth_m_gg = new TH1D("h_truth_m_gg", "Truth #gamma#gamma invariant mass from primary #pi^{0} #rightarrow 2#gamma;M_{#gamma#gamma}^{truth} [GeV];Pairs", 100, 0.0, 1.0);
  TH1D* h_truth_pair_e_asym = new TH1D("h_truth_pair_e_asym", "Truth #gamma#gamma energy asymmetry from primary #pi^{0} #rightarrow 2#gamma;(|E_{1}-E_{2}|)/(E_{1}+E_{2});Pairs", 100, -1.0, 1.0);
  TH1D* h_truth_gamma_e_in_acceptance = new TH1D("h_truth_gamma_e_in_acceptance", "Truth #gamma energy, both #gamma in CEMC acceptance;E_{#gamma}^{truth} [GeV];Photons", 200, 0.0, 20.0);
  TH1D* h_truth_m_gg_in_acceptance = new TH1D("h_truth_m_gg_in_acceptance", "Truth #gamma#gamma invariant mass, both #gamma in CEMC acceptance;M_{#gamma#gamma}^{truth} [GeV];Pairs", 100, 0.0, 1.0);
  TH1D* h_truth_pair_e_asym_in_acceptance = new TH1D("h_truth_pair_e_asym_in_acceptance", "Truth #gamma#gamma energy asymmetry, both #gamma in CEMC acceptance;(|E_{1}-E_{2}|)/(E_{1}+E_{2});Pairs", 100, -1.0, 1.0);
  TH1D* h_truth_summary_counts = new TH1D("h_truth_summary_counts", "Truth primary #pi^{0} summary;Category;Count", 7, 0.5, 7.5);
  TH1D* h_truth_summary_fractions = new TH1D("h_truth_summary_fractions", "Truth primary #pi^{0} fractions;Category;Fraction", 4, 0.5, 4.5);

  h_truth_summary_counts->GetXaxis()->SetBinLabel(1, "primary_pi0");
  h_truth_summary_counts->GetXaxis()->SetBinLabel(2, "pi0_to_2gamma");
  h_truth_summary_counts->GetXaxis()->SetBinLabel(3, "pi0_not_2gamma");
  h_truth_summary_counts->GetXaxis()->SetBinLabel(4, "both_gamma_in_acceptance");
  h_truth_summary_counts->GetXaxis()->SetBinLabel(5, "at_least_one_gamma_out_acceptance");
  h_truth_summary_counts->GetXaxis()->SetBinLabel(6, "missing_gamma_projection");
  h_truth_summary_counts->GetXaxis()->SetBinLabel(7, "direct_gamma_from_primary_pi0");
  h_truth_summary_fractions->GetXaxis()->SetBinLabel(1, "pi0_to_2gamma_over_primary_pi0");
  h_truth_summary_fractions->GetXaxis()->SetBinLabel(2, "both_in_acceptance_over_2gamma");
  h_truth_summary_fractions->GetXaxis()->SetBinLabel(3, "at_least_one_out_over_2gamma");
  h_truth_summary_fractions->GetXaxis()->SetBinLabel(4, "missing_projection_over_2gamma");

  TTree* truth_pi0_tree = new TTree("truth_pi0_tree", "Truth primary pi0 decay summary");
  int tree_event = 0;
  int tree_pi0_track_id = 0;
  int tree_n_direct_gamma = 0;
  int tree_is_pi0_to_2gamma = 0;
  int tree_gamma1_track_id = -999;
  int tree_gamma2_track_id = -999;
  int tree_gamma1_in_acceptance = 0;
  int tree_gamma2_in_acceptance = 0;
  int tree_both_gamma_in_acceptance = 0;
  int tree_at_least_one_gamma_out_acceptance = 0;
  int tree_missing_gamma_projection = 0;
  double tree_gamma1_e = TruthPi0Hist::kInvalidDouble;
  double tree_gamma2_e = TruthPi0Hist::kInvalidDouble;
  double tree_gamma1_px = TruthPi0Hist::kInvalidDouble;
  double tree_gamma1_py = TruthPi0Hist::kInvalidDouble;
  double tree_gamma1_pz = TruthPi0Hist::kInvalidDouble;
  double tree_gamma2_px = TruthPi0Hist::kInvalidDouble;
  double tree_gamma2_py = TruthPi0Hist::kInvalidDouble;
  double tree_gamma2_pz = TruthPi0Hist::kInvalidDouble;
  double tree_gamma1_eta1 = TruthPi0Hist::kInvalidDouble;
  double tree_gamma1_phi1 = TruthPi0Hist::kInvalidDouble;
  double tree_gamma2_eta1 = TruthPi0Hist::kInvalidDouble;
  double tree_gamma2_phi1 = TruthPi0Hist::kInvalidDouble;
  double tree_m_gg = TruthPi0Hist::kInvalidDouble;
  double tree_pair_e_asym = TruthPi0Hist::kInvalidDouble;

  truth_pi0_tree->Branch("event", &tree_event);
  truth_pi0_tree->Branch("pi0_track_id", &tree_pi0_track_id);
  truth_pi0_tree->Branch("n_direct_gamma", &tree_n_direct_gamma);
  truth_pi0_tree->Branch("is_pi0_to_2gamma", &tree_is_pi0_to_2gamma);
  truth_pi0_tree->Branch("gamma1_track_id", &tree_gamma1_track_id);
  truth_pi0_tree->Branch("gamma2_track_id", &tree_gamma2_track_id);
  truth_pi0_tree->Branch("gamma1_e", &tree_gamma1_e);
  truth_pi0_tree->Branch("gamma2_e", &tree_gamma2_e);
  truth_pi0_tree->Branch("gamma1_px", &tree_gamma1_px);
  truth_pi0_tree->Branch("gamma1_py", &tree_gamma1_py);
  truth_pi0_tree->Branch("gamma1_pz", &tree_gamma1_pz);
  truth_pi0_tree->Branch("gamma2_px", &tree_gamma2_px);
  truth_pi0_tree->Branch("gamma2_py", &tree_gamma2_py);
  truth_pi0_tree->Branch("gamma2_pz", &tree_gamma2_pz);
  truth_pi0_tree->Branch("gamma1_eta1", &tree_gamma1_eta1);
  truth_pi0_tree->Branch("gamma1_phi1", &tree_gamma1_phi1);
  truth_pi0_tree->Branch("gamma2_eta1", &tree_gamma2_eta1);
  truth_pi0_tree->Branch("gamma2_phi1", &tree_gamma2_phi1);
  truth_pi0_tree->Branch("gamma1_in_acceptance", &tree_gamma1_in_acceptance);
  truth_pi0_tree->Branch("gamma2_in_acceptance", &tree_gamma2_in_acceptance);
  truth_pi0_tree->Branch("both_gamma_in_acceptance", &tree_both_gamma_in_acceptance);
  truth_pi0_tree->Branch("at_least_one_gamma_out_acceptance", &tree_at_least_one_gamma_out_acceptance);
  truth_pi0_tree->Branch("missing_gamma_projection", &tree_missing_gamma_projection);
  truth_pi0_tree->Branch("m_gg", &tree_m_gg);
  truth_pi0_tree->Branch("pair_e_asym", &tree_pair_e_asym);

  Long64_t n_primary_pi0 = 0;
  Long64_t n_pi0_to_2gamma = 0;
  Long64_t n_pi0_not_2gamma = 0;
  Long64_t n_both_gamma_in_acceptance = 0;
  Long64_t n_at_least_one_gamma_out_acceptance = 0;
  Long64_t n_missing_gamma_projection = 0;
  Long64_t n_direct_gamma_from_primary_pi0 = 0;

  for (std::map<Pi0Key, Pi0>::const_iterator pi0_iter = pi0_by_key.begin(); pi0_iter != pi0_by_key.end(); ++pi0_iter)
  {
    const Pi0& pi0 = pi0_iter->second;
    const std::vector<Gamma>& gammas = pi0.gammas;
    const int n_direct_gamma = static_cast<int>(gammas.size());

    ++n_primary_pi0;
    n_direct_gamma_from_primary_pi0 += n_direct_gamma;
    h_truth_decay_ngamma->Fill(n_direct_gamma);

    tree_event = pi0.event;
    tree_pi0_track_id = pi0.track_id;
    tree_n_direct_gamma = n_direct_gamma;
    tree_is_pi0_to_2gamma = n_direct_gamma == 2 ? 1 : 0;
    tree_gamma1_track_id = -999;
    tree_gamma2_track_id = -999;
    tree_gamma1_in_acceptance = 0;
    tree_gamma2_in_acceptance = 0;
    tree_both_gamma_in_acceptance = 0;
    tree_at_least_one_gamma_out_acceptance = 0;
    tree_missing_gamma_projection = 0;
    tree_gamma1_e = TruthPi0Hist::kInvalidDouble;
    tree_gamma2_e = TruthPi0Hist::kInvalidDouble;
    tree_gamma1_px = TruthPi0Hist::kInvalidDouble;
    tree_gamma1_py = TruthPi0Hist::kInvalidDouble;
    tree_gamma1_pz = TruthPi0Hist::kInvalidDouble;
    tree_gamma2_px = TruthPi0Hist::kInvalidDouble;
    tree_gamma2_py = TruthPi0Hist::kInvalidDouble;
    tree_gamma2_pz = TruthPi0Hist::kInvalidDouble;
    tree_gamma1_eta1 = TruthPi0Hist::kInvalidDouble;
    tree_gamma1_phi1 = TruthPi0Hist::kInvalidDouble;
    tree_gamma2_eta1 = TruthPi0Hist::kInvalidDouble;
    tree_gamma2_phi1 = TruthPi0Hist::kInvalidDouble;
    tree_m_gg = TruthPi0Hist::kInvalidDouble;
    tree_pair_e_asym = TruthPi0Hist::kInvalidDouble;

    if (n_direct_gamma != 2)
    {
      ++n_pi0_not_2gamma;
      truth_pi0_tree->Fill();
      continue;
    }

    ++n_pi0_to_2gamma;

    const Gamma& first = gammas[0];
    const Gamma& second = gammas[1];
    const bool first_in_acceptance = TruthPi0Hist::in_acceptance(first, acceptance_eta_max);
    const bool second_in_acceptance = TruthPi0Hist::in_acceptance(second, acceptance_eta_max);
    const bool both_in_acceptance = first_in_acceptance && second_in_acceptance;
    const bool missing_projection = !first.has_projection || !second.has_projection;

    tree_gamma1_track_id = first.track_id;
    tree_gamma2_track_id = second.track_id;
    tree_gamma1_e = first.e;
    tree_gamma2_e = second.e;
    tree_gamma1_px = first.px;
    tree_gamma1_py = first.py;
    tree_gamma1_pz = first.pz;
    tree_gamma2_px = second.px;
    tree_gamma2_py = second.py;
    tree_gamma2_pz = second.pz;
    tree_gamma1_eta1 = first.eta1;
    tree_gamma1_phi1 = first.phi1;
    tree_gamma2_eta1 = second.eta1;
    tree_gamma2_phi1 = second.phi1;
    tree_gamma1_in_acceptance = first_in_acceptance ? 1 : 0;
    tree_gamma2_in_acceptance = second_in_acceptance ? 1 : 0;
    tree_both_gamma_in_acceptance = both_in_acceptance ? 1 : 0;
    tree_at_least_one_gamma_out_acceptance = both_in_acceptance ? 0 : 1;
    tree_missing_gamma_projection = missing_projection ? 1 : 0;
    tree_m_gg = TruthPi0Hist::invariant_mass(first, second);
    tree_pair_e_asym = TruthPi0Hist::energy_asymmetry(first, second);

    h_truth_gamma_e->Fill(first.e);
    h_truth_gamma_e->Fill(second.e);
    h_truth_m_gg->Fill(tree_m_gg);
    if (TruthPi0Hist::is_valid(tree_pair_e_asym))
    {
      h_truth_pair_e_asym->Fill(tree_pair_e_asym);
    }

    if (both_in_acceptance)
    {
      ++n_both_gamma_in_acceptance;
      h_truth_gamma_e_in_acceptance->Fill(first.e);
      h_truth_gamma_e_in_acceptance->Fill(second.e);
      h_truth_m_gg_in_acceptance->Fill(tree_m_gg);
      if (TruthPi0Hist::is_valid(tree_pair_e_asym))
      {
        h_truth_pair_e_asym_in_acceptance->Fill(tree_pair_e_asym);
      }
    }
    else
    {
      ++n_at_least_one_gamma_out_acceptance;
    }

    if (missing_projection)
    {
      ++n_missing_gamma_projection;
    }

    truth_pi0_tree->Fill();
  }

  const double frac_pi0_to_2gamma = n_primary_pi0 > 0 ? static_cast<double>(n_pi0_to_2gamma) / static_cast<double>(n_primary_pi0) : 0.0;
  const double frac_both_gamma_in_acceptance = n_pi0_to_2gamma > 0 ? static_cast<double>(n_both_gamma_in_acceptance) / static_cast<double>(n_pi0_to_2gamma) : 0.0;
  const double frac_at_least_one_gamma_out_acceptance = n_pi0_to_2gamma > 0 ? static_cast<double>(n_at_least_one_gamma_out_acceptance) / static_cast<double>(n_pi0_to_2gamma) : 0.0;
  const double frac_missing_gamma_projection = n_pi0_to_2gamma > 0 ? static_cast<double>(n_missing_gamma_projection) / static_cast<double>(n_pi0_to_2gamma) : 0.0;

  h_truth_summary_counts->SetBinContent(1, n_primary_pi0);
  h_truth_summary_counts->SetBinContent(2, n_pi0_to_2gamma);
  h_truth_summary_counts->SetBinContent(3, n_pi0_not_2gamma);
  h_truth_summary_counts->SetBinContent(4, n_both_gamma_in_acceptance);
  h_truth_summary_counts->SetBinContent(5, n_at_least_one_gamma_out_acceptance);
  h_truth_summary_counts->SetBinContent(6, n_missing_gamma_projection);
  h_truth_summary_counts->SetBinContent(7, n_direct_gamma_from_primary_pi0);
  h_truth_summary_fractions->SetBinContent(1, frac_pi0_to_2gamma);
  h_truth_summary_fractions->SetBinContent(2, frac_both_gamma_in_acceptance);
  h_truth_summary_fractions->SetBinContent(3, frac_at_least_one_gamma_out_acceptance);
  h_truth_summary_fractions->SetBinContent(4, frac_missing_gamma_projection);

  output->cd();
  h_truth_decay_ngamma->Write();
  h_truth_gamma_e->Write();
  h_truth_m_gg->Write();
  h_truth_pair_e_asym->Write();
  h_truth_gamma_e_in_acceptance->Write();
  h_truth_m_gg_in_acceptance->Write();
  h_truth_pair_e_asym_in_acceptance->Write();
  h_truth_summary_counts->Write();
  h_truth_summary_fractions->Write();
  truth_pi0_tree->Write();

  TParameter<double>("acceptance_eta_max", acceptance_eta_max).Write();
  TParameter<double>("n_primary_pi0", static_cast<double>(n_primary_pi0)).Write();
  TParameter<double>("n_pi0_to_2gamma", static_cast<double>(n_pi0_to_2gamma)).Write();
  TParameter<double>("n_pi0_not_2gamma", static_cast<double>(n_pi0_not_2gamma)).Write();
  TParameter<double>("n_both_gamma_in_acceptance", static_cast<double>(n_both_gamma_in_acceptance)).Write();
  TParameter<double>("n_at_least_one_gamma_out_acceptance", static_cast<double>(n_at_least_one_gamma_out_acceptance)).Write();
  TParameter<double>("n_missing_gamma_projection", static_cast<double>(n_missing_gamma_projection)).Write();
  TParameter<double>("frac_pi0_to_2gamma", frac_pi0_to_2gamma).Write();
  TParameter<double>("frac_both_gamma_in_acceptance", frac_both_gamma_in_acceptance).Write();
  TParameter<double>("frac_at_least_one_gamma_out_acceptance", frac_at_least_one_gamma_out_acceptance).Write();
  TParameter<double>("frac_missing_gamma_projection", frac_missing_gamma_projection).Write();

  output->Close();
  input->Close();

  std::cout << "MakeTruthPi0HistogramsFromEventDisplayTree - wrote " << output_file << std::endl;
  std::cout << "  n_primary_pi0 = " << n_primary_pi0 << std::endl;
  std::cout << "  n_pi0_to_2gamma = " << n_pi0_to_2gamma << " (" << frac_pi0_to_2gamma << ")" << std::endl;
  std::cout << "  n_both_gamma_in_acceptance = " << n_both_gamma_in_acceptance << " (" << frac_both_gamma_in_acceptance << " of pi0_to_2gamma)" << std::endl;
  std::cout << "  n_at_least_one_gamma_out_acceptance = " << n_at_least_one_gamma_out_acceptance << " (" << frac_at_least_one_gamma_out_acceptance << " of pi0_to_2gamma)" << std::endl;

  return 0;
}
