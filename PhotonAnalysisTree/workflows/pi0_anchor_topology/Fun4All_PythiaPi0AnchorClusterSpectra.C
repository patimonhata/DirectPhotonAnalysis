#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <g4calo/RawTowerBuilder.h>
#include <g4detectors/PHG4FullProjSpacalCellReco.h>

#include <TSystem.h>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "/sphenix/user/ryotaro/DirectPhotonAnalysis/PhotonAnalysisTree/install/include/PythiaPi0AnchorClusterSpectrum.h"

R__LOAD_LIBRARY(/sphenix/user/ryotaro/DirectPhotonAnalysis/PhotonAnalysisTree/install/lib64/libPhotonAnalysisTree.so)
R__LOAD_LIBRARY(libg4detectors.so)
R__LOAD_LIBRARY(libg4calo.so)

namespace
{
bool valid_suffix(const std::string& suffix)
{
  return suffix.size() > 5U &&
      suffix.compare(suffix.size() - 5U, 5U, ".root") == 0 &&
      suffix.find('/') == std::string::npos &&
      suffix.find('\\') == std::string::npos;
}
}

int Fun4All_PythiaPi0AnchorClusterSpectra(
    const std::string manifest_path,
    const long long manifest_begin,
    const long long manifest_end,
    const std::string output_file,
    const int n_bins = 100,
    const double et_max = 20.0,
    const double truth_eta_max = 0.7,
    const double anchor_cluster_eta_max = 0.7,
    const double partner_cluster_eta_max = -1.0,
    const double min_cluster_energy = 0.2,
    const double dominant_fraction_min = 0.5,
    const double anchor_pi0_fraction_min = 0.5,
    const double min_energy_contribution_fraction = 0.0,
    const double min_photon_energy_recovery = 0.5,
    const double max_abs_vertex_z = 60.0)
{
  if (manifest_path.empty() || output_file.empty() || manifest_begin < 0 ||
      manifest_end <= manifest_begin)
  {
    std::cerr << "Fun4All_PythiaPi0AnchorClusterSpectra - invalid argument"
              << std::endl;
    return EXIT_FAILURE;
  }

  std::ifstream manifest(manifest_path);
  if (!manifest)
  {
    std::cerr << "Fun4All_PythiaPi0AnchorClusterSpectra - cannot open manifest: "
              << manifest_path << std::endl;
    return EXIT_FAILURE;
  }
  std::vector<std::string> suffixes;
  std::set<std::string> unique_suffixes;
  std::string line;
  long long row = 0;
  while (std::getline(manifest, line) && row < manifest_end)
  {
    if (!line.empty() && line.back() == '\r')
    {
      line.pop_back();
    }
    if (row >= manifest_begin)
    {
      if (!valid_suffix(line) || !unique_suffixes.insert(line).second)
      {
        std::cerr
            << "Fun4All_PythiaPi0AnchorClusterSpectra - invalid or duplicate suffix at row "
            << row << ": " << line << std::endl;
        return EXIT_FAILURE;
      }
      suffixes.push_back(line);
    }
    ++row;
  }
  if (static_cast<long long>(suffixes.size()) !=
      manifest_end - manifest_begin)
  {
    std::cerr
        << "Fun4All_PythiaPi0AnchorClusterSpectra - incomplete manifest range"
        << std::endl;
    return EXIT_FAILURE;
  }

  Fun4AllServer* server = Fun4AllServer::instance();
  server->Verbosity(0);
  auto* calo_input = new Fun4AllDstInputManager("DST_CALO_CLUSTER");
  auto* mbd_input = new Fun4AllDstInputManager("DST_MBD_EPD");
  auto* truth_jet_input = new Fun4AllDstInputManager("DST_TRUTH_JET");
  auto* g4hits_input = new Fun4AllDstInputManager("G4HITS");
  for (const std::string& suffix : suffixes)
  {
    calo_input->AddFile("DST_CALO_CLUSTER_" + suffix);
    mbd_input->AddFile("DST_MBD_EPD_" + suffix);
    truth_jet_input->AddFile("DST_TRUTH_JET_" + suffix);
    g4hits_input->AddFile("G4Hits_" + suffix);
  }
  server->registerInputManager(calo_input);
  server->registerInputManager(mbd_input);
  server->registerInputManager(truth_jet_input);
  server->registerInputManager(g4hits_input);

  const char* calibration_root = std::getenv("CALIBRATIONROOT");
  if (!calibration_root)
  {
    std::cerr
        << "Fun4All_PythiaPi0AnchorClusterSpectra - CALIBRATIONROOT is not set"
        << std::endl;
    return EXIT_FAILURE;
  }
  auto* cemc_cells = new PHG4FullProjSpacalCellReco("CEMCCYLCELLRECO");
  cemc_cells->Detector("CEMC");
  cemc_cells->Verbosity(0);
  cemc_cells->get_light_collection_model().load_data_file(
      std::string(calibration_root) +
          "/CEMC/LightCollection/Prototype3Module.xml",
      "data_grid_light_guide_efficiency", "data_grid_fiber_trans");
  server->registerSubsystem(cemc_cells);

  auto* truth_towers = new RawTowerBuilder("CEMCTruthRawTowerBuilder");
  truth_towers->Detector("CEMC");
  truth_towers->set_sim_tower_node_prefix("SIM");
  truth_towers->set_towerinfo(RawTowerBuilder::ProcessTowerType::kRawTowerOnly);
  truth_towers->Verbosity(0);
  server->registerSubsystem(truth_towers);

  auto* accumulator = new PythiaPi0AnchorClusterSpectrum("PythiaPi0AnchorClusterSpectrum");
  accumulator->set_output_file_name(output_file);
  accumulator->set_manifest_path(manifest_path);
  accumulator->set_manifest_range(manifest_begin, manifest_end);
  accumulator->set_suffix_range(suffixes.front(), suffixes.back());
  accumulator->set_signal_embedding_id(1);
  accumulator->set_truth_node_name("G4TruthInfo");
  accumulator->set_hepmc_event_map_node_name("PHHepMCGenEventMap");
  accumulator->set_tower_node_name("TOWERINFO_CALIB_CEMC");
  accumulator->set_raw_truth_tower_node_name("TOWER_SIM_CEMC");
  accumulator->set_truth_cell_node_name("G4CELL_CEMC");
  accumulator->set_truth_hit_node_name("G4HIT_CEMC");
  accumulator->set_split_cluster_node_name("CLUSTERINFO_CEMC");
  accumulator->set_binning(n_bins, et_max);
  accumulator->set_truth_eta_max(truth_eta_max);
  accumulator->set_anchor_cluster_eta_max(anchor_cluster_eta_max);
  accumulator->set_partner_cluster_eta_max(partner_cluster_eta_max);
  accumulator->set_min_cluster_energy(min_cluster_energy);
  accumulator->set_dominant_fraction_min(dominant_fraction_min);
  accumulator->set_anchor_pi0_fraction_min(anchor_pi0_fraction_min);
  accumulator->set_min_energy_contribution_fraction(min_energy_contribution_fraction);
  accumulator->set_min_photon_energy_recovery(min_photon_energy_recovery);
  accumulator->set_max_abs_vertex_z(max_abs_vertex_z);
  accumulator->set_verbosity(1);
  server->registerSubsystem(accumulator);

  std::cout
      << "Fun4All_PythiaPi0AnchorClusterSpectra - range/files/output = ["
      << manifest_begin << ":" << manifest_end << "]/" << suffixes.size()
      << "/" << output_file << std::endl;
  const int run_status = server->run(0);
  const int end_status = server->End();
  delete server;
  const bool run_ok = run_status == Fun4AllReturnCodes::EVENT_OK ||
      run_status == Fun4AllReturnCodes::ABORTPROCESSING ||
      run_status == Fun4AllReturnCodes::ABORTEVENT;
  if (!run_ok || end_status != Fun4AllReturnCodes::EVENT_OK)
  {
    std::cerr << "Fun4All_PythiaPi0AnchorClusterSpectra - failed (run="
              << run_status << ", End=" << end_status << ")" << std::endl;
    gSystem->Exit(EXIT_FAILURE);
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
