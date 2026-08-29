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

R__ADD_INCLUDE_PATH(/sphenix/user/ryotaro/DirectPhotonAnalysis/PhotonAnalysisTree/install/include)

#include "/sphenix/user/ryotaro/DirectPhotonAnalysis/TopologyEventDisplay/install/include/TopologyEventDisplayDump.h"

R__LOAD_LIBRARY(/sphenix/user/ryotaro/DirectPhotonAnalysis/PhotonAnalysisTree/install/lib64/libPhotonAnalysisTree.so)
R__LOAD_LIBRARY(/sphenix/user/ryotaro/DirectPhotonAnalysis/TopologyEventDisplay/install/lib64/libTopologyEventDisplay.so)
R__LOAD_LIBRARY(libg4detectors.so)
R__LOAD_LIBRARY(libg4calo.so)

namespace
{
bool topology_display_valid_suffix(const std::string& suffix)
{
  return suffix.size() > 5U &&
      suffix.compare(suffix.size() - 5U, 5U, ".root") == 0 &&
      suffix.find('/') == std::string::npos &&
      suffix.find('\\') == std::string::npos;
}
}

int RunPythiaTopologyEventDump(
    const std::string manifest_path,
    const long long manifest_begin,
    const long long manifest_end,
    const std::string output_file,
    const int n_events = 0,
    const bool write_detail = true,
    const double anchor_cluster_eta_max = 0.7,
    const double partner_cluster_eta_max = -1.0,
    const double cemc_acceptance_eta_max = 1.1,
    const double min_cluster_energy = 0.2,
    const double anchor_pi0_fraction_min = 0.5,
    const double min_energy_contribution_fraction = 0.0,
    const double min_photon_energy_recovery = 0.5,
    const double min_direct_match_cluster_energy_coverage = 0.5,
    const double missing_diagnostic_max_delta_r = 0.15,
    const int first_event = 0,
    const double pre_cemc_interaction_radius = 90.0)
{
  if (manifest_path.empty() || output_file.empty() || manifest_begin < 0 ||
      manifest_end <= manifest_begin || n_events < 0 || first_event < 0)
  {
    std::cerr << "RunPythiaTopologyEventDump - invalid arguments" << std::endl;
    return EXIT_FAILURE;
  }
  std::ifstream manifest(manifest_path);
  if (!manifest)
  {
    std::cerr << "RunPythiaTopologyEventDump - cannot open " << manifest_path
              << std::endl;
    return EXIT_FAILURE;
  }
  std::vector<std::string> suffixes;
  std::set<std::string> unique;
  std::string line;
  long long row = 0;
  while (std::getline(manifest, line) && row < manifest_end)
  {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (row >= manifest_begin)
    {
      if (!topology_display_valid_suffix(line) || !unique.insert(line).second)
      {
        std::cerr << "RunPythiaTopologyEventDump - invalid/duplicate suffix at "
                  << row << ": " << line << std::endl;
        return EXIT_FAILURE;
      }
      suffixes.push_back(line);
    }
    ++row;
  }
  if (static_cast<long long>(suffixes.size()) != manifest_end - manifest_begin)
  {
    std::cerr << "RunPythiaTopologyEventDump - incomplete manifest range"
              << std::endl;
    return EXIT_FAILURE;
  }

  auto* server = Fun4AllServer::instance();
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
    std::cerr << "RunPythiaTopologyEventDump - CALIBRATIONROOT is not set"
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

  auto* display = new TopologyEventDisplayDump("TopologyEventDisplayDump");
  display->set_output_file_name(output_file);
  display->set_source_label("pythia");
  display->set_manifest_path(manifest_path);
  display->set_manifest_range(manifest_begin, manifest_end);
  display->set_sample_mode(photon_tree::Pi0SampleMode::pythia);
  display->set_signal_embedding_id(1);
  display->set_first_event(first_event);
  display->set_anchor_cluster_eta_max(anchor_cluster_eta_max);
  display->set_partner_cluster_eta_max(partner_cluster_eta_max);
  display->set_cemc_acceptance_eta_max(cemc_acceptance_eta_max);
  display->set_pre_cemc_interaction_radius(pre_cemc_interaction_radius);
  display->set_min_cluster_energy(min_cluster_energy);
  display->set_anchor_pi0_fraction_min(anchor_pi0_fraction_min);
  display->set_min_energy_contribution_fraction(min_energy_contribution_fraction);
  display->set_min_photon_energy_recovery(min_photon_energy_recovery);
  display->set_min_direct_match_cluster_energy_coverage(min_direct_match_cluster_energy_coverage);
  display->set_missing_diagnostic_max_delta_r(missing_diagnostic_max_delta_r);
  display->set_write_detail(write_detail);
  display->set_verbosity(1);
  server->registerSubsystem(display);

  std::cout << "RunPythiaTopologyEventDump - manifest range/output/detail = ["
            << manifest_begin << ":" << manifest_end << "]/" << output_file
            << "/" << write_detail << std::endl;
  const int run_status = server->run(n_events > 0 ? first_event + n_events : 0);
  const int end_status = server->End();
  delete server;
  const bool run_ok = run_status == Fun4AllReturnCodes::EVENT_OK ||
      run_status == Fun4AllReturnCodes::ABORTPROCESSING ||
      run_status == Fun4AllReturnCodes::ABORTEVENT;
  if (!run_ok || end_status != Fun4AllReturnCodes::EVENT_OK)
  {
    std::cerr << "RunPythiaTopologyEventDump - failed run/end = "
              << run_status << "/" << end_status << std::endl;
    gSystem->Exit(EXIT_FAILURE);
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
