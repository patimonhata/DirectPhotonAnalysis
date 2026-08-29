#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <g4calo/RawTowerBuilder.h>
#include <g4detectors/PHG4FullProjSpacalCellReco.h>

#include <TSystem.h>

#include <cstdlib>
#include <iostream>
#include <string>

R__ADD_INCLUDE_PATH(/sphenix/user/ryotaro/DirectPhotonAnalysis/PhotonAnalysisTree/install/include)

#include "/sphenix/user/ryotaro/DirectPhotonAnalysis/TopologyEventDisplay/install/include/TopologyEventDisplayDump.h"

R__LOAD_LIBRARY(/sphenix/user/ryotaro/DirectPhotonAnalysis/PhotonAnalysisTree/install/lib64/libPhotonAnalysisTree.so)
R__LOAD_LIBRARY(/sphenix/user/ryotaro/DirectPhotonAnalysis/TopologyEventDisplay/install/lib64/libTopologyEventDisplay.so)
R__LOAD_LIBRARY(libg4detectors.so)
R__LOAD_LIBRARY(libg4calo.so)

int RunSingleParticleTopologyEventDump(
    const std::string input_file,
    const std::string output_file,
    const int n_events = 0,
    const std::string cluster_node = "CLUSTERINFO_CEMC",
    const bool write_detail = true,
    const bool rebuild_truth_nodes = false,
    const double anchor_cluster_eta_max = 1.1,
    const double partner_cluster_eta_max = -1.0,
    const double cemc_acceptance_eta_max = 1.1,
    const double min_cluster_energy = 0.0,
    const double anchor_pi0_fraction_min = 0.5,
    const double min_energy_contribution_fraction = 0.0,
    const double min_photon_energy_recovery = 0.5,
    const double min_direct_match_cluster_energy_coverage = 0.5,
    const double missing_diagnostic_max_delta_r = 0.15,
    const int first_event = 0,
    const double pre_cemc_interaction_radius = 90.0)
{
  if (input_file.empty() || output_file.empty() || n_events < 0 || first_event < 0)
  {
    std::cerr << "RunSingleParticleTopologyEventDump - invalid arguments"
              << std::endl;
    return EXIT_FAILURE;
  }
  auto* server = Fun4AllServer::instance();
  server->Verbosity(0);
  auto* input = new Fun4AllDstInputManager("DSTIN");
  input->AddFile(input_file);
  server->registerInputManager(input);

  if (rebuild_truth_nodes)
  {
    const char* calibration_root = std::getenv("CALIBRATIONROOT");
    if (!calibration_root)
    {
      std::cerr << "RunSingleParticleTopologyEventDump - CALIBRATIONROOT is not set"
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
  }
  auto* display = new TopologyEventDisplayDump("TopologyEventDisplayDump");
  display->set_output_file_name(output_file);
  display->set_source_label(input_file);
  display->set_sample_mode(photon_tree::Pi0SampleMode::single_particle);
  display->set_cluster_node_name(cluster_node);
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

  std::cout << "RunSingleParticleTopologyEventDump - input/output/detail = "
            << input_file << "/" << output_file << "/" << write_detail
            << std::endl;
  const int run_status = server->run(n_events > 0 ? first_event + n_events : 0);
  const int end_status = server->End();
  delete server;
  const bool run_ok = run_status == Fun4AllReturnCodes::EVENT_OK ||
      run_status == Fun4AllReturnCodes::ABORTPROCESSING ||
      run_status == Fun4AllReturnCodes::ABORTEVENT;
  if (!run_ok || end_status != Fun4AllReturnCodes::EVENT_OK)
  {
    std::cerr << "RunSingleParticleTopologyEventDump - failed run/end = "
              << run_status << "/" << end_status << std::endl;
    gSystem->Exit(EXIT_FAILURE);
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
