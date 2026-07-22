#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>

#include <TSystem.h>

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "/sphenix/user/ryotaro/DirectPhotonAnalysis/PhotonAnalysisTree/install/include/PhotonAnalysisTree.h"

R__LOAD_LIBRARY(/sphenix/user/ryotaro/DirectPhotonAnalysis/Pi0Reconstruction/install/lib/libPi0Reconstruction.so)
R__LOAD_LIBRARY(/sphenix/user/ryotaro/DirectPhotonAnalysis/PhotonAnalysisTree/install/lib64/libPhotonAnalysisTree.so)

int Fun4All_PhotonAnalysisTree(
    const int process_id = 0,
    const int n_events = 0,
    const std::string input_directory =
        "/sphenix/user/ryotaro/Pi0DirectGammaSeparation/SinglePi0GunSimulation/output/newDST_pi0_5to15GeV_etapm1",
    const std::string output_directory =
        "/sphenix/user/ryotaro/DirectPhotonAnalysis/PhotonAnalysisTree/output/root",
    const int expected_primary_pdg = 111)
{
  if (process_id < 0 || n_events < 0 ||
      (expected_primary_pdg != 22 && expected_primary_pdg != 111))
  {
    std::cerr << "Fun4All_PhotonAnalysisTree - invalid argument" << std::endl;
    return EXIT_FAILURE;
  }

  std::ostringstream id;
  id << std::setw(6) << std::setfill('0') << process_id;
  const std::string id_string = id.str();
  const std::string input_file = input_directory +
      "/DST_single_pi0_reconstructedInfo_" + id_string + ".root";
  const std::string output_file = output_directory +
      "/photon_analysis_tree_" + id_string + ".root";

  Fun4AllServer* server = Fun4AllServer::instance();
  server->Verbosity(1);

  auto* input_manager = new Fun4AllDstInputManager("DSTIN");
  input_manager->AddFile(input_file);
  server->registerInputManager(input_manager);

  auto* tree_maker = new PhotonAnalysisTree("PhotonAnalysisTree");
  tree_maker->set_input_file_name(input_file);
  tree_maker->set_output_file_name(output_file);
  tree_maker->set_source_file_id(static_cast<unsigned int>(process_id));
  tree_maker->set_expected_primary_pdg(expected_primary_pdg);
  tree_maker->set_truth_node_name("G4TruthInfo");
  tree_maker->set_tower_node_name("TOWERINFO_CALIB_CEMC");
  tree_maker->set_tower_geom_node_name("TOWERGEOM_CEMC");
  tree_maker->set_split_cluster_node_name("CLUSTERINFO_CEMC");
  tree_maker->set_nosplit_cluster_node_name("CLUSTERINFO_CEMC_NO_SPLIT");
  tree_maker->set_acceptance_eta_max(1.1);
  tree_maker->set_min_cluster_energy(0.0);
  tree_maker->set_shower_shape_min_tower_energy(0.070);
  tree_maker->set_store_shower_shape_tower_patch(true);
  tree_maker->set_verbosity(1);
  server->registerSubsystem(tree_maker);

  std::cout << "Fun4All_PhotonAnalysisTree - input: " << input_file << '\n'
            << "Fun4All_PhotonAnalysisTree - output: " << output_file << std::endl;
  const int run_status = server->run(n_events);
  const int end_status = server->End();
  delete server;

  // Fun4AllServer::run(0) returns -1 after consuming the only DST input to EOF.
  // PhotonAnalysisTree itself never returns ABORTEVENT/SYNC_FAIL (-1), so this
  // value is an expected completion only for the explicit "run to EOF" mode.
  const bool run_ok = run_status == Fun4AllReturnCodes::EVENT_OK ||
      (n_events == 0 && run_status == Fun4AllReturnCodes::ABORTEVENT);
  if (!run_ok || end_status != Fun4AllReturnCodes::EVENT_OK)
  {
    std::cerr << "Fun4All_PhotonAnalysisTree - failed (run=" << run_status
              << ", End=" << end_status << ")" << std::endl;
    gSystem->Exit(EXIT_FAILURE);
    return EXIT_FAILURE;
  }
  if (run_status != Fun4AllReturnCodes::EVENT_OK)
  {
    std::cout << "Fun4All_PhotonAnalysisTree - completed at input EOF" << std::endl;
  }
  return EXIT_SUCCESS;
}
