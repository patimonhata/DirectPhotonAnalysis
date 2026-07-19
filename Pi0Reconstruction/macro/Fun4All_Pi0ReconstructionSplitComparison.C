#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <TSystem.h>

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "/sphenix/user/ryotaro/DirectPhotonAnalysis/Pi0Reconstruction/install/include/pi0reconstruction/Pi0Reconstruction.h"

R__LOAD_LIBRARY(/sphenix/user/ryotaro/DirectPhotonAnalysis/Pi0Reconstruction/install/lib/libPi0Reconstruction.so)

int Fun4All_Pi0ReconstructionSplitComparison(
    int processID = 0,
    const int n_events = 0,
    const std::string input_directory = "/sphenix/user/ryotaro/DirectPhotonAnalysis/SinglePi0GunSimulation/output/DST_pi0_5GeV_eta0_towerinfo",
    const std::string output_directory = "/sphenix/user/ryotaro/DirectPhotonAnalysis/Pi0Reconstruction/output/root")
{
  if (processID < 0) {
    std::cerr << "Fun4All_Pi0ReconstructionSplitComparison - processID must be non-negative" << std::endl;
    return EXIT_FAILURE;
  }

  std::ostringstream pid;
  pid << std::setw(6) << std::setfill('0') << processID;
  const std::string pid_str = pid.str();

  const std::string input_file = input_directory + "/DST_single_pi0_reconstructedInfo_" + pid_str + ".root";
  const std::string split_output_file = output_directory + "/pi0_reconstruction_SPLIT_" + pid_str + ".root";
  const std::string no_split_output_file = output_directory + "/pi0_reconstruction_NO_SPLIT_" + pid_str + ".root";

  Fun4AllServer *fun4all_server = Fun4AllServer::instance();
  fun4all_server->Verbosity(1);

  Fun4AllInputManager *input_manager = new Fun4AllDstInputManager("DSTIN");
  input_manager->Verbosity(1);
  input_manager->AddFile(input_file);
  fun4all_server->registerInputManager(input_manager);

  Pi0Reconstruction *split_reconstruction = new Pi0Reconstruction("Pi0Reconstruction_SPLIT");
  split_reconstruction->set_output_file_name(split_output_file);
  split_reconstruction->set_cluster_node_name("CLUSTERINFO_CEMC");
  split_reconstruction->set_process_id(static_cast<unsigned int>(processID));
  split_reconstruction->set_vertex_mode(Pi0Reconstruction::VertexMode::Origin);
  split_reconstruction->set_abort_on_missing_cluster_node(true);
  split_reconstruction->set_min_cluster_energy(0.0);
  split_reconstruction->set_mass_histogram_bins(100, 0.0, 1.0);
  fun4all_server->registerSubsystem(split_reconstruction);

  Pi0Reconstruction *no_split_reconstruction = new Pi0Reconstruction("Pi0Reconstruction_NO_SPLIT");
  no_split_reconstruction->set_output_file_name(no_split_output_file);
  no_split_reconstruction->set_cluster_node_name("CLUSTERINFO_CEMC_NO_SPLIT");
  no_split_reconstruction->set_process_id(static_cast<unsigned int>(processID));
  no_split_reconstruction->set_vertex_mode(Pi0Reconstruction::VertexMode::Origin);
  no_split_reconstruction->set_abort_on_missing_cluster_node(true);
  no_split_reconstruction->set_min_cluster_energy(0.0);
  no_split_reconstruction->set_mass_histogram_bins(100, 0.0, 1.0);
  fun4all_server->registerSubsystem(no_split_reconstruction);

  std::cout << "Fun4All_Pi0ReconstructionSplitComparison - input: " << input_file << std::endl;
  std::cout << "Fun4All_Pi0ReconstructionSplitComparison - SPLIT output: " << split_output_file << std::endl;
  std::cout << "Fun4All_Pi0ReconstructionSplitComparison - NO_SPLIT output: " << no_split_output_file << std::endl;
  std::cout << "Fun4All_Pi0ReconstructionSplitComparison - process ID: " << processID << std::endl;

  const int run_status = fun4all_server->run(n_events);
  const int end_status = fun4all_server->End();

  delete fun4all_server;

  if (run_status != 0 || end_status != 0) {
    std::cerr << "Fun4All_Pi0ReconstructionSplitComparison - Fun4All failed" << " (run=" << run_status << ", End=" << end_status << ")" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Fun4All_Pi0ReconstructionSplitComparison - all done" << std::endl;
  return EXIT_SUCCESS;
}
