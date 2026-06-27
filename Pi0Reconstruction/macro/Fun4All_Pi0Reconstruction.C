#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/recoConsts.h>

#include <TSystem.h>

#include <cstdlib>
#include <iostream>
#include <string>

#include "/sphenix/user/ryotaro/DirectPhotonAnalysis/Pi0Reconstruction/install/include/pi0reconstruction/Pi0Reconstruction.h"

R__LOAD_LIBRARY(/sphenix/user/ryotaro/DirectPhotonAnalysis/Pi0Reconstruction/install/lib/libPi0Reconstruction.so)

int Fun4All_Pi0Reconstruction(
    int processID=0,
    const int n_events = 0
    // const std::string input_file = "/sphenix/u/ryotaro/DirectPhotonAnalysis/SinglePi0GunSimulation/output/DST_pi0/DST_single_pi0_reconstructedInfo_000000_ClusterBuilder_100events.root",
    // const std::string output_file = "/sphenix/u/ryotaro/DirectPhotonAnalysis/Pi0Reconstruction/output/pi0_reconstruction.root")
)
{
  std::ostringstream pid;
  pid << std::setw(6) << std::setfill('0') << processID;
  std::string pid_str = pid.str();
  // const std::string input_file = Form("/sphenix/user/ryotaro/DirectPhotonAnalysis/SinglePi0GunSimulation/output/DST_pi0_6GeV_eta0/DST_single_pi0_reconstructedInfo_%s.root", pid_str.c_str());
  const std::string input_file = Form("/sphenix/user/ryotaro/DirectPhotonAnalysis/SinglePi0GunSimulation/output/DST_pi0_5GeV_eta0_towerinforyotaro/DST_single_pi0_reconstructedInfo_%s.root", pid_str.c_str());
  const std::string output_file = Form("/sphenix/user/ryotaro/DirectPhotonAnalysis/Pi0Reconstruction/output/root/pi0_reconstruction_%s.root", pid_str.c_str());
  
  

  Fun4AllServer *fun4all_server = Fun4AllServer::instance();
  fun4all_server->Verbosity(1);

  Fun4AllInputManager *input_manager = new Fun4AllDstInputManager("DSTIN");
  input_manager->Verbosity(1);
  input_manager->AddFile(input_file);
  fun4all_server->registerInputManager(input_manager);

  Pi0Reconstruction *pi0_reconstruction = new Pi0Reconstruction("Pi0Reconstruction");
  pi0_reconstruction->set_output_file_name(output_file);
  pi0_reconstruction->set_cluster_node_name("CLUSTERINFO_CEMC"); 
  // pi0_reconstruction->set_cluster_node_name("CLUSTERINFO_CEMC_NO_SPLIT"); 
  pi0_reconstruction->set_vertex_mode(Pi0Reconstruction::VertexMode::Origin);
  pi0_reconstruction->set_abort_on_missing_cluster_node(true);
  pi0_reconstruction->set_min_cluster_energy(0.00);
  pi0_reconstruction->set_mass_histogram_bins(100, 0.0, 1.0);
  fun4all_server->registerSubsystem(pi0_reconstruction);

  std::cout << "Fun4All_Pi0Reconstruction - input: " << input_file << std::endl;
  std::cout << "Fun4All_Pi0Reconstruction - output: " << output_file << std::endl;
  fun4all_server->run(n_events);
  fun4all_server->End();

  delete fun4all_server;
  std::cout << "Fun4All_Pi0Reconstruction - all done" << std::endl;

  return EXIT_SUCCESS;
}
