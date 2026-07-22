#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <TSystem.h>

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "/sphenix/user/ryotaro/DirectPhotonAnalysis/Pi0Reconstruction/install/include/pi0reconstruction/TowerClusterEnergyAudit.h"

R__LOAD_LIBRARY(/sphenix/user/ryotaro/DirectPhotonAnalysis/Pi0Reconstruction/install/lib/libPi0Reconstruction.so)

int Fun4All_TowerClusterEnergyAudit(
    int processID = 0,
    const int n_events = 0,
    const std::string tower_node = "TOWERINFO_CALIBryotaro_CEMC",
    const std::string split_cluster_node = "CLUSTERINFO_CEMC",
    const std::string no_split_cluster_node = "CLUSTERINFO_CEMC_NO_SPLIT",
    const std::string input_directory = "/sphenix/user/ryotaro/DirectPhotonAnalysis/SinglePi0GunSimulation/output/DST_pi0_5GeV_eta0_towerinforyotaro",
    const std::string output_tag = "")
{
  std::ostringstream pid;
  pid << std::setw(6) << std::setfill('0') << processID;
  const std::string pid_str = pid.str();

  const std::string input_file = input_directory + "/DST_single_pi0_reconstructedInfo_" + pid_str + ".root";
  const std::string output_suffix = output_tag.empty() ? pid_str : output_tag + "_" + pid_str;
  const std::string output_file = "/sphenix/user/ryotaro/DirectPhotonAnalysis/Pi0Reconstruction/output/tower_cluster_energy_audit_" + output_suffix + ".root";

  Fun4AllServer *fun4all_server = Fun4AllServer::instance();
  fun4all_server->Verbosity(1);

  Fun4AllInputManager *input_manager = new Fun4AllDstInputManager("DSTIN");
  input_manager->Verbosity(1);
  input_manager->AddFile(input_file);
  fun4all_server->registerInputManager(input_manager);

  TowerClusterEnergyAudit *audit = new TowerClusterEnergyAudit("TowerClusterEnergyAudit");
  audit->set_output_file_name(output_file);
  audit->set_tower_node_name(tower_node);
  audit->set_split_cluster_node_name(split_cluster_node);
  audit->set_no_split_cluster_node_name(no_split_cluster_node);
  audit->set_tower_energy_threshold(0.070);
  audit->set_cluster_energy_threshold(0.070);
  audit->set_abort_on_missing_nodes(true);
  fun4all_server->registerSubsystem(audit);

  std::cout << "Fun4All_TowerClusterEnergyAudit - input: " << input_file << std::endl;
  std::cout << "Fun4All_TowerClusterEnergyAudit - output: " << output_file << std::endl;
  std::cout << "Fun4All_TowerClusterEnergyAudit - tower node: " << tower_node << std::endl;
  std::cout << "Fun4All_TowerClusterEnergyAudit - split cluster node: " << split_cluster_node << std::endl;
  std::cout << "Fun4All_TowerClusterEnergyAudit - no-split cluster node: " << no_split_cluster_node << std::endl;

  fun4all_server->run(n_events);
  fun4all_server->End();

  delete fun4all_server;
  std::cout << "Fun4All_TowerClusterEnergyAudit - all done" << std::endl;

  return EXIT_SUCCESS;
}
