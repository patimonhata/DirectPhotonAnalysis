#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <TSystem.h>

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "/sphenix/user/ryotaro/DirectPhotonAnalysis/EventDisplay/install/include/pi0eventdisplaydump/Pi0EventDisplayDump.h"

R__LOAD_LIBRARY(/sphenix/user/ryotaro/DirectPhotonAnalysis/EventDisplay/install/lib/libPi0EventDisplayDump.so)

int Fun4All_Pi0EventDisplayDump(
    int processID = 0,
    const int n_events = 100,
    // const std::string cluster_node = "CLUSTERINFO_CEMC",
    const std::string cluster_node = "CLUSTERINFO_CEMC_NO_SPLIT",
    const std::string tower_node = "TOWERINFO_CALIB_CEMC"
    // const std::string tower_node = "TOWERINFO_CALIB_CEMC"
)
{
  std::ostringstream pid;
  pid << std::setw(6) << std::setfill('0') << processID;
  const std::string pid_str = pid.str();

  const std::string input_file = Form("/sphenix/user/ryotaro/DirectPhotonAnalysis/SinglePi0GunSimulation/output/DST_pi0_5GeV_eta0_towerinforyotaro/DST_single_pi0_reconstructedInfo_%s.root", pid_str.c_str());
  // const std::string output_file = Form("/sphenix/u/ryotaro/DirectPhotonAnalysis/EventDisplay/output/root/event_display_%s.root", pid_str.c_str());
  const std::string output_file = Form("/sphenix/user/ryotaro/DirectPhotonAnalysis/EventDisplay/output/root/event_display_%s.root", pid_str.c_str());

  Fun4AllServer* se = Fun4AllServer::instance();
  se->Verbosity(1);

  Fun4AllInputManager* input = new Fun4AllDstInputManager("DSTIN");
  input->Verbosity(1);
  input->AddFile(input_file);
  se->registerInputManager(input);

  Pi0EventDisplayDump* display = new Pi0EventDisplayDump("Pi0EventDisplayDump");
  display->set_output_file(output_file);
  display->set_truth_node("G4TruthInfo");
  display->set_cluster_node(cluster_node);
  display->set_tower_node(tower_node);
  display->set_tower_geom_node("TOWERGEOM_CEMC");
  display->set_cemc_hit_node("G4HIT_CEMC");
  display->set_cluster_energy_min(0.0);
  display->set_tower_energy_min(0.0);
  display->set_write_hits(true);
  display->set_verbosity(1);
  se->registerSubsystem(display);

  std::cout << "Fun4All_Pi0EventDisplayDump - input: " << input_file << std::endl;
  std::cout << "Fun4All_Pi0EventDisplayDump - output: " << output_file << std::endl;

  se->run(n_events);
  se->End();

  delete se;
  std::cout << "Fun4All_Pi0EventDisplayDump - all done" << std::endl;
  return EXIT_SUCCESS;
}
