#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllUtils.h>
#include <Calo_Calib.C>

#include "/sphenix/user/ryotaro/DisplayLegoPlot/install/include/displaylegoplot/DisplayLegoPlot.h"

R__LOAD_LIBRARY(libDisplayLegoPlot.so)

void Fun4All_DisplayLegoPlot(int process_id, int run, int n_events, bool save_tree=false){
  Fun4AllServer* fun4all_server = Fun4AllServer::instance();

  std::vector<Fun4AllInputManager*> input_managers; 

  std::string job_index = std::to_string( process_id );
  int job_index_len = 5;
  job_index.insert(0, job_index_len - job_index.size(), '0');

  std::string fileName = "DST_CALO_CLUSTER_pythia8_Jet30-0000000021-000000.root"; 
  input_managers.push_back(new Fun4AllDstInputManager(Form("DST_track")));
  input_managers.back() -> Verbosity(2);
  input_managers.back() -> AddFile( fileName );
  fun4all_server->registerInputManager(input_managers.back());
  
  pair<int, int> runseg = Fun4AllUtils::GetRunSegment(fileName);
  int runnumber = runseg.first;
  recoConsts *rc = recoConsts::instance();
  rc->set_StringFlag("CDB_GLOBALTAG", "MDC2"); /* From wiki: The latest version/tags of the CDB that is being updated continuously are "ProdA_2024" (for 2024 data), "ProdA_2023" (for 2023 data), and "MDC2" (for MC). */
  rc->set_uint64Flag("TIMESTAMP", runnumber);

  Process_Calo_Calib(); /* From wiki: Note you want to put this after you declare the input managers but before your analysis module. */

  DisplayLegoPlot* module_display_lego_plot = new DisplayLegoPlot("name", run, job_index, save_tree);

  // std::string output_directory = Form("output/%d", run);
  std::string output_directory = std::string("/sphenix/user/ryotaro/DisplayLegoPlot/output/") + std::to_string(run);;
  std::string final_output_file_name = Form("%d-%s.root", run, job_index.c_str());
  system(Form("if [ -f %s/completed/%s ]; then rm %s/completed/%s; fi;", output_directory.c_str(), final_output_file_name.c_str(), output_directory.c_str(), final_output_file_name.c_str()));  

  fun4all_server->registerSubsystem(module_display_lego_plot);
  fun4all_server->run(n_events);
  fun4all_server->End();
  
  system(Form("mv %s/%s %s/completed", output_directory.c_str(), final_output_file_name.c_str(), output_directory.c_str()));

  // delete module_carryover_monitor;
  // for (int i=0; i<8; i++) {
  //   delete input_managers[i]; /* This will cause segmentaion fault. Probably because input_managers will be deleted automatically when we delete fun4all_server.
  // }
  delete fun4all_server;


}
