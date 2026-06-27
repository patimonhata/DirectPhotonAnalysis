#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllUtils.h>

#include <string>
#include <vector>

#include "/sphenix/user/ryotaro/DirectPhotonAnalysis/EmcalEtViewer/install/include/cemctowerdumper/CemcClusterDumper.h"

R__LOAD_LIBRARY(libCemcClusterDumper.so)

void Fun4All_CemcClusterDumper(int process_id, int run, int n_events, bool save_tree=false, const std::string& cluster_node_name="CLUSTERINFO_CEMC"){
  Fun4AllServer* fun4all_server = Fun4AllServer::instance();

  std::vector<Fun4AllInputManager*> input_managers;

  std::string job_index = std::to_string(process_id);
  int job_index_len = 5;
  job_index.insert(0, job_index_len - job_index.size(), '0');

  std::string fileName = "/sphenix/user/ryotaro/DirectPhotonAnalysis/SinglePi0GunSimulation/output/DST_pi0_5GeV_eta0/DST_single_pi0_reconstructedInfo_000000.root";
  input_managers.push_back(new Fun4AllDstInputManager(Form("DST_cluster")));
  input_managers.back()->Verbosity(2);
  input_managers.back()->AddFile(fileName);
  fun4all_server->registerInputManager(input_managers.back());

  CemcClusterDumper* module_cemc_cluster_dumper = new CemcClusterDumper("CemcClusterDumper", run, job_index, save_tree);
  module_cemc_cluster_dumper->set_cluster_node_name(cluster_node_name);

  std::string output_directory = std::string("/sphenix/user/ryotaro/CemcClusterDumper/output/") + std::to_string(run);
  std::string final_output_file_name = Form("%d-%s.root", run, job_index.c_str());
  system(Form("mkdir -p %s/completed", output_directory.c_str()));
  system(Form("if [ -f %s/completed/%s ]; then rm %s/completed/%s; fi;", output_directory.c_str(), final_output_file_name.c_str(), output_directory.c_str(), final_output_file_name.c_str()));

  fun4all_server->registerSubsystem(module_cemc_cluster_dumper);
  fun4all_server->run(n_events);
  fun4all_server->End();

  system(Form("mv %s/%s %s/completed", output_directory.c_str(), final_output_file_name.c_str(), output_directory.c_str()));

  delete fun4all_server;
}
