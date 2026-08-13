#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>

#include <TSystem.h>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "/sphenix/user/ryotaro/DirectPhotonAnalysis/PhotonAnalysisTree/install/include/PythiaPi0ClusterMultiplicity.h"

R__LOAD_LIBRARY(/sphenix/user/ryotaro/DirectPhotonAnalysis/PhotonAnalysisTree/install/lib64/libPhotonAnalysisTree.so)

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

int Fun4All_PythiaPi0ClusterMultiplicity(
    const std::string manifest_path,
    const long long manifest_begin,
    const long long manifest_end,
    const std::string output_file,
    const int pt_bins = 100,
    const double pt_max = 20.0,
    const int multiplicity_max = 20,
    const int cluster_energy_bins = 100,
    const double cluster_energy_max = 20.0,
    const double truth_eta_max = 0.7,
    const double cluster_eta_max = 0.7)
{
  if (manifest_path.empty() || output_file.empty() || manifest_begin < 0 ||
      manifest_end <= manifest_begin)
  {
    std::cerr << "Fun4All_PythiaPi0ClusterMultiplicity - invalid argument"
              << std::endl;
    return EXIT_FAILURE;
  }
  std::ifstream manifest(manifest_path);
  if (!manifest)
  {
    std::cerr << "Fun4All_PythiaPi0ClusterMultiplicity - cannot open manifest: "
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
        std::cerr << "Fun4All_PythiaPi0ClusterMultiplicity - invalid or duplicate suffix at row "
                  << row << ": " << line << std::endl;
        return EXIT_FAILURE;
      }
      suffixes.push_back(line);
    }
    ++row;
  }
  if (static_cast<long long>(suffixes.size()) != manifest_end - manifest_begin)
  {
    std::cerr << "Fun4All_PythiaPi0ClusterMultiplicity - incomplete manifest range"
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

  auto* accumulator =
      new PythiaPi0ClusterMultiplicity("PythiaPi0ClusterMultiplicity");
  accumulator->set_output_file_name(output_file);
  accumulator->set_manifest_path(manifest_path);
  accumulator->set_manifest_range(manifest_begin, manifest_end);
  accumulator->set_suffix_range(suffixes.front(), suffixes.back());
  accumulator->set_signal_embedding_id(1);
  accumulator->set_truth_node_name("G4TruthInfo");
  accumulator->set_hepmc_event_map_node_name("PHHepMCGenEventMap");
  accumulator->set_tower_node_name("TOWERINFO_CALIB_CEMC");
  accumulator->set_raw_truth_tower_node_name("TOWER_CALIB_CEMC");
  accumulator->set_split_cluster_node_name("CLUSTERINFO_CEMC");
  accumulator->set_pt_binning(pt_bins, pt_max);
  accumulator->set_multiplicity_max(multiplicity_max);
  accumulator->set_cluster_energy_binning(
      cluster_energy_bins, cluster_energy_max);
  accumulator->set_truth_eta_max(truth_eta_max);
  accumulator->set_cluster_eta_max(cluster_eta_max);
  accumulator->set_verbosity(1);
  server->registerSubsystem(accumulator);

  std::cout << "Fun4All_PythiaPi0ClusterMultiplicity - range/files/output = ["
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
    std::cerr << "Fun4All_PythiaPi0ClusterMultiplicity - failed (run="
              << run_status << ", End=" << end_status << ")" << std::endl;
    gSystem->Exit(EXIT_FAILURE);
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
