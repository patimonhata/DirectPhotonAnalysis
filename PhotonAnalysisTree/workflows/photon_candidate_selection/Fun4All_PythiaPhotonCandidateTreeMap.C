#include <caloreco/RawClusterBuilderTopo.h>
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

#include "/sphenix/user/ryotaro/DirectPhotonAnalysis/PhotonAnalysisTree/install/include/PythiaPhotonCandidateTree.h"

R__LOAD_LIBRARY(/sphenix/user/ryotaro/DirectPhotonAnalysis/Pi0Reconstruction/install/lib/libPi0Reconstruction.so)
R__LOAD_LIBRARY(/sphenix/user/ryotaro/DirectPhotonAnalysis/PhotonAnalysisTree/install/lib64/libPhotonAnalysisTree.so)
R__LOAD_LIBRARY(libcalo_reco.so)
R__LOAD_LIBRARY(libg4detectors.so)
R__LOAD_LIBRARY(libg4calo.so)

namespace
{
bool valid_suffix(const std::string& suffix)
{
  return suffix.size() > 5U && suffix.compare(suffix.size() - 5U, 5U, ".root") == 0 &&
      suffix.find('/') == std::string::npos && suffix.find('\\') == std::string::npos;
}
}

int Fun4All_PythiaPhotonCandidateTreeMap(
    const std::string manifest_path,
    const long long manifest_begin,
    const long long manifest_end,
    const std::string output_file,
    const std::string sample_name,
    const unsigned int map_chunk_id,
    const int n_events = 0,
    const std::string model_file = "/sphenix/user/ryotaro/DirectPhotonAnalysis/PhotonAnalysisTree/models/model_ppg15_nominal_base_v3E_split_3to35_allplus3jet40_ppg12split_single_tmva.root")
{
  if (manifest_path.empty() || output_file.empty() || sample_name.empty() || model_file.empty() ||
      manifest_begin < 0 || manifest_end <= manifest_begin || n_events < 0)
  {
    std::cerr << "Fun4All_PythiaPhotonCandidateTreeMap - invalid argument" << std::endl;
    return EXIT_FAILURE;
  }

  std::ifstream manifest(manifest_path);
  if (!manifest)
  {
    std::cerr << "Fun4All_PythiaPhotonCandidateTreeMap - cannot open manifest: " << manifest_path << std::endl;
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
        std::cerr << "Fun4All_PythiaPhotonCandidateTreeMap - invalid or duplicate suffix at row " << row << ": " << line << std::endl;
        return EXIT_FAILURE;
      }
      suffixes.push_back(line);
    }
    ++row;
  }
  if (static_cast<long long>(suffixes.size()) != manifest_end - manifest_begin)
  {
    std::cerr << "Fun4All_PythiaPhotonCandidateTreeMap - incomplete manifest range" << std::endl;
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

  auto* em_tc = new RawClusterBuilderTopo("em_tc");
  em_tc->set_nodename("TOPOCLUSTER_ALLCALO");
  em_tc->set_enable_HCal(true);
  em_tc->set_enable_EMCal(true);
  em_tc->set_noise(0.006, 0.03, 0.09);
  em_tc->set_significance(4.0, 2.0, 1.0);
  em_tc->allow_corner_neighbor(true);
  em_tc->set_do_split(true);
  em_tc->set_minE_local_max(1.0, 2.0, 0.5);
  em_tc->set_R_shower(0.025);
  server->registerSubsystem(em_tc);

  const char* calibration_root = std::getenv("CALIBRATIONROOT");
  if (!calibration_root)
  {
    std::cerr << "Fun4All_PythiaPhotonCandidateTreeMap - CALIBRATIONROOT is not set" << std::endl;
    return EXIT_FAILURE;
  }
  auto* cemc_cells = new PHG4FullProjSpacalCellReco("CEMCCYLCELLRECO");
  cemc_cells->Detector("CEMC");
  cemc_cells->Verbosity(0);
  cemc_cells->get_light_collection_model().load_data_file(
      std::string(calibration_root) + "/CEMC/LightCollection/Prototype3Module.xml",
      "data_grid_light_guide_efficiency", "data_grid_fiber_trans");
  server->registerSubsystem(cemc_cells);

  auto* truth_towers = new RawTowerBuilder("CEMCTruthRawTowerBuilder");
  truth_towers->Detector("CEMC");
  truth_towers->set_sim_tower_node_prefix("SIM");
  truth_towers->set_towerinfo(RawTowerBuilder::ProcessTowerType::kRawTowerOnly);
  truth_towers->Verbosity(0);
  server->registerSubsystem(truth_towers);

  auto* producer = new PythiaPhotonCandidateTree("PythiaPhotonCandidateTree");
  producer->set_input_file_name(manifest_path);
  producer->set_output_file_name(output_file);
  producer->set_model_file_name(model_file);
  producer->set_sample_name(sample_name);
  producer->set_primary_input_manager(calo_input);
  producer->set_manifest_path(manifest_path);
  producer->set_manifest_range(manifest_begin, manifest_end);
  producer->set_suffix_range(suffixes.front(), suffixes.back());
  producer->set_map_chunk_id(map_chunk_id);
  producer->set_signal_embedding_id(1);
  producer->set_truth_jet_node_name("AntiKt_Truth_r04");
  producer->set_verbosity(1);
  server->registerSubsystem(producer);

  std::cout << "Fun4All_PythiaPhotonCandidateTreeMap - sample/range/files/output = "
            << sample_name << "/[" << manifest_begin << ":" << manifest_end << "]/" << suffixes.size() << "/" << output_file << std::endl;
  const int run_status = server->run(n_events);
  const int end_status = server->End();
  delete server;

  const bool run_ok = run_status == Fun4AllReturnCodes::EVENT_OK ||
      run_status == Fun4AllReturnCodes::ABORTPROCESSING || run_status == Fun4AllReturnCodes::ABORTEVENT;
  if (!run_ok || end_status != Fun4AllReturnCodes::EVENT_OK)
  {
    std::cerr << "Fun4All_PythiaPhotonCandidateTreeMap - failed (run=" << run_status << ", End=" << end_status << ")" << std::endl;
    gSystem->Exit(EXIT_FAILURE);
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
