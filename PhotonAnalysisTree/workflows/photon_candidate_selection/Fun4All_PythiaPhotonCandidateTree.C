#include <caloreco/RawClusterBuilderTopo.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <g4calo/RawTowerBuilder.h>
#include <g4detectors/PHG4FullProjSpacalCellReco.h>

#include <TSystem.h>

#include <cstdlib>
#include <exception>
#include <iostream>
#include <limits>
#include <string>

#include "/sphenix/user/ryotaro/DirectPhotonAnalysis/PhotonAnalysisTree/install/include/PythiaPhotonCandidateTree.h"

R__LOAD_LIBRARY(/sphenix/user/ryotaro/DirectPhotonAnalysis/Pi0Reconstruction/install/lib/libPi0Reconstruction.so)
R__LOAD_LIBRARY(/sphenix/user/ryotaro/DirectPhotonAnalysis/PhotonAnalysisTree/install/lib64/libPhotonAnalysisTree.so)
R__LOAD_LIBRARY(libcalo_reco.so)
R__LOAD_LIBRARY(libg4detectors.so)
R__LOAD_LIBRARY(libg4calo.so)

int Fun4All_PythiaPhotonCandidateTree(
    const std::string input_suffix = "pythia8_Jet5-0000000028-000000.root",
    const std::string sample_name = "jet5",
    const int n_events = 0,
    const std::string output_directory = "/sphenix/user/ryotaro/DirectPhotonAnalysis/PhotonAnalysisTree/output/photon_candidate_selection",
    const std::string model_file = "/sphenix/user/ryotaro/DirectPhotonAnalysis/PhotonAnalysisTree/models/model_ppg15_nominal_base_v3E_split_3to35_allplus3jet40_ppg12split_single_tmva.root")
{
  constexpr const char* extension = ".root";
  if (input_suffix.empty() || input_suffix.find('/') != std::string::npos ||
      input_suffix.size() <= 5U || input_suffix.compare(input_suffix.size() - 5U, 5U, extension) != 0 ||
      sample_name.empty() || n_events < 0 || output_directory.empty() || model_file.empty())
  {
    std::cerr << "Fun4All_PythiaPhotonCandidateTree - invalid argument" << std::endl;
    return EXIT_FAILURE;
  }

  const std::size_t segment_separator = input_suffix.rfind('-');
  const std::size_t extension_position = input_suffix.size() - 5U;
  if (segment_separator == std::string::npos || segment_separator + 1U >= extension_position)
  {
    std::cerr << "Fun4All_PythiaPhotonCandidateTree - suffix has no numeric segment: " << input_suffix << std::endl;
    return EXIT_FAILURE;
  }
  const std::string segment_text = input_suffix.substr(segment_separator + 1U, extension_position - segment_separator - 1U);
  std::size_t parsed = 0U;
  unsigned long segment_id = 0UL;
  try
  {
    segment_id = std::stoul(segment_text, &parsed);
  }
  catch (const std::exception&)
  {
    parsed = 0U;
  }
  if (parsed != segment_text.size() || segment_id > std::numeric_limits<unsigned int>::max())
  {
    std::cerr << "Fun4All_PythiaPhotonCandidateTree - invalid segment: " << input_suffix << std::endl;
    return EXIT_FAILURE;
  }

  const std::string calo_file = "DST_CALO_CLUSTER_" + input_suffix;
  const std::string mbd_file = "DST_MBD_EPD_" + input_suffix;
  const std::string jet_file = "DST_TRUTH_JET_" + input_suffix;
  const std::string hits_file = "G4Hits_" + input_suffix;
  const std::string output_tag = input_suffix.substr(0U, extension_position);
  const std::string output_file = output_directory + "/pythia_photon_candidate_tree_" + output_tag + ".root";

  Fun4AllServer* server = Fun4AllServer::instance();
  server->Verbosity(0);
  const auto register_input = [server](const std::string& name, const std::string& file) {
    auto* manager = new Fun4AllDstInputManager(name);
    manager->AddFile(file);
    server->registerInputManager(manager);
  };
  register_input("DST_CALO_CLUSTER", calo_file);
  register_input("DST_MBD_EPD", mbd_file);
  register_input("DST_TRUTH_JET", jet_file);
  register_input("G4HITS", hits_file);

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
    std::cerr << "Fun4All_PythiaPhotonCandidateTree - CALIBRATIONROOT is not set" << std::endl;
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
  producer->set_input_file_name(calo_file + ";" + mbd_file + ";" + jet_file + ";" + hits_file);
  producer->set_output_file_name(output_file);
  producer->set_model_file_name(model_file);
  producer->set_sample_name(sample_name);
  producer->set_source_file_id(static_cast<unsigned int>(segment_id));
  producer->set_signal_embedding_id(1);
  producer->set_truth_jet_node_name("AntiKt_Truth_r04");
  producer->set_verbosity(1);
  server->registerSubsystem(producer);

  std::cout << "Fun4All_PythiaPhotonCandidateTree - sample/input/output = "
            << sample_name << "/" << input_suffix << "/" << output_file << std::endl;
  const int run_status = server->run(n_events);
  const int end_status = server->End();
  delete server;

  const bool run_ok = run_status == Fun4AllReturnCodes::EVENT_OK ||
      (n_events == 0 && run_status == Fun4AllReturnCodes::ABORTPROCESSING);
  if (!run_ok || end_status != Fun4AllReturnCodes::EVENT_OK)
  {
    std::cerr << "Fun4All_PythiaPhotonCandidateTree - failed (run=" << run_status << ", End=" << end_status << ")" << std::endl;
    gSystem->Exit(EXIT_FAILURE);
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
