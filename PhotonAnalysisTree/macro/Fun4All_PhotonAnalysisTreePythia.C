#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>

#include <TSystem.h>

#include <cstdlib>
#include <exception>
#include <iostream>
#include <limits>
#include <string>

#include "/sphenix/user/ryotaro/DirectPhotonAnalysis/PhotonAnalysisTree/install/include/PythiaPhotonAnalysisTree.h"

R__LOAD_LIBRARY(/sphenix/user/ryotaro/DirectPhotonAnalysis/Pi0Reconstruction/install/lib/libPi0Reconstruction.so)
R__LOAD_LIBRARY(/sphenix/user/ryotaro/DirectPhotonAnalysis/PhotonAnalysisTree/install/lib64/libPhotonAnalysisTree.so)

int Fun4All_PhotonAnalysisTreePythia(
    const std::string input_suffix = "pythia8_Jet5-0000000028-000000.root",
    const int n_events = 0,
    const std::string output_directory = "/sphenix/user/ryotaro/DirectPhotonAnalysis/PhotonAnalysisTree/output/root")
{
  constexpr const char* root_extension = ".root";
  const std::string extension(root_extension);
  if (input_suffix.empty() ||
      input_suffix.find('/') != std::string::npos ||
      input_suffix.size() <= extension.size() ||
      input_suffix.compare(input_suffix.size() - extension.size(), extension.size(), extension) != 0 ||
      n_events < 0 || output_directory.empty()) {
    std::cerr << "Fun4All_PhotonAnalysisTreePythia - invalid argument" << std::endl;
    return EXIT_FAILURE;
  }

  const std::size_t segment_separator = input_suffix.rfind('-');
  const std::size_t extension_position = input_suffix.size() - extension.size();
  if (segment_separator == std::string::npos || segment_separator + 1 >= extension_position) {
    std::cerr << "Fun4All_PhotonAnalysisTreePythia - suffix has no numeric segment: " << input_suffix << std::endl;
    return EXIT_FAILURE;
  }
  const std::string segment_text = input_suffix.substr(segment_separator + 1, extension_position - segment_separator - 1);
  std::size_t parsed_characters = 0;
  unsigned long segment_id = 0;
  try {
    segment_id = std::stoul(segment_text, &parsed_characters);
  }
  catch (const std::exception&) {
    parsed_characters = 0;
  }
  if (parsed_characters != segment_text.size() || segment_id > std::numeric_limits<unsigned int>::max()) {
    std::cerr << "Fun4All_PhotonAnalysisTreePythia - invalid numeric segment in suffix: " << input_suffix << std::endl;
    return EXIT_FAILURE;
  }

  const std::string calo_file = "DST_CALO_CLUSTER_" + input_suffix;
  const std::string mbd_epd_file = "DST_MBD_EPD_" + input_suffix;
  const std::string truth_jet_file = "DST_TRUTH_JET_" + input_suffix;
  const std::string g4hits_file = "G4Hits_" + input_suffix;
  const std::string input_summary = calo_file + ";" + mbd_epd_file + ";" + truth_jet_file + ";" + g4hits_file;
  const std::string output_tag = input_suffix.substr(0, extension_position);
  const std::string output_file = output_directory + "/pythia_photon_analysis_tree_" + output_tag + ".root";

  Fun4AllServer* server = Fun4AllServer::instance();
  server->Verbosity(1);

  /* no-name function to register inputs */ 
  const auto register_input = [server](const std::string& manager_name, const std::string& file_name) {
    auto* input_manager = new Fun4AllDstInputManager(manager_name);
    input_manager->AddFile(file_name);
    server->registerInputManager(input_manager);
  };
  register_input("DST_CALO_CLUSTER", calo_file);
  register_input("DST_MBD_EPD", mbd_epd_file);
  register_input("DST_TRUTH_JET", truth_jet_file);
  register_input("G4HITS", g4hits_file);

  auto* tree_maker = new PythiaPhotonAnalysisTree("PythiaPhotonAnalysisTree");
  tree_maker->set_input_file_name(input_summary);
  tree_maker->set_output_file_name(output_file);
  tree_maker->set_source_file_id(static_cast<unsigned int>(segment_id));
  tree_maker->set_signal_embedding_id(1);
  tree_maker->set_truth_node_name("G4TruthInfo");
  tree_maker->set_hepmc_event_map_node_name("PHHepMCGenEventMap");
  tree_maker->set_tower_node_name("TOWERINFO_CALIB_CEMC");
  tree_maker->set_raw_truth_tower_node_name("TOWER_CALIB_CEMC");
  tree_maker->set_tower_geom_node_name("TOWERGEOM_CEMC");
  tree_maker->set_split_cluster_node_name("CLUSTERINFO_CEMC");
  tree_maker->set_nosplit_cluster_node_name("CLUSTERINFO_CEMC_NO_SPLIT");
  tree_maker->set_require_nosplit_cluster_node(false);
  tree_maker->set_min_cluster_energy(0.0);
  tree_maker->set_shower_shape_min_tower_energy(0.070);
  tree_maker->set_store_shower_shape_tower_patch(true);
  tree_maker->set_verbosity(1);
  server->registerSubsystem(tree_maker);

  std::cout << "Fun4All_PhotonAnalysisTreePythia - inputs:\n"
            << "  " << calo_file << '\n'
            << "  " << mbd_epd_file << '\n'
            << "  " << truth_jet_file << '\n'
            << "  " << g4hits_file << '\n'
            << "Fun4All_PhotonAnalysisTreePythia - output: " << output_file << std::endl;
  const int run_status = server->run(n_events);
  const int end_status = server->End();
  delete server;

  const bool run_ok = run_status == Fun4AllReturnCodes::EVENT_OK || (n_events == 0 && run_status == Fun4AllReturnCodes::ABORTEVENT);
  if (!run_ok || end_status != Fun4AllReturnCodes::EVENT_OK) {
    std::cerr << "Fun4All_PhotonAnalysisTreePythia - failed (run=" << run_status << ", End=" << end_status << ")" << std::endl;
    gSystem->Exit(EXIT_FAILURE);
    return EXIT_FAILURE;
  }
  if (run_status != Fun4AllReturnCodes::EVENT_OK) {
    std::cout << "Fun4All_PhotonAnalysisTreePythia - completed at input EOF" << std::endl;
  }
  return EXIT_SUCCESS;
}
