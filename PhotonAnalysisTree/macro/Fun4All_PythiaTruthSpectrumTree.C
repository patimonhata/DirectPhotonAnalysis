#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>

#include <TSystem.h>

#include <cstdlib>
#include <exception>
#include <iostream>
#include <limits>
#include <string>

#include "/sphenix/user/ryotaro/DirectPhotonAnalysis/PhotonAnalysisTree/install/include/PythiaTruthSpectrumTree.h"

R__LOAD_LIBRARY(/sphenix/user/ryotaro/DirectPhotonAnalysis/PhotonAnalysisTree/install/lib64/libPhotonAnalysisTree.so)

int Fun4All_PythiaTruthSpectrumTree(
    const std::string input_suffix = "pythia8_Detroit-0000000028-000000.root",
    const int n_events = 0,
    const std::string output_file = "/sphenix/user/ryotaro/DirectPhotonAnalysis/PhotonAnalysisTree/output/truth_root/pythia_truth_spectrum_tree.root")
{
  constexpr const char* root_extension = ".root";
  const std::string extension(root_extension);
  /* safe guard for input variables */
  if (input_suffix.empty() || input_suffix.find('/') != std::string::npos || input_suffix.size() <= extension.size() || input_suffix.compare(input_suffix.size() - extension.size(), extension.size(), extension) != 0 || n_events < 0 || output_file.empty()) {
    std::cerr << "Fun4All_PythiaTruthSpectrumTree - invalid argument" << std::endl;
    return EXIT_FAILURE;
  }

  const std::size_t segment_separator = input_suffix.rfind('-');
  const std::size_t extension_position = input_suffix.size() - extension.size();
  if (segment_separator == std::string::npos || segment_separator + 1U >= extension_position) {
    std::cerr << "Fun4All_PythiaTruthSpectrumTree - suffix has no numeric segment: " << input_suffix << std::endl;
    return EXIT_FAILURE;
  }
  const std::string segment_text = input_suffix.substr(segment_separator + 1U, extension_position - segment_separator - 1U);
  std::size_t parsed_characters = 0U;
  unsigned long segment_id = 0UL;
  try {
    segment_id = std::stoul(segment_text, &parsed_characters);
  } catch (const std::exception&) {
    parsed_characters = 0U;
  }
  if (parsed_characters != segment_text.size() || segment_id > std::numeric_limits<unsigned int>::max()) {
    std::cerr << "Fun4All_PythiaTruthSpectrumTree - invalid numeric segment: " << input_suffix << std::endl;
    return EXIT_FAILURE;
  }

  const std::string g4hits_file = "G4Hits_" + input_suffix;
  Fun4AllServer* server = Fun4AllServer::instance();
  server->Verbosity(0);

  auto* input_manager = new Fun4AllDstInputManager("G4HITS");
  input_manager->AddFile(g4hits_file);
  server->registerInputManager(input_manager);

  auto* tree_maker = new PythiaTruthSpectrumTree("PythiaTruthSpectrumTree");
  tree_maker->set_input_file_name(g4hits_file);
  tree_maker->set_output_file_name(output_file);
  tree_maker->set_source_file_id(static_cast<unsigned int>(segment_id));
  tree_maker->set_signal_embedding_id(1);
  tree_maker->set_hepmc_event_map_node_name("PHHepMCGenEventMap");
  tree_maker->set_truth_node_name("G4TruthInfo");
  tree_maker->set_verbosity(1);
  server->registerSubsystem(tree_maker);

  std::cout << "Fun4All_PythiaTruthSpectrumTree - input: " << g4hits_file << '\n'
            << "Fun4All_PythiaTruthSpectrumTree - output: " << output_file << std::endl;
  const int run_status = server->run(n_events);
  const int end_status = server->End();
  delete server;

  const bool run_ok =
      run_status == Fun4AllReturnCodes::EVENT_OK ||
      (n_events == 0 &&
       (run_status == Fun4AllReturnCodes::ABORTPROCESSING ||
        run_status == Fun4AllReturnCodes::ABORTEVENT));
  if (!run_ok || end_status != Fun4AllReturnCodes::EVENT_OK) {
    std::cerr << "Fun4All_PythiaTruthSpectrumTree - failed (run=" << run_status << ", End=" << end_status << ")" << std::endl;
    gSystem->Exit(EXIT_FAILURE);
    return EXIT_FAILURE;
  }
  if (run_status != Fun4AllReturnCodes::EVENT_OK)
  {
    std::cout << "Fun4All_PythiaTruthSpectrumTree - completed at input EOF" << std::endl;
  }
  return EXIT_SUCCESS;
}
