#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>

#include <TopoClusterHCalTree.h>

#include <TSystem.h>

#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

R__LOAD_LIBRARY(libg4dst.so)
R__LOAD_LIBRARY(libTopoClusterHCalStudy.so)

int Fun4All_TopoClusterHCalStudy(
    const int job_index = 0,
    const int n_events = 0,
    const std::string project_root = "/sphenix/user/ryotaro/DirectPhotonAnalysis",
    const std::string output_file_override = "")
{
  constexpr int files_per_sample = 500;
  if (job_index < 0 || job_index >= 2 * files_per_sample || n_events < 0)
  {
    std::cerr << "Fun4All_TopoClusterHCalStudy - JOB_INDEX must be in [0,999] "
              << "and N_EVENTS must be non-negative" << std::endl;
    return EXIT_FAILURE;
  }

  const bool is_pi0 = job_index >= files_per_sample;
  const int process_id = is_pi0 ? job_index - files_per_sample : job_index;
  const std::string sample = is_pi0 ? "pi0" : "gamma";

  std::ostringstream process_tag_stream;
  process_tag_stream << std::setw(6) << std::setfill('0') << process_id;
  const std::string process_tag = process_tag_stream.str();

  const std::string study_directory = project_root + "/TopoClusterHCalStudy";
  const std::string input_directory =
      project_root + "/SinglePi0GunSimulation/output/DST_" + sample +
      "_25to35GeV_etapm1_vertexpm60";
  const std::string input_file = input_directory + "/DST_single_" + sample +
      "_reconstructedInfo_" + process_tag + ".root";
  const std::string output_file = output_file_override.empty()
      ? study_directory + "/output/root/topocluster_hcal_" + sample + "_" + process_tag + ".root"
      : output_file_override;

  if (!std::filesystem::is_regular_file(input_file))
  {
    std::cerr << "Fun4All_TopoClusterHCalStudy - input does not exist: "
              << input_file << std::endl;
    return EXIT_FAILURE;
  }

  const std::filesystem::path output_path(output_file);
  std::error_code directory_error;
  std::filesystem::create_directories(output_path.parent_path(), directory_error);
  if (directory_error)
  {
    std::cerr << "Fun4All_TopoClusterHCalStudy - could not create output directory: "
              << directory_error.message() << std::endl;
    return EXIT_FAILURE;
  }

  auto *server = Fun4AllServer::instance();
  server->Verbosity(1);

  auto *input_manager = new Fun4AllDstInputManager("DSTIN");
  input_manager->AddFile(input_file);
  server->registerInputManager(input_manager);

  auto *tree_maker = new TopoClusterHCalTree("TopoClusterHCalTree");
  tree_maker->set_input_file_name(input_file);
  tree_maker->set_output_file_name(output_file);
  tree_maker->set_topocluster_node_name("TOPOCLUSTER_ALLCALO");
  tree_maker->set_sample_name(sample);
  tree_maker->set_sample_id(is_pi0 ? 1U : 0U);
  tree_maker->set_job_index(static_cast<unsigned int>(job_index));
  tree_maker->set_process_id(static_cast<unsigned int>(process_id));
  server->registerSubsystem(tree_maker);

  std::cout << "Fun4All_TopoClusterHCalStudy - sample: " << sample << '\n'
            << "Fun4All_TopoClusterHCalStudy - input: " << input_file << '\n'
            << "Fun4All_TopoClusterHCalStudy - output: " << output_file << std::endl;

  const int run_status = server->run(n_events);
  const int end_status = server->End();
  delete server;

  // run(0) consumes the sole input through EOF and reports ABORTEVENT (-1).
  const bool run_ok = run_status == Fun4AllReturnCodes::EVENT_OK ||
      (n_events == 0 && run_status == Fun4AllReturnCodes::ABORTEVENT);
  if (!run_ok || end_status != Fun4AllReturnCodes::EVENT_OK)
  {
    std::cerr << "Fun4All_TopoClusterHCalStudy - failed (run=" << run_status
              << ", End=" << end_status << ')' << std::endl;
    gSystem->Exit(EXIT_FAILURE);
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
