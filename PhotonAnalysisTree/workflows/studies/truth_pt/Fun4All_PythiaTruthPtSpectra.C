#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>

#include <TSystem.h>

#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "/sphenix/user/ryotaro/DirectPhotonAnalysis/PhotonAnalysisTree/install/include/PythiaTruthPtSpectrum.h"

R__LOAD_LIBRARY(/sphenix/user/ryotaro/DirectPhotonAnalysis/PhotonAnalysisTree/install/lib64/libPhotonAnalysisTree.so)

namespace
{
bool valid_suffix(const std::string& suffix)
{
  return suffix.size() > 5 && suffix.compare(suffix.size() - 5, 5, ".root") == 0 &&
      suffix.find('/') == std::string::npos && suffix.find('\\') == std::string::npos;
}
}

int Fun4All_PythiaTruthPtSpectra(
    const std::string manifest_path,
    const long long manifest_begin,
    const long long manifest_end,
    const std::string output_file,
    const int n_bins = 100,
    const double pt_max = 20.0,
    const double max_abs_eta = 0.7,
    const bool use_event_weight = false)
{
  if (manifest_path.empty() || output_file.empty() || manifest_begin < 0 || manifest_end <= manifest_begin ||
      n_bins <= 0 || !std::isfinite(pt_max) || pt_max <= 0.0 || !std::isfinite(max_abs_eta))
  {
    std::cerr << "Fun4All_PythiaTruthPtSpectra - invalid argument" << std::endl;
    return EXIT_FAILURE;
  }

  std::ifstream manifest(manifest_path);
  if (!manifest)
  {
    std::cerr << "Fun4All_PythiaTruthPtSpectra - cannot open manifest: " << manifest_path << std::endl;
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
        std::cerr << "Fun4All_PythiaTruthPtSpectra - invalid or duplicate suffix at row " << row << ": " << line << std::endl;
        return EXIT_FAILURE;
      }
      suffixes.push_back(line);
    }
    ++row;
  }
  if (static_cast<long long>(suffixes.size()) != manifest_end - manifest_begin)
  {
    std::cerr << "Fun4All_PythiaTruthPtSpectra - incomplete manifest range" << std::endl;
    return EXIT_FAILURE;
  }

  constexpr const char* input_file_prefix = "G4Hits_";
  Fun4AllServer* server = Fun4AllServer::instance();
  server->Verbosity(0);
  auto* input_manager = new Fun4AllDstInputManager("G4HITS");
  for (const std::string& suffix : suffixes)
  {
    input_manager->AddFile(std::string(input_file_prefix) + suffix);
  }
  server->registerInputManager(input_manager);

  auto* accumulator = new PythiaTruthPtSpectrum("PythiaTruthPtSpectrum");
  accumulator->set_output_file_name(output_file);
  accumulator->set_manifest_path(manifest_path);
  accumulator->set_manifest_range(manifest_begin, manifest_end);
  accumulator->set_suffix_range(suffixes.front(), suffixes.back());
  accumulator->set_input_file_prefix(input_file_prefix);
  accumulator->set_signal_embedding_id(1);
  accumulator->set_hepmc_event_map_node_name("PHHepMCGenEventMap");
  accumulator->set_truth_node_name("G4TruthInfo");
  accumulator->set_binning(n_bins, pt_max);
  accumulator->set_max_abs_eta(max_abs_eta);
  accumulator->set_use_event_weight(use_event_weight);
  accumulator->set_verbosity(1);
  server->registerSubsystem(accumulator);

  std::cout << "Fun4All_PythiaTruthPtSpectra - range/files/output = [" << manifest_begin << ":" << manifest_end << "]/"
            << suffixes.size() << "/" << output_file << std::endl;
  const int run_status = server->run(0);
  const int end_status = server->End();
  delete server;
  const bool run_ok = run_status == Fun4AllReturnCodes::EVENT_OK || run_status == Fun4AllReturnCodes::ABORTPROCESSING ||
      run_status == Fun4AllReturnCodes::ABORTEVENT;
  if (!run_ok || end_status != Fun4AllReturnCodes::EVENT_OK)
  {
    std::cerr << "Fun4All_PythiaTruthPtSpectra - failed (run=" << run_status << ", End=" << end_status << ")" << std::endl;
    gSystem->Exit(EXIT_FAILURE);
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
