#include <uspin/SpinDBContentv1.h>
#include <uspin/SpinDBOutput.h>

#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <system_error>
#include <vector>
#include <unistd.h>

#if defined(__CLING__)
R__LOAD_LIBRARY(libuspin.so)
#endif

int select_spin_qa_runs(
    const char* golden_path = "PhotonAnalysisTree/input/physics_data/Full_ppGoldenRunList_Version3.list",
    const char* output_path = "PhotonAnalysisTree/input/physics_data/run2pp_spin_qa.list")
{
  std::ifstream golden(golden_path);
  if (!golden)  {
    std::cerr << "Cannot open golden run list: " << golden_path << '\n';
    return EXIT_FAILURE;
  }

  std::error_code error;
  if (std::filesystem::equivalent(golden_path, output_path, error)) {
    std::cerr << "The golden run list and output must be different files.\n";
    return EXIT_FAILURE;
  }

  std::set<int> runs;
  std::string line;
  std::size_t line_number = 0;
  while (std::getline(golden, line)) {
    ++line_number;
    std::istringstream record(line);
    record >> std::ws;
    if (record.eof()) continue;

    int run = 0;
    char extra = '\0';
    if (!(record >> run) || run <= 0 || (record >> extra) || !runs.insert(run).second) {
      std::cerr << "Invalid or duplicate run number in " << golden_path << ':' << line_number << '\n';
      return EXIT_FAILURE;
    }
  }
  if (!golden.eof() || golden.bad() || runs.empty()) {
    std::cerr << "Cannot read a nonempty golden run list: " << golden_path << '\n';
    return EXIT_FAILURE;
  }

  SpinDBOutput spin_db("phnxrc");
  if (spin_db.StoreDBContent(*runs.begin(), *runs.rbegin()) != 1) {
    std::cerr << "Cannot load default Spin DB rows for runs " << *runs.begin() << " through " << *runs.rbegin() << '\n';
    return EXIT_FAILURE;
  }

  SpinDBContentv1 spin_entry;
  SpinDBContent* content = &spin_entry;
  std::vector<int> selected;
  std::size_t missing = 0;
  std::size_t flagged = 0;
  for (const int run : runs) {
    if (spin_db.GetDBContentStore(content, run) != 1) {
      ++missing;
    }
    else if (content->GetBadRunFlag() == 0) {
      selected.push_back(run);
    }
    else {
      ++flagged;
    }
  }

  std::string temporary_path = std::string(output_path) + ".XXXXXX";
  const int descriptor = mkstemp(temporary_path.data());
  if (descriptor < 0) {
    std::perror(temporary_path.c_str());
    return EXIT_FAILURE;
  }

  std::FILE* output = fdopen(descriptor, "w");
  if (!output) {
    std::perror(temporary_path.c_str());
    close(descriptor);
    std::remove(temporary_path.c_str());
    return EXIT_FAILURE;
  }

  for (const int run : selected) std::fprintf(output, "%d\n", run);
  const bool write_ok = std::ferror(output) == 0;
  const bool close_ok = std::fclose(output) == 0;
  if (!write_ok || !close_ok || std::rename(temporary_path.c_str(), output_path) != 0) {
    std::cerr << "Cannot write Spin QA run list: " << output_path << '\n';
    std::remove(temporary_path.c_str());
    return EXIT_FAILURE;
  }

  std::cout << "Golden runs: " << runs.size() << ", accepted: " << selected.size()
            << ", bad-run flag nonzero: " << flagged << ", missing default row: " << missing << '\n'
            << "Wrote " << output_path << '\n';
  return EXIT_SUCCESS;
}
