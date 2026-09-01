#include "check_pythia_photon_candidate_tree.C"

#include <string>

namespace
{
template <class T>
bool map_bind(TTree* tree, const char* name, T* address)
{
  return tree->GetBranch(name) && tree->SetBranchAddress(name, address) >= 0;
}
}

int check_pythia_photon_candidate_map(const char* path)
{
  const int tree_status = check_pythia_photon_candidate_tree(path);
  if (tree_status != 0)
  {
    return tree_status;
  }

  TFile file(path, "READ");
  auto* events = file.Get<TTree>("event_tree");
  auto* metadata = file.Get<TTree>("metadata");
  int schema_version = -1;
  std::string* input_manifest = nullptr;
  std::string* first_input_suffix = nullptr;
  std::string* last_input_suffix = nullptr;
  std::string* sample_name = nullptr;
  long long manifest_begin = -1;
  long long manifest_end = -1;
  long long input_file_count = -1;
  unsigned int map_chunk_id = 0U;
  unsigned long long n_events_written = 0ULL;
  unsigned long long n_clusters_region_a = 0ULL;
  unsigned long long n_clusters_region_b = 0ULL;
  unsigned long long n_clusters_region_c = 0ULL;
  unsigned long long n_clusters_region_d = 0ULL;
  unsigned long long n_clusters_final_photon = 0ULL;

  bool ok = true;
  ok &= map_bind(metadata, "schema_version", &schema_version);
  ok &= map_bind(metadata, "input_manifest", &input_manifest);
  ok &= map_bind(metadata, "manifest_begin", &manifest_begin);
  ok &= map_bind(metadata, "manifest_end", &manifest_end);
  ok &= map_bind(metadata, "input_file_count", &input_file_count);
  ok &= map_bind(metadata, "first_input_suffix", &first_input_suffix);
  ok &= map_bind(metadata, "last_input_suffix", &last_input_suffix);
  ok &= map_bind(metadata, "map_chunk_id", &map_chunk_id);
  ok &= map_bind(metadata, "sample_name", &sample_name);
  ok &= map_bind(metadata, "n_events_written", &n_events_written);
  ok &= map_bind(metadata, "n_clusters_region_a", &n_clusters_region_a);
  ok &= map_bind(metadata, "n_clusters_region_b", &n_clusters_region_b);
  ok &= map_bind(metadata, "n_clusters_region_c", &n_clusters_region_c);
  ok &= map_bind(metadata, "n_clusters_region_d", &n_clusters_region_d);
  ok &= map_bind(metadata, "n_clusters_final_photon", &n_clusters_final_photon);
  if (!ok || metadata->GetEntry(0) <= 0 || schema_version != 3 || !input_manifest || input_manifest->empty() ||
      manifest_begin < 0 || manifest_end <= manifest_begin || input_file_count != manifest_end - manifest_begin ||
      !first_input_suffix || first_input_suffix->empty() || !last_input_suffix || last_input_suffix->empty() ||
      !sample_name || sample_name->empty() || n_events_written != static_cast<unsigned long long>(events->GetEntries()))
  {
    std::cerr << "Invalid photon-candidate map metadata" << std::endl;
    return 4;
  }

  unsigned long long event_uid = 0ULL;
  unsigned int region_a_count = 0U;
  unsigned int region_b_count = 0U;
  unsigned int region_c_count = 0U;
  unsigned int region_d_count = 0U;
  unsigned int final_photon_count = 0U;
  std::vector<unsigned char>* region_a = nullptr;
  std::vector<unsigned char>* region_b = nullptr;
  std::vector<unsigned char>* region_c = nullptr;
  std::vector<unsigned char>* region_d = nullptr;
  std::vector<unsigned char>* final_photon = nullptr;
  ok = true;
  ok &= map_bind(events, "event_uid", &event_uid);
  ok &= map_bind(events, "region_a_count", &region_a_count);
  ok &= map_bind(events, "region_b_count", &region_b_count);
  ok &= map_bind(events, "region_c_count", &region_c_count);
  ok &= map_bind(events, "region_d_count", &region_d_count);
  ok &= map_bind(events, "final_photon_count", &final_photon_count);
  ok &= map_bind(events, "split_cluster_pass_region_a", &region_a);
  ok &= map_bind(events, "split_cluster_pass_region_b", &region_b);
  ok &= map_bind(events, "split_cluster_pass_region_c", &region_c);
  ok &= map_bind(events, "split_cluster_pass_region_d", &region_d);
  ok &= map_bind(events, "split_cluster_pass_final_photon", &final_photon);
  if (!ok)
  {
    std::cerr << "Missing photon-candidate map count branch" << std::endl;
    return 5;
  }

  std::set<unsigned long long> unique_event_uids;
  unsigned long long total_a = 0ULL;
  unsigned long long total_b = 0ULL;
  unsigned long long total_c = 0ULL;
  unsigned long long total_d = 0ULL;
  unsigned long long total_final = 0ULL;
  for (Long64_t entry = 0; entry < events->GetEntries(); ++entry)
  {
    events->GetEntry(entry);
    const auto count_pass = [](const std::vector<unsigned char>& flags) {
      return static_cast<unsigned int>(std::count_if(flags.begin(), flags.end(), [](unsigned char value) { return value != 0U; }));
    };
    if (!unique_event_uids.insert(event_uid).second ||
        region_a_count != count_pass(*region_a) || region_b_count != count_pass(*region_b) ||
        region_c_count != count_pass(*region_c) || region_d_count != count_pass(*region_d) ||
        final_photon_count != count_pass(*final_photon))
    {
      std::cerr << "Invalid event UID or ABCD count at entry " << entry << std::endl;
      return 6;
    }
    total_a += region_a_count;
    total_b += region_b_count;
    total_c += region_c_count;
    total_d += region_d_count;
    total_final += final_photon_count;
  }
  if (total_a != n_clusters_region_a || total_b != n_clusters_region_b || total_c != n_clusters_region_c ||
      total_d != n_clusters_region_d || total_final != n_clusters_final_photon)
  {
    std::cerr << "Map metadata ABCD totals do not match event_tree" << std::endl;
    return 7;
  }

  std::cout << "check_pythia_photon_candidate_map - chunk/range/files/A/B/C/D/final = "
            << map_chunk_id << "/[" << manifest_begin << ":" << manifest_end << "]/" << input_file_count << "/"
            << total_a << "/" << total_b << "/" << total_c << "/" << total_d << "/" << total_final << std::endl;
  return 0;
}
