#ifndef PHOTON_CANDIDATE_SETTINGS_H
#define PHOTON_CANDIDATE_SETTINGS_H
#include <TFile.h>
#include <TTree.h>
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace photon_candidate_settings
{
// Stable order for diagnostic output metadata; map metadata retains named scalars.
inline bool read(const std::string& path, std::vector<double>& values)
{
  TFile file(path.c_str(), "READ");
  auto* tree = file.Get<TTree>("metadata");
  if (file.IsZombie() || !tree || tree->GetEntries() != 1) return false;
  const char* names[] = {"min_cluster_energy", "pi0_partner_min_energy", "eta_partner_min_energy", "pi0_mass_min", "pi0_mass_max", "eta_mass_min", "eta_mass_max", "missing_energy_min", "missing_energy_max", "min_photon_energy_recovery", "partner_diagnostic_min_cluster_energy"};
  values.assign(sizeof(names) / sizeof(*names), -1.0);
  for (std::size_t i = 0; i < values.size(); ++i)
    if (!tree->GetBranch(names[i]) || tree->SetBranchAddress(names[i], &values[i]) < 0) return false;
  if (tree->GetEntry(0) <= 0) return false;
  return std::all_of(values.begin(), values.end(), [](double v) { return std::isfinite(v) && v >= 0.0; }) &&
      values[3] < values[4] && values[5] < values[6] && values[7] < values[8] && values[9] == 0.5 && values[10] == 0.0;
}
inline bool same(const std::vector<double>& a, const std::vector<double>& b, std::size_t begin = 0)
{
  if (a.size() != 11 || b.size() != 11) return false;
  for (std::size_t i = begin; i < a.size(); ++i)
    if (!std::isfinite(a[i]) || !std::isfinite(b[i]) || std::abs(a[i] - b[i]) > 1e-11 * std::max({1.0, std::abs(a[i]), std::abs(b[i])})) return false;
  return true;
}
inline bool validate_partials(TTree& tree, const char* branch)
{
  int source_schema = -1, topology_version = -1;
  std::vector<double>* settings = nullptr;
  if (!tree.GetBranch(branch) || tree.SetBranchAddress(branch, &settings) < 0 ||
      tree.SetBranchAddress("source_map_schema_version", &source_schema) < 0 ||
      tree.SetBranchAddress("pi0_topology_algorithm_version", &topology_version) < 0) return false;
  std::vector<double> reference;
  bool ok = tree.GetEntries() > 0;
  for (Long64_t entry = 0; ok && entry < tree.GetEntries(); ++entry)
  {
    ok = tree.GetEntry(entry) > 0 && settings && source_schema == 5 && topology_version == 9;
    if (!ok) break;
    if (entry == 0) reference = *settings;
    ok = same(reference, *settings);
  }
  tree.ResetBranchAddresses();
  return ok;
}

}
#endif
