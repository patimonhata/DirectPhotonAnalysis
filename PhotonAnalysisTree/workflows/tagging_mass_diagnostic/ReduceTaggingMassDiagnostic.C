#include "TaggingMassDiagnostic.h"

int ReduceTaggingMassDiagnostic(const std::string family, const std::string map_root, const std::string output_base,
    const bool require_complete, const std::string sample, const int shard_index = 0, const int shard_count = 1, const Long64_t max_events = 0)
{
  using namespace tagging_mass;
  const auto definitions = samples(family);
  const auto definition = std::find_if(definitions.begin(), definitions.end(), [&](const auto& s) { return s.name == sample; });
  if (definition == definitions.end() || map_root.empty() || output_base.empty() || shard_count <= 0 || shard_index < 0 || shard_index >= shard_count || max_events < 0) return 1;
  const std::string output_path = output_base + "/tagging_mass_diagnostic.root";
  if (!gSystem->AccessPathName(output_path.c_str())) { std::cerr << "Output already exists: " << output_path << std::endl; return 2; }
  Metadata metadata;
  metadata.family = family; metadata.sample = sample; metadata.require_complete = require_complete;
  metadata.shard_index = shard_index; metadata.shard_count = shard_count; metadata.max_events = max_events;
  Input input;
  double tag_threshold = -1;
  if (!inspect(*definition, map_root + "/" + sample + "/map_*.root", require_complete, input, metadata.release, metadata.model,
               tag_threshold, metadata.selection_settings, metadata.region_settings)) { std::cerr << "Map inspection failed" << std::endl; return 3; }
  metadata.map_count = input.maps; metadata.sumw = input.sumw;
  const auto& settings = metadata.selection_settings;
  TChain tree("event_tree");
  if (tree.Add(input.pattern.c_str()) <= 0) return 3;
  metadata.total_entries = tree.GetEntries();
  metadata.entry_begin = metadata.total_entries * shard_index / shard_count;
  metadata.entry_end = metadata.total_entries * (shard_index + 1) / shard_count;
  if (max_events > 0) metadata.entry_end = std::min(metadata.entry_end, metadata.entry_begin + max_events);
  if (metadata.entry_end <= metadata.entry_begin) return 3;
  tree.SetBranchStatus("*", false); tree.SetCacheSize(64LL * 1024 * 1024);
  unsigned char weight_valid = 0, stitch_valid = 0, stitch_pass = 0;
  unsigned int ncluster = 0;
  double weight_numerator = 0;
  std::vector<unsigned int>* id = nullptr;
  std::vector<double> *energy = nullptr, *et = nullptr, *eta = nullptr, *phi = nullptr;
  std::vector<unsigned char> *region_a = nullptr, *anchor_valid = nullptr, *prompt = nullptr, *pi0_tag = nullptr, *eta_tag = nullptr;
  std::vector<float> *fraction = nullptr, *main_fraction = nullptr, *truth_mass = nullptr, *partner_e = nullptr;
  std::vector<int> *topology = nullptr, *partner_id = nullptr, *truth_status = nullptr;
#define B(name, var) if (!bind_active(tree, name, &var)) { std::cerr << "Missing/incompatible branch: " << name << std::endl; return 4; }
  B("event_weight_valid", weight_valid); B("sample_stitching_valid", stitch_valid); B("sample_stitching_pass", stitch_pass);
  B("weight_numerator_pb", weight_numerator); B("split_ncluster", ncluster); B("split_cluster_id", id);
  B("split_cluster_e", energy); B("split_cluster_et", et); B("split_cluster_eta", eta); B("split_cluster_phi", phi);
  B("split_cluster_pass_region_a", region_a); B("split_cluster_pi0_anchor_valid", anchor_valid);
  B("split_cluster_pi0_anchor_main_fraction", main_fraction); B("split_cluster_pi0_anchor_topology", topology);
  B("split_cluster_pi0_anchor_truth_partner_mass", truth_mass); B("split_cluster_pi0_anchor_truth_partner_cluster_e", partner_e);
  B("split_cluster_pi0_anchor_truth_partner_cluster_id", partner_id);
  B("split_cluster_pi0_anchor_truth_partner_tag_status", truth_status);
  B("split_cluster_truth_prompt_cluster", prompt); B("split_cluster_truth_dominant_fraction", fraction);
  B("split_cluster_pi0_tag", pi0_tag); B("split_cluster_eta_tag", eta_tag);
#undef B
  tree.StopCacheLearningPhase();
  std::array<std::unique_ptr<Histograms>, 6> groups;
  for (std::size_t i = 0; i < groups.size(); ++i) groups[i] = std::make_unique<Histograms>(kSelectionKeys[i]);
  Flow flow;
  for (Long64_t entry = metadata.entry_begin; entry < metadata.entry_end; ++entry)
  {
    if (tree.GetEntry(entry) <= 0) return 5;
#define S(v) if (!v || v->size() != ncluster) return 5
    S(id); S(energy); S(et); S(eta); S(phi); S(region_a); S(anchor_valid); S(main_fraction); S(topology); S(truth_mass);
    S(partner_e); S(partner_id); S(truth_status); S(prompt); S(fraction); S(pi0_tag); S(eta_tag);
#undef S
    if (!weight_valid || !stitch_valid || !stitch_pass) continue;
    const double weight = weight_numerator / input.sumw;
    if (!std::isfinite(weight)) return 5;
    std::set<unsigned int> ids;
    for (unsigned int i = 0; i < ncluster; ++i)
      if (!ids.insert((*id)[i]).second || !std::isfinite((*energy)[i]) || !std::isfinite((*eta)[i]) || !std::isfinite((*phi)[i]) ||
          !std::isfinite((*et)[i]) || !((*energy)[i] > settings[0])) return 5;
    for (unsigned int i = 0; i < ncluster; ++i)
    {
      if (!(*region_a)[i]) continue;
      const double candidate_et = (*et)[i];
      const bool truth_pi0 = (*anchor_valid)[i] && (*main_fraction)[i] > 0.5F;
      const bool truth_prompt = (*prompt)[i] && (*fraction)[i] > 0.5F;
      if (truth_pi0 && truth_prompt) { std::cerr << "Overlapping truth origins" << std::endl; return 5; }
      if (truth_pi0 && (*topology)[i] == 1)
      {
        flow.fill(0, candidate_et, weight);
        const double mass = (*truth_mass)[i], pe = (*partner_e)[i];
        if ((*partner_id)[i] < 0 || static_cast<unsigned int>((*partner_id)[i]) == (*id)[i] || !std::isfinite(mass) || mass < 0 || !std::isfinite(pe) || pe <= 0)
          flow.fill(3, candidate_et, weight);
        else
        {
          flow.fill(1, candidate_et, weight);
          groups[0]->fill(mass, candidate_et, weight, settings[3], settings[4]);
          // Use the pre-rounding map decision, not the float-stored energy at the boundary.
          if ((*truth_status)[i] == 4)
          {
            flow.fill(2, candidate_et, weight);
            groups[1]->fill(mass, candidate_et, weight, settings[3], settings[4]);
          }
        }
      }
      if (!truth_prompt) continue;
      flow.fill(4, candidate_et, weight);
      bool found[2] = {false, false};
      for (unsigned int j = 0; j < ncluster; ++j)
      {
        if (i == j || !((*energy)[j] > std::min(settings[1], settings[2]))) continue;
        const double mass = pair_mass((*energy)[i], (*eta)[i], (*phi)[i], (*energy)[j], (*eta)[j], (*phi)[j]);
        if (!std::isfinite(mass) || mass < 0) return 5;
        for (int meson = 0; meson < 2; ++meson)
        {
          if (!((*energy)[j] > settings[1 + meson])) continue;
          const double low = settings[3 + 2 * meson], high = settings[4 + 2 * meson];
          groups[2 + 2 * meson]->fill(mass, candidate_et, weight, low, high);
          if (mass > low && mass < high)
          {
            found[meson] = true;
            groups[3 + 2 * meson]->fill(mass, candidate_et, weight, low, high);
          }
        }
      }
      if (found[0] != bool((*pi0_tag)[i]) || found[1] != bool((*eta_tag)[i]))
      {
        std::cerr << "Recomputed veto mismatch at entry " << entry << ", cluster " << (*id)[i] << std::endl;
        return 6;
      }
      if (found[0]) flow.fill(5, candidate_et, weight);
      if (found[1]) flow.fill(6, candidate_et, weight);
      if (found[0] || found[1]) flow.fill(7, candidate_et, weight);
    }
  }
  if (gSystem->mkdir(output_base.c_str(), true) != 0 && gSystem->AccessPathName(output_base.c_str())) return 7;
  const std::string temporary = output_path + ".tmp." + std::to_string(gSystem->GetPid());
  TFile output(temporary.c_str(), "CREATE");
  if (output.IsZombie()) return 7;
  for (std::size_t i = 0; i < groups.size(); ++i) { output.mkdir(kSelectionKeys[i])->cd(); groups[i]->write(); }
  output.cd(); flow.write(); metadata.write(); output.Close();
  if (output.TestBit(TFile::kWriteError) || !gSystem->AccessPathName(output_path.c_str()) || gSystem->Rename(temporary.c_str(), output_path.c_str()) != 0) return 7;
  std::cout << "Reduced " << sample << " shard " << shard_index << "/" << shard_count << ", " << metadata.entry_end - metadata.entry_begin << " events: " << output_path << std::endl;
  return 0;
}
