#ifndef PHOTON_CANDIDATE_REDUCE_IMPLEMENTATION
#define PHOTON_CANDIDATE_REDUCE_IMPLEMENTATION

#include "CandidateSelection.h"

namespace photon_candidate_analysis
{
int reduce(
    const std::string family,
    const std::string map_root,
    const std::string output_base,
    const bool require_complete,
    const int n_bins,
    const double et_max,
    const std::string sample_name,
    const int shard_index)
{
  using namespace candidate_composition;
  const std::vector<SampleDefinition> family_samples = sample_definitions(family);
  const auto selected_sample = std::find_if(family_samples.begin(), family_samples.end(), [&](const SampleDefinition& sample) { return sample.name == sample_name; });
  const int shard_count = required_shard_count(sample_name);
  if (selected_sample == family_samples.end() || shard_index < 0 || shard_index >= shard_count || map_root.empty() ||
      n_bins <= 0 || !std::isfinite(et_max) || et_max <= 0.0) return 1;
  std::string normalized_map_root = map_root;
  while (normalized_map_root.size() > 1U && normalized_map_root.back() == '/') normalized_map_root.pop_back();
  const std::string configuration_root = normalized_map_root.substr(normalized_map_root.find_last_of('/') + 1U) == "maps"
      ? normalized_map_root.substr(0, normalized_map_root.find_last_of('/')) : normalized_map_root;
  const std::string configuration = configuration_root.substr(configuration_root.find_last_of('/') + 1U);
  const std::string default_output_suffix = "/partials/" + family + "/" + sample_name + "/shard_" + std::to_string(shard_index) + "/photon_candidate_selection";
  const std::string resolved_output_base = output_base.empty()
      ? "/sphenix/user/ryotaro/DirectPhotonAnalysis/PhotonAnalysisTree/output/photon_candidate_selection/productions/topology_v12_eta_v2_decay_v1/" +
            configuration + default_output_suffix
      : output_base;
  if (!make_output_directory(resolved_output_base)) return 2;

  std::array<std::unique_ptr<Histograms>, kSelectionCount> composition_histograms;
  std::array<std::unique_ptr<Spectra>, kSelectionCount> topology_histograms;
  std::array<pi0_bdt::Histograms, kSelectionCount> bdt_histograms;
  std::array<std::unique_ptr<Spectra>, kSelectionCount> eta_histograms;
  std::array<pi0_bdt::Histograms, kSelectionCount> eta_bdt_histograms;
  for (std::size_t selection_index = 0; selection_index < kSelectionCount; ++selection_index)
  {
    composition_histograms[selection_index] = std::make_unique<Histograms>(n_bins, et_max);
    topology_histograms[selection_index] = std::make_unique<Spectra>(n_bins, et_max);
    eta_histograms[selection_index] = std::make_unique<Spectra>(n_bins, et_max);
  }
  std::vector<std::string> sample_names;
  std::vector<unsigned long long> sample_map_counts;
  std::vector<double> sample_sum_generator_weights;
  std::string analysis_release, model_sha256;
  double min_cluster_energy = -1.0;
  double partner_diagnostic_min_cluster_energy = -1.0;
  double pi0_partner_min_energy = -1.0;
  double eta_partner_min_energy = -1.0;
  double pi0_mass_min = -1.0;
  double pi0_mass_max = -1.0;
  double eta_mass_min = -1.0;
  double eta_mass_max = -1.0;
  double missing_energy_min = -1.0;
  double missing_energy_max = -1.0;
  double min_photon_energy_recovery = -1.0;
  int pi0_topology_algorithm_version = -1;
  int eta_topology_algorithm_version = -1;
  decay_origin::Selections diagnostic_histograms;
  for (auto& h : diagnostic_histograms) h = std::make_unique<decay_origin::Histograms>(n_bins, et_max);
  other_origin::Selections other_histograms;
  for (auto& h : other_histograms) h = std::make_unique<other_origin::Histograms>(n_bins, et_max);
  eta_decay::Selections eta_decay_histograms;
  for (auto& hist : eta_decay_histograms) hist = std::make_unique<eta_decay::Histograms>(n_bins, et_max);
  std::array<Summary, kSelectionCount> composition_summaries;
  unsigned long long events_written = 0, expected_events_written = 0, events_stitch_pass = 0;
  unsigned long long region_a_clusters = 0, region_a_prompt_clusters = 0, region_a_anchor_clusters = 0;
  std::array<unsigned long long, kSelectionCount> topology_selected_clusters = {};
  std::array<unsigned long long, kSelectionCount> topology_selected_anchor_clusters = {};

  const SampleDefinition& sample = *selected_sample;
  unsigned long long shard_map_begin = 0, shard_map_end = 0, total_map_count = 0;
  std::vector<MapMetadata> maps;
  if (!collect_maps(map_root + "/" + sample.name + "/map_*.root", sample, require_complete, maps)) return 3;
  if (maps.empty()) return 3;
  analysis_release = maps.front().analysis_release;
  model_sha256 = maps.front().model_sha256;
  min_cluster_energy = maps.front().min_cluster_energy;
  partner_diagnostic_min_cluster_energy = maps.front().partner_diagnostic_min_cluster_energy;
  pi0_partner_min_energy = maps.front().pi0_partner_min_energy;
  eta_partner_min_energy = maps.front().eta_partner_min_energy;
  pi0_mass_min = maps.front().pi0_mass_min;
  pi0_mass_max = maps.front().pi0_mass_max;
  eta_mass_min = maps.front().eta_mass_min;
  eta_mass_max = maps.front().eta_mass_max;
  missing_energy_min = maps.front().missing_energy_min;
  missing_energy_max = maps.front().missing_energy_max;
  min_photon_energy_recovery = maps.front().min_photon_energy_recovery;
  set_missing_labels(missing_energy_min, missing_energy_max);
  pi0_topology_algorithm_version = maps.front().pi0_topology_algorithm_version;
  eta_topology_algorithm_version = maps.front().eta_topology_algorithm_version;
  double sample_sumw = 0.0;
  for (const auto& map : maps) sample_sumw += map.sum_generator_weight_processed;
  if (!std::isfinite(sample_sumw) || sample_sumw <= 0.0) return 4;
  total_map_count = maps.size();
  shard_map_begin = total_map_count * static_cast<unsigned long long>(shard_index) / static_cast<unsigned long long>(shard_count);
  shard_map_end = total_map_count * static_cast<unsigned long long>(shard_index + 1) / static_cast<unsigned long long>(shard_count);
  if (shard_map_begin >= shard_map_end) return 4;
  for (unsigned long long map_index = shard_map_begin; map_index < shard_map_end; ++map_index)
    expected_events_written += maps[map_index].events_written;
  sample_names.push_back(sample.name);
  sample_map_counts.push_back(maps.size());
  sample_sum_generator_weights.push_back(sample_sumw);

  for (unsigned long long map_index = shard_map_begin; map_index < shard_map_end; ++map_index)
  {
    const MapMetadata& map = maps[map_index];
    TFile file(map.path.c_str(), "READ");
    auto* tree = file.Get<TTree>("event_tree");
    if (file.IsZombie() || !tree) return 5;
    decay_origin::Reader diagnostic;
    unsigned int ncluster = 0;
    unsigned char weight_valid = 0, stitch_valid = 0, stitch_pass = 0;
    double weight_numerator = 0.0;
    std::vector<double>* et = nullptr;
    std::vector<float>* bdt_score = nullptr;
    std::vector<unsigned char>* bdt_valid = nullptr;
    std::vector<unsigned char>* pass_nonisolated = nullptr;
    std::vector<unsigned char>* pass_kinematics = nullptr;
    std::vector<unsigned char>* pass_preselection = nullptr;
    std::vector<unsigned char>* pass_tight = nullptr;
    std::vector<unsigned char>* pass_isolated = nullptr;
    std::vector<unsigned char>* pass_region_a = nullptr;
    std::vector<unsigned char>* pi0_tag = nullptr;
    std::vector<unsigned char>* eta_tag = nullptr;
    std::vector<unsigned char>* pass_final_photon = nullptr;
    std::vector<unsigned char>* truth_valid = nullptr;
    std::vector<unsigned char>* prompt_flag = nullptr;
    std::vector<float>* dominant_fraction = nullptr;
    std::vector<unsigned char>* pi0_valid = nullptr;
    std::vector<float>* pi0_fraction = nullptr;
    std::vector<int>* topology = nullptr;
    std::vector<int>* missing_category = nullptr;
    std::vector<unsigned char>* eta_valid = nullptr;
    std::vector<float>* eta_main_fraction = nullptr;
    std::vector<int>* eta_topology = nullptr;
    std::vector<int>* eta_missing_category = nullptr;
    std::vector<unsigned int>* offset = nullptr;
    std::vector<int>* g4_pdg = nullptr;
    std::vector<int>* embedding = nullptr;
    std::vector<float>* fraction = nullptr;
    std::vector<int>* source = nullptr;
    std::vector<int>* parent = nullptr;
    std::vector<int>* photon_category = nullptr;
    std::vector<int>* eta_parent_barcode = nullptr;
    std::vector<int>* eta_parent_track_id = nullptr;
    std::vector<int>* eta_decay_mode = nullptr;
    std::vector<int>* eta_split_status = nullptr;
    std::vector<int>* pi0_parent_barcode = nullptr;
    std::vector<int>* pi0_parent_track_id = nullptr;

    tree->SetBranchStatus("*", false);
    tree->SetCacheSize(kEventTreeCacheSize);
    bool ok = bind_active(tree, "split_ncluster", &ncluster) &&
        bind_active(tree, "event_weight_valid", &weight_valid) &&
        bind_active(tree, "sample_stitching_valid", &stitch_valid) &&
        bind_active(tree, "sample_stitching_pass", &stitch_pass) &&
        bind_active(tree, "weight_numerator_pb", &weight_numerator) &&
        bind_active(tree, "split_cluster_bdt_score", &bdt_score) &&
        bind_active(tree, "split_cluster_bdt_valid", &bdt_valid) &&
        bind_active(tree, "split_cluster_pass_nonisolated", &pass_nonisolated) &&
        bind_active(tree, "split_cluster_et", &et) &&
        bind_active(tree, "split_cluster_pass_kinematics", &pass_kinematics) &&
        bind_active(tree, "split_cluster_pass_preselection", &pass_preselection) &&
        bind_active(tree, "split_cluster_pass_tight", &pass_tight) &&
        bind_active(tree, "split_cluster_pass_isolated", &pass_isolated) &&
        bind_active(tree, "split_cluster_pass_region_a", &pass_region_a) &&
        bind_active(tree, "split_cluster_pi0_tag", &pi0_tag) &&
        bind_active(tree, "split_cluster_eta_tag", &eta_tag) &&
        bind_active(tree, "split_cluster_pass_final_photon", &pass_final_photon) &&
        bind_active(tree, "split_cluster_truth_valid", &truth_valid) &&
        bind_active(tree, "split_cluster_truth_prompt_cluster", &prompt_flag) &&
        bind_active(tree, "split_cluster_truth_dominant_fraction", &dominant_fraction) &&
        bind_active(tree, "split_cluster_pi0_anchor_valid", &pi0_valid) &&
        bind_active(tree, "split_cluster_pi0_anchor_main_fraction", &pi0_fraction) &&
        bind_active(tree, "split_cluster_pi0_anchor_topology", &topology) &&
        bind_active(tree, "split_cluster_pi0_anchor_missing_category", &missing_category) &&
        bind_active(tree, "split_cluster_eta_anchor_valid", &eta_valid) &&
        bind_active(tree, "split_cluster_eta_anchor_main_fraction", &eta_main_fraction) &&
        bind_active(tree, "split_cluster_eta_anchor_topology", &eta_topology) &&
        bind_active(tree, "split_cluster_eta_anchor_missing_category", &eta_missing_category) &&
        bind_active(tree, "split_cluster_truth_contributor_offset", &offset) &&
        bind_active(tree, "split_cluster_truth_contributor_g4_pdg_id", &g4_pdg) &&
        bind_active(tree, "split_cluster_truth_contributor_embedding_id", &embedding) &&
        bind_active(tree, "split_cluster_truth_contributor_fraction", &fraction) &&
        bind_active(tree, "split_cluster_truth_contributor_photon_source", &source) &&
        bind_active(tree, "split_cluster_truth_contributor_classification_parent_pdg", &parent);
    ok &= diagnostic.bind(tree);
    ok &= bind_active(tree, "split_cluster_truth_contributor_photon_category", &photon_category);
    ok &= bind_active(tree, "split_cluster_truth_contributor_eta_parent_barcode", &eta_parent_barcode);
    ok &= bind_active(tree, "split_cluster_truth_contributor_eta_parent_track_id", &eta_parent_track_id);
    ok &= bind_active(tree, "split_cluster_truth_contributor_eta_decay_mode", &eta_decay_mode);
    ok &= bind_active(tree, "split_cluster_truth_contributor_eta_split_status", &eta_split_status);
    ok &= bind_active(tree, "split_cluster_truth_contributor_pi0_parent_barcode", &pi0_parent_barcode);
    ok &= bind_active(tree, "split_cluster_truth_contributor_pi0_parent_track_id", &pi0_parent_track_id);
    if (!ok) return 5;
    tree->StopCacheLearningPhase();
    events_written += static_cast<unsigned long long>(tree->GetEntries());

    for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry)
    {
      tree->GetEntry(entry);
      if (!bdt_score || !bdt_valid || !pass_nonisolated || bdt_score->size() != ncluster || bdt_valid->size() != ncluster || pass_nonisolated->size() != ncluster) return 6;
      if (!et || !pass_kinematics || !pass_preselection || !pass_tight || !pass_isolated || !pass_region_a || !pi0_tag || !eta_tag ||
          !pass_final_photon || !truth_valid || !prompt_flag || !dominant_fraction || !pi0_valid || !pi0_fraction || !topology || !missing_category ||
          !offset || !g4_pdg || !embedding || !fraction || !source || !parent || et->size() != ncluster || pass_kinematics->size() != ncluster ||
          pass_preselection->size() != ncluster || pass_tight->size() != ncluster || pass_isolated->size() != ncluster ||
          pass_region_a->size() != ncluster || pi0_tag->size() != ncluster || eta_tag->size() != ncluster || pass_final_photon->size() != ncluster ||
          truth_valid->size() != ncluster || prompt_flag->size() != ncluster || dominant_fraction->size() != ncluster ||
          pi0_valid->size() != ncluster || pi0_fraction->size() != ncluster || topology->size() != ncluster || missing_category->size() != ncluster ||
          !valid_contributors(ncluster, *offset, *g4_pdg, *embedding, *fraction, *source, *parent)) return 6;
      if (!eta_valid || !eta_main_fraction || !eta_topology || !eta_missing_category || eta_valid->size() != ncluster ||
          eta_main_fraction->size() != ncluster || eta_topology->size() != ncluster || eta_missing_category->size() != ncluster) return 6;
      if (!photon_category || photon_category->size() != fraction->size() || !diagnostic.consistent(ncluster, *offset, *fraction)) return 6;
      if (!eta_parent_barcode || eta_parent_barcode->size() != fraction->size()) return 6;
      if (!eta_parent_track_id || eta_parent_track_id->size() != fraction->size()) return 6;
      if (!eta_decay_mode || eta_decay_mode->size() != fraction->size()) return 6;
      if (!eta_split_status || eta_split_status->size() != fraction->size()) return 6;
      if (!pi0_parent_barcode || pi0_parent_barcode->size() != fraction->size()) return 6;
      if (!pi0_parent_track_id || pi0_parent_track_id->size() != fraction->size()) return 6;
      if (!weight_valid || !stitch_valid || !stitch_pass) continue;
      ++events_stitch_pass;
      const double weight = weight_numerator / sample_sumw;
      if (!std::isfinite(weight)) return 6;
      for (std::size_t cluster = 0; cluster < ncluster; ++cluster)
      {
        const double cluster_et = (*et)[cluster];
        if (!std::isfinite(cluster_et)) return 6;
        const auto eta_origin = eta_origin::summarize((*offset)[cluster], (*offset)[cluster + 1], *embedding, *fraction,
            *eta_parent_barcode, *eta_parent_track_id, *eta_decay_mode, *eta_split_status, *pi0_parent_barcode, *pi0_parent_track_id);
        const double eta_contribution = eta_origin.fraction;
        if (!std::isfinite(eta_contribution) || eta_contribution < -1e-6 || eta_contribution > 1.0 + 1e-5) return 6;
        const bool is_prompt = (*prompt_flag)[cluster] && (*dominant_fraction)[cluster] > kMajorityThreshold;
        const bool is_pi0 = (*pi0_valid)[cluster] && (*pi0_fraction)[cluster] > kMajorityThreshold;
        const bool is_eta = !is_pi0 && eta_contribution > kMajorityThreshold;
        const int majority_count = static_cast<int>(is_prompt) + static_cast<int>(is_pi0) + static_cast<int>(is_eta);
        const bool half_boundary = ((*prompt_flag)[cluster] && (*dominant_fraction)[cluster] == kMajorityThreshold) ||
            ((*pi0_valid)[cluster] && (*pi0_fraction)[cluster] == kMajorityThreshold) || eta_contribution == kMajorityThreshold;
        const std::size_t category = origin_category(is_prompt, is_pi0, is_eta, (*topology)[cluster]);
        if (category == category_count) return 6;
        const bool computed_final_photon = (*pass_region_a)[cluster] && !(*pi0_tag)[cluster] && !(*eta_tag)[cluster];
        if (static_cast<bool>((*pass_final_photon)[cluster]) != computed_final_photon) return 6;
        if ((*pass_region_a)[cluster])
        {
          ++region_a_clusters;
          if ((*prompt_flag)[cluster]) ++region_a_prompt_clusters;
          if (category >= pi0_separated && category <= pi0_other) ++region_a_anchor_clusters;
        }
        std::array<bool, kSelectionCount> selected_by = {};
        for (std::size_t selection_index = 0; selection_index < kSelectionCount; ++selection_index)
          selected_by[selection_index] = pass_selection(selection_index, (*pass_kinematics)[cluster], (*pass_preselection)[cluster], (*pass_tight)[cluster],
              (*pass_isolated)[cluster], (*pass_region_a)[cluster], (*pi0_tag)[cluster], (*eta_tag)[cluster]);
        for (std::size_t topology_selection = 0; topology_selection < kSelectionCount; ++topology_selection)
        {
          if (!selected_by[topology_selection]) continue;
          ++topology_selected_clusters[topology_selection];
          Spectra& spectra = *topology_histograms[topology_selection];
          if ((*prompt_flag)[cluster]) spectra.fill(0, cluster_et, weight);
          auto& eta_spectra = *eta_histograms[topology_selection];
          if ((*prompt_flag)[cluster]) eta_spectra.fill(0, cluster_et, weight);
          if (category == eta && (*eta_valid)[cluster] && (*eta_main_fraction)[cluster] > kMajorityThreshold)
          {
            if (!eta_topology::fill(eta_spectra, (*eta_topology)[cluster], (*eta_missing_category)[cluster], cluster_et, weight) ||
                !eta_bdt_histograms[topology_selection].fill(eta_topology::bdt_category((*eta_topology)[cluster]), (*bdt_score)[cluster],
                    (*bdt_valid)[cluster], (*pass_isolated)[cluster], (*pass_nonisolated)[cluster], cluster_et, weight)) return 6;
          }
          if (category < pi0_separated || category > pi0_other) continue;
          if (!bdt_histograms[topology_selection].fill(static_cast<int>(category - pi0_separated), (*bdt_score)[cluster], (*bdt_valid)[cluster],
              (*pass_isolated)[cluster], (*pass_nonisolated)[cluster], cluster_et, weight)) return 6;
          ++topology_selected_anchor_clusters[topology_selection];
          spectra.fill(1, cluster_et, weight);
          const int topology_value = (*topology)[cluster];
          if (topology_value == 1) spectra.fill(2, cluster_et, weight);
          else if (topology_value == 2) spectra.fill(3, cluster_et, weight);
          else if (topology_value == 4) spectra.fill(4, cluster_et, weight);
          else if (topology_value == 0) spectra.fill(12, cluster_et, weight);
          else if (topology_value == 3)
          {
            spectra.fill(5, cluster_et, weight);
            const int missing_value = (*missing_category)[cluster];
            if (missing_value == 1) spectra.fill(6, cluster_et, weight);
            else if (missing_value == 4) spectra.fill(7, cluster_et, weight);
            else if (missing_value == 2) spectra.fill(8, cluster_et, weight);
            else if (missing_value == 5) spectra.fill(9, cluster_et, weight);
            else if (missing_value == 6) spectra.fill(10, cluster_et, weight);
            else if (missing_value == 3) spectra.fill(11, cluster_et, weight);
            else return 6;
          }
          else return 6;
        }
        if (std::none_of(selected_by.begin(), selected_by.end(), [](bool selected) { return selected; })) continue;
        for (std::size_t selection_index = 0; selection_index < kSelectionCount; ++selection_index)
        {
          if (!selected_by[selection_index]) continue;
          Summary& summary = composition_summaries[selection_index];
          Histograms& spectra = *composition_histograms[selection_index];
          ++summary.selected;
          spectra.fill(denominator, cluster_et, weight);
          if (!(*truth_valid)[cluster]) ++summary.invalid_truth;
          if (half_boundary) ++summary.half_boundary;
          if (majority_count > 1) ++summary.overlap;
          spectra.fill(category, cluster_et, weight);
          if (category == prompt) ++summary.prompt;
          else if (category >= pi0_separated && category <= pi0_other) ++summary.pi0;
          else if (category == eta)
          {
            ++summary.eta;
            eta_decay_histograms[selection_index]->fill(eta_origin.mode, eta_origin.unresolved, cluster_et, weight);
          }
          else
          {
            ++summary.other;
            diagnostic_histograms[selection_index]->fill(decay_origin::classify(diagnostic, cluster), cluster_et, weight, (*bdt_score)[cluster], (*bdt_valid)[cluster]);
            const auto detail = other_origin::classify((*truth_valid)[cluster], (*offset)[cluster], (*offset)[cluster + 1],
                *g4_pdg, *embedding, *fraction, *parent, *photon_category, majority_count > 1);
            other_histograms[selection_index]->fill(detail, cluster_et, weight, (*bdt_score)[cluster], (*bdt_valid)[cluster]);
          }
        }
      }
    }
  }

  if (events_written != expected_events_written || sample_names.empty()) return 7;
  for (std::size_t selection_index = 0; selection_index < kSelectionCount; ++selection_index)
  {
    const Summary& summary = composition_summaries[selection_index];
    if (!consistent_pi0_bdt(*composition_histograms[selection_index], *topology_histograms[selection_index], bdt_histograms[selection_index]) ||
        !eta_topology::valid(*composition_histograms[selection_index], *eta_histograms[selection_index], eta_bdt_histograms[selection_index])) return 7;
    if (summary.selected != summary.prompt + summary.pi0 + summary.eta + summary.other ||
        !candidate_composition::valid_partition(composition_histograms[selection_index]->counts) ||
        !candidate_composition::valid_partition(composition_histograms[selection_index]->weighted) ||
        !::valid_partition(topology_histograms[selection_index]->counts) || !::valid_partition(topology_histograms[selection_index]->weighted_pb)) return 7;
  }
  for (std::size_t i = 0; i < kSelectionCount; ++i)
    if (!eta_decay_histograms[i]->valid(*composition_histograms[i]->counts[eta], *composition_histograms[i]->weighted[eta],
                                      *eta_histograms[i]->counts[1], *eta_histograms[i]->weighted_pb[1])) return 7;
  for (std::size_t i = 0; i < kSelectionCount; ++i)
    if (!other_histograms[i]->valid(*composition_histograms[i]->counts[other], *composition_histograms[i]->weighted[other])) return 7;
  for (std::size_t i = 0; i < kSelectionCount; ++i)
    if (!diagnostic_histograms[i]->valid(*composition_histograms[i]->counts[other], *composition_histograms[i]->weighted[other])) return 7;
  TFile output((resolved_output_base + ".root").c_str(), "RECREATE");
  if (output.IsZombie()) return 8;
  auto* diagnostic_root = output.mkdir("decay_origin");
  if (!diagnostic_root) return 8;
  for (std::size_t i = 0; i < kSelectionCount; ++i)
  {
    auto* directory = diagnostic_root->mkdir(kSelectionKeys[i]);
    if (!directory || !diagnostic_histograms[i]->write(*directory)) return 8;
  }
  auto* other_root = output.mkdir("other_origin");
  if (!other_root) return 8;
  for (std::size_t i = 0; i < kSelectionCount; ++i)
  {
    auto* directory = other_root->mkdir(kSelectionKeys[i]);
    if (!directory || !other_histograms[i]->write(*directory)) return 8;
  }
  auto* decay_root = output.mkdir("eta_decay");
  if (!decay_root) return 8;
  for (std::size_t i = 0; i < kSelectionCount; ++i)
  {
    auto* directory = decay_root->mkdir(kSelectionKeys[i]);
    if (!directory || !eta_decay_histograms[i]->write(*directory)) return 8;
  }
  TDirectory* composition_root_directory = output.mkdir("composition");
  if (!composition_root_directory) return 8;
  for (std::size_t selection_index = 0; selection_index < kSelectionCount; ++selection_index)
  {
    TDirectory* composition_directory = composition_root_directory->mkdir(kSelectionKeys[selection_index]);
    if (!composition_directory) return 8;
    composition_directory->cd();
    for (std::size_t index = 0; index < category_count; ++index)
    {
      composition_histograms[selection_index]->counts[index]->Write();
      composition_histograms[selection_index]->weighted[index]->Write();
    }
  }
  output.cd();
  TDirectory* topology_root_directory = output.mkdir("anchor_topology");
  if (!topology_root_directory) return 8;
  for (std::size_t topology_selection = 0; topology_selection < kSelectionCount; ++topology_selection)
  {
    TDirectory* topology_directory = topology_root_directory->mkdir(kSelectionKeys[topology_selection]);
    if (!topology_directory) return 8;
    topology_directory->cd();
    for (std::size_t index = 0; index < kSpectrumCount; ++index)
    {
      topology_histograms[topology_selection]->counts[index]->Write();
      topology_histograms[topology_selection]->weighted_pb[index]->Write();
    }
  }
  output.cd();
  auto* bdt_directory = output.mkdir("pi0_bdt");
  if (!bdt_directory) return 8;
  for (std::size_t i = 0; i < kSelectionCount; ++i)
  {
    auto* selection_directory = bdt_directory->mkdir(kSelectionKeys[i]);
    if (!selection_directory || !bdt_histograms[i].write(*selection_directory)) return 8;
  }
  output.cd();
  if (!eta_topology::write_partial(output, eta_histograms, eta_bdt_histograms)) return 8;
  output.cd();
  int schema_version = 10, source_schema_version = 8, signal_embedding_id = kSignalEmbeddingId;
  int decay_origin_algorithm_version = maps.front().decay_origin_algorithm_version;
  double decay_origin_max_proper_length_cm = maps.front().decay_origin_max_proper_length_cm;
  double decay_origin_max_flight_cm = maps.front().decay_origin_max_flight_cm;
  int other_origin_algorithm_version = other_origin::kAlgorithmVersion;
  double majority_threshold = kMajorityThreshold;
  std::string majority_comparison = "strictly_greater_than";
  std::string eta_definition = "max_same_eta_parent_descendant_fraction_after_single_pi0_majority_priority";
  std::string other_definition = "no_unique_prompt_pi0_eta_strict_majority_including_exact_0p5";
  std::string weight_definition = "cross_section_pb_times_generator_weight_divided_by_full_sample_sum_generator_weight_processed";
  std::string metadata_family = family;
  bool metadata_require_complete = require_complete;
  int metadata_n_bins = n_bins;
  double metadata_et_max = et_max;
  double sample_cross_section_pb = sample.cross_section_pb;
  TTree metadata("metadata", "Photon-candidate selection partial metadata");
  metadata.Branch("schema_version", &schema_version);
  metadata.Branch("source_map_schema_version", &source_schema_version);
  metadata.Branch("other_origin_algorithm_version", &other_origin_algorithm_version);
  metadata.Branch("decay_origin_algorithm_version", &decay_origin_algorithm_version);
  metadata.Branch("decay_origin_max_proper_length_cm", &decay_origin_max_proper_length_cm);
  metadata.Branch("decay_origin_max_flight_cm", &decay_origin_max_flight_cm);
  metadata.Branch("family", &metadata_family);
  metadata.Branch("map_root", &normalized_map_root);
  metadata.Branch("require_complete", &metadata_require_complete);
  metadata.Branch("n_bins", &metadata_n_bins);
  metadata.Branch("et_max", &metadata_et_max);
  metadata.Branch("min_cluster_energy", &min_cluster_energy);
  metadata.Branch("partner_diagnostic_min_cluster_energy", &partner_diagnostic_min_cluster_energy);
  metadata.Branch("pi0_partner_min_energy", &pi0_partner_min_energy);
  metadata.Branch("eta_partner_min_energy", &eta_partner_min_energy);
  metadata.Branch("pi0_mass_min", &pi0_mass_min);
  metadata.Branch("pi0_mass_max", &pi0_mass_max);
  metadata.Branch("eta_mass_min", &eta_mass_min);
  metadata.Branch("eta_mass_max", &eta_mass_max);
  metadata.Branch("missing_energy_min", &missing_energy_min);
  metadata.Branch("missing_energy_max", &missing_energy_max);
  metadata.Branch("min_photon_energy_recovery", &min_photon_energy_recovery);
  metadata.Branch("pi0_topology_algorithm_version", &pi0_topology_algorithm_version);
  metadata.Branch("eta_topology_algorithm_version", &eta_topology_algorithm_version);
  metadata.Branch("signal_embedding_id", &signal_embedding_id);
  metadata.Branch("majority_threshold", &majority_threshold);
  metadata.Branch("majority_comparison", &majority_comparison);
  metadata.Branch("eta_definition", &eta_definition);
  metadata.Branch("other_definition", &other_definition);
  metadata.Branch("weight_definition", &weight_definition);
  metadata.Branch("analysis_release", &analysis_release);
  metadata.Branch("model_sha256", &model_sha256);
  metadata.Branch("sample_names", &sample_names);
  metadata.Branch("sample_map_counts", &sample_map_counts);
  metadata.Branch("sample_sum_generator_weights", &sample_sum_generator_weights);
  metadata.Branch("sample_cross_section_pb", &sample_cross_section_pb);
  metadata.Branch("events_written", &events_written);
  metadata.Branch("events_stitch_pass", &events_stitch_pass);
  metadata.Branch("region_a_clusters", &region_a_clusters);
  metadata.Branch("region_a_prompt_clusters", &region_a_prompt_clusters);
  metadata.Branch("region_a_anchor_clusters", &region_a_anchor_clusters);
  metadata.Fill();
  metadata.Write();
  TTree composition_summary("composition_summary", "Per-selection candidate-composition partial summary");
  unsigned int composition_selection_index = 0;
  std::string composition_selection_key;
  unsigned long long composition_selected = 0, composition_prompt = 0, composition_pi0 = 0, composition_eta = 0, composition_other = 0;
  unsigned long long composition_overlap = 0, composition_half_boundary = 0, composition_invalid_truth = 0;
  composition_summary.Branch("selection_index", &composition_selection_index);
  composition_summary.Branch("selection_key", &composition_selection_key);
  composition_summary.Branch("selected_clusters", &composition_selected);
  composition_summary.Branch("prompt_clusters", &composition_prompt);
  composition_summary.Branch("pi0_clusters", &composition_pi0);
  composition_summary.Branch("eta_clusters", &composition_eta);
  composition_summary.Branch("other_clusters", &composition_other);
  composition_summary.Branch("overlap_clusters", &composition_overlap);
  composition_summary.Branch("half_boundary_clusters", &composition_half_boundary);
  composition_summary.Branch("invalid_truth_clusters", &composition_invalid_truth);
  for (composition_selection_index = 0; composition_selection_index < kSelectionCount; ++composition_selection_index)
  {
    composition_selection_key = kSelectionKeys[composition_selection_index];
    const Summary& summary = composition_summaries[composition_selection_index];
    composition_selected = summary.selected;
    composition_prompt = summary.prompt;
    composition_pi0 = summary.pi0;
    composition_eta = summary.eta;
    composition_other = summary.other;
    composition_overlap = summary.overlap;
    composition_half_boundary = summary.half_boundary;
    composition_invalid_truth = summary.invalid_truth;
    composition_summary.Fill();
  }
  composition_summary.Write();
  TTree topology_summary("topology_summary", "Per-selection anchor-topology partial summary");
  unsigned int topology_selection_index = 0;
  std::string topology_selection_key;
  unsigned long long topology_selection_selected = 0, topology_selection_anchor = 0;
  topology_summary.Branch("selection_index", &topology_selection_index);
  topology_summary.Branch("selection_key", &topology_selection_key);
  topology_summary.Branch("selected_clusters", &topology_selection_selected);
  topology_summary.Branch("selected_anchor_clusters", &topology_selection_anchor);
  for (topology_selection_index = 0; topology_selection_index < kSelectionCount; ++topology_selection_index)
  {
    topology_selection_key = kSelectionKeys[topology_selection_index];
    topology_selection_selected = topology_selected_clusters[topology_selection_index];
    topology_selection_anchor = topology_selected_anchor_clusters[topology_selection_index];
    topology_summary.Fill();
  }
  topology_summary.Write();
  std::string shard_sample_name = sample_name;
  int metadata_shard_index = shard_index, metadata_shard_count = shard_count;
  TTree shard_metadata("shard_metadata", "Photon-candidate selection partial shard metadata");
  shard_metadata.Branch("sample_name", &shard_sample_name);
  shard_metadata.Branch("shard_index", &metadata_shard_index);
  shard_metadata.Branch("shard_count", &metadata_shard_count);
  shard_metadata.Branch("map_index_begin", &shard_map_begin);
  shard_metadata.Branch("map_index_end", &shard_map_end);
  shard_metadata.Branch("total_map_count", &total_map_count);
  shard_metadata.Fill();
  shard_metadata.Write();
  output.Close();
  std::cout << "ReducePythiaPhotonCandidateSelection - family/sample/selections/output = "
            << family << "/" << sample_name << "/" << kSelectionCount << "/" << resolved_output_base;
  std::cout << " (partial " << sample_name << " shard " << shard_index << "/" << shard_count << ", maps " << shard_map_begin << ":" << shard_map_end << ")";
  std::cout << std::endl;
  const bool have_overlap = std::any_of(composition_summaries.begin(), composition_summaries.end(), [](const Summary& summary) { return summary.overlap != 0ULL; });
  return have_overlap ? 9 : 0;
}
}

#endif
