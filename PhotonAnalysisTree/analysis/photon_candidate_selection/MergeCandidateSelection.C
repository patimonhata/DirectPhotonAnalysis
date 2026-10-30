#ifndef PHOTON_CANDIDATE_MERGE_IMPLEMENTATION
#define PHOTON_CANDIDATE_MERGE_IMPLEMENTATION

#include "CandidateMerge.h"
#include "EtaTopologyPlots.h"
#include "OtherOriginPlots.h"
#include "DecayOriginPlots.h"
#include "CommonCompositionPlots.h"

namespace photon_candidate_analysis
{
int merge(
    const std::string family,
    const std::string partial_root,
    const std::string output_base)
{
  using namespace candidate_composition;
  using namespace candidate_composition_merge;
  const std::vector<SampleDefinition> samples = sample_definitions(family);
  if (samples.empty() || partial_root.empty() || output_base.empty()) return 1;
  std::string normalized_partial_root = partial_root;
  while (normalized_partial_root.size() > 1U && normalized_partial_root.back() == '/') normalized_partial_root.pop_back();
  if (gSystem->AccessPathName(output_base.c_str()) && gSystem->mkdir(output_base.c_str(), true) != 0) return 2;

  std::array<std::unique_ptr<Histograms>, kSelectionCount> composition_histograms;
  std::array<Summary, kSelectionCount> composition_summaries;
  std::array<std::unique_ptr<Spectra>, kSelectionCount> topology_histograms;
  std::array<pi0_bdt::Histograms, kSelectionCount> bdt_histograms;
  std::array<std::unique_ptr<Spectra>, kSelectionCount> eta_histograms;
  std::array<pi0_bdt::Histograms, kSelectionCount> eta_bdt_histograms;
  std::vector<unsigned long long> sample_events_written, sample_events_stitch_pass, sample_region_a_clusters;
  std::vector<unsigned long long> sample_region_a_prompt_clusters, sample_region_a_anchor_clusters;
  std::array<std::vector<unsigned long long>, kSelectionCount> sample_selected_clusters, sample_selected_anchor_clusters;
  std::vector<double> sample_cross_sections_pb;
  Metadata combined;
  eta_decay::Selections eta_decay_histograms;
  decay_origin::Selections diagnostic_histograms;
  other_origin::Selections other_histograms;
  common_composition::Selections common_histograms;
  std::array<std::unique_ptr<Histograms>, kSelectionCount> legacy_histograms;
  bool have_reference = false;
  for (const SampleDefinition& sample : samples)
  {
    const int shard_count = required_shard_count(sample.name);
    unsigned long long next_map_index = 0, total_map_count = 0;
    double sample_sumw = 0.0;
    unsigned long long sample_event_count = 0, sample_stitch_count = 0, sample_region_a_count = 0;
    unsigned long long sample_region_a_prompt_count = 0, sample_region_a_anchor_count = 0;
    std::array<unsigned long long, kSelectionCount> sample_selected = {}, sample_selected_anchor = {};
    for (int shard_index = 0; shard_index < shard_count; ++shard_index)
    {
      const std::string path = normalized_partial_root + "/" + sample.name + "/shard_" + std::to_string(shard_index) + "/photon_candidate_selection.root";
      TFile file(path.c_str(), "READ");
      Metadata metadata;
      ShardMetadata shard;
      CompositionSummary composition_summary;
      TopologySummary topology_summary;
      if (file.IsZombie() || !read_metadata(file, metadata) || !read_shard_metadata(file, shard) || !read_composition_summary(file, composition_summary) || !read_topology_summary(file, topology_summary))
      {
        std::cerr << "Missing or invalid composition partial: " << path << std::endl;
        return 3;
      }
      if (metadata.schema_version != 11 || metadata.source_map_schema_version != 9 || metadata.family != family ||
          metadata.sample_names.size() != 1U || metadata.sample_names.front() != sample.name ||
          metadata.sample_map_counts.size() != 1U || metadata.sample_sum_generator_weights.size() != 1U || !same_double(metadata.sample_cross_section_pb, sample.cross_section_pb) ||
          metadata.events_stitch_pass > metadata.events_written || metadata.region_a_prompt_clusters > metadata.region_a_clusters ||
          metadata.region_a_anchor_clusters > metadata.region_a_clusters ||
          shard.sample_name != sample.name || shard.shard_index != shard_index || shard.shard_count != shard_count ||
          shard.total_map_count != metadata.sample_map_counts.front() || shard.total_map_count == 0 ||
          shard.map_index_begin != shard.total_map_count * static_cast<unsigned long long>(shard_index) / static_cast<unsigned long long>(shard_count) ||
          shard.map_index_end != shard.total_map_count * static_cast<unsigned long long>(shard_index + 1) / static_cast<unsigned long long>(shard_count) ||
          shard.map_index_begin != next_map_index || shard.map_index_begin >= shard.map_index_end)
      {
        std::cerr << "Unexpected partial metadata: " << path << std::endl;
        return 3;
      }
      if (!have_reference)
      {
        combined = metadata;
        set_missing_labels(combined.missing_energy_min, combined.missing_energy_max);
        combined.sample_names.clear();
        combined.sample_map_counts.clear();
        combined.sample_sum_generator_weights.clear();
        combined.events_written = 0;
        combined.events_stitch_pass = 0;
        combined.region_a_clusters = 0;
        combined.region_a_prompt_clusters = 0;
        combined.region_a_anchor_clusters = 0;
        have_reference = true;
      }
      else if (!compatible(metadata, combined))
      {
        std::cerr << "Incompatible composition partial metadata: " << path << std::endl;
        return 4;
      }
      if (shard_index == 0)
      {
        total_map_count = shard.total_map_count;
        sample_sumw = metadata.sample_sum_generator_weights.front();
        if (!std::isfinite(sample_sumw) || sample_sumw <= 0.0)
        {
          std::cerr << "Invalid full-sample normalization metadata: " << path << std::endl;
          return 4;
        }
      }
      else if (shard.total_map_count != total_map_count || !same_double(metadata.sample_sum_generator_weights.front(), sample_sumw))
      {
        std::cerr << "Inconsistent shard normalization metadata: " << path << std::endl;
        return 4;
      }

      for (std::size_t composition_selection = 0; composition_selection < kSelectionCount; ++composition_selection)
      {
        if (!common_histograms[composition_selection]) common_histograms[composition_selection] = std::make_unique<common_composition::Histograms>(metadata.n_bins, metadata.et_max);
        auto* common_directory = file.GetDirectory((std::string("common_composition/") + kSelectionKeys[composition_selection]).c_str());
        if (!common_directory || !common_histograms[composition_selection]->add(*common_directory)) return 5;
        if (!composition_histograms[composition_selection])
        {
          composition_histograms[composition_selection] = std::make_unique<Histograms>(metadata.n_bins, metadata.et_max);
          for (std::size_t index = 0; index < category_count; ++index)
          {
            composition_histograms[composition_selection]->counts[index]->Reset("ICES");
            composition_histograms[composition_selection]->weighted[index]->Reset("ICES");
          }
        }
        for (std::size_t index = 0; index < category_count; ++index)
        {
          const std::string prefix = std::string("composition/") + kSelectionKeys[composition_selection] + "/h_candidate_" + candidate_composition::kKeys[index] + "_et_";
          TH1D* count = file.Get<TH1D>((prefix + "count").c_str());
          TH1D* weighted_pb = file.Get<TH1D>((prefix + "pb").c_str());
          if (!valid_histogram(count, metadata) || !valid_histogram(weighted_pb, metadata) ||
              !composition_histograms[composition_selection]->counts[index]->Add(count) ||
              !composition_histograms[composition_selection]->weighted[index]->Add(weighted_pb))
          {
            std::cerr << "Missing or incompatible candidate-composition histogram in " << path << std::endl;
            return 5;
          }
        }
        auto* bdt_directory = file.GetDirectory((std::string("pi0_bdt/") + kSelectionKeys[composition_selection]).c_str());
        if (!bdt_directory || !bdt_histograms[composition_selection].add(*bdt_directory)) return 5;
        auto* eta_bdt_directory = file.GetDirectory((std::string("eta_bdt/") + kSelectionKeys[composition_selection]).c_str());
        if (!eta_bdt_directory || !eta_bdt_histograms[composition_selection].add(*eta_bdt_directory)) return 5;
        if (!eta_histograms[composition_selection]) eta_histograms[composition_selection] = std::make_unique<Spectra>(metadata.n_bins, metadata.et_max);
        for (std::size_t index = 0; index < kSpectrumCount; ++index)
        {
          const std::string prefix = std::string("eta_anchor_topology/") + kSelectionKeys[composition_selection] + "/h_region_a_" + eta_topology::spectrum_key(index) + "_et_";
          auto* count = file.Get<TH1D>((prefix + "count").c_str());
          auto* weighted = file.Get<TH1D>((prefix + "pb").c_str());
          if (!valid_histogram(count, metadata) || !valid_histogram(weighted, metadata) ||
              !eta_histograms[composition_selection]->counts[index]->Add(count) || !eta_histograms[composition_selection]->weighted_pb[index]->Add(weighted)) return 5;
        }
        if (!eta_decay_histograms[composition_selection]) eta_decay_histograms[composition_selection] = std::make_unique<eta_decay::Histograms>(metadata.n_bins, metadata.et_max);
        auto* decay_directory = file.GetDirectory((std::string("eta_decay/") + kSelectionKeys[composition_selection]).c_str());
        if (!decay_directory || !eta_decay_histograms[composition_selection]->add(*decay_directory)) return 5;
        if (!other_histograms[composition_selection]) other_histograms[composition_selection] = std::make_unique<other_origin::Histograms>(metadata.n_bins, metadata.et_max);
        auto* other_directory = file.GetDirectory((std::string("other_origin/") + kSelectionKeys[composition_selection]).c_str());
        if (!other_directory || !other_histograms[composition_selection]->add(*other_directory)) return 5;
        if (!diagnostic_histograms[composition_selection]) diagnostic_histograms[composition_selection] = std::make_unique<decay_origin::Histograms>(metadata.n_bins, metadata.et_max);
        auto* diagnostic_directory = file.GetDirectory((std::string("decay_origin/") + kSelectionKeys[composition_selection]).c_str());
        if (!diagnostic_directory || !diagnostic_histograms[composition_selection]->add(*diagnostic_directory)) return 5;
        add_summary(composition_summaries[composition_selection], composition_summary.selections[composition_selection]);
      }
      for (std::size_t topology_selection = 0; topology_selection < kSelectionCount; ++topology_selection)
      {
        if (!topology_histograms[topology_selection])
        {
          topology_histograms[topology_selection] = std::make_unique<Spectra>(metadata.n_bins, metadata.et_max);
          for (std::size_t index = 0; index < kSpectrumCount; ++index)
          {
            topology_histograms[topology_selection]->counts[index]->Reset("ICES");
            topology_histograms[topology_selection]->weighted_pb[index]->Reset("ICES");
            topology_histograms[topology_selection]->counts[index]->SetDirectory(nullptr);
            topology_histograms[topology_selection]->weighted_pb[index]->SetDirectory(nullptr);
          }
        }
        for (std::size_t index = 0; index < kSpectrumCount; ++index)
        {
          const std::string prefix = std::string("anchor_topology/") + kSelectionKeys[topology_selection] + "/h_region_a_" + ::kKeys[index] + "_et_";
          TH1D* count = file.Get<TH1D>((prefix + "count").c_str());
          TH1D* weighted_pb = file.Get<TH1D>((prefix + "pb").c_str());
          if (!valid_histogram(count, metadata) || !valid_histogram(weighted_pb, metadata) ||
              !topology_histograms[topology_selection]->counts[index]->Add(count) ||
              !topology_histograms[topology_selection]->weighted_pb[index]->Add(weighted_pb))
          {
            std::cerr << "Missing or incompatible anchor-topology histogram in " << path << std::endl;
            return 5;
          }
        }
        sample_selected[topology_selection] += topology_summary.selected_clusters[topology_selection];
        sample_selected_anchor[topology_selection] += topology_summary.selected_anchor_clusters[topology_selection];
      }
      sample_event_count += metadata.events_written;
      sample_stitch_count += metadata.events_stitch_pass;
      sample_region_a_count += metadata.region_a_clusters;
      sample_region_a_prompt_count += metadata.region_a_prompt_clusters;
      sample_region_a_anchor_count += metadata.region_a_anchor_clusters;
      add_metadata_counters(combined, metadata);
      next_map_index = shard.map_index_end;
    }
    if (next_map_index != total_map_count) return 6;
    combined.sample_names.push_back(sample.name);
    combined.sample_map_counts.push_back(total_map_count);
    combined.sample_sum_generator_weights.push_back(sample_sumw);
    sample_cross_sections_pb.push_back(sample.cross_section_pb);
    sample_events_written.push_back(sample_event_count);
    sample_events_stitch_pass.push_back(sample_stitch_count);
    sample_region_a_clusters.push_back(sample_region_a_count);
    sample_region_a_prompt_clusters.push_back(sample_region_a_prompt_count);
    sample_region_a_anchor_clusters.push_back(sample_region_a_anchor_count);
    for (std::size_t topology_selection = 0; topology_selection < kSelectionCount; ++topology_selection)
    {
      sample_selected_clusters[topology_selection].push_back(sample_selected[topology_selection]);
      sample_selected_anchor_clusters[topology_selection].push_back(sample_selected_anchor[topology_selection]);
    }
  }

  if (!have_reference) return 6;
  plot_caption = {family, "", combined.min_cluster_energy, combined.pi0_partner_min_energy, combined.eta_partner_min_energy,
                  combined.pi0_mass_min, combined.pi0_mass_max, combined.eta_mass_min, combined.eta_mass_max};
  for (std::size_t selection_index = 0; selection_index < kSelectionCount; ++selection_index)
  {
    const Summary& summary = composition_summaries[selection_index];
    if (!composition_histograms[selection_index] || summary.selected != summary.prompt + summary.pi0 + summary.eta + summary.other ||
        !candidate_composition::valid_partition(composition_histograms[selection_index]->counts) ||
        !candidate_composition::valid_partition(composition_histograms[selection_index]->weighted) ||
        !topology_histograms[selection_index] || !::valid_partition(topology_histograms[selection_index]->counts) ||
        !::valid_partition(topology_histograms[selection_index]->weighted_pb)) return 6;
    for (std::size_t index = 0; index < category_count; ++index)
    {
      if (!valid_subset(*composition_histograms[selection_index]->counts[index], *composition_histograms[0]->counts[index],
                        *composition_histograms[selection_index]->weighted[index], *composition_histograms[0]->weighted[index])) return 6;
      const auto& previous = *composition_histograms[kPreviousCutSelections[selection_index]];
      if (!valid_subset(*composition_histograms[selection_index]->counts[index], *previous.counts[index],
                        *composition_histograms[selection_index]->weighted[index], *previous.weighted[index])) return 6;
    }
    for (std::size_t index = 0; index < kSpectrumCount; ++index)
    {
      if (!valid_subset(*topology_histograms[selection_index]->counts[index], *topology_histograms[0]->counts[index],
                        *topology_histograms[selection_index]->weighted_pb[index], *topology_histograms[0]->weighted_pb[index])) return 6;
    }
  }

  for (std::size_t i = 0; i < kSelectionCount; ++i)
  {
    if (!consistent_pi0_bdt(*composition_histograms[i], *topology_histograms[i], bdt_histograms[i]) ||
        !eta_topology::valid(*composition_histograms[i], *eta_histograms[i], eta_bdt_histograms[i])) return 6;
    for (std::size_t index = 0; index < kSpectrumCount; ++index)
      if (!valid_subset(*eta_histograms[i]->counts[index], *eta_histograms[0]->counts[index],
                        *eta_histograms[i]->weighted_pb[index], *eta_histograms[0]->weighted_pb[index])) return 6;
    const std::string directory = output_base + "/" + kSelectionKeys[i];
    if (!make_output_directory(directory + "/pi0_bdt_score")) return 7;
    plot_caption.selection = kSelectionLabels[i];
    if (!bdt_histograms[i].export_files(directory, kSelectionLabels[i], i == 3 || i >= 4, caption_lines(plot_caption))) return 7;
    TFile bdt_output((directory + "/pi0_bdt_score.root").c_str(), "UPDATE");
    if (bdt_output.IsZombie()) return 7;
    TNamed model("model_sha256", combined.model_sha256.c_str()); model.Write();
    TNamed configuration("configuration", output_base.c_str()); configuration.Write();
    TNamed selection_definition("selection_definition", kSelectionDefinitions[i]); selection_definition.Write();
    TTree bdt_metadata("metadata", "Pi0 BDT source and selection metadata");
    int bdt_schema_version = 1;
    std::string bdt_selection_key = kSelectionKeys[i];
    bdt_metadata.Branch("schema_version", &bdt_schema_version);
    bdt_metadata.Branch("source_partial_schema_version", &combined.schema_version);
    bdt_metadata.Branch("source_map_schema_version", &combined.source_map_schema_version);
    bdt_metadata.Branch("common_composition_algorithm_version", &combined.common_composition_algorithm_version);
    bdt_metadata.Branch("selection_key", &bdt_selection_key);
    bdt_metadata.Branch("family", &combined.family);
    bdt_metadata.Branch("map_root", &combined.map_root);
    bdt_metadata.Branch("require_complete", &combined.require_complete);
    bdt_metadata.Branch("model_sha256", &combined.model_sha256);
    bdt_metadata.Branch("min_cluster_energy", &combined.min_cluster_energy);
    bdt_metadata.Branch("pi0_partner_min_energy", &combined.pi0_partner_min_energy);
    bdt_metadata.Branch("eta_partner_min_energy", &combined.eta_partner_min_energy);
    bdt_metadata.Branch("pi0_mass_min", &combined.pi0_mass_min);
    bdt_metadata.Branch("pi0_mass_max", &combined.pi0_mass_max);
    bdt_metadata.Branch("eta_mass_min", &combined.eta_mass_min);
    bdt_metadata.Branch("eta_mass_max", &combined.eta_mass_max);
    bdt_metadata.Branch("pi0_topology_algorithm_version", &combined.pi0_topology_algorithm_version);
    bdt_metadata.Branch("sample_names", &combined.sample_names);
    bdt_metadata.Branch("sample_map_counts", &combined.sample_map_counts);
    bdt_metadata.Branch("sample_sum_generator_weights", &combined.sample_sum_generator_weights);
    bdt_metadata.Branch("sample_cross_sections_pb", &sample_cross_sections_pb);
    bdt_metadata.Fill();
    bdt_metadata.Write();
  }

  TFile composition_output((output_base + "/candidate_composition.root").c_str(), "RECREATE");
  if (composition_output.IsZombie()) return 7;
  for (std::size_t composition_selection = 0; composition_selection < kSelectionCount; ++composition_selection)
  {
    plot_caption.selection = kSelectionLabels[composition_selection];
    Histograms& histograms = *composition_histograms[composition_selection];
    std::array<std::unique_ptr<TH1D>, category_count> fractions;
    for (std::size_t index = 1; index < category_count; ++index)
    {
      const std::string name = index == prompt ? "h_photon_candidate_purity" : std::string("h_candidate_") + candidate_composition::kKeys[index] + "_fraction";
      fractions[index] = fraction_histogram(*histograms.weighted[index], *histograms.weighted[denominator], name);
    }
    const std::string selection_output_base = output_base + "/" + kSelectionKeys[composition_selection] + "/photon_candidate_composition";
    if (!make_output_directory(selection_output_base)) return 7;
    auto origin_density = make_candidate_origin_density(histograms, std::string("h_") + kSelectionKeys[composition_selection] + "_candidate_origin_");
    draw_candidate_origin_spectrum(origin_density, selection_output_base + "_truth_origin.pdf", kSelectionKeys[composition_selection]);
    draw_stack(fractions, selection_output_base + "_category_fraction_stack.pdf", false);
    draw_stack(fractions, selection_output_base + "_category_fraction_stack_detailed.pdf", true);
    auto missing_fractions = composition_missing_fractions(histograms, *topology_histograms[composition_selection]);
    if (missing_fractions.size() != kMissingSpectrumIndices.size())
    {
      std::cerr << "Cannot split composition missing: its population differs from anchor topology for " << kSelectionKeys[composition_selection] << std::endl;
      return 7;
    }
    draw_superdetailed_stack(fractions, missing_fractions, selection_output_base + "_category_fraction_stack_superdetailed.pdf");
    std::array<std::unique_ptr<SurvivalCurve>, category_count> survival;
    for (std::size_t index = 0; index < category_count; ++index)
    {
      const std::string key = index == denominator ? "all" : candidate_composition::kKeys[index];
      const std::string label = index == denominator ? "All candidates" : candidate_composition::kLabels[index];
      const int color = index == denominator ? kBlack : candidate_composition::kColors[index];
      survival[index] = make_survival_curve(*histograms.weighted[index], *composition_histograms[0]->weighted[index],
          "h_candidate_" + key + "_survival_fraction_relative_to_kinematic", label, color);
    }
    std::vector<const TH1D*> selected_pi0_inputs, kinematic_pi0_inputs;
    for (std::size_t index = pi0_separated; index <= pi0_unknown; ++index)
    {
      selected_pi0_inputs.push_back(histograms.weighted[index].get());
      kinematic_pi0_inputs.push_back(composition_histograms[0]->weighted[index].get());
    }
    auto selected_pi0 = sum_histograms(selected_pi0_inputs, "h_selected_pi0_for_survival");
    auto kinematic_pi0 = sum_histograms(kinematic_pi0_inputs, "h_kinematic_pi0_for_survival");
    if (!selected_pi0 || !kinematic_pi0) return 7;
    auto pi0_survival = make_survival_curve(*selected_pi0, *kinematic_pi0,
        "h_candidate_pi0_survival_fraction_relative_to_kinematic", "#pi^{0}", candidate_composition::kColors[pi0_separated]);
    const std::vector<const SurvivalCurve*> summary_survival = {
        survival[prompt].get(), pi0_survival.get(), survival[eta].get(), survival[other].get()};
    const std::vector<const SurvivalCurve*> detailed_survival = {
        survival[prompt].get(), survival[pi0_separated].get(), survival[pi0_merged].get(),
        survival[pi0_single_contaminated].get(), survival[pi0_missing].get(), survival[pi0_other].get(), survival[pi0_not_applicable].get(), survival[pi0_unknown].get(), survival[eta].get(), survival[other].get()};
    draw_survival_curves(summary_survival, *histograms.weighted[denominator],
        selection_output_base + "_survival_fraction_relative_to_kinematic.pdf",
        std::string("c_") + kSelectionKeys[composition_selection] + "_candidate_survival", false);
    draw_survival_curves(detailed_survival, *histograms.weighted[denominator],
        selection_output_base + "_survival_fraction_relative_to_kinematic_detailed.pdf",
        std::string("c_") + kSelectionKeys[composition_selection] + "_candidate_survival_detailed", true);
    const auto& previous = *composition_histograms[kPreviousCutSelections[composition_selection]];
    const std::string previous_label = kPreviousCutLabels[composition_selection];
    std::array<std::unique_ptr<SurvivalCurve>, category_count> previous_survival;
    for (std::size_t index = 0; index < category_count; ++index)
    {
      const std::string key = index == denominator ? "all" : candidate_composition::kKeys[index];
      const std::string label = index == denominator ? "All candidates" : candidate_composition::kLabels[index];
      const int color = index == denominator ? kBlack : candidate_composition::kColors[index];
      previous_survival[index] = make_survival_curve(*histograms.weighted[index], *previous.weighted[index],
          "h_candidate_" + key + "_survival_fraction_relative_to_previous_cut", label, color, previous_label);
    }
    std::vector<const TH1D*> previous_pi0_inputs;
    for (std::size_t index = pi0_separated; index <= pi0_unknown; ++index) previous_pi0_inputs.push_back(previous.weighted[index].get());
    auto previous_pi0 = sum_histograms(previous_pi0_inputs, "h_previous_pi0_for_survival");
    if (!previous_pi0) return 7;
    auto previous_pi0_survival = make_survival_curve(*selected_pi0, *previous_pi0,
        "h_candidate_pi0_survival_fraction_relative_to_previous_cut", "#pi^{0}", candidate_composition::kColors[pi0_separated], previous_label);
    const std::vector<const SurvivalCurve*> previous_summary_survival = {
        previous_survival[prompt].get(), previous_pi0_survival.get(), previous_survival[eta].get(), previous_survival[other].get()};
    draw_survival_curves(previous_summary_survival, *histograms.weighted[denominator],
        selection_output_base + "_survival_fraction_relative_to_previous_cut.pdf",
        std::string("c_") + kSelectionKeys[composition_selection] + "_candidate_survival_previous_cut", false, previous_label);
    TDirectory* directory = composition_output.mkdir(kSelectionKeys[composition_selection]);
    if (!directory) return 7;
    directory->cd();
    for (std::size_t index = 0; index < category_count; ++index)
    {
      histograms.counts[index]->Write();
      histograms.weighted[index]->Write();
      if (index > 0) fractions[index]->Write();
      survival[index]->histogram->Write();
      previous_survival[index]->histogram->Write();
    }
    for (auto& fraction : missing_fractions) fraction->Write();
    pi0_survival->histogram->Write();
    previous_pi0_survival->histogram->Write();
    composition_output.cd();
  }
  std::vector<std::string> composition_selection_keys(kSelectionKeys.begin(), kSelectionKeys.end());
  std::vector<std::string> composition_selection_labels(kSelectionLabels.begin(), kSelectionLabels.end());
  std::vector<std::string> composition_selection_definitions(kSelectionDefinitions.begin(), kSelectionDefinitions.end());
  std::vector<std::string> previous_cut_denominator_selections;
  for (std::size_t index : kPreviousCutSelections) previous_cut_denominator_selections.emplace_back(kSelectionKeys[index]);
  std::string previous_cut_survival_fraction_definition = "weighted_selected_category_yield_divided_by_weighted_previous_cut_category_yield_per_cluster_et_bin";
  int composition_schema_version = 8;
  int composition_source_partial_schema_version = combined.schema_version;
  std::string survival_fraction_denominator_selection = "kinematic";
  std::string survival_fraction_definition = "weighted_selected_category_yield_divided_by_weighted_kinematic_category_yield_per_cluster_et_bin";
  std::string survival_fraction_uncertainty = "weighted_subset_covariance";
  std::string survival_fraction_zero_denominator = "stored_as_zero_and_not_drawn";
  TTree metadata("metadata", "Photon-candidate composition merge metadata");
  metadata.Branch("schema_version", &composition_schema_version);
  metadata.Branch("source_partial_schema_version", &composition_source_partial_schema_version);
  metadata.Branch("source_map_schema_version", &combined.source_map_schema_version);
  metadata.Branch("common_composition_algorithm_version", &combined.common_composition_algorithm_version);
  metadata.Branch("family", &combined.family);
  metadata.Branch("map_root", &combined.map_root);
  metadata.Branch("require_complete", &combined.require_complete);
  metadata.Branch("n_bins", &combined.n_bins);
  metadata.Branch("et_max", &combined.et_max);
  metadata.Branch("selection_keys", &composition_selection_keys);
  metadata.Branch("selection_labels", &composition_selection_labels);
  metadata.Branch("selection_definitions", &composition_selection_definitions);
  metadata.Branch("min_cluster_energy", &combined.min_cluster_energy);
  metadata.Branch("partner_diagnostic_min_cluster_energy", &combined.partner_diagnostic_min_cluster_energy);
  metadata.Branch("pi0_partner_min_energy", &combined.pi0_partner_min_energy);
  metadata.Branch("eta_partner_min_energy", &combined.eta_partner_min_energy);
  metadata.Branch("pi0_mass_min", &combined.pi0_mass_min);
  metadata.Branch("pi0_mass_max", &combined.pi0_mass_max);
  metadata.Branch("eta_mass_min", &combined.eta_mass_min);
  metadata.Branch("eta_mass_max", &combined.eta_mass_max);
  metadata.Branch("missing_energy_min", &combined.missing_energy_min);
  metadata.Branch("missing_energy_max", &combined.missing_energy_max);
  metadata.Branch("min_photon_energy_recovery", &combined.min_photon_energy_recovery);
  metadata.Branch("pi0_topology_algorithm_version", &combined.pi0_topology_algorithm_version);
  metadata.Branch("signal_embedding_id", &combined.signal_embedding_id);
  metadata.Branch("majority_threshold", &combined.majority_threshold);
  metadata.Branch("majority_comparison", &combined.majority_comparison);
  metadata.Branch("eta_definition", &combined.eta_definition);
  metadata.Branch("other_definition", &combined.other_definition);
  metadata.Branch("weight_definition", &combined.weight_definition);
  metadata.Branch("survival_fraction_denominator_selection", &survival_fraction_denominator_selection);
  metadata.Branch("survival_fraction_definition", &survival_fraction_definition);
  metadata.Branch("previous_cut_denominator_selections", &previous_cut_denominator_selections);
  metadata.Branch("previous_cut_survival_fraction_definition", &previous_cut_survival_fraction_definition);
  metadata.Branch("survival_fraction_uncertainty", &survival_fraction_uncertainty);
  metadata.Branch("survival_fraction_zero_denominator", &survival_fraction_zero_denominator);
  metadata.Branch("analysis_release", &combined.analysis_release);
  metadata.Branch("model_sha256", &combined.model_sha256);
  metadata.Branch("sample_names", &combined.sample_names);
  metadata.Branch("sample_map_counts", &combined.sample_map_counts);
  metadata.Branch("sample_sum_generator_weights", &combined.sample_sum_generator_weights);
  metadata.Fill();
  metadata.Write();
  TTree composition_summary_output("selection_summary", "Per-selection candidate-composition merge summary");
  unsigned int composition_selection_index = 0;
  std::string composition_selection_key, composition_selection_label;
  unsigned long long composition_selected = 0, composition_prompt = 0, composition_pi0 = 0, composition_eta = 0, composition_other = 0;
  unsigned long long composition_overlap = 0, composition_half_boundary = 0, composition_invalid_truth = 0;
  composition_summary_output.Branch("selection_index", &composition_selection_index);
  composition_summary_output.Branch("selection_key", &composition_selection_key);
  composition_summary_output.Branch("selection_label", &composition_selection_label);
  composition_summary_output.Branch("selected_clusters", &composition_selected);
  composition_summary_output.Branch("prompt_clusters", &composition_prompt);
  composition_summary_output.Branch("pi0_clusters", &composition_pi0);
  composition_summary_output.Branch("eta_clusters", &composition_eta);
  composition_summary_output.Branch("other_clusters", &composition_other);
  composition_summary_output.Branch("overlap_clusters", &composition_overlap);
  composition_summary_output.Branch("half_boundary_clusters", &composition_half_boundary);
  composition_summary_output.Branch("invalid_truth_clusters", &composition_invalid_truth);
  for (composition_selection_index = 0; composition_selection_index < kSelectionCount; ++composition_selection_index)
  {
    composition_selection_key = kSelectionKeys[composition_selection_index];
    composition_selection_label = kSelectionLabels[composition_selection_index];
    const Summary& summary = composition_summaries[composition_selection_index];
    composition_selected = summary.selected;
    composition_prompt = summary.prompt;
    composition_pi0 = summary.pi0;
    composition_eta = summary.eta;
    composition_other = summary.other;
    composition_overlap = summary.overlap;
    composition_half_boundary = summary.half_boundary;
    composition_invalid_truth = summary.invalid_truth;
    composition_summary_output.Fill();
  }
  composition_summary_output.Write();
  composition_output.Close();

  if (gSystem->AccessPathName(output_base.c_str()) && gSystem->mkdir(output_base.c_str(), true) != 0) return 7;
  TFile topology_output((output_base + "/pi0_anchor_topology.root").c_str(), "RECREATE");
  if (topology_output.IsZombie()) return 7;
  for (std::size_t topology_selection = 0; topology_selection < kSelectionCount; ++topology_selection)
  {
    plot_caption.selection = kSelectionLabels[topology_selection];
    Spectra& spectra = *topology_histograms[topology_selection];
    const std::string selection_output_base = output_base + "/" + kSelectionKeys[topology_selection] + "/region_a_pi0_anchor_topology";
    if (!make_output_directory(selection_output_base)) return 7;
    std::array<std::unique_ptr<TH1D>, kSpectrumCount> density;
    for (std::size_t index = 0; index < kSpectrumCount; ++index)
    {
      density[index].reset(static_cast<TH1D*>(spectra.weighted_pb[index]->Clone((std::string("h_") + kSelectionKeys[topology_selection] + "_" + ::kKeys[index] + "_et_pb_per_gev").c_str())));
      density[index]->SetDirectory(nullptr);
      density[index]->Scale(1.0, "width");
      density[index]->SetLineColor(::kColors[index]);
      style_axes(density[index].get(), "Weighted Counts [a.u.]");
      spectra.counts[index]->SetLineColor(::kColors[index]);
      spectra.weighted_pb[index]->SetLineColor(::kColors[index]);
    }
    std::vector<std::size_t> detailed_spectrum_indices(kSpectrumCount);
    std::iota(detailed_spectrum_indices.begin(), detailed_spectrum_indices.end(), std::size_t{0});
    const std::vector<std::size_t> summary_spectrum_indices(kSummarySpectrum.begin(), kSummarySpectrum.end());
    const std::vector<std::size_t> detailed_category_indices(kDetailedCategories.begin(), kDetailedCategories.end());
    const std::vector<std::size_t> summary_category_indices(kSummaryCategories.begin(), kSummaryCategories.end());
    auto detailed_fractions = make_fractions(spectra.weighted_pb, detailed_category_indices, std::string("h_") + kSelectionKeys[topology_selection] + "_detailed_");
    auto summary_fractions = make_fractions(spectra.weighted_pb, summary_category_indices, std::string("h_") + kSelectionKeys[topology_selection] + "_summary_");
    draw_spectrum(density, detailed_spectrum_indices, selection_output_base + "_detailed.pdf", true);
    draw_spectrum(density, summary_spectrum_indices, selection_output_base + ".pdf", false);
    draw_fraction_lines(detailed_fractions, detailed_category_indices, selection_output_base + "_category_fractions_detailed.pdf", true);
    draw_fraction_lines(summary_fractions, summary_category_indices, selection_output_base + "_category_fractions.pdf", false);
    draw_fraction_stack(detailed_fractions, detailed_category_indices, selection_output_base + "_category_fraction_stack_detailed.pdf", true);
    draw_fraction_stack(summary_fractions, summary_category_indices, selection_output_base + "_category_fraction_stack.pdf", false);
    std::array<std::unique_ptr<SurvivalCurve>, kSpectrumCount> topology_survival;
    for (std::size_t index = 1; index < kSpectrumCount; ++index)
    {
      const int color = index == 1 ? kBlack : ::kColors[index];
      const std::string label = index == 1 ? "All #pi^{0}-main anchors" : ::kLabels[index];
      topology_survival[index] = make_survival_curve(*spectra.weighted_pb[index], *topology_histograms[0]->weighted_pb[index],
          std::string("h_") + kSelectionKeys[topology_selection] + "_" + ::kKeys[index] + "_survival_fraction_relative_to_kinematic",
          label, color);
    }
    std::vector<const SurvivalCurve*> topology_summary_survival = {topology_survival[1].get()};
    for (std::size_t index : summary_category_indices) topology_summary_survival.push_back(topology_survival[index].get());
    std::vector<const SurvivalCurve*> topology_detailed_survival = {topology_survival[1].get()};
    for (std::size_t index : detailed_category_indices) topology_detailed_survival.push_back(topology_survival[index].get());
    draw_survival_curves(topology_summary_survival, *spectra.weighted_pb[1],
        selection_output_base + "_survival_fraction_relative_to_kinematic.pdf",
        std::string("c_") + kSelectionKeys[topology_selection] + "_pi0_anchor_survival", false);
    draw_survival_curves(topology_detailed_survival, *spectra.weighted_pb[1],
        selection_output_base + "_survival_fraction_relative_to_kinematic_detailed.pdf",
        std::string("c_") + kSelectionKeys[topology_selection] + "_pi0_anchor_survival_detailed", true);

    TDirectory* directory = topology_output.mkdir(kSelectionKeys[topology_selection]);
    if (!directory) return 7;
    directory->cd();
    for (std::size_t index = 0; index < kSpectrumCount; ++index)
    {
      spectra.counts[index]->Write();
      spectra.weighted_pb[index]->Write();
      density[index]->Write();
      if (index > 0) topology_survival[index]->histogram->Write();
    }
    for (auto& histogram : detailed_fractions) histogram->Write();
    for (auto& histogram : summary_fractions) histogram->Write();
    topology_output.cd();
  }

  int topology_schema_version = 6;
  int topology_source_schema_version = combined.source_map_schema_version;
  int topology_source_partial_schema_version = combined.schema_version;
  std::vector<std::string> topology_selection_keys(kSelectionKeys.begin(), kSelectionKeys.end());
  std::vector<std::string> topology_selection_labels(kSelectionLabels.begin(), kSelectionLabels.end());
  std::vector<std::string> topology_selection_definitions(kSelectionDefinitions.begin(), kSelectionDefinitions.end());
  std::string topology_prompt_definition = "selected_and_truth_prompt_cluster";
  std::string topology_anchor_definition = "composition_pi0_valid_anchor_strict_majority_no_origin_overlap";
  std::string topology_source = "stored_map_topology_evaluated_with_strict_min_cluster_energy_from_metadata";
  unsigned int topology_sample_count = combined.sample_names.size();
  TTree topology_metadata("metadata", "Pi0-anchor topology merge metadata");
  topology_metadata.Branch("anchor_definition", &topology_anchor_definition);
  topology_metadata.Branch("schema_version", &topology_schema_version);
  topology_metadata.Branch("source_partial_schema_version", &topology_source_partial_schema_version);
  topology_metadata.Branch("source_map_schema_version", &topology_source_schema_version);
  topology_metadata.Branch("family", &combined.family);
  topology_metadata.Branch("map_root", &combined.map_root);
  topology_metadata.Branch("require_complete", &combined.require_complete);
  topology_metadata.Branch("n_bins", &combined.n_bins);
  topology_metadata.Branch("et_max", &combined.et_max);
  topology_metadata.Branch("selection_keys", &topology_selection_keys);
  topology_metadata.Branch("selection_labels", &topology_selection_labels);
  topology_metadata.Branch("selection_definitions", &topology_selection_definitions);
  topology_metadata.Branch("prompt_definition", &topology_prompt_definition);
  topology_metadata.Branch("topology_source", &topology_source);
  topology_metadata.Branch("analysis_release", &combined.analysis_release);
  topology_metadata.Branch("model_sha256", &combined.model_sha256);
  topology_metadata.Branch("weight_definition", &combined.weight_definition);
  topology_metadata.Branch("survival_fraction_denominator_selection", &survival_fraction_denominator_selection);
  topology_metadata.Branch("survival_fraction_definition", &survival_fraction_definition);
  topology_metadata.Branch("survival_fraction_uncertainty", &survival_fraction_uncertainty);
  topology_metadata.Branch("survival_fraction_zero_denominator", &survival_fraction_zero_denominator);
  topology_metadata.Branch("sample_count", &topology_sample_count);
  topology_metadata.Branch("sample_names", &combined.sample_names);
  topology_metadata.Branch("min_cluster_energy", &combined.min_cluster_energy);
  topology_metadata.Branch("partner_diagnostic_min_cluster_energy", &combined.partner_diagnostic_min_cluster_energy);
  topology_metadata.Branch("pi0_partner_min_energy", &combined.pi0_partner_min_energy);
  topology_metadata.Branch("eta_partner_min_energy", &combined.eta_partner_min_energy);
  topology_metadata.Branch("pi0_mass_min", &combined.pi0_mass_min);
  topology_metadata.Branch("pi0_mass_max", &combined.pi0_mass_max);
  topology_metadata.Branch("eta_mass_min", &combined.eta_mass_min);
  topology_metadata.Branch("eta_mass_max", &combined.eta_mass_max);
  topology_metadata.Branch("missing_energy_min", &combined.missing_energy_min);
  topology_metadata.Branch("missing_energy_max", &combined.missing_energy_max);
  topology_metadata.Branch("min_photon_energy_recovery", &combined.min_photon_energy_recovery);
  topology_metadata.Branch("pi0_topology_algorithm_version", &combined.pi0_topology_algorithm_version);
  topology_metadata.Branch("sample_map_counts", &combined.sample_map_counts);
  topology_metadata.Branch("sample_events_written", &sample_events_written);
  topology_metadata.Branch("sample_events_stitch_pass", &sample_events_stitch_pass);
  topology_metadata.Branch("sample_region_a_clusters", &sample_region_a_clusters);
  topology_metadata.Branch("sample_region_a_prompt_clusters", &sample_region_a_prompt_clusters);
  topology_metadata.Branch("sample_region_a_anchor_clusters", &sample_region_a_anchor_clusters);
  topology_metadata.Branch("sample_cross_sections_pb", &sample_cross_sections_pb);
  topology_metadata.Branch("sample_sum_generator_weights", &combined.sample_sum_generator_weights);
  topology_metadata.Fill();
  topology_metadata.Write();

  TTree selection_summary("selection_summary", "Per-selection merge summary");
  unsigned int topology_selection_index = 0;
  std::string topology_selection_key, topology_selection_label;
  std::vector<unsigned long long> selected_clusters, selected_anchor_clusters;
  selection_summary.Branch("selection_index", &topology_selection_index);
  selection_summary.Branch("selection_key", &topology_selection_key);
  selection_summary.Branch("selection_label", &topology_selection_label);
  selection_summary.Branch("sample_selected_clusters", &selected_clusters);
  selection_summary.Branch("sample_selected_anchor_clusters", &selected_anchor_clusters);
  for (topology_selection_index = 0; topology_selection_index < kSelectionCount; ++topology_selection_index)
  {
    topology_selection_key = kSelectionKeys[topology_selection_index];
    topology_selection_label = kSelectionLabels[topology_selection_index];
    selected_clusters = sample_selected_clusters[topology_selection_index];
    selected_anchor_clusters = sample_selected_anchor_clusters[topology_selection_index];
    selection_summary.Fill();
  }
  selection_summary.Write();
  topology_output.Close();

  if (!eta_topology::export_outputs(output_base, combined, sample_cross_sections_pb, composition_histograms, eta_histograms, eta_bdt_histograms, &eta_decay_histograms)) return 7;

  for (std::size_t i = 0; i < kSelectionCount; ++i)
  {
    if (!common_histograms[i]->valid() || !common_composition::valid_populations(*common_histograms[i], *composition_histograms[i], *topology_histograms[i], *eta_histograms[i])) return 7;
    legacy_histograms[i] = std::make_unique<Histograms>(combined.n_bins, combined.et_max);
    for (int old = 0; old < 3; ++old)
    {
      const auto category = old == 0 ? prompt : old == 1 ? pi0_other : eta;
      legacy_histograms[i]->counts[category]->Add(common_histograms[i]->legacy_counts[old].get());
      legacy_histograms[i]->weighted[category]->Add(common_histograms[i]->legacy_weighted[old].get());
    }
    legacy_histograms[i]->counts[other]->Add(common_histograms[i]->legacy_counts[3].get());
    legacy_histograms[i]->weighted[other]->Add(common_histograms[i]->legacy_weighted[3].get());
    legacy_histograms[i]->counts[denominator]->Add(common_histograms[i]->counts[0].get());
    legacy_histograms[i]->weighted[denominator]->Add(common_histograms[i]->weighted[0].get());
  }
  if (!common_composition::export_outputs(output_base, combined, common_histograms)) return 7;
  if (!other_origin::export_outputs(output_base, combined, legacy_histograms, other_histograms)) return 7;
  if (!decay_origin::export_outputs(output_base, combined, legacy_histograms, diagnostic_histograms)) return 7;

  std::cout << "MergePythiaPhotonCandidateSelection - family/selections/samples/output = "
            << family << "/" << kSelectionCount << "/" << combined.sample_names.size() << "/" << output_base << std::endl;
  const bool have_overlap = std::any_of(composition_summaries.begin(), composition_summaries.end(), [](const Summary& summary) { return summary.overlap != 0ULL; });
  return have_overlap ? 9 : 0;
}
}

#endif
