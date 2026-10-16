#ifndef ETA_TOPOLOGY_PLOTS_H
#define ETA_TOPOLOGY_PLOTS_H

namespace eta_topology
{
using candidate_composition_merge::SurvivalCurve;
using candidate_composition_merge::make_survival_curve;
using candidate_composition_merge::draw_survival_curves;
using candidate_composition_merge::sum_histograms;

inline void write_metadata(TFile& output, const candidate_composition_merge::Metadata& source, const std::vector<double>& cross_sections)
{
  output.cd();
  // Local copies keep branch addresses valid until Fill completes.
  auto metadata = source;
  auto sample_cross_sections = cross_sections;
  int schema_version = 1, parent_pdg = 221;
  double majority_threshold = 0.5;
  std::string population = "composition_eta_subset_direct_eta_to_two_photons_same_parent_fraction_gt_0p5";
  std::string residual = "legacy_composition_eta_minus_selected_eta_two_photon_anchors";
  std::string denominator = "kinematic";
  std::string uncertainty = "weighted_subset_covariance";
  TTree tree("metadata", "Eta topology source and population metadata");
  tree.Branch("schema_version", &schema_version);
  tree.Branch("parent_pdg", &parent_pdg);
  tree.Branch("majority_threshold", &majority_threshold);
  tree.Branch("population", &population);
  tree.Branch("residual_eta_definition", &residual);
  tree.Branch("survival_fraction_denominator_selection", &denominator);
  tree.Branch("fraction_uncertainty", &uncertainty);
  tree.Branch("source_map_schema_version", &metadata.source_map_schema_version);
  tree.Branch("source_partial_schema_version", &metadata.schema_version);
#define ETA_METADATA(member) tree.Branch(#member, &metadata.member)
  ETA_METADATA(family); ETA_METADATA(map_root); ETA_METADATA(require_complete); ETA_METADATA(model_sha256); ETA_METADATA(analysis_release);
  ETA_METADATA(eta_topology_algorithm_version); ETA_METADATA(pi0_topology_algorithm_version);
  ETA_METADATA(min_cluster_energy); ETA_METADATA(eta_partner_min_energy); ETA_METADATA(pi0_partner_min_energy);
  ETA_METADATA(eta_mass_min); ETA_METADATA(eta_mass_max); ETA_METADATA(pi0_mass_min); ETA_METADATA(pi0_mass_max);
  ETA_METADATA(missing_energy_min); ETA_METADATA(missing_energy_max); ETA_METADATA(min_photon_energy_recovery);
  ETA_METADATA(sample_names); ETA_METADATA(sample_map_counts); ETA_METADATA(sample_sum_generator_weights);
#undef ETA_METADATA
  tree.Branch("sample_cross_sections_pb", &sample_cross_sections);
  tree.Fill(); tree.Write();
}

// Subtract a strict subset, including its Sumw2, rather than treating the two histograms as independent.
inline std::unique_ptr<TH1D> residual_histogram(const TH1D& total, const TH1D& subset, const std::string& name)
{
  auto result = std::unique_ptr<TH1D>(static_cast<TH1D*>(total.Clone(name.c_str())));
  result->SetDirectory(nullptr);
  for (int bin = 0; bin <= total.GetNbinsX() + 1; ++bin)
  {
    result->SetBinContent(bin, total.GetBinContent(bin) - subset.GetBinContent(bin));
    result->SetBinError(bin, std::sqrt(std::max(0.0, std::pow(total.GetBinError(bin), 2) - std::pow(subset.GetBinError(bin), 2))));
  }
  return result;
}

struct CompositionComponents
{
  std::vector<std::unique_ptr<TH1D>> counts, weighted;
  std::vector<std::string> labels;
  std::vector<int> colors;
};
inline CompositionComponents composition_components(const candidate_composition::Histograms& composition, const Spectra& topology, bool missing_detail)
{
  using namespace candidate_composition;
  CompositionComponents result;
  const auto append = [&](const TH1D& count, const TH1D& weighted, const std::string& key, const std::string& label, int color) {
    result.counts.emplace_back(static_cast<TH1D*>(count.Clone(("h_candidate_" + key + "_et_count").c_str())));
    result.weighted.emplace_back(static_cast<TH1D*>(weighted.Clone(("h_candidate_" + key + "_et_pb").c_str())));
    result.counts.back()->SetDirectory(nullptr); result.weighted.back()->SetDirectory(nullptr);
    result.labels.push_back(label); result.colors.push_back(color);
  };
  append(*composition.counts[prompt], *composition.weighted[prompt], "prompt", "Prompt-#gamma", candidate_composition::kColors[prompt]);
  std::vector<const TH1D*> pi0_counts, pi0_weights;
  for (std::size_t index = pi0_separated; index <= pi0_other; ++index)
  {
    pi0_counts.push_back(composition.counts[index].get()); pi0_weights.push_back(composition.weighted[index].get());
  }
  auto pi0_count = sum_histograms(pi0_counts, "h_pi0_count");
  auto pi0_weight = sum_histograms(pi0_weights, "h_pi0_pb");
  append(*pi0_count, *pi0_weight, "pi0", "#pi^{0}", candidate_composition::kColors[pi0_separated]);
  const std::vector<std::size_t> indices = missing_detail ? std::vector<std::size_t>(kDetailedCategories.begin(), kDetailedCategories.end())
                                                       : std::vector<std::size_t>(kSummaryCategories.begin(), kSummaryCategories.end());
  for (std::size_t index : indices)
    append(*topology.counts[index], *topology.weighted_pb[index], std::string("eta_") + spectrum_key(index), "#eta#rightarrow#gamma#gamma: " + ::kLabels[index], ::kColors[index]);
  auto residual_count = residual_histogram(*composition.counts[eta], *topology.counts[1], "h_eta_residual_count");
  auto residual_weight = residual_histogram(*composition.weighted[eta], *topology.weighted_pb[1], "h_eta_residual_pb");
  append(*residual_count, *residual_weight, "eta_residual", "Other #eta contributions", candidate_composition::kColors[eta]);
  append(*composition.counts[other], *composition.weighted[other], "other", "Other origin", candidate_composition::kColors[other]);
  return result;
}

inline bool export_outputs(const std::string& output_base, const candidate_composition_merge::Metadata& metadata,
                           const std::vector<double>& cross_sections,
                           const std::array<std::unique_ptr<candidate_composition::Histograms>, kSelectionCount>& composition,
                           const std::array<std::unique_ptr<Spectra>, kSelectionCount>& topology_histograms,
                           const std::array<pi0_bdt::Histograms, kSelectionCount>& bdt)
{
  const auto original_labels = ::kLabels;
  struct RestoreLabels { decltype(::kLabels) labels; ~RestoreLabels() { ::kLabels = labels; } } restore{original_labels};
  ::kLabels[1] = "#eta#rightarrow#gamma#gamma-main anchor";
  TFile topology_output((output_base + "/eta_anchor_topology.root").c_str(), "RECREATE");
  TFile composition_output((output_base + "/candidate_composition_eta_topology.root").c_str(), "RECREATE");
  if (topology_output.IsZombie() || composition_output.IsZombie()) return false;
  for (std::size_t topology_selection = 0; topology_selection < kSelectionCount; ++topology_selection)
  {
    plot_caption.selection = kSelectionLabels[topology_selection];
    Spectra& spectra = *topology_histograms[topology_selection];
    const std::string selection_output_base = output_base + "/" + kSelectionKeys[topology_selection] + "/region_a_eta_anchor_topology";
    if (!make_output_directory(selection_output_base)) return false;
    std::array<std::unique_ptr<TH1D>, kSpectrumCount> density;
    for (std::size_t index = 0; index < kSpectrumCount; ++index)
    {
      density[index].reset(static_cast<TH1D*>(spectra.weighted_pb[index]->Clone((std::string("h_") + kSelectionKeys[topology_selection] + "_" + spectrum_key(index) + "_et_pb_per_gev").c_str())));
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
      const std::string label = index == 1 ? "All #eta#rightarrow#gamma#gamma-main anchors" : ::kLabels[index];
      topology_survival[index] = make_survival_curve(*spectra.weighted_pb[index], *topology_histograms[0]->weighted_pb[index],
          std::string("h_") + kSelectionKeys[topology_selection] + "_" + spectrum_key(index) + "_survival_fraction_relative_to_kinematic",
          label, color);
    }
    std::vector<const SurvivalCurve*> topology_summary_survival = {topology_survival[1].get()};
    for (std::size_t index : summary_category_indices) topology_summary_survival.push_back(topology_survival[index].get());
    std::vector<const SurvivalCurve*> topology_detailed_survival = {topology_survival[1].get()};
    for (std::size_t index : detailed_category_indices) topology_detailed_survival.push_back(topology_survival[index].get());
    draw_survival_curves(topology_summary_survival, *spectra.weighted_pb[1],
        selection_output_base + "_survival_fraction_relative_to_kinematic.pdf",
        std::string("c_") + kSelectionKeys[topology_selection] + "_eta_anchor_survival", false);
    draw_survival_curves(topology_detailed_survival, *spectra.weighted_pb[1],
        selection_output_base + "_survival_fraction_relative_to_kinematic_detailed.pdf",
        std::string("c_") + kSelectionKeys[topology_selection] + "_eta_anchor_survival_detailed", true);

    TDirectory* directory = topology_output.mkdir(kSelectionKeys[topology_selection]);
    if (!directory) return false;
    directory->cd();
    for (std::size_t index = 0; index < kSpectrumCount; ++index)
    {
      spectra.counts[index]->Write((std::string("h_region_a_") + spectrum_key(index) + "_et_count").c_str());
      spectra.weighted_pb[index]->Write((std::string("h_region_a_") + spectrum_key(index) + "_et_pb").c_str());
      density[index]->Write();
      if (index > 0) topology_survival[index]->histogram->Write();
    }
    for (auto& histogram : detailed_fractions) histogram->Write();
    for (auto& histogram : summary_fractions) histogram->Write();
    topology_output.cd();
  }

  write_metadata(topology_output, metadata, cross_sections);
  write_metadata(composition_output, metadata, cross_sections);
  for (std::size_t selection = 0; selection < kSelectionCount; ++selection)
  {
    const std::string directory = output_base + "/" + kSelectionKeys[selection];
    plot_caption.selection = kSelectionLabels[selection];
    auto caption = caption_lines(plot_caption);
    caption[2] += "; #eta#rightarrow#gamma#gamma-main anchors";
    if (!bdt[selection].export_files(directory, kSelectionLabels[selection], selection == 3 || selection >= 4, caption, "eta")) return false;
    { TFile scores((directory + "/eta_bdt_score.root").c_str(), "UPDATE"); write_metadata(scores, metadata, cross_sections); }
    auto* root_directory = composition_output.mkdir(kSelectionKeys[selection]);
    if (!root_directory) return false;
    for (bool detailed : {false, true})
    {
      const std::string detail = detailed ? "superdetailed" : "detailed";
      auto components = composition_components(*composition[selection], *topology_histograms[selection], detailed);
      auto kinematic = composition_components(*composition[0], *topology_histograms[0], detailed);
      std::vector<std::unique_ptr<TH1D>> fractions;
      std::vector<TH1D*> fraction_pointers;
      std::vector<std::unique_ptr<SurvivalCurve>> survival;
      std::vector<const SurvivalCurve*> survival_pointers;
      for (std::size_t index = 0; index < components.weighted.size(); ++index)
      {
        const std::string key = components.weighted[index]->GetName();
        fractions.push_back(candidate_composition::fraction_histogram(*components.weighted[index], *composition[selection]->weighted[candidate_composition::denominator], key + "_fraction"));
        fraction_pointers.push_back(fractions.back().get());
        survival.push_back(make_survival_curve(*components.weighted[index], *kinematic.weighted[index], key + "_survival_fraction_relative_to_kinematic", components.labels[index], components.colors[index]));
        survival_pointers.push_back(survival.back().get());
      }
      const std::string base = directory + "/photon_candidate_composition_eta_topology";
      candidate_composition::draw_composition_components(fraction_pointers, components.labels, components.colors,
          base + "_category_fraction_stack_" + detail + ".pdf", "superdetailed");
      if (!detailed) draw_survival_curves(survival_pointers, *composition[selection]->weighted[candidate_composition::denominator],
          base + "_survival_fraction_relative_to_kinematic_detailed.pdf", "c_eta_composition_survival", true);
      auto* detail_directory = root_directory->mkdir(detail.c_str());
      if (!detail_directory) return false;
      detail_directory->cd();
      composition[selection]->counts[candidate_composition::denominator]->Write();
      composition[selection]->weighted[candidate_composition::denominator]->Write();
      for (std::size_t index = 0; index < components.weighted.size(); ++index)
      {
        components.counts[index]->Write(); components.weighted[index]->Write(); fractions[index]->Write(); survival[index]->histogram->Write();
      }
    }
  }
  return true;
}
}
#endif
