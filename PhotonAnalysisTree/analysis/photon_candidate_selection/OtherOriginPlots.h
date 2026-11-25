#ifndef OTHER_ORIGIN_PLOTS_H
#define OTHER_ORIGIN_PLOTS_H
#include "CandidateMerge.h"
#include <fstream>
#include <iomanip>
namespace other_origin
{
struct Components
{
  std::vector<std::unique_ptr<TH1D>> counts, weighted;
  std::vector<std::string> labels;
  std::vector<int> colors;
};
inline Components components(const candidate_composition::Histograms& composition, const Histograms& other, bool detailed, bool all)
{
  Components result;
  const auto append = [&](const TH1D& n, const TH1D& w, const std::string& key, const std::string& label, int color) {
    result.counts.emplace_back(static_cast<TH1D*>(n.Clone(("h_other_composition_" + key + "_et_count").c_str())));
    result.weighted.emplace_back(static_cast<TH1D*>(w.Clone(("h_other_composition_" + key + "_et_pb").c_str())));
    result.counts.back()->SetDirectory(nullptr); result.weighted.back()->SetDirectory(nullptr);
    result.labels.push_back(label); result.colors.push_back(color);
  };
  if (all)
  {
    using namespace candidate_composition;
    append(*composition.counts[prompt], *composition.weighted[prompt], "prompt", "Prompt-#gamma", candidate_composition::kColors[prompt]);
    auto n = std::unique_ptr<TH1D>(static_cast<TH1D*>(composition.counts[pi0_separated]->Clone("h_other_all_pi0_n")));
    auto w = std::unique_ptr<TH1D>(static_cast<TH1D*>(composition.weighted[pi0_separated]->Clone("h_other_all_pi0_w")));
    n->SetDirectory(nullptr); w->SetDirectory(nullptr);
    for (std::size_t i = pi0_separated + 1; i <= pi0_other; ++i) { n->Add(composition.counts[i].get()); w->Add(composition.weighted[i].get()); }
    append(*n, *w, "pi0", "#pi^{0}", candidate_composition::kColors[pi0_separated]);
    append(*composition.counts[eta], *composition.weighted[eta], "eta", "#eta", candidate_composition::kColors[eta]);
  }
  const std::size_t start = result.counts.size();
  for (std::size_t i = 0; i < category_count; ++i)
  {
    if (detailed) append(*other.counts[i], *other.weighted[i], keys[i], labels[i], colors[i]);
    else
    {
      const int group = summary_category(i);
      if (result.counts.size() <= start + group)
        append(*other.counts[i], *other.weighted[i], summary_keys[group], summary_labels[group], summary_colors[group]);
      else { result.counts[start + group]->Add(other.counts[i].get()); result.weighted[start + group]->Add(other.weighted[i].get()); }
    }
  }
  return result;
}
inline void metadata(TFile& file, const candidate_composition_merge::Metadata& source)
{
  file.cd(); auto value = source;
  int schema = 1, algorithm = kAlgorithmVersion; double threshold = .5;
  std::string population = "existing_composition_other_only";
  std::string charged = "single_signal_g4_contributor_abs_pdg_11_13_211_321_2212_fraction_gt_0p5";
  std::string photon = "single_signal_g4_primary_photon_fraction_gt_0p5_hepmc_physical_parent_meson_or_baryon_except_pi0_eta";
  std::string purity = "majority_is_not_a_clean_single_photon_training_label;other_primary_gamma_fraction_does_not_include_g4_primary_pi0_descendants";
  TTree tree("metadata", "Other-origin diagnostic definitions and source");
  tree.Branch("schema_version", &schema); tree.Branch("other_origin_algorithm_version", &algorithm);
  tree.Branch("majority_threshold", &threshold); tree.Branch("population", &population);
  tree.Branch("charged_definition", &charged); tree.Branch("decay_photon_definition", &photon); tree.Branch("training_label_note", &purity);
#define OTHER_METADATA(member) tree.Branch(#member, &value.member)
  OTHER_METADATA(family); OTHER_METADATA(map_root); OTHER_METADATA(require_complete); OTHER_METADATA(analysis_release); OTHER_METADATA(model_sha256);
  OTHER_METADATA(source_map_schema_version); OTHER_METADATA(pi0_topology_algorithm_version); OTHER_METADATA(eta_topology_algorithm_version);
  OTHER_METADATA(sample_names); OTHER_METADATA(sample_map_counts); OTHER_METADATA(sample_sum_generator_weights);
  OTHER_METADATA(min_cluster_energy); OTHER_METADATA(pi0_partner_min_energy); OTHER_METADATA(eta_partner_min_energy);
#undef OTHER_METADATA
  tree.Branch("source_partial_schema_version", &value.schema_version);
  value.mbd.create_branches(tree);
  tree.Fill(); tree.Write();
}
inline bool export_outputs(const std::string& base, const candidate_composition_merge::Metadata& source,
                           const std::array<std::unique_ptr<candidate_composition::Histograms>, kSelectionCount>& composition,
                           const Selections& histograms)
{
  using namespace candidate_composition_merge;
  TFile file((base + "/other_origin.root").c_str(), "RECREATE"); if (file.IsZombie()) return false;
  metadata(file, source);
  std::ofstream csv(base + "/other_origin_summary.csv"); if (!csv) return false;
  csv << "selection,bin,et_low,et_high,category,raw_count,weighted_pb,other_raw_count,other_weighted_pb,all_raw_count,all_weighted_pb,fraction_of_other,fraction_of_all\n" << std::setprecision(17);
  for (std::size_t selection = 0; selection < kSelectionCount; ++selection)
  {
    const auto& c = *composition[selection]; const auto& h = *histograms[selection];
    if (!h.valid(*c.counts[candidate_composition::other], *c.weighted[candidate_composition::other])) return false;
    auto* directory = file.mkdir(kSelectionKeys[selection]); if (!directory || !h.write(*directory)) return false;
    plot_caption.selection = kSelectionLabels[selection];
    const std::string output = base + "/" + kSelectionKeys[selection];
    for (bool detailed : {false, true})
    {
      const std::string suffix = detailed ? "_detailed" : "";
      for (bool all : {false, true})
      {
        auto parts = components(c, h, detailed, all);
        const auto& denominator = *c.weighted[all ? candidate_composition::denominator : candidate_composition::other];
        std::vector<std::unique_ptr<TH1D>> fractions; std::vector<TH1D*> pointers;
        for (std::size_t i = 0; i < parts.weighted.size(); ++i)
        {
          fractions.push_back(candidate_composition::fraction_histogram(*parts.weighted[i], denominator,
              std::string(parts.weighted[i]->GetName()) + (all ? "_fraction_all" : "_fraction_other")));
          pointers.push_back(fractions.back().get());
        }
        const std::string prefix = all ? "photon_candidate_composition_other_origin" : "photon_candidate_other_origin";
        candidate_composition::draw_composition_components(pointers, parts.labels, parts.colors,
            output + "/" + prefix + "_category_fraction_stack" + suffix + ".pdf", "superdetailed", all ? "Fraction of selected candidates" : "Fraction within Other");
        directory->cd(); for (const auto& histogram : fractions) histogram->Write();
      }
      auto selected = components(c, h, detailed, false);
      auto kinematic = components(*composition[0], *histograms[0], detailed, false);
      std::vector<std::unique_ptr<SurvivalCurve>> survival; std::vector<const SurvivalCurve*> pointers;
      for (std::size_t i = 0; i < selected.weighted.size(); ++i)
      {
        survival.push_back(make_survival_curve(*selected.weighted[i], *kinematic.weighted[i],
            std::string(selected.weighted[i]->GetName()) + "_survival_fraction_relative_to_kinematic", selected.labels[i], selected.colors[i]));
        pointers.push_back(survival.back().get());
      }
      draw_survival_curves(pointers, *c.weighted[candidate_composition::other],
          output + "/photon_candidate_other_origin_survival_fraction_relative_to_kinematic" + suffix + ".pdf", "c_other_origin_survival", true);
      directory->cd(); for (const auto& curve : survival) curve->histogram->Write();
    }
    const auto& total = *c.counts[candidate_composition::denominator];
    for (int bin = -1; bin <= total.GetNbinsX() + 1; ++bin)
    {
      const auto value = [bin](const TH1D& histogram) { return bin < 0 ? histogram.Integral(0, histogram.GetNbinsX() + 1) : histogram.GetBinContent(bin); };
      const double other_w = value(*c.weighted[candidate_composition::other]), all_w = value(*c.weighted[candidate_composition::denominator]);
      for (std::size_t i = 0; i < category_count; ++i)
        csv << kSelectionKeys[selection] << ',' << bin << ',' << (bin < 0 ? 0 : total.GetXaxis()->GetBinLowEdge(bin)) << ','
            << (bin < 0 ? total.GetXaxis()->GetXmax() : total.GetXaxis()->GetBinUpEdge(bin)) << ',' << keys[i] << ','
            << value(*h.counts[i]) << ',' << value(*h.weighted[i]) << ',' << value(*c.counts[candidate_composition::other]) << ',' << other_w << ','
            << value(total) << ',' << all_w << ',' << (other_w != 0 ? value(*h.weighted[i]) / other_w : 0) << ','
            << (all_w != 0 ? value(*h.weighted[i]) / all_w : 0) << '\n';
    }
  }
  return bool(csv);
}
}
#endif
