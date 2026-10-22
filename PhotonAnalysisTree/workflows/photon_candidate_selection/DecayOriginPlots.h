#ifndef DECAY_ORIGIN_PLOTS_H
#define DECAY_ORIGIN_PLOTS_H
namespace decay_origin
{
inline bool export_outputs(const std::string& base, const candidate_composition_merge::Metadata& source,
                           const std::array<std::unique_ptr<candidate_composition::Histograms>, kSelectionCount>& composition,
                           const Selections& selections)
{
  using namespace candidate_composition_merge;
  TFile file((base + "/decay_origin.root").c_str(), "RECREATE"); if (file.IsZombie()) return false;
  auto saved = source; TTree metadata("metadata", "Supplemental decay-origin diagnostic policy");
  std::string population = "existing_composition_other_only";
  std::string note = "pi0/eta diphoton majority is not a topology or clean single-photon label;zero flight is in first log bin;unknown flight/photon share in underflow";
  metadata.Branch("population", &population); metadata.Branch("diagnostic_note", &note);
#define META(member) metadata.Branch(#member, &saved.member)
  META(family); META(map_root); META(model_sha256); META(source_map_schema_version); META(decay_origin_algorithm_version);
  META(decay_origin_max_proper_length_cm); META(decay_origin_max_flight_cm); META(sample_names); META(sample_map_counts); META(sample_sum_generator_weights);
#undef META
  metadata.Fill(); metadata.Write();
  std::ofstream csv(base + "/decay_origin_summary.csv"); if (!csv) return false;
  csv << "selection,bin,et_low,et_high,category,raw_count,weighted_pb,other_raw_count,other_weighted_pb,all_raw_count,all_weighted_pb,fraction_of_other,fraction_of_all\n" << std::setprecision(17);
  for (std::size_t selection = 0; selection < kSelectionCount; ++selection)
  {
    const auto& c = *composition[selection]; const auto& h = *selections[selection];
    const auto& other_n = *c.counts[candidate_composition::other]; const auto& other_w = *c.weighted[candidate_composition::other];
    const auto& all_n = *c.counts[candidate_composition::denominator]; const auto& all_w = *c.weighted[candidate_composition::denominator];
    if (!h.valid(other_n, other_w)) return false;
    auto* directory = file.mkdir(kSelectionKeys[selection]); if (!directory || !h.write(*directory)) return false;
    plot_caption.selection = kSelectionLabels[selection];
    std::vector<std::unique_ptr<TH1D>> fractions; std::vector<TH1D*> fraction_pointers;
    std::vector<std::unique_ptr<SurvivalCurve>> survival; std::vector<const SurvivalCurve*> survival_pointers;
    std::vector<std::string> plot_labels; std::vector<int> plot_colors;
    for (std::size_t i = 0; i < category_count; ++i)
    {
      fractions.push_back(candidate_composition::fraction_histogram(*h.weighted[i], other_w, std::string(h.weighted[i]->GetName()) + "_fraction_of_other"));
      fraction_pointers.push_back(fractions.back().get()); plot_labels.push_back(labels[i]); plot_colors.push_back(colors[i]);
      survival.push_back(make_survival_curve(*h.weighted[i], *selections[0]->weighted[i], std::string(h.weighted[i]->GetName()) + "_survival", labels[i], colors[i]));
      survival_pointers.push_back(survival.back().get());
    }
    const std::string output = base + "/" + kSelectionKeys[selection];
    candidate_composition::draw_composition_components(fraction_pointers, plot_labels, plot_colors,
        output + "/photon_candidate_other_decay_origin_category_fraction_stack.pdf", "superdetailed", "Fraction within Other");
    draw_survival_curves(survival_pointers, other_w, output + "/photon_candidate_other_decay_origin_survival_fraction_relative_to_kinematic.pdf", "c_decay_origin_survival", true);
    directory->cd(); for (const auto& f : fractions) f->Write(); for (const auto& curve : survival) curve->histogram->Write();
    for (int bin = -1; bin <= all_n.GetNbinsX() + 1; ++bin)
    {
      const auto value = [bin](const TH1D& h) { return bin < 0 ? h.Integral(0, h.GetNbinsX() + 1) : h.GetBinContent(bin); };
      for (std::size_t i = 0; i < category_count; ++i)
        csv << kSelectionKeys[selection] << ',' << bin << ',' << (bin < 0 ? 0 : all_n.GetXaxis()->GetBinLowEdge(bin)) << ','
            << (bin < 0 ? all_n.GetXaxis()->GetXmax() : all_n.GetXaxis()->GetBinUpEdge(bin)) << ',' << keys[i] << ','
            << value(*h.counts[i]) << ',' << value(*h.weighted[i]) << ',' << value(other_n) << ',' << value(other_w) << ','
            << value(all_n) << ',' << value(all_w) << ',' << (value(other_w) != 0 ? value(*h.weighted[i]) / value(other_w) : 0) << ','
            << (value(all_w) != 0 ? value(*h.weighted[i]) / value(all_w) : 0) << '\n';
    }
  }
  return bool(csv);
}
}
#endif
