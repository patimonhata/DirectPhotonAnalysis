#ifndef COMMON_COMPOSITION_PLOTS_H
#define COMMON_COMPOSITION_PLOTS_H
#include "CandidateMerge.h"
namespace common_composition
{
inline bool export_outputs(const std::string& base, const candidate_composition_merge::Metadata& source, const Selections& selections)
{
  using namespace candidate_composition_merge;
  TFile file((base + "/common_composition.root").c_str(), "RECREATE"); if (file.IsZombie()) return false;
  auto saved = source; TTree metadata("metadata", "Boundary independent composition and matching quality, all origins");
  std::string policy = "units_before_downstream_hadron_parent;strict_majority;original_shower_deposit_denominator;no_lifetime_or_flight_cut";
  std::string matching = "split:complete_and_conserved;topology:complete_or_partial_coverage_ge_0.5";
  metadata.Branch("policy", &policy); metadata.Branch("matching", &matching);
#define M(member) metadata.Branch(#member, &saved.member)
  M(common_composition_algorithm_version); M(source_map_schema_version); M(schema_version); M(model_sha256); M(family); M(map_root);
  M(pi0_topology_algorithm_version); M(eta_topology_algorithm_version); M(sample_names); M(sample_map_counts); M(sample_sum_generator_weights);
#undef M
  metadata.Fill(); metadata.Write();
  std::ofstream csv(base + "/common_composition_migration.csv"), detail_csv(base + "/common_composition_summary.csv");
  if (!csv || !detail_csv) return false;
  csv << "selection,bin,et_low,et_high,legacy_origin,common_detail,count,weighted_pb,sumw2\n" << std::setprecision(17);
  detail_csv << "selection,bin,category,count,weighted_pb,sumw2,fraction_of_all\n" << std::setprecision(17);
  for (std::size_t s = 0; s < kSelectionCount; ++s)
  {
    const auto& h = *selections[s]; if (!h.valid()) return false;
    auto* d = file.mkdir(kSelectionKeys[s]); if (!d || !h.write(*d)) return false;
    plot_caption.selection = kSelectionLabels[s]; std::vector<std::unique_ptr<TH1D>> fractions;
    std::vector<TH1D*> pointers; std::vector<std::string> plot_labels; std::vector<int> colors;
    std::vector<std::unique_ptr<SurvivalCurve>> survival, previous_survival;
    std::vector<const SurvivalCurve*> curves, previous_curves;
    constexpr std::array<int, 12> palette = {kRed+1,kAzure+7,kMagenta+1,kCyan+2,kGreen+2,kGray+1,kViolet+1,kYellow+2,kOrange+7,kPink+7,kSpring+5,kBlue+1};
    for (std::size_t i = 1; i < keys.size(); ++i)
    {
      fractions.push_back(candidate_composition::fraction_histogram(*h.weighted[i], *h.weighted[0], std::string(keys[i]) + "_fraction"));
      pointers.push_back(fractions.back().get()); plot_labels.emplace_back(labels[i]); colors.push_back(palette[(i-1)%palette.size()]);
      survival.push_back(make_survival_curve(*h.weighted[i], *selections[0]->weighted[i], std::string(keys[i]) + "_survival_kinematic", labels[i], colors.back()));
      previous_survival.push_back(make_survival_curve(*h.weighted[i], *selections[kPreviousCutSelections[s]]->weighted[i], std::string(keys[i]) + "_survival_previous", labels[i], colors.back(), kPreviousCutLabels[s]));
      curves.push_back(survival.back().get()); previous_curves.push_back(previous_survival.back().get());
      d->cd(); fractions.back()->Write(); survival.back()->histogram->Write(); previous_survival.back()->histogram->Write();
    }
    const auto output = base + "/" + kSelectionKeys[s] + "/photon_candidate_common_composition";
    candidate_composition::draw_composition_components(pointers, plot_labels, colors, output + "_category_fraction_stack.pdf", "superdetailed");
    draw_survival_curves(curves, *h.weighted[0], output + "_survival_fraction_relative_to_kinematic.pdf", "c_common_survival", true);
    draw_survival_curves(previous_curves, *h.weighted[0], output + "_survival_fraction_relative_to_previous_cut.pdf", "c_common_previous", true);
    for (int b = 0; b <= h.counts[0]->GetNbinsX() + 1; ++b)
    {
      for (std::size_t i = 1; i < keys.size(); ++i)
      {
        const auto* axis = h.counts[0]->GetXaxis();
        detail_csv << kSelectionKeys[s] << ',' << b << ',' << keys[i] << ',' << h.counts[i]->GetBinContent(b) << ',' << h.weighted[i]->GetBinContent(b) << ','
            << std::pow(h.weighted[i]->GetBinError(b),2) << ',' << fractions[i-1]->GetBinContent(b) << '\n';
        for (int old = 0; old < 4; ++old)
          csv << kSelectionKeys[s] << ',' << b << ',' << axis->GetBinLowEdge(b) << ',' << axis->GetBinUpEdge(b) << ',' << kCandidateOriginKeys[old] << ','
              << keys[i] << ',' << h.migration_count.GetBinContent(b, old+1, i) << ',' << h.migration_weighted.GetBinContent(b, old+1, i) << ','
              << std::pow(h.migration_weighted.GetBinError(h.migration_weighted.GetBin(b, old+1, i)),2) << '\n';
      }
    }
  }
  return bool(csv) && bool(detail_csv);
}
}
#endif
