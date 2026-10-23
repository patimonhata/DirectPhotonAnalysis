#include "ReduceHybridThresholdDiagnostic.C"

#include <fstream>
#include <iomanip>

namespace
{
double hybrid_subset_error(const TH1D& numerator, const TH1D& denominator, int bin)
{
  const double den = denominator.GetBinContent(bin);
  if (den == 0) return 0;
  const double num = numerator.GetBinContent(bin);
  const double fraction = num / den;
  const double num_w2 = std::pow(numerator.GetBinError(bin), 2);
  const double complement_w2 = std::max(0.0, std::pow(denominator.GetBinError(bin), 2) - num_w2);
  return std::sqrt(std::max(0.0, (std::pow(1.0 - fraction, 2) * num_w2 + std::pow(fraction, 2) * complement_w2) / std::pow(den, 2)));
}

double ratio_error(double numerator, double numerator_w2, double denominator, double denominator_w2, double covariance)
{
  if (denominator == 0) return 0;
  const double variance = numerator_w2 / std::pow(denominator, 2) + std::pow(numerator, 2) * denominator_w2 / std::pow(denominator, 4) -
      2.0 * numerator * covariance / std::pow(denominator, 3);
  return std::sqrt(std::max(0.0, variance));
}

double purity_variance(double signal, double signal_w2, double background, double background_w2)
{
  const double total = signal + background;
  if (total == 0) return 0;
  return (std::pow(background, 2) * signal_w2 + std::pow(signal, 2) * background_w2) / std::pow(total, 4);
}
}

int SummarizeHybridThresholdDiagnostic(const std::string input_file, const std::string output_prefix)
{
  TFile input(input_file.c_str(), "READ");
  auto* metadata = input.Get<TTree>("metadata");
  if (input.IsZombie() || !metadata || metadata->GetEntries() != 16 || output_prefix.empty()) return 1;
  if (!photon_candidate_settings::validate_partials(*metadata, "low_selection_settings")) return 1;
  if (!photon_candidate_settings::validate_partials(*metadata, "high_selection_settings")) return 1;
  auto get = [&](const std::string& wp, const std::string& key, const std::string& suffix) {
    return input.Get<TH1D>(("working_points/" + wp + "/h_" + wp + "_" + key + "_et_" + suffix).c_str());
  };
  for (const char* wp : kWorkingPoints)
    for (const char* key : kHybridFlowKeys)
      if (!get(wp, key, "count") || !get(wp, key, "pb")) return 2;

  unsigned long long uid_mismatches = 0, missing_anchors = 0, region_mismatches = 0;
  metadata->SetBranchAddress("uid_mismatches", &uid_mismatches);
  metadata->SetBranchAddress("missing_anchors", &missing_anchors);
  metadata->SetBranchAddress("region_mismatches", &region_mismatches);
  for (Long64_t entry = 0; entry < metadata->GetEntries(); ++entry)
  {
    metadata->GetEntry(entry);
    if (uid_mismatches != 0 || missing_anchors != 0 || region_mismatches != 0) return 2;
  }

  std::ofstream flow_output(output_prefix + "_flow.tsv");
  std::ofstream metric_output(output_prefix + "_metrics.tsv");
  std::ofstream comparison_output(output_prefix + "_comparisons.tsv");
  if (!flow_output || !metric_output || !comparison_output) return 3;
  flow_output << std::setprecision(12) << "et_low\tet_high\tworking_point\tkey\traw_count\tweighted_pb\tweighted_pb_error\n";
  for (int bin = 1; bin <= kEtBins; ++bin)
  {
    for (const char* wp : kWorkingPoints)
    {
      auto* axis = get(wp, "all_before", "pb");
      for (const char* key : kHybridFlowKeys)
      {
        auto* count = get(wp, key, "count");
        auto* pb = get(wp, key, "pb");
        flow_output << axis->GetXaxis()->GetBinLowEdge(bin) << '\t' << axis->GetXaxis()->GetBinUpEdge(bin) << '\t' << wp << '\t' << key << '\t'
                    << count->GetBinContent(bin) << '\t' << pb->GetBinContent(bin) << '\t' << pb->GetBinError(bin) << '\n';
      }
    }
  }

  metric_output << std::setprecision(12) << "et_low\tet_high\tworking_point\tmetric\tvalue\tstat_error\tnumerator_pb\tdenominator_pb\n";
  auto write_ratio = [&](int bin, const char* wp, const std::string& metric, const TH1D& numerator, const TH1D& denominator) {
    const double den = denominator.GetBinContent(bin);
    const double num = numerator.GetBinContent(bin);
    metric_output << denominator.GetXaxis()->GetBinLowEdge(bin) << '\t' << denominator.GetXaxis()->GetBinUpEdge(bin) << '\t' << wp << '\t' << metric << '\t'
                  << (den != 0 ? num / den : 0) << '\t' << hybrid_subset_error(numerator, denominator, bin) << '\t' << num << '\t' << den << '\n';
  };
  for (int bin = 1; bin <= kEtBins; ++bin)
  {
    for (const char* wp : kWorkingPoints)
    {
      write_ratio(bin, wp, "purity_before", *get(wp, "prompt_before", "pb"), *get(wp, "all_before", "pb"));
      write_ratio(bin, wp, "purity_after", *get(wp, "prompt_after", "pb"), *get(wp, "all_after", "pb"));
      write_ratio(bin, wp, "all_survival", *get(wp, "all_after", "pb"), *get(wp, "all_before", "pb"));
      write_ratio(bin, wp, "prompt_survival", *get(wp, "prompt_after", "pb"), *get(wp, "prompt_before", "pb"));
      write_ratio(bin, wp, "background_survival", *get(wp, "background_after", "pb"), *get(wp, "background_before", "pb"));
      write_ratio(bin, wp, "truth_pi0_tag_fraction", *get(wp, "truth_pi0_pi0_veto", "pb"), *get(wp, "truth_pi0_before", "pb"));
      write_ratio(bin, wp, "truth_pi0_union_survival", *get(wp, "truth_pi0_after", "pb"), *get(wp, "truth_pi0_before", "pb"));
    }
  }

  comparison_output << std::setprecision(12)
                    << "et_low\tet_high\treference\talternative\tmetric\tvalue\tstat_error\treference_value\talternative_value\tcovariance_w2\n";
  const std::array<std::pair<int, int>, 3> comparisons = {{{0, 2}, {1, 2}, {0, 1}}};
  for (int bin = 1; bin <= kEtBins; ++bin)
  {
    for (const auto& comparison : comparisons)
    {
      const int reference = comparison.first, alternative = comparison.second;
      const std::string pair_name = std::string(kWorkingPoints[std::min(reference, alternative)]) + "__" +
          kWorkingPoints[std::max(reference, alternative)];
      for (const char* population : {"all", "prompt", "background"})
      {
        auto* ref = get(kWorkingPoints[reference], std::string(population) + "_after", "pb");
        auto* alt = get(kWorkingPoints[alternative], std::string(population) + "_after", "pb");
        auto* overlap = input.Get<TH1D>(("overlap/h_" + std::string(population) + "_after_overlap_" + pair_name + "_et_pb").c_str());
        if (!overlap) return 2;
        const double ref_value = ref->GetBinContent(bin), alt_value = alt->GetBinContent(bin);
        const double covariance = std::pow(overlap->GetBinError(bin), 2);
        comparison_output << ref->GetXaxis()->GetBinLowEdge(bin) << '\t' << ref->GetXaxis()->GetBinUpEdge(bin) << '\t'
                          << kWorkingPoints[reference] << '\t' << kWorkingPoints[alternative] << '\t' << population << "_yield_ratio\t"
                          << (ref_value != 0 ? alt_value / ref_value : 0) << '\t'
                          << ratio_error(alt_value, std::pow(alt->GetBinError(bin), 2), ref_value, std::pow(ref->GetBinError(bin), 2), covariance) << '\t'
                          << ref_value << '\t' << alt_value << '\t' << covariance << '\n';
      }

      auto* ref_signal = get(kWorkingPoints[reference], "prompt_after", "pb");
      auto* ref_background = get(kWorkingPoints[reference], "background_after", "pb");
      auto* alt_signal = get(kWorkingPoints[alternative], "prompt_after", "pb");
      auto* alt_background = get(kWorkingPoints[alternative], "background_after", "pb");
      auto* signal_overlap = input.Get<TH1D>(("overlap/h_prompt_after_overlap_" + pair_name + "_et_pb").c_str());
      auto* background_overlap = input.Get<TH1D>(("overlap/h_background_after_overlap_" + pair_name + "_et_pb").c_str());
      const double rs = ref_signal->GetBinContent(bin), rb = ref_background->GetBinContent(bin);
      const double as = alt_signal->GetBinContent(bin), ab = alt_background->GetBinContent(bin);
      const double ref_total = rs + rb, alt_total = as + ab;
      if (ref_total == 0 || alt_total == 0) continue;
      const double ref_purity = rs / ref_total, alt_purity = as / alt_total;
      const double purity_covariance = rb / std::pow(ref_total, 2) * ab / std::pow(alt_total, 2) * std::pow(signal_overlap->GetBinError(bin), 2) +
          rs / std::pow(ref_total, 2) * as / std::pow(alt_total, 2) * std::pow(background_overlap->GetBinError(bin), 2);
      const double difference_variance = purity_variance(rs, std::pow(ref_signal->GetBinError(bin), 2), rb, std::pow(ref_background->GetBinError(bin), 2)) +
          purity_variance(as, std::pow(alt_signal->GetBinError(bin), 2), ab, std::pow(alt_background->GetBinError(bin), 2)) - 2.0 * purity_covariance;
      comparison_output << ref_signal->GetXaxis()->GetBinLowEdge(bin) << '\t' << ref_signal->GetXaxis()->GetBinUpEdge(bin) << '\t'
                        << kWorkingPoints[reference] << '\t' << kWorkingPoints[alternative] << "\tpurity_difference\t"
                        << alt_purity - ref_purity << '\t' << std::sqrt(std::max(0.0, difference_variance)) << '\t'
                        << ref_purity << '\t' << alt_purity << '\t' << purity_covariance << '\n';
    }
  }

  const int bin = get("both_0p5", "all_before", "pb")->FindBin(5.5);
  for (const char* wp : kWorkingPoints)
  {
    auto* prompt_after = get(wp, "prompt_after", "pb");
    auto* all_after = get(wp, "all_after", "pb");
    auto* prompt_before = get(wp, "prompt_before", "pb");
    auto* background_after = get(wp, "background_after", "pb");
    auto* background_before = get(wp, "background_before", "pb");
    std::cout << "ET=[5,6) " << wp << ": purity=" << prompt_after->GetBinContent(bin) / all_after->GetBinContent(bin)
              << ", prompt/background survival=" << prompt_after->GetBinContent(bin) / prompt_before->GetBinContent(bin) << "/"
              << background_after->GetBinContent(bin) / background_before->GetBinContent(bin) << std::endl;
  }
  return 0;
}
