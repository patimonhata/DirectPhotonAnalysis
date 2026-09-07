#include "ReduceTaggingPurityDiagnostic.C"

#include <fstream>
#include <iomanip>

namespace
{
double subset_error(const TH1D& numerator, const TH1D& denominator, int bin)
{
  const double den = denominator.GetBinContent(bin);
  if (den == 0) return 0;
  const double num = numerator.GetBinContent(bin);
  const double fraction = num / den;
  const double num_w2 = std::pow(numerator.GetBinError(bin), 2);
  const double complement_w2 = std::max(0.0, std::pow(denominator.GetBinError(bin), 2) - num_w2);
  return std::sqrt(std::max(0.0, (std::pow(1.0 - fraction, 2) * num_w2 + std::pow(fraction, 2) * complement_w2) / std::pow(den, 2)));
}

std::unique_ptr<TH1D> add_histograms(const TH1D& left, const TH1D& right, const std::string& name)
{
  auto result = std::unique_ptr<TH1D>(static_cast<TH1D*>(left.Clone(name.c_str())));
  result->SetDirectory(nullptr);
  result->Add(&right);
  return result;
}
}

int SummarizeTaggingPurityDiagnostic(const std::string input_file, const std::string output_prefix)
{
  TFile input(input_file.c_str(), "READ");
  auto* metadata = input.Get<TTree>("metadata");
  if (input.IsZombie() || !metadata || metadata->GetEntries() != 16 || output_prefix.empty()) return 1;
  if (!photon_candidate_settings::validate_partials(*metadata, "selection_settings")) return 1;
  auto get = [&](const std::string& key, const std::string& suffix) {
    return input.Get<TH1D>(("flow/h_" + key + "_et_" + suffix).c_str());
  };
  for (const char* key : kFlowKeys)
    if (!get(key, "count") || !get(key, "pb")) return 2;

  std::ofstream flow_output(output_prefix + "_flow.tsv");
  std::ofstream metric_output(output_prefix + "_metrics.tsv");
  if (!flow_output || !metric_output) return 3;
  flow_output << std::setprecision(12) << "et_low\tet_high\tkey\traw_count\tweighted_pb\tweighted_pb_error\n";
  for (int bin = 1; bin <= kEtBins; ++bin)
  {
    auto* axis = get("all_before", "pb");
    for (const char* key : kFlowKeys)
    {
      auto* count = get(key, "count");
      auto* pb = get(key, "pb");
      flow_output << axis->GetXaxis()->GetBinLowEdge(bin) << '\t' << axis->GetXaxis()->GetBinUpEdge(bin) << '\t' << key << '\t'
                  << count->GetBinContent(bin) << '\t' << pb->GetBinContent(bin) << '\t' << pb->GetBinError(bin) << '\n';
    }
  }

  metric_output << std::setprecision(12) << "et_low\tet_high\tmetric\tvalue\tstat_error\tnumerator_pb\tdenominator_pb\n";
  auto write_ratio = [&](int bin, const std::string& key, const TH1D& numerator, const TH1D& denominator) {
    const double den = denominator.GetBinContent(bin);
    const double num = numerator.GetBinContent(bin);
    metric_output << denominator.GetXaxis()->GetBinLowEdge(bin) << '\t' << denominator.GetXaxis()->GetBinUpEdge(bin) << '\t' << key << '\t'
                  << (den != 0 ? num / den : 0) << '\t' << subset_error(numerator, denominator, bin) << '\t' << num << '\t' << den << '\n';
  };
  auto* all_before = get("all_before", "pb");
  auto* all_after = get("all_after", "pb");
  auto* prompt_before = get("prompt_before", "pb");
  auto* prompt_after = get("prompt_after", "pb");
  auto* background_before = get("background_before", "pb");
  auto* background_after = get("background_after", "pb");
  auto* truth_pi0_before = get("truth_pi0_before", "pb");
  auto* truth_pi0_after = get("truth_pi0_after", "pb");
  auto prompt_pi0_veto = add_histograms(*get("prompt_pi0_only_veto", "pb"), *get("prompt_both_veto", "pb"), "prompt_pi0_veto");
  auto prompt_eta_veto = add_histograms(*get("prompt_eta_only_veto", "pb"), *get("prompt_both_veto", "pb"), "prompt_eta_veto");
  for (int bin = 1; bin <= kEtBins; ++bin)
  {
    write_ratio(bin, "purity_before", *prompt_before, *all_before);
    write_ratio(bin, "purity_after", *prompt_after, *all_after);
    write_ratio(bin, "all_survival", *all_after, *all_before);
    write_ratio(bin, "prompt_survival", *prompt_after, *prompt_before);
    write_ratio(bin, "background_survival", *background_after, *background_before);
    write_ratio(bin, "truth_pi0_survival", *truth_pi0_after, *truth_pi0_before);
    write_ratio(bin, "prompt_pi0_veto_fraction", *prompt_pi0_veto, *prompt_before);
    write_ratio(bin, "prompt_eta_veto_fraction", *prompt_eta_veto, *prompt_before);
  }
  flow_output.close();
  metric_output.close();
  const int bin = all_before->FindBin(5.5);
  std::cout << "ET=[5,6): purity before/after = " << prompt_before->GetBinContent(bin) / all_before->GetBinContent(bin) << "/"
            << prompt_after->GetBinContent(bin) / all_after->GetBinContent(bin) << ", prompt/background survival = "
            << prompt_after->GetBinContent(bin) / prompt_before->GetBinContent(bin) << "/"
            << background_after->GetBinContent(bin) / background_before->GetBinContent(bin) << std::endl;
  return 0;
}
