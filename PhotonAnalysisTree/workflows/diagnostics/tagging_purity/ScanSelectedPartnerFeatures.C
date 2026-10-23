#include "ReduceTaggingPurityDiagnostic.C"

#include <fstream>
#include <iomanip>
#include <limits>

int ScanSelectedPartnerFeatures(const std::string input_file, const std::string output_file, const double minimum_true_efficiency = 0.9,
    const std::string roc_output_file = "")
{
  if (input_file.empty() || output_file.empty() || !(minimum_true_efficiency > 0 && minimum_true_efficiency <= 1)) return 1;
  TFile input(input_file.c_str(), "READ");
  if (input.IsZombie()) return 1;
  std::ofstream output(output_file);
  const bool write_roc = !roc_output_file.empty();
  std::ofstream roc_output;
  if (write_roc) roc_output.open(roc_output_file);
  if (!output || (write_roc && !roc_output)) return 2;
  output << std::setprecision(12)
         << "feature\tdirection\tcut\ttrue_partner_efficiency\tprompt_combinatorial_acceptance\tprompt_combinatorial_rejection"
            "\ttrue_partner_pb\tprompt_combinatorial_pb\n";
  if (write_roc)
    roc_output << std::setprecision(12)
               << "feature\tdirection\tcut\ttrue_partner_efficiency\tprompt_combinatorial_acceptance\tprompt_combinatorial_rejection\n";

  for (const char* feature : kFeatureKeys)
  {
    auto* truth = input.Get<TH1D>(("selected_partner_features/truth_pi0_selected_truth_partner/h_truth_pi0_selected_truth_partner_" +
        std::string(feature) + "_pb").c_str());
    auto* prompt = input.Get<TH1D>(("selected_partner_features/prompt_combinatorial/h_prompt_combinatorial_" + std::string(feature) + "_pb").c_str());
    if (!truth || !prompt || truth->GetNbinsX() != prompt->GetNbinsX()) return 3;
    const int bins = truth->GetNbinsX();
    const double truth_total = truth->Integral(0, bins + 1);
    const double prompt_total = prompt->Integral(0, bins + 1);
    if (!(truth_total > 0) || !(prompt_total > 0)) continue;

    std::string best_direction;
    double best_cut = 0, best_truth_efficiency = 0, best_prompt_acceptance = std::numeric_limits<double>::infinity();
    for (int bin = 0; bin <= bins + 1; ++bin)
    {
      const double truth_efficiency = truth->Integral(0, bin) / truth_total;
      const double prompt_acceptance = prompt->Integral(0, bin) / prompt_total;
      const double cut = bin == 0 ? -std::numeric_limits<double>::infinity() :
          bin == bins + 1 ? std::numeric_limits<double>::infinity() : truth->GetXaxis()->GetBinUpEdge(bin);
      if (write_roc)
        roc_output << feature << "\tless_equal\t" << cut << '\t' << truth_efficiency << '\t' << prompt_acceptance << '\t' << 1.0 - prompt_acceptance << '\n';
      if (truth_efficiency >= minimum_true_efficiency && prompt_acceptance < best_prompt_acceptance)
      {
        best_direction = "less_equal";
        best_cut = cut;
        best_truth_efficiency = truth_efficiency;
        best_prompt_acceptance = prompt_acceptance;
      }
    }
    for (int bin = 0; bin <= bins + 1; ++bin)
    {
      const double truth_efficiency = truth->Integral(bin, bins + 1) / truth_total;
      const double prompt_acceptance = prompt->Integral(bin, bins + 1) / prompt_total;
      const double cut = bin == 0 ? -std::numeric_limits<double>::infinity() :
          bin == bins + 1 ? std::numeric_limits<double>::infinity() : truth->GetXaxis()->GetBinLowEdge(bin);
      if (write_roc)
        roc_output << feature << "\tgreater_equal\t" << cut << '\t' << truth_efficiency << '\t' << prompt_acceptance << '\t'
                   << 1.0 - prompt_acceptance << '\n';
      if (truth_efficiency >= minimum_true_efficiency && prompt_acceptance < best_prompt_acceptance)
      {
        best_direction = "greater_equal";
        best_cut = cut;
        best_truth_efficiency = truth_efficiency;
        best_prompt_acceptance = prompt_acceptance;
      }
    }
    output << feature << '\t' << best_direction << '\t' << best_cut << '\t' << best_truth_efficiency << '\t' << best_prompt_acceptance << '\t'
           << 1.0 - best_prompt_acceptance << '\t' << truth_total << '\t' << prompt_total << '\n';
  }
  return 0;
}
