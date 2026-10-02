#include "Pi0BdtHistograms.h"
#include <iostream>

// Redraw only the two BDT PDFs in each selection directory; retain ROOT/CSV metadata and all other plots.
int ReplotPi0BdtScores(const std::string output_base)
{
  const std::array<const char*, 6> selections = {"kinematic", "preselection", "preselection_tight", "preselection_isolation", "region_a", "region_a_tagging_veto"};
  std::array<std::unique_ptr<pi0_bdt::Histograms>, 6> histograms;
  std::array<std::string, 6> labels;
  for (std::size_t i = 0; i < selections.size(); ++i)
  {
    const std::string path = output_base + "/" + selections[i] + "/pi0_bdt_score.root";
    TFile input(path.c_str(), "READ");
    auto* selection = input.Get<TNamed>("selection");
    histograms[i] = std::make_unique<pi0_bdt::Histograms>();
    if (input.IsZombie() || !selection || !histograms[i]->add(input) || !histograms[i]->valid())
    {
      std::cerr << "Missing or incompatible BDT histograms: " << path << std::endl;
      return 1;
    }
    labels[i] = selection->GetTitle();
  }
  for (std::size_t i = 0; i < selections.size(); ++i)
    if (!histograms[i]->export_plots(output_base + "/" + selections[i], labels[i], i >= 3)) return 2;
  return 0;
}
