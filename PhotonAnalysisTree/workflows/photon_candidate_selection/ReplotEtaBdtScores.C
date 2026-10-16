#include "ReplotPi0BdtScores.C"

int ReplotEtaBdtScores(const std::string output_base)
{
  const std::array<const char*, 6> selections = {"kinematic", "preselection", "preselection_tight", "preselection_isolation", "region_a", "region_a_tagging_veto"};
  for (std::size_t i = 0; i < selections.size(); ++i)
  {
    const std::string directory = output_base + "/" + selections[i];
    TFile input((directory + "/eta_bdt_score.root").c_str(), "READ");
    auto* selection = input.Get<TNamed>("selection");
    pi0_bdt::Histograms scores;
    if (input.IsZombie() || !selection || !scores.add(input) || !scores.valid()) return 1;
    const std::string label = selection->GetTitle();
    const auto caption = pi0_bdt::read_plot_caption(input, label);
    if (!scores.export_plots(directory, label, i == 3 || i >= 4, caption, "eta")) return 2;
  }
  return 0;
}
