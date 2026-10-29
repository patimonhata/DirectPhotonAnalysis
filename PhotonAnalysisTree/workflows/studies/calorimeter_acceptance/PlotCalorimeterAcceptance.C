#include "../../../analysis/calorimeter_acceptance/PlotAcceptance.C"

void PlotCalorimeterAcceptance(const std::string output)
{
  try { calo_acceptance_plot::plot(output); }
  catch (const std::exception& error) { std::cerr << error.what() << std::endl; gSystem->Exit(1); }
}
