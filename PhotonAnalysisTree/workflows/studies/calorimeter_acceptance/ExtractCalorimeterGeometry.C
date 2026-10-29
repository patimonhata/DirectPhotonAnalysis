#include "../../../analysis/calorimeter_acceptance/ExtractTowerGeometry.C"

void ExtractCalorimeterGeometry(const std::string output)
{
  try { calo_acceptance::extract(output); }
  catch (const std::exception& error) { std::cerr << error.what() << std::endl; gSystem->Exit(1); }
}
