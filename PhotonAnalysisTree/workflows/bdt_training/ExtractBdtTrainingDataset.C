#include "../../analysis/bdt_training/BdtTrainingClusters.h"

int ExtractBdtTrainingDataset(const std::string jet_map_root, const std::string photonjet_map_root, const std::string output_path)
{
  std::vector<bdt_training::SampleInput> inputs;
  for (const auto& sample : sample_definitions("jet")) inputs.push_back({jet_map_root, sample.name});
  for (const auto& sample : sample_definitions("photonjet")) inputs.push_back({photonjet_map_root, sample.name});
  return bdt_training::write_clusters(inputs, output_path);
}
