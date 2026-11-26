#include "../../analysis/bdt_training/BdtTrainingClusters.h"

int ExtractBdtTrainingClusters(const std::string map_root, const std::string sample_name, const std::string output_path, int shard_index = 0, int shard_count = 1)
{
  return bdt_training::write_clusters({{map_root, sample_name, shard_index, shard_count}}, output_path);
}
