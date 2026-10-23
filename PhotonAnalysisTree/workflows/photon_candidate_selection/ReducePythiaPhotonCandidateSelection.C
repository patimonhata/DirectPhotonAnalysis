#include "../../analysis/photon_candidate_selection/ReduceCandidateSelection.C"

int ReducePythiaPhotonCandidateSelection(
    const std::string family,
    const std::string map_root,
    const std::string output_base,
    const bool require_complete,
    const int n_bins,
    const double et_max,
    const std::string sample_name,
    const int shard_index)
{
  return photon_candidate_analysis::reduce(family, map_root, output_base, require_complete, n_bins, et_max, sample_name, shard_index);
}
