#include "../../analysis/photon_candidate_selection/MergeCandidateSelection.C"

int MergePythiaPhotonCandidateSelection(
    const std::string family,
    const std::string partial_root,
    const std::string output_base)
{
  return photon_candidate_analysis::merge(family, partial_root, output_base);
}
