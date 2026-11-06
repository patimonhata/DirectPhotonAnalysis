#!/usr/bin/env bash
set -euo pipefail

usage()
{
  echo "Usage: $0 FAMILY PARTIAL_ROOT OUTPUT_BASE" >&2
}

if (( $# != 3 )); then
  usage
  exit 2
fi

workflow_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
module_dir=$(cd "$workflow_dir/../.." && pwd)
source "$module_dir/workflows/common/runtime.sh"
photon_dispatch "$module_dir" "workflows/photon_candidate_selection/run_merge.sh" "$@"
family=$1
partial_root=$2
output_base=$3

if [[ "$family" != jet && "$family" != photonjet ]]; then
  echo "FAMILY must be jet or photonjet: $family" >&2
  exit 2
fi
if [[ ! -d "$partial_root" ]]; then
  echo "PARTIAL_ROOT is not a directory: $partial_root" >&2
  exit 2
fi
if [[ -z "$output_base" || "$partial_root" == *\"* || "$partial_root" == *\\* ||
      "$output_base" == *\"* || "$output_base" == *\\* ]]; then
  echo "PARTIAL_ROOT and OUTPUT_BASE must be non-empty paths without quotes or backslashes" >&2
  exit 2
fi
partial_root=$(cd "$partial_root" && pwd)

photon_setup_runtime "$module_dir"
photon_check_output "$output_base/candidate_composition.root" --stage merge
photon_check_inputs reduce "$partial_root/*/shard_*/photon_candidate_selection.root"
root -l -b -q "$workflow_dir/MergePythiaPhotonCandidateSelection.C(\"$family\",\"$partial_root\",\"$output_base\")"

while IFS= read -r -d '' result; do
  photon_record_execution merge "$result" "$@"
done < <(find "$output_base" -type f -name '*.root' -print0)

echo "Merged selection output base: $output_base"
