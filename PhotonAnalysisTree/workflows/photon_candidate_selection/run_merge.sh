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

set +u
source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.565
set -u
root -l -b -q "$workflow_dir/MergePythiaPhotonCandidateSelection.C(\"$family\",\"$partial_root\",\"$output_base\")"

echo "Merged selection output base: $output_base"
