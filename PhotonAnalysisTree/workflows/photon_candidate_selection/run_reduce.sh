#!/usr/bin/env bash
set -euo pipefail

usage()
{
  echo "Usage: $0 FAMILY MAP_ROOT OUTPUT_BASE REQUIRE_COMPLETE [N_BINS] [ET_MAX_GEV]" >&2
}

if (( $# < 4 || $# > 6 )); then
  usage
  exit 2
fi

workflow_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
family=$1
map_root=$2
output_base=$3
require_complete=$4
n_bins=${5:-200}
et_max=${6:-40.0}

if [[ "$family" != jet && "$family" != photonjet ]]; then
  echo "FAMILY must be jet or photonjet: $family" >&2
  exit 2
fi
if [[ "$require_complete" != true && "$require_complete" != false ]]; then
  echo "REQUIRE_COMPLETE must be true or false: $require_complete" >&2
  exit 2
fi
if ! [[ "$n_bins" =~ ^[1-9][0-9]*$ ]]; then
  echo "N_BINS must be a positive integer: $n_bins" >&2
  exit 2
fi
if ! [[ "$et_max" =~ ^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$ ]] || [[ "$et_max" =~ ^0*([.]0*)?$ ]]; then
  echo "ET_MAX_GEV must be a positive finite number: $et_max" >&2
  exit 2
fi
if [[ ! -d "$map_root" ]]; then
  echo "MAP_ROOT is not a directory: $map_root" >&2
  exit 2
fi
if [[ -z "$output_base" || "$map_root" == *\"* || "$map_root" == *\\* || "$output_base" == *\"* || "$output_base" == *\\* ]]; then
  echo "MAP_ROOT and OUTPUT_BASE must be non-empty paths without quotes or backslashes" >&2
  exit 2
fi
map_root=$(cd "$map_root" && pwd)

set +u
source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.565
set -u
root -l -b -q \
  "$workflow_dir/ReducePythiaRegionAAnchorTopology.C+(\"$family\",\"$map_root\",\"$output_base\",$require_complete,$n_bins,$et_max)"

echo "Photon-candidate reduce output base: $output_base"
