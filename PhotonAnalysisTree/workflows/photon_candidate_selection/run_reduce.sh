#!/usr/bin/env bash
set -euo pipefail

usage()
{
  echo "Usage: $0 FAMILY MAP_ROOT OUTPUT_BASE SELECTION REQUIRE_COMPLETE N_BINS ET_MAX_GEV SAMPLE_NAME SHARD_INDEX" >&2
}

if (( $# != 9 )); then
  usage
  exit 2
fi

workflow_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
family=$1
map_root=$2
output_base=$3
selection=$4
require_complete=$5
n_bins=$6
et_max=$7
sample_name=$8
shard_index=$9

if [[ "$family" != jet && "$family" != photonjet ]]; then
  echo "FAMILY must be jet or photonjet: $family" >&2
  exit 2
fi
if [[ -z "$sample_name" ]]; then
  echo "SAMPLE_NAME must not be empty" >&2
  exit 2
fi
case "$family:$sample_name" in
  jet:jet3|jet:jet5|jet:jet8|jet:jet12|jet:jet20|jet:jet30|jet:jet40|photonjet:photonjet3|photonjet:photonjet5|photonjet:photonjet10|photonjet:photonjet20) ;;
  *)
    echo "SAMPLE_NAME does not belong to FAMILY: $family/$sample_name" >&2
    exit 2
    ;;
esac
case "$selection" in
  kinematic|preselection|preselection_tight|preselection_isolation|region_a|region_a_tagging_veto|final_photon) ;;
  *)
    echo "Unsupported SELECTION: $selection" >&2
    exit 2
    ;;
esac
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
if ! [[ "$shard_index" =~ ^[0-9]+$ ]] || { [[ "$sample_name" == jet12 ]] && (( shard_index > 9 )); } || { [[ "$sample_name" != jet12 ]] && (( shard_index != 0 )); }; then
  echo "SHARD_INDEX must be 0-9 for jet12 and 0 for every other sample: $sample_name/$shard_index" >&2
  exit 2
fi
if [[ ! -d "$map_root" ]]; then
  echo "MAP_ROOT is not a directory: $map_root" >&2
  exit 2
fi
if [[ -z "$output_base" || "$map_root" == *\"* || "$map_root" == *\\* || "$output_base" == *\"* || "$output_base" == *\\* || "$sample_name" == *\"* || "$sample_name" == *\\* ]]; then
  echo "MAP_ROOT and OUTPUT_BASE must be non-empty paths without quotes or backslashes" >&2
  exit 2
fi
map_root=$(cd "$map_root" && pwd)

set +u
source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.565
set -u
root -l -b -q "$workflow_dir/ReducePythiaPhotonCandidateSelection.C+(\"$family\",\"$map_root\",\"$output_base\",\"$selection\",$require_complete,$n_bins,$et_max,\"$sample_name\",$shard_index)"

echo "Photon-candidate selection partial output base: $output_base"
