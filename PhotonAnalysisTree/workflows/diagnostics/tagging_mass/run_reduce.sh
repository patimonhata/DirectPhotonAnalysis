#!/usr/bin/env bash
set -euo pipefail

if (( $# != 7 )); then
  echo "Usage: $0 FAMILY MAP_ROOT OUTPUT_BASE REQUIRE_COMPLETE SAMPLE_NAME SHARD_INDEX SHARD_COUNT" >&2
  exit 2
fi
workflow_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
family=$1
map_root=$2
output_base=$3
require_complete=$4
sample_name=$5
shard_index=$6
shard_count=$7
if [[ "$family" != jet && "$family" != photonjet ]]; then echo "FAMILY must be jet or photonjet" >&2; exit 2; fi
if [[ "$require_complete" != true && "$require_complete" != false ]]; then echo "REQUIRE_COMPLETE must be true or false" >&2; exit 2; fi
if ! [[ "$shard_index" =~ ^[0-9]+$ && "$shard_count" =~ ^[1-9][0-9]*$ ]] || (( shard_index >= shard_count )); then
  echo "Require 0 <= SHARD_INDEX < SHARD_COUNT" >&2
  exit 2
fi
case "$family:$sample_name" in
  jet:jet3|jet:jet5|jet:jet8|jet:jet12|jet:jet20|jet:jet30|jet:jet40|photonjet:photonjet3|photonjet:photonjet5|photonjet:photonjet10|photonjet:photonjet20) ;;
  *) echo "SAMPLE_NAME does not belong to FAMILY: $family/$sample_name" >&2; exit 2 ;;
esac
if [[ ! -d "$map_root" || -z "$output_base" || "$map_root$output_base$sample_name" == *\"* || "$map_root$output_base$sample_name" == *\\* ]]; then
  echo "Invalid input/output path" >&2
  exit 2
fi
map_root=$(cd "$map_root" && pwd)
set +u
source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.565
set -u
root -l -b -q -e ".L $workflow_dir/ReduceTaggingMassDiagnostic.C" \
  -e "gSystem->Exit(ReduceTaggingMassDiagnostic(\"$family\",\"$map_root\",\"$output_base\",$require_complete,\"$sample_name\",$shard_index,$shard_count,0))"
