#!/usr/bin/env bash
set -euo pipefail
if (( $# != 8 )); then
  echo "Usage: $0 LOW_MAP_ROOT HIGH_MAP_ROOT OUTPUT_FILE SAMPLE SHARD_INDEX SHARD_COUNT REQUIRE_COMPLETE MAX_EVENTS" >&2
  exit 2
fi
workflow_dir=$(cd "$(dirname "$0")" && pwd)
low_map_root=$1
high_map_root=$2
output_file=$3
sample=$4
shard_index=$5
shard_count=$6
require_complete=$7
max_events=$8
if [[ ! -d "$low_map_root" || ! -d "$high_map_root" || -e "$output_file" || "$low_map_root$high_map_root$output_file$sample" == *\"* ||
      "$low_map_root$high_map_root$output_file$sample" == *\\* ]]; then
  echo "Invalid input/output path or output already exists" >&2
  exit 2
fi
if [[ "$require_complete" != true && "$require_complete" != false ]]; then echo "REQUIRE_COMPLETE must be true or false" >&2; exit 2; fi
if ! [[ "$shard_index" =~ ^[0-9]+$ && "$shard_count" =~ ^[1-9][0-9]*$ && "$max_events" =~ ^[0-9]+$ ]]; then
  echo "Invalid numeric argument" >&2
  exit 2
fi
low_map_root=$(cd "$low_map_root" && pwd)
high_map_root=$(cd "$high_map_root" && pwd)
mkdir -p "$(dirname "$output_file")"
temporary="$output_file.tmp.$$"
trap 'rm -f "$temporary"' EXIT
set +u
source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.565
set -u
root -l -b -q -e ".L $workflow_dir/ReduceHybridThresholdDiagnostic.C" \
  -e "gSystem->Exit(ReduceHybridThresholdDiagnostic(\"$low_map_root\",\"$high_map_root\",\"$temporary\",\"$sample\",$shard_index,$shard_count,$require_complete,$max_events))"
mv "$temporary" "$output_file"
trap - EXIT
