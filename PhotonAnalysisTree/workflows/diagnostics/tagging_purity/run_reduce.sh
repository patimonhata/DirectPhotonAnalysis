#!/usr/bin/env bash
set -euo pipefail
if (( $# != 8 )); then
  echo "Usage: $0 MAP_ROOT OUTPUT_FILE EXPECTED_THRESHOLD SAMPLE SHARD_INDEX SHARD_COUNT REQUIRE_COMPLETE MAX_EVENTS" >&2
  exit 2
fi
workflow_dir=$(cd "$(dirname "$0")" && pwd)
map_root=$1
output_file=$2
threshold=$3
sample=$4
shard_index=$5
shard_count=$6
require_complete=$7
max_events=$8
if [[ ! -d "$map_root" || -e "$output_file" || "$map_root$output_file$sample" == *\"* || "$map_root$output_file$sample" == *\\* ]]; then
  echo "Invalid input/output path or output already exists" >&2
  exit 2
fi
if [[ "$require_complete" != true && "$require_complete" != false ]]; then echo "REQUIRE_COMPLETE must be true or false" >&2; exit 2; fi
if ! [[ "$threshold" =~ ^[0-9]+([.][0-9]+)?$ && "$shard_index" =~ ^[0-9]+$ && "$shard_count" =~ ^[1-9][0-9]*$ && "$max_events" =~ ^[0-9]+$ ]]; then
  echo "Invalid numeric argument" >&2
  exit 2
fi
map_root=$(cd "$map_root" && pwd)
mkdir -p "$(dirname "$output_file")"
temporary="$output_file.tmp.$$"
trap 'rm -f "$temporary"' EXIT
set +u
source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.565
set -u
root -l -b -q -e ".L $workflow_dir/ReduceTaggingPurityDiagnostic.C" \
  -e "gSystem->Exit(ReduceTaggingPurityDiagnostic(\"$map_root\",\"$temporary\",$threshold,\"$sample\",$shard_index,$shard_count,$require_complete,$max_events))"
mv "$temporary" "$output_file"
trap - EXIT
