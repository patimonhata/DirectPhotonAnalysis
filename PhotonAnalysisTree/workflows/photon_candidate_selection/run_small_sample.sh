#!/usr/bin/env bash
set -euo pipefail

usage()
{
  echo "Usage: $0 SAMPLE_NAME N_SEGMENTS [FILES_PER_MAP] [OUTPUT_ROOT] [N_EVENTS_PER_MAP]" >&2
}

if (( $# < 2 || $# > 5 )); then
  usage
  exit 2
fi

workflow_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
module_dir=$(cd "$workflow_dir/../.." && pwd)
sample_name=$1
segment_count=$2
files_per_map=${3:-10}
output_root=${4:-$module_dir/output/qa/photon_candidate_selection/${sample_name}_${segment_count}segments}
n_events=${5:-0}

for value in "$segment_count" "$files_per_map" "$n_events"; do
  if ! [[ "$value" =~ ^[0-9]+$ ]]; then
    echo "N_SEGMENTS, FILES_PER_MAP, and N_EVENTS_PER_MAP must be non-negative integers: $value" >&2
    exit 2
  fi
done
if (( segment_count == 0 || files_per_map == 0 )); then
  echo "N_SEGMENTS and FILES_PER_MAP must be positive" >&2
  exit 2
fi

case "$sample_name" in
  jet3|jet5|jet8|jet12|jet20|jet30|jet40)
    family=jet
    ;;
  photonjet3|photonjet5|photonjet10|photonjet20)
    family=photonjet
    ;;
  *)
    echo "Unsupported SAMPLE_NAME: $sample_name" >&2
    exit 2
    ;;
esac

input_manifest="$module_dir/input/$sample_name/segments.list"
if [[ ! -r "$input_manifest" ]]; then
  echo "Input manifest is not readable: $input_manifest" >&2
  exit 2
fi
manifest_rows=$(wc -l < "$input_manifest")
if (( segment_count > manifest_rows )); then
  echo "Requested $segment_count segments, but $input_manifest has only $manifest_rows rows" >&2
  exit 2
fi

mkdir -p -- "$output_root"
output_root=$(cd "$output_root" && pwd)
map_root="$output_root/maps"
map_output="$map_root/$sample_name"
output_base="$output_root/plots/region_a_pi0_anchor_topology"
map_count=$(((segment_count + files_per_map - 1) / files_per_map))

echo "Small-sample QA: sample=$sample_name segments=$segment_count maps=$map_count events_per_map=$n_events"
for ((job_index = 0; job_index < map_count; ++job_index)); do
  "$workflow_dir/run_map.sh" "$job_index" 0 "$segment_count" "$files_per_map" "$input_manifest" "$sample_name" "$map_output" "$n_events"
done

set +u
source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.565
set -u
root -l -b -q \
  "$workflow_dir/ReducePythiaRegionAAnchorTopology.C(\"$family\",\"$map_root\",\"$output_base\",false,200,40.0)"

echo "Small-sample QA maps: $map_output"
echo "Small-sample QA plots: $output_root/plots"
