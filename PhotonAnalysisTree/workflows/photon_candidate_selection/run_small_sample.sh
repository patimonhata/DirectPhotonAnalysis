#!/usr/bin/env bash
set -euo pipefail

usage()
{
  echo "Usage: $0 SAMPLE_NAME N_SEGMENTS [FILES_PER_MAP] [MIN_CLUSTER_ENERGY_GEV] [OUTPUT_ROOT] [N_EVENTS_PER_MAP] [TAGGING_PARTNER_MIN_ENERGY_GEV]" >&2
}

if (( $# < 2 || $# > 7 )); then
  usage
  exit 2
fi

workflow_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
module_dir=$(cd "$workflow_dir/../.." && pwd)
sample_name=$1
segment_count=$2
files_per_map=${3:-10}
min_cluster_energy=${4:-0.1}
output_root=${5:-}
n_events=${6:-0}
tagging_partner_min_energy=${7:-$min_cluster_energy}

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
if ! [[ "$min_cluster_energy" =~ ^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$ ]]; then
  echo "MIN_CLUSTER_ENERGY_GEV must be a non-negative finite number: $min_cluster_energy" >&2
  exit 2
fi
if ! [[ "$tagging_partner_min_energy" =~ ^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$ ]]; then
  echo "TAGGING_PARTNER_MIN_ENERGY_GEV must be a non-negative finite number: $tagging_partner_min_energy" >&2
  exit 2
fi
printf -v canonical_threshold "%.12g" "$min_cluster_energy"
threshold_tag=${canonical_threshold//./p}
threshold_tag=${threshold_tag//+/}
threshold_tag=${threshold_tag//-/m}
if [[ -z "$output_root" ]]; then
  output_root="$module_dir/output/qa/photon_candidate_selection/cluster_e_gt_${threshold_tag}/${sample_name}_${segment_count}segments"
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
partial_root="$output_root/reduce/$sample_name"
map_count=$(((segment_count + files_per_map - 1) / files_per_map))

echo "Small-sample QA: sample=$sample_name segments=$segment_count maps=$map_count events_per_map=$n_events min_cluster_energy=$min_cluster_energy tagging_partner_min_energy=$tagging_partner_min_energy"
for ((job_index = 0; job_index < map_count; ++job_index)); do
  "$workflow_dir/run_map.sh" "$job_index" 0 "$segment_count" "$files_per_map" "$input_manifest" "$sample_name" "$map_output" "$n_events" "$min_cluster_energy" "$tagging_partner_min_energy"
done

shard_count=1
if [[ "$sample_name" == jet12 ]]; then shard_count=10; fi
if (( map_count < shard_count )); then
  echo "Jet12 QA requires at least 10 maps so every fixed shard is non-empty: maps=$map_count" >&2
  exit 2
fi
for ((shard_index = 0; shard_index < shard_count; ++shard_index)); do
  output_base="$partial_root/shard_${shard_index}/photon_candidate_selection"
  "$workflow_dir/run_reduce.sh" "$family" "$map_root" "$output_base" false 200 40.0 "$sample_name" "$shard_index"
done

echo "Small-sample QA maps: $map_output"
echo "Small-sample QA partials: $partial_root"
