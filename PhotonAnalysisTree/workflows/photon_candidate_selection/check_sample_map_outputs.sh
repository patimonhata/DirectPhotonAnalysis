#!/usr/bin/env bash
set -euo pipefail

workflow_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
module_dir=$(cd "$workflow_dir/../.." && pwd)
usage="usage: workflows/photon_candidate_selection/check_sample_map_outputs.sh SAMPLE_NAME [FILES_PER_JOB] [DEEP_VALIDATION] [OUTPUT_DIRECTORY] [MIN_CLUSTER_ENERGY_GEV] [PI0_PARTNER_MIN_ENERGY_GEV] [ETA_PARTNER_MIN_ENERGY_GEV] [PI0_MASS_MIN] [PI0_MASS_MAX] [ETA_MASS_MIN] [ETA_MASS_MAX] [MISSING_ENERGY_MIN] [MISSING_ENERGY_MAX]"
sample_name=${1:?$usage}
files_per_job=${2:-10}
deep_validation=${3:-false}
output_directory=${4:-}

min_cluster_energy=${5:-0.1}
tagging_partner_min_energy=${6:-$min_cluster_energy}
eta_partner_min_energy=${7:-$tagging_partner_min_energy}
pi0_mass_min=${8:-0.10}
pi0_mass_max=${9:-0.20}
eta_mass_min=${10:-0.45}
eta_mass_max=${11:-0.65}
missing_energy_min=${12:-0.2}
missing_energy_max=${13:-0.5}
for value in "$eta_partner_min_energy" "$pi0_mass_min" "$pi0_mass_max" "$eta_mass_min" "$eta_mass_max" "$missing_energy_min" "$missing_energy_max"; do
  if ! [[ "$value" =~ ^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$ ]]; then
    echo "Invalid numeric selection setting: $value" >&2
    exit 2
  fi
done
case "$sample_name" in
  photonjet3|photonjet5|photonjet10|photonjet20|jet3|jet5|jet8|jet12|jet20|jet30|jet40) ;;
  *)
    echo "Unsupported SAMPLE_NAME: $sample_name" >&2
    exit 2
    ;;
esac
if ! [[ "$files_per_job" =~ ^[1-9][0-9]*$ ]]; then
  echo "FILES_PER_JOB must be a positive integer: $files_per_job" >&2
  exit 2
fi
if [[ "$deep_validation" != true && "$deep_validation" != false ]]; then
  echo "DEEP_VALIDATION must be true or false" >&2
  exit 2
fi

input_manifest="$module_dir/input/$sample_name/segments.list"
if ! [[ "$min_cluster_energy" =~ ^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$ ]]; then
  echo "MIN_CLUSTER_ENERGY_GEV must be a non-negative finite number: $min_cluster_energy" >&2
  exit 2
fi
if ! [[ "$tagging_partner_min_energy" =~ ^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$ ]]; then
  echo "TAGGING_PARTNER_MIN_ENERGY_GEV must be a non-negative finite number: $tagging_partner_min_energy" >&2
  exit 2
fi
if [[ -z "$output_directory" ]]; then
  output_directory="$module_dir/output/intermediate_files/photon_candidate_selection/cluster_e_gt_0p1/$sample_name"
fi
if [[ ! -r "$input_manifest" ]]; then
  echo "Input manifest is not readable: $input_manifest" >&2
  exit 2
fi
total_files=$(wc -l < "$input_manifest")
expected_maps=$(((total_files + files_per_job - 1) / files_per_job))
missing=0
invalid=0

if [[ "$deep_validation" == true ]]; then
  set +u
  source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.565
  set -u
fi
for ((chunk = 0; chunk < expected_maps; ++chunk)); do
  printf -v chunk_tag "%06d" "$chunk"
  path="$output_directory/map_${chunk_tag}.root"
  if [[ ! -s "$path" ]]; then
    if (( missing < 20 )); then
      echo "MISSING $path"
    fi
    ((missing += 1))
    continue
  fi
  if [[ "$deep_validation" == true ]] && ! root -l -b -q "$workflow_dir/check_pythia_photon_candidate_map.C(\"$path\",${min_cluster_energy},${tagging_partner_min_energy},${eta_partner_min_energy},${pi0_mass_min},${pi0_mass_max},${eta_mass_min},${eta_mass_max},${missing_energy_min},${missing_energy_max})"; then
    echo "INVALID $path"
    ((invalid += 1))
  fi
done
if (( missing > 20 )); then
  echo "MISSING ... and $((missing - 20)) more files"
fi

echo "check_sample_map_outputs - sample/topology-min-E/tagging-min-E/files/files_per_map/expected/missing/invalid = " \
  "$sample_name/$min_cluster_energy/$tagging_partner_min_energy/$total_files/$files_per_job/$expected_maps/$missing/$invalid"
if (( missing != 0 || invalid != 0 )); then
  exit 1
fi
