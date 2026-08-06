#!/bin/bash
set -eo pipefail

study_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
root_dir="$study_dir/output/root"
merge_dir="$study_dir/output/merge"

mkdir -p "$merge_dir"

gamma_files=("$root_dir"/topocluster_hcal_gamma_*.root)
pi0_files=("$root_dir"/topocluster_hcal_pi0_*.root)
if ((${#gamma_files[@]} != 500 || ${#pi0_files[@]} != 500)); then
  echo "Expected 500 gamma and 500 pi0 files; found ${#gamma_files[@]} and ${#pi0_files[@]}" >&2
  exit 2
fi

source /opt/sphenix/core/bin/sphenix_setup.sh -n "${SPHENIX_RELEASE:-ana}"

gamma_output="$merge_dir/topocluster_hcal_gamma_merged.root"
pi0_output="$merge_dir/topocluster_hcal_pi0_merged.root"
combined_output="$merge_dir/topocluster_hcal_gamma_pi0_merged.root"
for output in "$gamma_output" "$pi0_output" "$combined_output"; do
  if [[ -e "$output" ]]; then
    echo "Refusing to overwrite existing merge output: $output" >&2
    exit 3
  fi
done

hadd "$gamma_output" "${gamma_files[@]}"
hadd "$pi0_output" "${pi0_files[@]}"
hadd "$combined_output" "$gamma_output" "$pi0_output"
