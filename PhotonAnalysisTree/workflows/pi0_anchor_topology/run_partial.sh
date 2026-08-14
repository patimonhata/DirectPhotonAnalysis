#!/bin/bash
set -eo pipefail

workflow_dir=$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")
module_dir=$(cd "$workflow_dir/../.." && pwd)
usage="usage: workflows/pi0_anchor_topology/run_partial.sh JOB_INDEX CHUNK_OFFSET TOTAL_FILES FILES_PER_JOB INPUT_MANIFEST OUTPUT_DIRECTORY [N_BINS] [ET_MAX] [TRUTH_ETA_MAX] [ANCHOR_CLUSTER_ETA_MAX] [PARTNER_CLUSTER_ETA_MAX] [MIN_CLUSTER_ENERGY] [DOMINANT_FRACTION_MIN] [ANCHOR_PI0_FRACTION_MIN] [MIN_ENERGY_CONTRIBUTION_FRACTION]"
job_index=${1:?$usage}
chunk_offset=${2:?$usage}
total_files=${3:?$usage}
files_per_job=${4:?$usage}
input_manifest=${5:?$usage}
output_directory=${6:?$usage}
n_bins=${7:-100}
et_max=${8:-20.0}
truth_eta_max=${9:-0.7}
anchor_cluster_eta_max=${10:-0.7}
partner_cluster_eta_max=${11:--1.0}
min_cluster_energy=${12:-0.2}
dominant_fraction_min=${13:-0.5}
anchor_pi0_fraction_min=${14:-0.5}
min_energy_contribution_fraction=${15:-0.3}

for value in "$job_index" "$chunk_offset" "$total_files" "$files_per_job" "$n_bins"; do
  if ! [[ "$value" =~ ^[0-9]+$ ]]; then
    echo "Chunk, file counts, and N_BINS must be non-negative integers: $value" >&2
    exit 2
  fi
done
if (( total_files == 0 || files_per_job == 0 || n_bins == 0 )); then
  echo "TOTAL_FILES, FILES_PER_JOB, and N_BINS must be positive" >&2
  exit 2
fi
if [[ ! -r "$input_manifest" ]]; then
  echo "Input manifest is not readable: $input_manifest" >&2
  exit 2
fi

chunk_index=$((job_index + chunk_offset))
manifest_begin=$((chunk_index * files_per_job))
if (( manifest_begin >= total_files )); then
  echo "Chunk $chunk_index begins beyond TOTAL_FILES=$total_files" >&2
  exit 2
fi
manifest_end=$((manifest_begin + files_per_job))
if (( manifest_end > total_files )); then
  manifest_end=$total_files
fi
printf -v chunk_tag "%06d" "$chunk_index"
final_output="$output_directory/partial_${chunk_tag}.root"
if [[ -e "$final_output" || -L "$final_output" ]]; then
  echo "Refusing to overwrite existing partial output: $final_output" >&2
  exit 3
fi
mkdir -p -- "$output_directory"
temporary_output=$(mktemp "$output_directory/.pi0-anchor-topology-partial.XXXXXX")
cleanup()
{
  if [[ -n "${temporary_output:-}" && -e "$temporary_output" ]]; then
    rm -f -- "$temporary_output"
  fi
}
trap cleanup EXIT

source /opt/sphenix/core/bin/sphenix_setup.sh -n ana
export LD_LIBRARY_PATH="$module_dir/install/lib64:/sphenix/user/ryotaro/DirectPhotonAnalysis/Pi0Reconstruction/install/lib:${LD_LIBRARY_PATH:-}"
root -l -b -q   "$workflow_dir/Fun4All_PythiaPi0AnchorClusterSpectra.C(\"${input_manifest}\",${manifest_begin},${manifest_end},\"${temporary_output}\",${n_bins},${et_max},${truth_eta_max},${anchor_cluster_eta_max},${partner_cluster_eta_max},${min_cluster_energy},${dominant_fraction_min},${anchor_pi0_fraction_min},${min_energy_contribution_fraction})"
root -l -b -q   "$workflow_dir/check_pythia_pi0_anchor_cluster_partial.C(\"${temporary_output}\")"

mv -- "$temporary_output" "$final_output"
temporary_output=
trap - EXIT
echo "Pythia pi0 anchor-topology partial ROOT file: $final_output"
