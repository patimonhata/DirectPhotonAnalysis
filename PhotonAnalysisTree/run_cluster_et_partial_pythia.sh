#!/bin/bash
set -eo pipefail

module_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
usage="usage: run_cluster_et_partial_pythia.sh JOB_INDEX CHUNK_OFFSET TOTAL_FILES FILES_PER_JOB INPUT_MANIFEST OUTPUT_DIRECTORY [N_BINS] [ET_MAX] [TRUTH_ETA_MAX] [CLUSTER_ETA_MAX] [MIN_CLUSTER_ENERGY] [DOMINANT_FRACTION_MIN] [PI0_CONTRIBUTOR_FRACTION_MIN] [SEPARATED_DR] [MERGED_DR] [RESPONSE_MIN] [RESPONSE_MAX]"
job_index=${1:?$usage}
chunk_offset=${2:?$usage}
total_files=${3:?$usage}
files_per_job=${4:?$usage}
input_manifest=${5:?$usage}
output_directory=${6:?$usage}
n_bins=${7:-100}
et_max=${8:-20.0}
truth_eta_max=${9:-0.7}
cluster_eta_max=${10:-0.7}
min_cluster_energy=${11:-0.2}
dominant_fraction_min=${12:-0.5}
pi0_contributor_fraction_min=${13:-0.5}
separated_dr=${14:-0.03}
merged_dr=${15:-0.06}
response_min=${16:-0.5}
response_max=${17:-1.5}

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
temporary_output=$(mktemp "$output_directory/.cluster-et-partial.XXXXXX")
cleanup()
{
  if [[ -n "${temporary_output:-}" && -e "$temporary_output" ]]; then
    rm -f -- "$temporary_output"
  fi
}
trap cleanup EXIT

source /opt/sphenix/core/bin/sphenix_setup.sh -n ana
export LD_LIBRARY_PATH="$module_dir/install/lib64:/sphenix/user/ryotaro/DirectPhotonAnalysis/Pi0Reconstruction/install/lib:${LD_LIBRARY_PATH:-}"
root -l -b -q \
  "$module_dir/macro/Fun4All_PythiaClusterEtSpectra.C(\"${input_manifest}\",${manifest_begin},${manifest_end},\"${temporary_output}\",${n_bins},${et_max},${truth_eta_max},${cluster_eta_max},${min_cluster_energy},${dominant_fraction_min},${pi0_contributor_fraction_min},${separated_dr},${merged_dr},${response_min},${response_max})"
root -l -b -q \
  "$module_dir/macro/check_pythia_cluster_et_partial.C(\"${temporary_output}\")"

mv -- "$temporary_output" "$final_output"
temporary_output=
trap - EXIT
echo "Pythia cluster ET partial ROOT file: $final_output"
