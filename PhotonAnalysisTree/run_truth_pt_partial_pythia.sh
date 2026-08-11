#!/bin/bash
set -eo pipefail

module_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
usage="usage: run_truth_pt_partial_pythia.sh JOB_INDEX CHUNK_OFFSET TOTAL_FILES FILES_PER_JOB INPUT_MANIFEST TREE_INPUT_DIRECTORY OUTPUT_DIRECTORY [N_BINS] [PT_MAX] [MAX_ABS_ETA] [USE_EVENT_WEIGHT]"
job_index=${1:?$usage}
chunk_offset=${2:?$usage}
total_files=${3:?$usage}
files_per_job=${4:?$usage}
input_manifest=${5:?$usage}
tree_input_directory=${6:?$usage}
output_directory=${7:?$usage}
n_bins=${8:-100}
pt_max=${9:-20.0}
max_abs_eta=${10:-0.7}
use_event_weight=${11:-false}

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
chunk_index=$((job_index + chunk_offset))
if [[ "$use_event_weight" != "true" && "$use_event_weight" != "false" &&
      "$use_event_weight" != "0" && "$use_event_weight" != "1" ]]; then
  echo "USE_EVENT_WEIGHT must be true, false, 0, or 1" >&2
  exit 2
fi
if [[ ! -r "$input_manifest" ]]; then
  echo "Input manifest is not readable: $input_manifest" >&2
  exit 2
fi
if [[ ! -d "$tree_input_directory" ]]; then
  echo "Tree input directory does not exist: $tree_input_directory" >&2
  exit 2
fi

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
if [[ ! -w "$output_directory" ]]; then
  echo "Output directory is not writable: $output_directory" >&2
  exit 3
fi

temporary_output=$(mktemp "$output_directory/.truth-pt-partial.XXXXXX")
cleanup()
{
  if [[ -n "${temporary_output:-}" && -e "$temporary_output" ]]; then
    rm -f -- "$temporary_output"
  fi
}
trap cleanup EXIT

source /opt/sphenix/core/bin/sphenix_setup.sh -n ana
root -l -b -q \
  "$module_dir/macro/AccumulatePythiaTruthPtSpectra.C(\"${input_manifest}\",\"${tree_input_directory}\",${manifest_begin},${manifest_end},\"${temporary_output}\",${n_bins},${pt_max},${max_abs_eta},${use_event_weight})"
root -l -b -q \
  "$module_dir/macro/check_pythia_truth_pt_partial.C(\"${temporary_output}\")"

mv -- "$temporary_output" "$final_output"
temporary_output=
trap - EXIT
echo "Pythia truth pT partial ROOT file: $final_output"
