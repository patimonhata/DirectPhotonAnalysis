#!/bin/bash
set -eo pipefail

study_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
process_id=${1:?usage: run_job.sh PROCESS_ID SAMPLE INPUT_DIRECTORY [N_EVENTS]}
sample=${2:?usage: run_job.sh PROCESS_ID SAMPLE INPUT_DIRECTORY [N_EVENTS]}
input_directory=${3:?usage: run_job.sh PROCESS_ID SAMPLE INPUT_DIRECTORY [N_EVENTS]}
n_events=${4:-0}

if ! [[ "$process_id" =~ ^[0-9]+$ ]]; then
  echo "PROCESS_ID must be a non-negative integer" >&2
  exit 2
fi
if [[ "$sample" != "gamma" && "$sample" != "pi0" ]]; then
  echo "SAMPLE must be gamma or pi0" >&2
  exit 2
fi
if ! [[ "$n_events" =~ ^[0-9]+$ ]]; then
  echo "N_EVENTS must be a non-negative integer" >&2
  exit 2
fi
if [[ ! -d "$input_directory" ]]; then
  echo "INPUT_DIRECTORY does not exist: $input_directory" >&2
  exit 2
fi

printf -v process_tag '%06d' "$process_id"
output_file="$study_dir/output/root/topocluster_hcal_${sample}_${process_tag}.root"
if [[ -e "$output_file" ]]; then
  echo "Refusing to overwrite existing output: $output_file" >&2
  exit 3
fi

source /opt/sphenix/core/bin/sphenix_setup.sh -n "${SPHENIX_RELEASE:-ana}"
export LD_LIBRARY_PATH="$study_dir/install/lib:$study_dir/install/lib64:${LD_LIBRARY_PATH:-}"
export ROOT_INCLUDE_PATH="$study_dir/install/include:${ROOT_INCLUDE_PATH:-}"

root -l -b -q \
  "$study_dir/macro/Fun4All_TopoClusterHCalStudy.C(${process_id},${n_events},\"${sample}\",\"${input_directory}\",\"${study_dir}/output/root\")"
