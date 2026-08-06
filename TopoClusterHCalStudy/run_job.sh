#!/bin/bash
set -eo pipefail

study_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
job_index=${1:?usage: run_job.sh JOB_INDEX [N_EVENTS]}
n_events=${2:-0}

if ! [[ "$job_index" =~ ^[0-9]+$ ]] || ((job_index >= 1000)); then
  echo "JOB_INDEX must be an integer in [0,999]" >&2
  exit 2
fi
if ! [[ "$n_events" =~ ^[0-9]+$ ]]; then
  echo "N_EVENTS must be a non-negative integer" >&2
  exit 2
fi

if ((job_index < 500)); then
  sample=gamma
  process_id=$job_index
else
  sample=pi0
  process_id=$((job_index - 500))
fi
printf -v process_tag '%06d' "$process_id"

output_file="$study_dir/output/root/topocluster_hcal_${sample}_${process_tag}.root"
if [[ -e "$output_file" ]]; then
  echo "Refusing to overwrite existing output: $output_file" >&2
  exit 3
fi

source /opt/sphenix/core/bin/sphenix_setup.sh -n "${SPHENIX_RELEASE:-ana}"
export LD_LIBRARY_PATH="$study_dir/install/lib64:$study_dir/install/lib:${LD_LIBRARY_PATH:-}"
export ROOT_INCLUDE_PATH="$study_dir/install/include:${ROOT_INCLUDE_PATH:-}"

root -l -b -q \
  "$study_dir/macro/Fun4All_TopoClusterHCalStudy.C(${job_index},${n_events})"
