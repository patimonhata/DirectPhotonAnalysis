#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  run_add_bdt.sh INPUT.root OUTPUT.root [MODEL.root]

Adds the PPG15/PPG12 base_v3E SPLIT photon-ID BDT score to event_tree.
The input file is never modified.
EOF
}

if [[ $# -lt 2 || $# -gt 3 ]]; then
  usage >&2
  exit 2
fi

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
input_file=$1
output_file=$2
model_file=${3:-/sphenix/user/shuhangli/ppg12/FunWithxgboost/binned_models/model_base_v3E_split_single_tmva.root}

if [[ ! -r "$input_file" ]]; then
  echo "Input file is not readable: $input_file" >&2
  exit 3
fi
if [[ ! -r "$model_file" ]]; then
  echo "Model file is not readable: $model_file" >&2
  exit 4
fi
if [[ $(realpath -m -- "$input_file") == $(realpath -m -- "$output_file") ]]; then
  echo "Input and output must be different files." >&2
  exit 5
fi
if [[ -e "$output_file" ]]; then
  echo "Output already exists; remove it or choose a new path: $output_file" >&2
  exit 6
fi

mkdir -p -- "$(dirname -- "$output_file")"
root_call="${script_dir}/macros/add_ppg15_bdt_to_pi0_tree.C(\"${input_file}\",\"${output_file}\",\"${model_file}\")"
root -l -b -q "$root_call"

if [[ ! -s "$output_file" ]]; then
  echo "Output file was not created: $output_file" >&2
  exit 7
fi

check_call="${script_dir}/macros/check_bdt_output.C(\"${output_file}\")"
root -l -b -q "$check_call"
