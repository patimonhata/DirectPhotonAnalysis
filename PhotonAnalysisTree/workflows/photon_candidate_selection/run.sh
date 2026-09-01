#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 || $# -gt 4 ]]; then
  echo "Usage: $0 INPUT_SUFFIX SAMPLE_NAME [N_EVENTS] [OUTPUT_DIRECTORY]" >&2
  exit 2
fi

workflow_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
input_suffix=$1
sample_name=$2
n_events=${3:-0}
output_directory=${4:-/sphenix/user/ryotaro/DirectPhotonAnalysis/PhotonAnalysisTree/output/intermediate_files/photon_candidate_selection/${sample_name}}
model_file=/sphenix/user/ryotaro/DirectPhotonAnalysis/PhotonAnalysisTree/models/model_ppg15_nominal_base_v3E_split_3to35_allplus3jet40_ppg12split_single_tmva.root

set +u
source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.565
set -u
root -l -b -q "${workflow_dir}/Fun4All_PythiaPhotonCandidateTree.C(\"${input_suffix}\",\"${sample_name}\",${n_events},\"${output_directory}\",\"${model_file}\")"
