#!/bin/bash
set -eo pipefail

workflow_dir=$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")
module_dir=$(cd "$workflow_dir/../.." && pwd)
process_id=${1:?usage: workflows/single_particle/run_tree.sh PROCESS_ID [N_EVENTS]}
n_events=${2:-0}
if ! [[ "$process_id" =~ ^[0-9]+$ ]] || ! [[ "$n_events" =~ ^[0-9]+$ ]]; then
  echo "PROCESS_ID and N_EVENTS must be non-negative integers" >&2
  exit 2
fi

printf -v process_tag '%06d' "$process_id"
output_file="$module_dir/output/root/photon_analysis_tree_${process_tag}.root"
if [[ -e "$output_file" ]]; then
  echo "Refusing to overwrite existing tree output: $output_file" >&2
  exit 3
fi

source /opt/sphenix/core/bin/sphenix_setup.sh -n ana
export LD_LIBRARY_PATH="$module_dir/install/lib64:/sphenix/user/ryotaro/DirectPhotonAnalysis/Pi0Reconstruction/install/lib:${LD_LIBRARY_PATH:-}"

root -l -b -q "$workflow_dir/Fun4All_PhotonAnalysisTree.C(${process_id},${n_events})"
root -l -b -q -e ".L $workflow_dir/check_tree.C" -e "gSystem->Exit(check_tree(\"$output_file\",true));"
