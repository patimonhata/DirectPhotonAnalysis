#!/bin/bash
set -eo pipefail

workflow_dir=$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")
module_dir=$(cd "$workflow_dir/../.." && pwd)
source "$module_dir/workflows/common/runtime.sh"
photon_dispatch "$module_dir" "workflows/single_particle/run_tree.sh" "$@"
process_id=${1:?usage: workflows/single_particle/run_tree.sh PROCESS_ID [N_EVENTS]}
n_events=${2:-0}
if ! [[ "$process_id" =~ ^[0-9]+$ ]] || ! [[ "$n_events" =~ ^[0-9]+$ ]]; then
  echo "PROCESS_ID and N_EVENTS must be non-negative integers" >&2
  exit 2
fi

printf -v process_tag '%06d' "$process_id"
output_file="$module_dir/output/single_particle/productions/historical/default/maps/photon_analysis_tree_${process_tag}.root"
if [[ -e "$output_file" ]]; then
  echo "Refusing to overwrite existing tree output: $output_file" >&2
  exit 3
fi

photon_setup_runtime "$module_dir"

root -l -b -q "$workflow_dir/Fun4All_PhotonAnalysisTree.C(${process_id},${n_events})"
root -l -b -q -e ".L $module_dir/validation/single_particle/check_tree.C" -e "gSystem->Exit(check_tree(\"$output_file\",true));"
