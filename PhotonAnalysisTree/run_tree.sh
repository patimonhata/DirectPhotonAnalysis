#!/bin/bash
set -eo pipefail

module_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd) # the absolute path to the directory in which this script is placed.
process_id=${1:?usage: run_tree.sh PROCESS_ID [N_EVENTS]}
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

# SingleGun:
# root -l -b -q "$module_dir/macro/Fun4All_PhotonAnalysisTree.C(${process_id},${n_events})"
# Pythia Jet5:
root -l -b -q "$module_dir/macro/Fun4All_PhotonAnalysisTreePythia.C(${process_id},${n_events})"
