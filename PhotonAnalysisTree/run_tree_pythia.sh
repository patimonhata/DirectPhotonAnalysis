#!/bin/bash
set -eo pipefail

module_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
process_id=${1:?usage: run_tree_pythia.sh PROCESS_ID [N_EVENTS]}
n_events=${2:-0}
if ! [[ "$process_id" =~ ^[0-9]+$ ]] || ! [[ "$n_events" =~ ^[0-9]+$ ]]; then
  echo "PROCESS_ID and N_EVENTS must be non-negative integers" >&2
  exit 2
fi

printf -v process_tag '%06d' "$process_id"
output_file="$module_dir/output/root/pythia_photon_analysis_tree_${process_tag}.root"
if [[ -e "$output_file" ]]; then
  echo "Refusing to overwrite existing Pythia tree output: $output_file" >&2
  exit 3
fi

source /opt/sphenix/core/bin/sphenix_setup.sh -n ana
export LD_LIBRARY_PATH="$module_dir/install/lib64:/sphenix/user/ryotaro/DirectPhotonAnalysis/Pi0Reconstruction/install/lib:${LD_LIBRARY_PATH:-}"

# The four matching stream files must be visible by basename in the working
# directory (for example as catalog-created symlinks).
root -l -b -q "$module_dir/macro/Fun4All_PhotonAnalysisTreePythia.C(${process_id},${n_events})"
