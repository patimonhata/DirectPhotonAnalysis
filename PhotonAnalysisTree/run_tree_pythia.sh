#!/bin/bash
set -eo pipefail

module_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
input_suffix=${1:?usage: run_tree_pythia.sh INPUT_SUFFIX [N_EVENTS] [OUTPUT_DIRECTORY]}
n_events=${2:-0}
output_directory=${3:-$module_dir/output/root}
if ! [[ "$input_suffix" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*\.root$ ]]; then
  echo "INPUT_SUFFIX must be a basename ending in .root: $input_suffix" >&2
  exit 2
fi
if ! [[ "$input_suffix" =~ -[0-9]+\.root$ ]]; then
  echo "INPUT_SUFFIX must end in a numeric segment followed by .root: $input_suffix" >&2
  exit 2
fi
if ! [[ "$n_events" =~ ^[0-9]+$ ]]; then
  echo "N_EVENTS must be a non-negative integer" >&2
  exit 2
fi

output_tag=${input_suffix%.root}
output_file="$output_directory/pythia_photon_analysis_tree_${output_tag}.root"
if [[ -e "$output_file" ]]; then
  echo "Refusing to overwrite existing Pythia tree output: $output_file" >&2
  exit 3
fi

source /opt/sphenix/core/bin/sphenix_setup.sh -n ana
export LD_LIBRARY_PATH="$module_dir/install/lib64:/sphenix/user/ryotaro/DirectPhotonAnalysis/Pi0Reconstruction/install/lib:${LD_LIBRARY_PATH:-}"

# Fun4All resolves each DST through its basename. The macro adds the four
# stream-specific prefixes to this common input suffix.
root -l -b -q "$module_dir/macro/Fun4All_PhotonAnalysisTreePythia.C(\"${input_suffix}\",${n_events},\"${output_directory}\")"
