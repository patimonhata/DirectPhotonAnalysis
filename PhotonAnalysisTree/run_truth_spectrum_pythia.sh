#!/bin/bash
set -eo pipefail

module_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
usage="usage: run_truth_spectrum_pythia.sh INPUT_SUFFIX [N_EVENTS] [OUTPUT_DIRECTORY]"
input_suffix=${1:?$usage}
n_events=${2:-0}
output_directory=${3:-$module_dir/output/truth_root}
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
final_output="$output_directory/pythia_truth_spectrum_tree_${output_tag}.root"
if [[ -e "$final_output" || -L "$final_output" ]]; then
  echo "Refusing to overwrite existing truth spectrum output: $final_output" >&2
  exit 3
fi
mkdir -p -- "$output_directory"
if [[ ! -w "$output_directory" ]]; then
  echo "Output directory is not writable: $output_directory" >&2
  exit 3
fi

temporary_output=$(mktemp "$output_directory/.pythia-truth-spectrum.XXXXXX")
cleanup()
{
  if [[ -n "${temporary_output:-}" && -e "$temporary_output" ]]; then
    rm -f -- "$temporary_output"
  fi
}
trap cleanup EXIT

source /opt/sphenix/core/bin/sphenix_setup.sh -n ana
export LD_LIBRARY_PATH="$module_dir/install/lib64:/sphenix/user/ryotaro/DirectPhotonAnalysis/Pi0Reconstruction/install/lib:${LD_LIBRARY_PATH:-}"

root -l -b -q \
  "$module_dir/macro/Fun4All_PythiaTruthSpectrumTree.C(\"${input_suffix}\",${n_events},\"${temporary_output}\")"
root -l -b -q \
  "$module_dir/macro/check_pythia_truth_spectrum_tree.C(\"${temporary_output}\")"

mv -- "$temporary_output" "$final_output"
temporary_output=
trap - EXIT
echo "Pythia truth spectrum ROOT file: $final_output"
