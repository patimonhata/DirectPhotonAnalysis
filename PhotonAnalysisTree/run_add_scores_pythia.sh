#!/bin/bash
set -eo pipefail

module_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
usage="usage: run_add_scores_pythia.sh INPUT FINAL_OUTPUT SPLIT_ONNX_MODEL [SPLIT_BDT_MODEL] [SPLIT_BDT_PPG15V1_MODEL]"
input_file=${1:?$usage}
final_output=${2:?$usage}
split_onnx_model=${3:?$usage}

split_bdt_model=${4:-/sphenix/user/shuhangli/ppg12/FunWithxgboost/binned_models/model_base_v3E_split_single_tmva.root}
split_bdt_ppg15v1_model=${5:-/sphenix/user/jaein213/photon/BDT/PPG15PhotonAN/bdt_training/outputs/nominal_base_v3E_split_3to35_allplus3jet40_ppg12split_xgb/model_ppg15_nominal_base_v3E_split_3to35_allplus3jet40_ppg12split_single_tmva.root}

if [[ ! -f "$input_file" ]]; then
  echo "Input file does not exist: $input_file" >&2
  exit 1
fi
if [[ -e "$final_output" || -L "$final_output" ]]; then
  echo "Refusing to overwrite final output: $final_output" >&2
  exit 1
fi
if [[ "$(realpath -m -- "$input_file")" == "$(realpath -m -- "$final_output")" ]]; then
  echo "Input and final output paths must differ" >&2
  exit 1
fi
for model in "$split_bdt_model" "$split_onnx_model" "$split_bdt_ppg15v1_model"; do
  if [[ ! -f "$model" ]]; then
    echo "Model file does not exist: $model" >&2
    exit 1
  fi
done

final_directory=$(dirname -- "$final_output")
mkdir -p -- "$final_directory"
if [[ ! -w "$final_directory" ]]; then
  echo "Final output directory is not writable: $final_directory" >&2
  exit 1
fi

temporary_file=$(mktemp "$final_directory/.pythia-scored.XXXXXX")
cleanup()
{
  if [[ -n "${temporary_file:-}" && -e "$temporary_file" ]]; then
    rm -f -- "$temporary_file"
  fi
}
trap cleanup EXIT
cp --reflink=auto -- "$input_file" "$temporary_file"
chmod u+w -- "$temporary_file"

source /opt/sphenix/core/bin/sphenix_setup.sh -n ana
export LD_LIBRARY_PATH="$module_dir/install/lib64:${LD_LIBRARY_PATH:-}"

# Validate the unscored schema before making any changes to the temporary copy.
root -l -b -q "$module_dir/macro/check_pythia_tree.C(\"$temporary_file\")"

root -l -b -q \
  -e ".L $module_dir/macro/add_split_bdt.C" \
  -e "gSystem->Exit(add_split_bdt(\"$temporary_file\",\"$split_bdt_model\"));"
root -l -b -q \
  -e ".L $module_dir/macro/add_split_bdt_ppg15v1.C" \
  -e "gSystem->Exit(add_split_bdt_ppg15v1(\"$temporary_file\",\"$split_bdt_ppg15v1_model\"));"
root -l -b -q \
  -e ".L $module_dir/macro/add_split_gamma_onnx.C" \
  -e "gSystem->Exit(add_split_gamma_onnx(\"$temporary_file\",\"$split_onnx_model\"));"

# Preserve schema-v4 truth consistency and validate all split score vectors.
root -l -b -q "$module_dir/macro/check_pythia_tree.C(\"$temporary_file\")"
root -l -b -q "$module_dir/macro/check_pythia_scored_tree.C(\"$temporary_file\")"

mv -- "$temporary_file" "$final_output"
temporary_file=
trap - EXIT
echo "Scored Pythia ROOT file: $final_output"
