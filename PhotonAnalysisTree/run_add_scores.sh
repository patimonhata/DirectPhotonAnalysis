#!/bin/bash
set -eo pipefail

module_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd) # the absolute path to the directory in which this script is placed.
usage="usage: run_add_scores.sh INPUT FINAL_OUTPUT [SPLIT_BDT_MODEL] [NOSPLIT_BDT_MODEL] [SPLIT_ONNX_MODEL] [NOSPLIT_ONNX_MODEL] [SPLIT_BDT_PPG15V1_MODEL]"
input_file=${1:?$usage}
final_output=${2:?$usage}

if [[ ! -f "$input_file" ]]; then
  echo "Input file does not exist: $input_file" >&2
  exit 1
fi
if [[ ! -w "$input_file" ]]; then
  echo "Input file is not writable: $input_file" >&2
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
final_directory=$(dirname -- "$final_output")
mkdir -p -- "$final_directory"
if [[ ! -w "$final_directory" ]]; then
  echo "Final output directory is not writable: $final_directory" >&2
  exit 1
fi

split_bdt_model=${3:-/sphenix/user/shuhangli/ppg12/FunWithxgboost/binned_models/model_base_v3E_split_single_tmva.root}
nosplit_bdt_model=${4:-/sphenix/user/shuhangli/ppg12/FunWithxgboost/binned_models/model_base_v3E_nosplit_single_tmva.root}
split_onnx_model=${5:-$module_dir/models/best_model.onnx}
nosplit_onnx_model=${6:-$module_dir/models/best_model.onnx}
split_bdt_ppg15v1_model=${7:-/sphenix/user/jaein213/photon/BDT/PPG15PhotonAN/bdt_training/outputs/nominal_base_v3E_split_3to35_allplus3jet40_ppg12split_xgb/model_ppg15_nominal_base_v3E_split_3to35_allplus3jet40_ppg12split_single_tmva.root}

source /opt/sphenix/core/bin/sphenix_setup.sh -n ana
export LD_LIBRARY_PATH="$module_dir/install/lib64:${LD_LIBRARY_PATH:-}"

#for single gun
root -l -b -q \
  -e ".L $module_dir/macro/add_split_bdt.C" \
  -e "gSystem->Exit(add_split_bdt(\"$input_file\",\"$split_bdt_model\"));"
root -l -b -q \
  -e ".L $module_dir/macro/add_split_bdt_ppg15v1.C" \
  -e "gSystem->Exit(add_split_bdt_ppg15v1(\"$input_file\",\"$split_bdt_ppg15v1_model\"));"
root -l -b -q \
  -e ".L $module_dir/macro/add_nosplit_bdt.C" \
  -e "gSystem->Exit(add_nosplit_bdt(\"$input_file\",\"$nosplit_bdt_model\"));"
root -l -b -q \
  -e ".L $module_dir/macro/add_split_gamma_onnx.C" \
  -e "gSystem->Exit(add_split_gamma_onnx(\"$input_file\",\"$split_onnx_model\"));"
  # -e "gSystem->Exit(add_split_gamma_onnx(\"$input_file\",\"$split_onnx_model\"));"
root -l -b -q \
  -e ".L $module_dir/macro/add_nosplit_gamma_onnx.C" \
  -e "gSystem->Exit(add_nosplit_gamma_onnx(\"$input_file\",\"$nosplit_onnx_model\"));"
root -l -b -q \
  -e ".L $module_dir/macro/check_scored_tree.C" \
  -e "gSystem->Exit(check_scored_tree(\"$input_file\"));"

mv -- "$input_file" "$final_output"
echo "Scored ROOT file: $final_output"

# for pythia
# root -l -b -q \
#    -e ".L $module_dir/macro/add_split_gamma_onnx.C" \
#    -e "gSystem->Exit(add_split_gamma_onnx(\"$input_file\",\"$split_onnx_model\"));"
