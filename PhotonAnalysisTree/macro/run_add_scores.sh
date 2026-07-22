#!/bin/bash
set -eo pipefail

module_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
input_file=${1:?usage: run_add_scores.sh INPUT BDT_OUTPUT FINAL_OUTPUT [BDT_MODEL] [ONNX_MODEL]}
bdt_output=${2:?usage: run_add_scores.sh INPUT BDT_OUTPUT FINAL_OUTPUT [BDT_MODEL] [ONNX_MODEL]}
final_output=${3:?usage: run_add_scores.sh INPUT BDT_OUTPUT FINAL_OUTPUT [BDT_MODEL] [ONNX_MODEL]}
bdt_model=${4:-/sphenix/user/shuhangli/ppg12/FunWithxgboost/binned_models/model_base_v3E_split_single_tmva.root}
onnx_model=${5:-$module_dir/models/best_model.onnx}

source /opt/sphenix/core/bin/sphenix_setup.sh -n ana
export LD_LIBRARY_PATH="$module_dir/install/lib64:${LD_LIBRARY_PATH:-}"

root -l -b -q \
  -e ".L $module_dir/macro/add_split_bdt.C" \
  -e "gSystem->Exit(add_split_bdt(\"$input_file\",\"$bdt_output\",\"$bdt_model\"));"
root -l -b -q \
  -e ".L $module_dir/macro/add_nosplit_gamma_onnx.C" \
  -e "gSystem->Exit(add_nosplit_gamma_onnx(\"$bdt_output\",\"$final_output\",\"$onnx_model\"));"
root -l -b -q \
  -e ".L $module_dir/macro/check_scored_tree.C" \
  -e "gSystem->Exit(check_scored_tree(\"$final_output\"));"
