#!/usr/bin/env bash
set -euo pipefail

workflow_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
module_dir=$(cd "$workflow_dir/../.." && pwd)
usage="usage: workflows/photon_candidate_selection/run_map.sh JOB_INDEX CHUNK_OFFSET TOTAL_FILES FILES_PER_JOB INPUT_MANIFEST SAMPLE_NAME OUTPUT_DIRECTORY [N_EVENTS] [MIN_CLUSTER_ENERGY_GEV]"
job_index=${1:?$usage}
chunk_offset=${2:?$usage}
total_files=${3:?$usage}
files_per_job=${4:?$usage}
input_manifest=${5:?$usage}
sample_name=${6:?$usage}
output_directory=${7:?$usage}
n_events=${8:-0}
min_cluster_energy=${9:-0.1}
for value in "$job_index" "$chunk_offset" "$total_files" "$files_per_job" "$n_events"; do
  if ! [[ "$value" =~ ^[0-9]+$ ]]; then
    echo "Chunk, file counts, and N_EVENTS must be non-negative integers: $value" >&2
    exit 2
  fi
done
if (( total_files == 0 || files_per_job == 0 )); then
  echo "TOTAL_FILES and FILES_PER_JOB must be positive" >&2
  exit 2
fi
if ! [[ "$min_cluster_energy" =~ ^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$ ]]; then
  echo "MIN_CLUSTER_ENERGY_GEV must be a non-negative finite number: $min_cluster_energy" >&2
  exit 2
fi
case "$sample_name" in
  photonjet3|photonjet5|photonjet10|photonjet20|jet3|jet5|jet8|jet12|jet20|jet30|jet40) ;;
  *)
    echo "Unsupported SAMPLE_NAME: $sample_name" >&2
    exit 2
    ;;
esac
if [[ ! -r "$input_manifest" ]]; then
  echo "Input manifest is not readable: $input_manifest" >&2
  exit 2
fi
manifest_rows=$(wc -l < "$input_manifest")
if (( manifest_rows < total_files )); then
  echo "Manifest has $manifest_rows rows, fewer than TOTAL_FILES=$total_files" >&2
  exit 2
fi

chunk_index=$((job_index + chunk_offset))
manifest_begin=$((chunk_index * files_per_job))
if (( manifest_begin >= total_files )); then
  echo "Chunk $chunk_index begins beyond TOTAL_FILES=$total_files" >&2
  exit 2
fi
manifest_end=$((manifest_begin + files_per_job))
if (( manifest_end > total_files )); then
  manifest_end=$total_files
fi
printf -v chunk_tag "%06d" "$chunk_index"
final_output="$output_directory/map_${chunk_tag}.root"
if [[ -e "$final_output" || -L "$final_output" ]]; then
  echo "Refusing to overwrite existing map output: $final_output" >&2
  exit 3
fi
mkdir -p -- "$output_directory"
temporary_output=$(mktemp --suffix=.root "$output_directory/.photon-candidate-map.XXXXXX")
cleanup()
{
  if [[ -n "${temporary_output:-}" && -e "$temporary_output" ]]; then
    rm -f -- "$temporary_output"
  fi
}
trap cleanup EXIT

model_file="$module_dir/models/model_ppg15_nominal_base_v3E_split_3to35_allplus3jet40_ppg12split_single_tmva.root"
set +u
source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.565
set -u
export LD_LIBRARY_PATH="$module_dir/install/lib64:/sphenix/user/ryotaro/DirectPhotonAnalysis/Pi0Reconstruction/install/lib:${LD_LIBRARY_PATH:-}"
root -l -b -q \
  "$workflow_dir/Fun4All_PythiaPhotonCandidateTreeMap.C(\"${input_manifest}\",${manifest_begin},${manifest_end},\"${temporary_output}\",\"${sample_name}\",${chunk_index},${n_events},${min_cluster_energy},\"${model_file}\")"
root -l -b -q \
  "$workflow_dir/check_pythia_photon_candidate_map.C(\"${temporary_output}\",${min_cluster_energy})"

mv -- "$temporary_output" "$final_output"
temporary_output=
trap - EXIT
echo "Pythia photon-candidate map ROOT file: $final_output"
