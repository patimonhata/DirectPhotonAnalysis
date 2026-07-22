#!/bin/bash
set -eo pipefail

module_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
scored_dir=${1:-$module_dir/output/root}
merged_file=${2:-$module_dir/output/merged/all.root}
manifest_file=${3:-$module_dir/output/merged/manifest.json}
bdt_model=${4:-/sphenix/user/shuhangli/ppg12/FunWithxgboost/binned_models/model_base_v3E_split_single_tmva.root}
onnx_model=${5:-$module_dir/models/best_model.onnx}
onnx_metadata="$module_dir/models/onnx_metadata.json"

for command in jq sha256sum; do
  if ! command -v "$command" >/dev/null 2>&1; then
    echo "Required command not found: $command" >&2
    exit 2
  fi
done
source /opt/sphenix/core/bin/sphenix_setup.sh -n ana
if ! command -v root >/dev/null 2>&1; then
  echo "Required command not found after sPHENIX setup: root" >&2
  exit 2
fi
for path in "$scored_dir" "$merged_file" "$bdt_model" "$onnx_model" "$onnx_metadata"; do
  if [[ ! -e "$path" ]]; then
    echo "Required input does not exist: $path" >&2
    exit 3
  fi
done
if [[ -e "$manifest_file" && ${FORCE:-0} != 1 ]]; then
  echo "Refusing to overwrite manifest: $manifest_file (set FORCE=1 to replace it)" >&2
  exit 4
fi

mapfile -t scored_files < <(
  find "$scored_dir" -maxdepth 1 -type f -regextype posix-extended \
    -regex '.*/photon_analysis_tree_[0-9]{6}_scored\.root' | sort
)
scored_count=${#scored_files[@]}
if (( scored_count == 0 )); then
  echo "No scored ROOT files found in $scored_dir" >&2
  exit 5
fi

root_output=$(root -l -b -q -e "
  TFile file(\"$merged_file\", \"READ\");
  auto* events = file.Get<TTree>(\"event_tree\");
  auto* metadata = file.Get<TTree>(\"metadata\");
  if (file.IsZombie() || !events || !metadata) gSystem->Exit(20);
  std::string* input_file = nullptr;
  UInt_t source_file_id = 0;
  ULong64_t processed = 0, written = 0, invalid_truth = 0, invalid_detector = 0;
  if (metadata->SetBranchAddress(\"input_file\", &input_file) < 0 ||
      metadata->SetBranchAddress(\"source_file_id\", &source_file_id) < 0 ||
      metadata->SetBranchAddress(\"n_events_processed\", &processed) < 0 ||
      metadata->SetBranchAddress(\"n_events_written\", &written) < 0 ||
      metadata->SetBranchAddress(\"n_events_invalid_truth\", &invalid_truth) < 0 ||
      metadata->SetBranchAddress(\"n_events_invalid_detector\", &invalid_detector) < 0) gSystem->Exit(21);
  std::string first_input, last_input;
  UInt_t first_id = 0, last_id = 0;
  ULong64_t sum_processed = 0, sum_written = 0, sum_invalid_truth = 0, sum_invalid_detector = 0;
  for (Long64_t entry = 0; entry < metadata->GetEntries(); ++entry) {
    if (metadata->GetEntry(entry) <= 0 || !input_file) gSystem->Exit(22);
    if (entry == 0) { first_input = *input_file; first_id = source_file_id; }
    if (entry + 1 == metadata->GetEntries()) { last_input = *input_file; last_id = source_file_id; }
    sum_processed += processed;
    sum_written += written;
    sum_invalid_truth += invalid_truth;
    sum_invalid_detector += invalid_detector;
  }
  std::cout << \"MANIFEST_DATA\\t\" << events->GetEntries()
            << \"\\t\" << metadata->GetEntries()
            << \"\\t\" << file.GetListOfKeys()->GetSize()
            << \"\\t\" << first_id << \"\\t\" << last_id
            << \"\\t\" << sum_processed << \"\\t\" << sum_written
            << \"\\t\" << sum_invalid_truth << \"\\t\" << sum_invalid_detector
            << \"\\t\" << first_input << \"\\t\" << last_input << std::endl;
  gSystem->Exit(0);")
manifest_line=$(printf '%s\n' "$root_output" | grep '^MANIFEST_DATA' || true)
if [[ -z "$manifest_line" ]]; then
  echo "Could not read merged ROOT metadata" >&2
  exit 6
fi
IFS=$'\t' read -r _ event_count metadata_entries root_key_count first_id last_id \
  events_processed events_written invalid_truth invalid_detector first_input last_input <<< "$manifest_line"
if (( metadata_entries != scored_count )); then
  echo "metadata entries ($metadata_entries) != scored file count ($scored_count)" >&2
  exit 7
fi
if (( root_key_count != 2 || event_count != events_written )); then
  echo "Merged ROOT validation failed: keys=$root_key_count events=$event_count written=$events_written" >&2
  exit 8
fi
if [[ $(dirname "$first_input") != $(dirname "$last_input") ]]; then
  echo "First and last source files are in different directories" >&2
  exit 9
fi

onnx_sha256=$(sha256sum "$onnx_model" | awk '{print $1}')
declared_onnx_sha256=$(jq -r '.onnx_sha256' "$onnx_metadata")
if [[ "$onnx_sha256" != "$declared_onnx_sha256" ]]; then
  echo "ONNX SHA-256 does not match models/onnx_metadata.json" >&2
  exit 10
fi
bdt_sha256=$(sha256sum "$bdt_model" | awk '{print $1}')
root_sha256=$(sha256sum "$merged_file" | awk '{print $1}')
onnx_metadata_sha256=$(sha256sum "$onnx_metadata" | awk '{print $1}')
created_at=$(date -u '+%Y-%m-%dT%H:%M:%SZ')
git_commit=$(git -C "$module_dir" rev-parse HEAD 2>/dev/null || printf 'unknown')
if [[ -n $(git -C "$module_dir" status --porcelain --untracked-files=normal -- . 2>/dev/null) ]]; then
  git_dirty=true
else
  git_dirty=false
fi
source_directory=$(dirname "$first_input")
source_dataset=$(basename "$source_directory")
dataset_name="${source_dataset}_scored"
mkdir -p "$(dirname "$manifest_file")"
temporary_manifest=$(mktemp "$(dirname "$manifest_file")/.manifest.XXXXXX")
trap 'rm -f "$temporary_manifest"' EXIT

jq -n \
  --arg dataset "$dataset_name" \
  --arg created_at "$created_at" \
  --arg merged_file "$merged_file" \
  --arg merged_sha256 "$root_sha256" \
  --arg scored_directory "$scored_dir" \
  --arg source_directory "$source_directory" \
  --arg first_input "$first_input" \
  --arg last_input "$last_input" \
  --arg bdt_model "$bdt_model" \
  --arg bdt_sha256 "$bdt_sha256" \
  --arg onnx_model "$onnx_model" \
  --arg onnx_sha256 "$onnx_sha256" \
  --arg onnx_metadata "$onnx_metadata" \
  --arg onnx_metadata_sha256 "$onnx_metadata_sha256" \
  --arg git_commit "$git_commit" \
  --argjson git_dirty "$git_dirty" \
  --argjson scored_count "$scored_count" \
  --argjson event_count "$event_count" \
  --argjson metadata_entries "$metadata_entries" \
  --argjson first_id "$first_id" \
  --argjson last_id "$last_id" \
  --argjson events_processed "$events_processed" \
  --argjson events_written "$events_written" \
  --argjson invalid_truth "$invalid_truth" \
  --argjson invalid_detector "$invalid_detector" \
  --slurpfile onnx_deployment "$onnx_metadata" \
  '{
    manifest_version: 1,
    dataset: $dataset,
    created_at_utc: $created_at,
    files: {
      merged_root: {path: $merged_file, sha256: $merged_sha256},
      scored_inputs: {
        directory: $scored_directory,
        filename_pattern: "photon_analysis_tree_[0-9]{6}_scored.root",
        count: $scored_count
      },
      root_layout: ["event_tree", "metadata"]
    },
    source_dataset: {
      directory: $source_directory,
      first_file: $first_input,
      last_file: $last_input,
      source_file_id: {first: $first_id, last: $last_id}
    },
    counts: {
      event_tree_entries: $event_count,
      metadata_entries: $metadata_entries,
      events_processed: $events_processed,
      events_written: $events_written,
      events_invalid_truth: $invalid_truth,
      events_invalid_detector: $invalid_detector
    },
    scores: {
      split_bdt: {
        model: {path: $bdt_model, sha256: $bdt_sha256, key: "myBDT"},
        score_branch: "split_cluster_bdt_base_v3E_score",
        valid_branch: "split_cluster_bdt_base_v3E_valid",
        feature_order: ["ET", "weta_cogx", "wphi_cogx", "vertex_z", "cluster_eta", "e11_over_e33", "et1", "et2", "et3", "et4", "e32_over_e35"],
        domain_warning: "Documented model performance bins begin at cluster ET=6 GeV; lower-ET scores are extrapolations."
      },
      nosplit_gamma_onnx: {
        model: {path: $onnx_model, sha256: $onnx_sha256},
        deployment_metadata: {path: $onnx_metadata, sha256: $onnx_metadata_sha256, content: $onnx_deployment[0]},
        score_branch: "nosplit_cluster_p_gamma",
        valid_branch: "nosplit_cluster_p_gamma_valid",
        global_features: ["cluster_eta", "log1p(cluster_energy/cosh(cluster_eta))", "log1p(cluster_ntower)"],
        point_features: ["tower_eta-cluster_eta", "wrapped(tower_phi-cluster_phi)", "log1p(tower_energy)", "tower_energy/cluster_energy"],
        domain_warning: "The model was trained only on events with exactly one no-split cluster; multi-cluster scores require separate validation."
      }
    },
    code: {git_commit: $git_commit, working_tree_dirty: $git_dirty}
  }' > "$temporary_manifest"

mv "$temporary_manifest" "$manifest_file"
trap - EXIT
echo "Wrote dataset manifest: $manifest_file"
