#!/usr/bin/env bash
set -euo pipefail
if (( $# != 1 )); then
  echo "Usage: $0 OUTPUT_DIRECTORY" >&2
  exit 2
fi
workflow_dir=$(cd "$(dirname "$0")" && pwd)
module_dir=$(cd "$workflow_dir/../../.." && pwd)
source "$module_dir/workflows/common/runtime.sh"
photon_dispatch "$module_dir" "workflows/diagnostics/tagging_purity/run_hybrid_merge.sh" "$@"
output_directory=$1
partial_root=/sphenix/user/ryotaro/DirectPhotonAnalysis/PhotonAnalysisTree/output/tagging_purity_diagnostic/productions/historical/hybrid_threshold/partials/jet
mkdir -p "$output_directory"
output_root=$output_directory/hybrid_threshold_diagnostic.root
if [[ -e "$output_root" ]]; then echo "Output already exists: $output_root" >&2; exit 2; fi
inputs=()
for sample in jet3 jet5 jet8 jet20 jet30 jet40; do
  inputs+=("$partial_root/$sample/shard_0/hybrid_threshold_diagnostic.root")
done
for shard in {0..9}; do
  inputs+=("$partial_root/jet12/shard_$shard/hybrid_threshold_diagnostic.root")
done
for input in "${inputs[@]}"; do
  [[ -s "$input" ]] || { echo "Missing partial: $input" >&2; exit 2; }
done
temporary="$output_root.tmp.$$"
trap 'rm -f "$temporary"' EXIT
photon_setup_runtime "$module_dir"
hadd -f "$temporary" "${inputs[@]}"
root -l -b -q -e ".L $workflow_dir/SummarizeHybridThresholdDiagnostic.C" \
  -e "gSystem->Exit(SummarizeHybridThresholdDiagnostic(\"$temporary\",\"$output_directory/hybrid_threshold_diagnostic\"))"
mv "$temporary" "$output_root"
trap - EXIT
