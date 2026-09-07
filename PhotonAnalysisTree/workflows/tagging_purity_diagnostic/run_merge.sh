#!/usr/bin/env bash
set -euo pipefail
if (( $# != 2 )); then
  echo "Usage: $0 CONFIGURATION OUTPUT_DIRECTORY" >&2
  exit 2
fi
workflow_dir=$(cd "$(dirname "$0")" && pwd)
configuration=$1
output_directory=$2
partial_root=/sphenix/user/ryotaro/DirectPhotonAnalysis/PhotonAnalysisTree/output/plots/tagging_purity_diagnostic/reduce/$configuration/jet/partial
mkdir -p "$output_directory"
output_root=$output_directory/tagging_purity_diagnostic.root
if [[ -e "$output_root" ]]; then echo "Output already exists: $output_root" >&2; exit 2; fi
inputs=()
for sample in jet3 jet5 jet8 jet20 jet30 jet40; do
  inputs+=("$partial_root/$sample/shard_0/tagging_purity_diagnostic.root")
done
for shard in {0..9}; do
  inputs+=("$partial_root/jet12/shard_$shard/tagging_purity_diagnostic.root")
done
for input in "${inputs[@]}"; do
  [[ -s "$input" ]] || { echo "Missing partial: $input" >&2; exit 2; }
done
temporary="$output_root.tmp.$$"
trap 'rm -f "$temporary"' EXIT
set +u
source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.565
set -u
hadd -f "$temporary" "${inputs[@]}"
root -l -b -q -e ".L $workflow_dir/SummarizeTaggingPurityDiagnostic.C" \
  -e "gSystem->Exit(SummarizeTaggingPurityDiagnostic(\"$temporary\",\"$output_directory/tagging_purity_diagnostic\"))"
mv "$temporary" "$output_root"
trap - EXIT
