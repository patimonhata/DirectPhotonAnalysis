#!/usr/bin/env bash
set -euo pipefail
workflow_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
module_dir=$(cd "$workflow_dir/../.." && pwd)
usage="Usage: $0 merge CONFIGURATION NEW_OUTPUT_DIRECTORY | reduce CONFIGURATION SAMPLE SHARD_INDEX NEW_OUTPUT_BASE"
mode=${1:?$usage}
configuration=${2:?$usage}
case "$configuration" in
  ClusterE0p12_Pi0Partner0p15|ClusterE0p12_Pi0Partner0p15_Pi0Window0p08-0p18|ClusterE0p12_Pi0Partner0p5|ClusterE0p12_Pi0Partner0p5_Pi0Window0p08-0p18) ;;
  *) echo "Unknown retained v11 configuration: $configuration" >&2; exit 2 ;;
esac
reference="$module_dir/output/photon_candidate_selection/productions/topology_v11/$configuration"
case "$configuration" in
  ClusterE0p12_Pi0Partner0p15*) commit=5c1a09ba54b091e2f3f4419d8106e4027cc9160f ;; # partial schema 5
  *) commit=f9813670f083bbcf6e3d8116a7d87aae1f868601 ;; # partial schema 6
esac
mkdir -p "$module_dir/build"
cache=$(mktemp -d "$module_dir/build/.reference-v11.XXXXXX")
trap 'rm -rf -- "$cache"' EXIT
git -C "$module_dir/.." archive "$commit" PhotonAnalysisTree/workflows/photon_candidate_selection PhotonAnalysisTree/macro/Utilities | tar -x -C "$cache"
old_workflow="$cache/PhotonAnalysisTree/workflows/photon_candidate_selection"
if [[ "$mode" == merge && $# == 3 ]]; then
  output=$3
  [[ ! -e "$output" ]] || { echo "Output already exists: $output" >&2; exit 3; }
  macro=MergePythiaPhotonCandidateSelection
  call="${macro}(\"jet\",\"$reference/partials/jet\",\"$output\")"
elif [[ "$mode" == reduce && $# == 5 ]]; then
  sample=$3; shard=$4; output=$5
  case "$sample:$shard" in jet3:0|jet5:0|jet8:0|jet20:0|jet30:0|jet40:0|jet12:[0-9]) ;; *) echo "Invalid sample/shard" >&2; exit 2 ;; esac
  [[ ! -e "$output.root" ]] || { echo "Output already exists: $output.root" >&2; exit 3; }
  macro=ReducePythiaPhotonCandidateSelection
  call="${macro}(\"jet\",\"$reference/maps\",\"$output\",true,40,40.0,\"$sample\",$shard)"
else
  echo "$usage" >&2; exit 2
fi
if [[ "$output" == *\"* || "$output" == *\\* ]]; then echo "Output cannot contain quotes or backslashes" >&2; exit 2; fi
set +u
source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.565
set -u
root -l -b -q -e ".L $old_workflow/$macro.C" -e "gSystem->Exit($call);"
