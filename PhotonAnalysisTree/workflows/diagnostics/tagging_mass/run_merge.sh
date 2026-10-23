#!/usr/bin/env bash
set -euo pipefail

if (( $# != 3 )); then
  echo "Usage: $0 FAMILY PARTIAL_ROOT OUTPUT_BASE" >&2
  exit 2
fi
workflow_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
family=$1
partial_root=$2
output_base=$3
if [[ "$family" != jet && "$family" != photonjet ]]; then echo "FAMILY must be jet or photonjet" >&2; exit 2; fi
if [[ ! -d "$partial_root" || -z "$output_base" || "$partial_root$output_base" == *\"* || "$partial_root$output_base" == *\\* ]]; then
  echo "Invalid partial/output path" >&2
  exit 2
fi
partial_root=$(cd "$partial_root" && pwd)
set +u
source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.565
set -u
root -l -b -q -e ".L $workflow_dir/MergeTaggingMassDiagnostic.C" \
  -e "gSystem->Exit(MergeTaggingMassDiagnostic(\"$family\",\"$partial_root\",\"$output_base\"))"
