#!/usr/bin/env bash
set -euo pipefail

if (( $# < 5 || $# > 7 )); then
  echo "Usage: $0 FAMILY MAP_ROOT OUTPUT_BASE REQUIRE_COMPLETE TAGGING_PARTNER_MIN_ENERGY_GEV [SAMPLE_FILTER] [MAX_EVENTS_PER_SAMPLE]" >&2
  exit 2
fi
workflow_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
family=$1
map_root=$2
output_base=$3
require_complete=$4
tagging_threshold=$5
sample_filter=${6:-}
max_events=${7:-0}
if [[ "$family" != jet && "$family" != photonjet ]]; then echo "FAMILY must be jet or photonjet" >&2; exit 2; fi
if [[ "$require_complete" != true && "$require_complete" != false ]]; then echo "REQUIRE_COMPLETE must be true or false" >&2; exit 2; fi
if ! [[ "$tagging_threshold" =~ ^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$ ]]; then echo "Invalid tagging threshold" >&2; exit 2; fi
if ! [[ "$max_events" =~ ^[0-9]+$ ]]; then echo "MAX_EVENTS_PER_SAMPLE must be a non-negative integer" >&2; exit 2; fi
case "$family:$sample_filter" in
  jet:|jet:jet3|jet:jet5|jet:jet8|jet:jet12|jet:jet20|jet:jet30|jet:jet40|photonjet:|photonjet:photonjet3|photonjet:photonjet5|photonjet:photonjet10|photonjet:photonjet20) ;;
  *) echo "SAMPLE_FILTER does not belong to FAMILY" >&2; exit 2 ;;
esac
if [[ ! -d "$map_root" || -z "$output_base" || "$map_root$output_base$sample_filter" == *\"* || "$map_root$output_base$sample_filter" == *\\* ]]; then
  echo "Invalid input/output path" >&2
  exit 2
fi
map_root=$(cd "$map_root" && pwd)
set +u
source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.565
set -u
root -l -b -q -e ".L $workflow_dir/PlotDisplacedPartnerMassDiagnostic.C" -e "gSystem->Exit(PlotDisplacedPartnerMassDiagnostic(\"$family\",\"$map_root\",\"$output_base\",$require_complete,$tagging_threshold,\"$sample_filter\",$max_events))"
