#!/usr/bin/env bash
set -euo pipefail
# Run after sourcing ana.565 and building the library; no installation required.
module_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
build_dir=${1:-$module_dir/build/candidate}
test_dir=$(mktemp -d /tmp/photon-candidate-tests.XXXXXX)
trap 'rm -rf "$test_dir"' EXIT
export LD_LIBRARY_PATH="$build_dir:$module_dir/../Pi0Reconstruction/install/lib:${LD_LIBRARY_PATH:-}"
g++ -std=c++17 -I"$module_dir/src" $(root-config --cflags) "$module_dir/tests/test_pi0_anchor_classification.cc" \
  -L"$build_dir" -lPhotonAnalysisTree -o "$test_dir/classification"
"$test_dir/classification"
root -l -b -q "$module_dir/tests/test_candidate_metadata.C"
root -l -b -q "$module_dir/tests/test_candidate_plot_content.C"
for script in "$module_dir/workflows/photon_candidate_selection/"*.sh "$module_dir/workflows/pi0_anchor_topology/"*.sh; do
  bash -n "$script"
done
