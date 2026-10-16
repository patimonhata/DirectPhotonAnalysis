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
g++ -std=c++17 -I"$module_dir/src" -I"$OFFLINE_MAIN/include" "$module_dir/tests/test_hepmc_diphoton_decay.cc" \
  -L"$OFFLINE_MAIN/lib" -L"$OFFLINE_MAIN/lib64" -lHepMC -o "$test_dir/hepmc_decay"
"$test_dir/hepmc_decay"
g++ -std=c++17 -I"$module_dir/src" -I"$OFFLINE_MAIN/include" -I"$module_dir/../Pi0Reconstruction/install/include/pi0reconstruction" $(root-config --cflags) "$module_dir/tests/test_eta_map_branches.cc" \
  -L"$build_dir" -lPhotonAnalysisTree $(root-config --libs) -o "$test_dir/eta_map"
"$test_dir/eta_map" "$test_dir/eta_map.root"
root -l -b -q "$module_dir/tests/test_candidate_metadata.C"
root -l -b -q "$module_dir/tests/test_candidate_plot_content.C"
python3 "$module_dir/tests/test_pi0_bdt.py"
root -l -b -q "$module_dir/tests/test_eta_topology.C"
python3 "$module_dir/tests/test_eta_pipeline.py"
for script in "$module_dir/workflows/photon_candidate_selection/"*.sh "$module_dir/workflows/pi0_anchor_topology/"*.sh; do
  bash -n "$script"
done
