#!/usr/bin/env bash
set -euo pipefail
# Run after sourcing ana.565 and building the library; no installation required.
module_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
prefix=${PHOTON_INSTALL_PREFIX:-$module_dir/install/dev/ana.565}
build_dir=${1:-$module_dir/build/dev/ana.565/photon}
export ROOT_INCLUDE_PATH="$prefix/include:$prefix/dependencies/Pi0Reconstruction/include/pi0reconstruction:${ROOT_INCLUDE_PATH:-}"
test_dir=$(mktemp -d /tmp/photon-candidate-tests.XXXXXX)
trap 'rm -rf "$test_dir"' EXIT
export LD_LIBRARY_PATH="$build_dir:$prefix/dependencies/Pi0Reconstruction/lib:${LD_LIBRARY_PATH:-}"
g++ -std=c++17 -I"$module_dir/src" $(root-config --cflags) "$module_dir/tests/test_pi0_anchor_classification.cc" \
  -L"$build_dir" -lPhotonAnalysisTree -o "$test_dir/classification"
"$test_dir/classification"
g++ -std=c++17 -I"$module_dir/src" "$module_dir/tests/test_common_composition.cc" -o "$test_dir/common_composition"
"$test_dir/common_composition"
root -l -b -q "$module_dir/tests/test_common_decay_truth.C"
g++ -std=c++17 -I"$module_dir/src" -I"$OFFLINE_MAIN/include" "$module_dir/tests/test_hepmc_diphoton_decay.cc" \
  -L"$OFFLINE_MAIN/lib" -L"$OFFLINE_MAIN/lib64" -lHepMC -o "$test_dir/hepmc_decay"
"$test_dir/hepmc_decay"
g++ -std=c++17 -I"$module_dir/src" -I"$OFFLINE_MAIN/include" -I"$prefix/dependencies/Pi0Reconstruction/include/pi0reconstruction" $(root-config --cflags) "$module_dir/tests/test_eta_map_branches.cc" \
  -L"$build_dir" -lPhotonAnalysisTree $(root-config --libs) -o "$test_dir/eta_map"
"$test_dir/eta_map" "$test_dir/eta_map.root"
root -l -b -q "$module_dir/tests/test_candidate_metadata.C"
root -l -b -q "$module_dir/tests/test_candidate_plot_content.C"
python3 "$module_dir/tests/test_pi0_bdt.py"
root -l -b -q "$module_dir/tests/test_eta_topology.C"
root -l -b -q "$module_dir/tests/test_other_origin.C"
root -l -b -q "$module_dir/tests/test_decay_origin.C"
root -l -b -q "$module_dir/tests/test_decay_origin_reduce.C"
root -l -b -q "$module_dir/tests/test_eta_cascade.C"
root -l -b -q "$module_dir/tests/test_eta_decay_flight.C"
python3 "$module_dir/tests/test_eta_pipeline.py"
python3 "$module_dir/tests/test_tagging_mass.py"
python3 "$module_dir/tests/test_software.py"
while IFS= read -r -d '' script; do
  bash -n "$script"
done < <(find "$module_dir/workflows" "$module_dir/validation" -name '*.sh' -print0)
