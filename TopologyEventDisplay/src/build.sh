#!/bin/bash
set -eo pipefail

module_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
project_dir=$(cd "$module_dir/.." && pwd)
build_dir="${TOPOLOGY_BUILD_DIR:-$project_dir/build/current}"
install_prefix="${TOPOLOGY_INSTALL_PREFIX:-$project_dir/install}"
photon_prefix="${PHOTON_ANALYSIS_TREE_PREFIX:-$project_dir/../PhotonAnalysisTree/install/dev/ana.565}"
source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.565

cmake -S "$module_dir" -B "$build_dir" \
  -DCMAKE_INSTALL_PREFIX="$install_prefix" \
  -DPHOTON_ANALYSIS_TREE_PREFIX="$photon_prefix" \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build "$build_dir" --parallel "${BUILD_JOBS:-4}"
cmake --install "$build_dir"
