#!/bin/bash
set -eo pipefail

module_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
project_dir=$(cd "$module_dir/.." && pwd)
build_dir="$project_dir/build/current"
source /opt/sphenix/core/bin/sphenix_setup.sh -n ana

cmake -S "$module_dir" -B "$build_dir" \
  -DCMAKE_INSTALL_PREFIX="$project_dir/install" \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build "$build_dir" --parallel "${BUILD_JOBS:-4}"
cmake --install "$build_dir"
