#!/bin/bash
set -eo pipefail

module_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source /opt/sphenix/core/bin/sphenix_setup.sh -n ana

cmake -S "$module_dir" -B "$module_dir/build" \
  -DCMAKE_INSTALL_PREFIX="$module_dir/install" \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build "$module_dir/build" --parallel "${BUILD_JOBS:-4}"
cmake --install "$module_dir/build"
