#!/bin/bash
set -eo pipefail

source /opt/sphenix/core/bin/sphenix_setup.sh -n "${SPHENIX_RELEASE:-ana}"

source_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
study_dir=$(cd "$source_dir/.." && pwd)
build_dir="$study_dir/build/${SPHENIX_RELEASE:-ana}"
mkdir -p "$study_dir/output/root" "$study_dir/output/merge" "$study_dir/output/condor"

cmake -S "$source_dir" -B "$build_dir" \
  -DCMAKE_INSTALL_PREFIX="$study_dir/install" \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build "$build_dir" --parallel "${BUILD_JOBS:-4}"
cmake --install "$build_dir"
