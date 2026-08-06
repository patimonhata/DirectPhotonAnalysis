#!/bin/bash
set -eo pipefail

release=${SPHENIX_RELEASE:-ana}
source /opt/sphenix/core/bin/sphenix_setup.sh -n "$release"

source_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
study_dir=$(cd "$source_dir/.." && pwd)
build_dir="$study_dir/build/$release"

mkdir -p "$study_dir/output/root" "$study_dir/output/merge" "$study_dir/output/condor"
mkdir -p "$build_dir" "$study_dir/install"

cd "$build_dir"
"$source_dir/autogen.sh" --prefix="$study_dir/install"
make -j"${BUILD_JOBS:-4}"
make install
