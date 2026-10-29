#!/bin/sh
set -eu
study_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
exec python3 "$study_dir/run.py" "$@"
