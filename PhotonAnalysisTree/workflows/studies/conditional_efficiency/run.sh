#!/usr/bin/env bash
set -euo pipefail
workflow_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
exec python3 -B "$workflow_dir/../../common/replay_study.py" "$(basename "$workflow_dir")" "$@"
