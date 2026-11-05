#!/usr/bin/env bash
set -euo pipefail
if (( $# < 4 || $# > 5 )); then
  echo 'Usage: run_common_matching.sh FROZEN_SOFTWARE_ROOT INPUT_MANIFEST MANIFEST_ROW NEW_OUTPUT_DIRECTORY [N_EVENTS=30]' >&2
  exit 2
fi
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
module_dir=$(cd "$script_dir/../.." && pwd)
export PHOTON_SOFTWARE_ROOT=$1
manifest=$(readlink -f "$2")
row=$3
output=$(readlink -m "$4")
events=${5:-30}
[[ $row =~ ^[0-9]+$ && $events =~ ^[0-9]+$ ]] && (( events > 0 )) || exit 2
[[ -r "$manifest" && ! -e "$output" ]] || { echo 'Require readable manifest and a new output directory' >&2; exit 2; }
# ROOT receives quoted C++ string arguments; reject characters requiring escaping.
for path in "$manifest" "$output" "$script_dir"; do
  [[ "$path" != *'"'* && "$path" != *'\'* && "$path" != *$'\n'* ]] || exit 2
done
source "$module_dir/workflows/common/runtime.sh"
photon_setup_runtime "$module_dir"
root -l -b -q "$script_dir/DiagnoseCommonMatching.C(\"$manifest\",$row,\"$output\",$events)"
# CSV streams are closed by server destruction; record exact frozen software and diagnostic code.
python3 - "$PHOTON_INSTALL_PREFIX" "$output" "$script_dir" <<'PY'
import hashlib, json, pathlib, sys
software, output, scripts = map(pathlib.Path, sys.argv[1:])
m = json.loads((software / 'manifest.json').read_text())
record = {'software_root': str(software), 'software_git_commit': m.get('source', {}).get('git_commit'),
          'software_manifest_sha256': hashlib.sha256((software / 'manifest.json').read_bytes()).hexdigest(),
          'diagnostic_sha256': {name: hashlib.sha256((scripts / name).read_bytes()).hexdigest()
                                for name in ['DiagnoseCommonMatching.C', 'run_common_matching.sh']}}
(output / 'diagnostic_provenance.json').write_text(json.dumps(record, indent=2) + '\n')
PY
