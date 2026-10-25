#!/usr/bin/env bash
# Shared by execution entrypoints; source this file, then call photon_setup_runtime.
photon_dispatch()
{
  local module=$1 entry=$2
  shift 2
  if [[ -z ${PHOTON_SOFTWARE_ROOT:-} && -z ${PHOTON_PRODUCTION_MANIFEST:-} ]]; then return; fi
  local selected
  selected=$(python3 "$module/tools/software.py" resolve) || return
  local fields
  mapfile -t fields <<< "$selected"
  local frozen_module=${fields[2]}
  if [[ $(readlink -f "$module") != $(readlink -f "$frozen_module") ]]; then
    exec bash "$frozen_module/$entry" "$@"
  fi
}

photon_setup_runtime()
{
  local module=$1 selected
  selected=$(python3 "$module/tools/software.py" resolve) || return
  local fields
  mapfile -t fields <<< "$selected"
  export PHOTON_RELEASE=${fields[0]} PHOTON_INSTALL_PREFIX=${fields[1]} PHOTON_SOURCE_ROOT=${fields[2]}
  set +u
  source /opt/sphenix/core/bin/sphenix_setup.sh -n "$PHOTON_RELEASE"
  set -u
  export LD_LIBRARY_PATH="$PHOTON_INSTALL_PREFIX/lib64:$PHOTON_INSTALL_PREFIX/dependencies/Pi0Reconstruction/lib:${LD_LIBRARY_PATH:-}"
  export ROOT_INCLUDE_PATH="$PHOTON_INSTALL_PREFIX/include:$PHOTON_INSTALL_PREFIX/dependencies/Pi0Reconstruction/include/pi0reconstruction:${ROOT_INCLUDE_PATH:-}"
}

photon_check_output()
{
  python3 "$PHOTON_SOURCE_ROOT/tools/software.py" check-output "$1"
}

photon_check_inputs()
{
  python3 "$PHOTON_SOURCE_ROOT/tools/software.py" check-inputs "$1" "$2"
}

photon_record_execution()
{
  python3 "$PHOTON_SOURCE_ROOT/tools/software.py" record "$@"
}
