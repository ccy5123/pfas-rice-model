#!/bin/bash
# =============================================================================
# SessionStart hook -- Claude Code on the web
# -----------------------------------------------------------------------------
# Provisions a fresh web container so the model, tests and the live HYDRUS-1D
# soil-coupling mode work out of the box:
#   1. the Python stack (numpy/scipy/.../rdkit/phydrus from requirements.txt) + pytest
#   2. the optional HYDRUS-1D FORTRAN engine (gfortran build of the vendored source)
#
# The HYDRUS build is BEST-EFFORT and must never block the session: the analytic
# `cwo_profile='flooded'` exposure path works without the engine, so a build
# failure just leaves the live "Run HYDRUS-1D" mode unavailable.
#
# Synchronous (no async): guarantees deps are present before the agent runs, so
# tests/linters never race an unfinished install. Switch to async for faster
# startup if you prefer (see the session-start-hook skill).
# =============================================================================
set -euo pipefail

# Web only -- local sessions already have the developer's environment.
if [ "${CLAUDE_CODE_REMOTE:-}" != "true" ]; then
  exit 0
fi

cd "${CLAUDE_PROJECT_DIR:-$(pwd)}"

# Use sudo only when not already root.
SUDO=""
[ "$(id -u)" -ne 0 ] && SUDO="sudo"

# 1. Python stack (model + Streamlit app + tests). requirements.txt pulls in
#    rdkit (SMILES input) and phydrus (HYDRUS driver); pip caches make re-runs
#    fast once the container is cached.
python -m pip install -q -r requirements.txt
python -m pip install -q pytest

# 2. HYDRUS-1D engine -- OPTIONAL. Best-effort; explicit early returns so a
#    failure is caught by the `|| echo` below and does NOT abort the session.
build_hydrus() {
  local src="external/hydrus_source/source"
  [ -f "$src/HYDRUS.FOR" ] || { echo "session-start: vendored HYDRUS source absent -> skip build"; return 0; }
  [ -x "$src/hydrus" ]      && { echo "session-start: HYDRUS engine already built"; return 0; }
  if ! command -v gfortran >/dev/null 2>&1; then
    export DEBIAN_FRONTEND=noninteractive
    $SUDO apt-get update -qq        || return 1
    $SUDO apt-get install -y -qq gfortran || return 1
  fi
  cp external/hydrus_source/makefile "$src/" || return 1
  ( cd "$src" && make )                        || return 1
  echo "session-start: HYDRUS engine built"
}
build_hydrus || echo "session-start: HYDRUS build skipped/failed (optional; cwo_profile='flooded' still works)" >&2

echo "session-start: setup complete"
