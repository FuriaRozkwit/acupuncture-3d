#!/usr/bin/env bash
# Uruchomienie atlasu w WSL2 z oknem graficznym przez WSLg.
# Powod: Smart App Control na Windows blokuje niepodpisane biblioteki VTK.
#
#   wsl -d Ubuntu-24.04 -u root bash /mnt/c/.../tools/run_wsl.sh [opcje acu3d]
set -euo pipefail

VENV="${VENV:-/opt/acu3d-venv}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if [ ! -x "$VENV/bin/python" ]; then
    echo "brak srodowiska w $VENV - uruchom najpierw tools/setup_wsl.sh" >&2
    exit 2
fi
cd "$HERE"
exec "$VENV/bin/python" -m acu3d "$@"
