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

# Auto-konfiguracja i samonaprawa gniazd WSLg dla usera root
if [ ! -d "/run/user/0" ]; then
    mkdir -p /run/user/0
    chmod 700 /run/user/0
fi
if [ ! -e "/run/user/0/wayland-0" ] && [ -e "/mnt/wslg/runtime-dir/wayland-0" ]; then
    ln -sf /mnt/wslg/runtime-dir/wayland-0 /run/user/0/wayland-0
    ln -sf /mnt/wslg/runtime-dir/wayland-0.lock /run/user/0/wayland-0.lock 2>/dev/null || true
fi
if [ ! -d "/tmp/.X11-unix" ] && [ -d "/mnt/wslg/.X11-unix" ]; then
    ln -sf /mnt/wslg/.X11-unix /tmp/.X11-unix
fi

export DISPLAY="${DISPLAY:-:0}"
export WAYLAND_DISPLAY="${WAYLAND_DISPLAY:-wayland-0}"
export XDG_RUNTIME_DIR="${XDG_RUNTIME_DIR:-/run/user/0}"

cd "$HERE"
exec "$VENV/bin/python" -m acu3d "$@"
