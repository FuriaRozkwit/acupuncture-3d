#!/usr/bin/env bash
# Przegladarka z trybem edycji punktow (gizmo + panel + profile) w WSL2 przez WSLg.
# To samo przygotowanie srodowiska co tools/run_wsl.sh, ale uruchamia acu3d.edit -
# run_wsl.sh ma na sztywno "python -m acu3d".
#
#   wsl -d Ubuntu-24.04 -u root bash tools/run_edit_wsl.sh [--fabryczna | --profil NAZWA]
set -euo pipefail

VENV="${VENV:-/opt/acu3d-venv}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if [ ! -x "$VENV/bin/python" ]; then
    echo "brak srodowiska w $VENV - uruchom najpierw tools/setup_wsl.sh" >&2
    exit 2
fi

# gniazda WSLg dla roota (jak w run_wsl.sh)
mkdir -p /run/user/0 && chmod 700 /run/user/0
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
exec "$VENV/bin/python" -m acu3d.edit "$@"
