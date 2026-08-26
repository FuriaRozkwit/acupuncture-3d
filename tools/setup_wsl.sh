#!/usr/bin/env bash
# Srodowisko uruchomieniowe w WSL2 (Ubuntu) - obejscie blokady Smart App Control.
#
# Windows z wlaczonym Smart App Control blokuje niepodpisane biblioteki VTK
# (dziennik CodeIntegrity, Id 3077: "did not meet the Enterprise signing level
# requirements"). Binaria linuksowe w WSL2 tej polityce nie podlegaja, a WSLg
# daje okno graficzne, wiec aplikacja dziala bez oslabiania zabezpieczen systemu.
#
# Uruchomienie z Windows:
#   wsl -d Ubuntu-24.04 -u root bash /mnt/c/Users/<user>/repos/acupuncture-3d/tools/setup_wsl.sh
set -euo pipefail

VENV="${VENV:-/opt/acu3d-venv}"

echo "== pakiety systemowe (OpenGL i X dla VTK) =="
export DEBIAN_FRONTEND=noninteractive
apt-get update -qq
apt-get install -y -qq --no-install-recommends \
    python3-venv python3-dev \
    libgl1 libglx-mesa0 libegl1 \
    libxrender1 libxext6 libsm6 libxt6 libx11-6 libxcursor1 libxinerama1 \
    libglib2.0-0t64 libgomp1 >/dev/null

echo "== srodowisko Pythona w $VENV =="
python3 -m venv "$VENV"
"$VENV/bin/pip" install --quiet --upgrade pip
"$VENV/bin/pip" install --quiet numpy "pyvista>=0.48"

echo "== kontrola =="
"$VENV/bin/python" - <<'PY'
import pyvista as pv
import numpy as np
from vtkmodules.vtkCommonCore import vtkVersion
print(f"pyvista {pv.__version__}, vtk {vtkVersion.GetVTKVersion()}, numpy {np.__version__}")
mesh = pv.Sphere()
print(f"siatka testowa: {mesh.n_cells} trojkatow")
PY

echo
echo "gotowe. Uruchomienie aplikacji:"
echo "  wsl -d Ubuntu-24.04 -u root bash -lc '$VENV/bin/python -m acu3d' "
echo "  (z katalogu projektu; patrz tools/run_wsl.sh)"
