"""Skora z BodyParts3D dopasowana do anatomii Z-Anatomy -> assets/ready/skin.vtp.

Z-Anatomy nie zawiera powloki ciala, wiec skore bierzemy z BodyParts3D - tego samego
modelu, z ktorego Z-Anatomy sie wywodzi. Oba zbiory sa w innych ukladach odniesienia,
wiec przesuniecie liczymy, a nie zgadujemy: porownujemy polozenie TEJ SAMEJ kosci
(kosc krzyzowa i zuchwa) w obu zbiorach. Skala wychodzi 1.0, wiec wystarczy translacja.

    python tools/prepare_skin.py [--target-polys 90000]
"""
from __future__ import annotations

import argparse
import warnings
import zipfile
from pathlib import Path

import numpy as np
import pyvista as pv

warnings.simplefilter("ignore")
ROOT = Path(__file__).resolve().parent.parent
DL = ROOT / "assets" / "download"
OBJ = DL / "isa_BP3D_4.0_obj_99"
REF = ROOT / "assets" / "mesh" / "ref"
READY = ROOT / "assets" / "ready"

SKIN = "FJ2810.obj"                       # FMA7163 skin
PAIRS = [("sacrum.ply", "FJ3393.obj"), ("mandible.ply", "FJ3289.obj")]


def rot180(pts: np.ndarray) -> np.ndarray:
    """Z ukladu Z-Anatomy/BP3D (+X lewa, +Y tyl) do acu3d (+X prawa, +Y przod)."""
    out = pts.copy()
    out[:, 0] *= -1.0
    out[:, 1] *= -1.0
    return out


def bp_mesh(name: str) -> pv.PolyData:
    path = OBJ / name
    if not path.exists():                 # wypakuj pojedynczy element z archiwum
        zf = DL / "bp3d_isa.zip"
        if not zf.exists():
            raise SystemExit(f"brak {zf} - uruchom najpierw tools/fetch_assets.py")
        with zipfile.ZipFile(zf) as z:
            inner = [n for n in z.namelist() if n.endswith(name)][0]
            z.extract(inner, DL)
    m = pv.read(path).triangulate()
    m.points /= 1000.0                    # BodyParts3D zapisuje w milimetrach
    m.points = rot180(np.asarray(m.points))
    return m


def bbox_center(mesh: pv.PolyData) -> np.ndarray:
    b = np.array(mesh.bounds).reshape(3, 2)
    return b.mean(axis=1)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target-polys", type=int, default=90000)
    args = ap.parse_args()

    offsets = []
    for zname, bname in PAIRS:
        zpath = REF / zname
        if not zpath.exists():
            print(f"  pomijam {zname} - brak pliku wzorcowego")
            continue
        z = pv.read(zpath)
        z.points = rot180(np.asarray(z.points))
        d = bbox_center(z) - bbox_center(bp_mesh(bname))
        offsets.append(d)
        print(f"  {zname[:-4]:10s} przesuniecie {np.round(d, 4)}")
    if not offsets:
        raise SystemExit("brak kosci wzorcowych - uruchom blender_export_named.py")
    off = np.mean(offsets, axis=0)
    off[0] = 0.0                          # os posrodkowa: przesuniecie w bok = 0
    print(f"  przyjete przesuniecie: {np.round(off, 4)} m")

    skin = bp_mesh(SKIN)
    skin.points += off
    before = skin.n_cells
    if skin.n_cells > args.target_polys:
        skin = skin.decimate_pro(1 - args.target_polys / skin.n_cells,
                                 feature_angle=45, preserve_topology=True)
    READY.mkdir(parents=True, exist_ok=True)
    skin.save(READY / "skin.vtp", binary=True)
    print(f"  skora: {before} -> {skin.n_cells} trojkatow, "
          f"zakres {[round(b, 3) for b in skin.bounds]}")
    print(f"  zapisano {READY / 'skin.vtp'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
