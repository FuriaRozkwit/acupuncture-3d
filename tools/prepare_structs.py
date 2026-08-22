"""Drobne struktury anatomiczne -> assets/ready/struct/*.vtp + data/structures.json.

Kosci srodrecza i srodstopia, paliczki kazdego palca, kosci nadgarstka i stepu,
elementy malzowiny, gałka oczna. To sa punkty odniesienia, ktorych uzywa podrecznik
przy punktach na palcach, uchu i wokol oka - dotad kotwiczonych "na oko".

    python tools/prepare_structs.py
"""
from __future__ import annotations

import json
import warnings
from pathlib import Path

import numpy as np
import pyvista as pv

warnings.simplefilter("ignore")
ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "assets" / "mesh" / "struct"
OUT = ROOT / "assets" / "ready" / "struct"


def principal_axis(points: np.ndarray):
    c = points.mean(axis=0)
    _, s, vt = np.linalg.svd(points - c, full_matrices=False)
    axis = vt[0]
    proj = (points - c) @ axis
    if proj.max() + proj.min() < 0:            # os w strone "dalsza" dla powtarzalnosci
        axis, proj = -axis, -proj
    return c, axis, float(proj.min()), float(proj.max())


def main() -> int:
    if not RAW.exists():
        raise SystemExit(f"brak {RAW} - uruchom najpierw blender_export_named.py")
    OUT.mkdir(parents=True, exist_ok=True)
    entries = {}
    for f in sorted(RAW.glob("*.ply")):
        m = pv.read(f).triangulate()
        p = np.asarray(m.points)
        p[:, 0] *= -1.0                        # uklad Z-Anatomy -> acu3d
        p[:, 1] *= -1.0
        m.points = p
        m.save(OUT / f"{f.stem}.vtp", binary=True)
        c, axis, lo, hi = principal_axis(p)
        entries[f.stem] = {
            "file": f"struct/{f.stem}.vtp",
            "polys": int(m.n_cells),
            "centroid": [round(float(x), 5) for x in c],
            "axis": [round(float(x), 5) for x in axis],
            "extent": [round(lo, 5), round(hi, 5)],
            "bounds": [round(float(b), 5) for b in m.bounds],
        }
    (ROOT / "data" / "structures.json").write_text(
        json.dumps({"source": "Z-Anatomy (CC BY-SA 4.0)",
                    "coordinate_system": "+X = prawa strona, +Y = przod, +Z = gora",
                    "structures": entries}, ensure_ascii=False, indent=1),
        encoding="utf-8")
    print(f"struktur: {len(entries)} -> data/structures.json")
    for k in list(entries)[:3]:
        print(f"  {k}: {entries[k]['polys']} trojkatow")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
