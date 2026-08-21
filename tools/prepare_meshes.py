"""Przygotowanie siatek Z-Anatomy do uzycia w przegladarce.

Wejscie:  assets/mesh/        (surowy eksport z Blendera, PLY w ukladzie Z-Anatomy)
Wyjscie:  assets/ready/       (siatki po obrocie i decymacji, format .vtp)
          data/muscles.json   (katalog: nazwa, strona, os dzialania, objetosc, plik)

Uklad Z-Anatomy: +X = strona lewa, +Y = tyl, +Z = gora, metry.
Uklad acu3d:     +X = strona prawa, +Y = przod, +Z = gora, metry.
Roznica to obrot o 180 stopni wokol osi Z (x -> -x, y -> -y) - bez odbicia,
wiec strony ciala pozostaja stronami, a nie swoim lustrem.

    python tools/prepare_meshes.py [--max-polys 4000] [--skeleton-polys 220000]
"""
from __future__ import annotations

import argparse
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pyvista as pv

warnings.simplefilter("ignore")
ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "assets" / "mesh"
OUT = ROOT / "assets" / "ready"

# czesci, ktorych nie chcemy w atlasie: kaletki, powiezi, sciegna jako osobne bryly
SKIP = ("bursa", "fascia", "aponeurosis", "sheath", "retinaculum", "raphe",
        "trochlea", "ligament", "septum", "system.g", "tendinous arch")


def rotate_z180(mesh: pv.PolyData) -> pv.PolyData:
    mesh.points[:, 0] *= -1.0
    mesh.points[:, 1] *= -1.0
    return mesh


def principal_axis(points: np.ndarray):
    """Os glowna chmury punktow - kierunek, wzdluz ktorego miesien sie kurczy."""
    c = points.mean(axis=0)
    u, s, vt = np.linalg.svd(points - c, full_matrices=False)
    axis = vt[0]
    proj = (points - c) @ axis
    if abs(axis[2]) > 1e-9 and axis[2] < 0:      # os zawsze "w gore" dla powtarzalnosci
        axis, proj = -axis, -proj
    return c, axis, float(proj.min()), float(proj.max()), float(s[1] / max(s[0], 1e-9))


def decimate(mesh: pv.PolyData, target: int) -> pv.PolyData:
    if mesh.n_cells <= target:
        return mesh
    return mesh.decimate_pro(1.0 - target / mesh.n_cells, feature_angle=45.0, preserve_topology=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-polys", type=int, default=4000, help="limit trojkatow na miesien")
    ap.add_argument("--skeleton-polys", type=int, default=220000)
    args = ap.parse_args()

    cat_path = RAW / "catalog.json"
    if not cat_path.exists():
        print(f"brak {cat_path} - najpierw uruchom blender_export.py", file=sys.stderr)
        return 2
    cat = json.loads(cat_path.read_text(encoding="utf-8"))
    (OUT / "muscles").mkdir(parents=True, exist_ok=True)

    muscles, skipped = [], 0
    total_in = total_out = 0
    for i, part in enumerate(cat["parts"], 1):
        low = part["name"].lower()
        # ".g" to znaczniki grup odstawione poza ciało, nie anatomia
        if any(s in low for s in SKIP) or part["raw"].endswith(".g")                 or abs(part["centroid"][0]) > 0.45:
            skipped += 1
            continue
        mesh = pv.read(RAW / part["file"])
        total_in += mesh.n_cells
        mesh = rotate_z180(mesh)
        mesh = decimate(mesh.triangulate(), args.max_polys)
        total_out += mesh.n_cells
        c, axis, lo, hi, flat = principal_axis(np.asarray(mesh.points))
        fn = Path(part["file"]).with_suffix(".vtp").name
        mesh.save(OUT / "muscles" / fn, binary=True)
        muscles.append({
            "id": Path(fn).stem,
            "name": part["name"],
            "side": part["side"],
            "file": f"muscles/{fn}",
            "polys": int(mesh.n_cells),
            "centroid": [round(float(x), 5) for x in c],
            "axis": [round(float(x), 5) for x in axis],
            "extent": [round(lo, 5), round(hi, 5)],
            "flatness": round(flat, 3),
            "volume_l": round(float(abs(mesh.volume)) * 1000, 4),
        })
        if i % 120 == 0:
            print(f"  ...{i}/{len(cat['parts'])}")

    sk = pv.read(RAW / cat["groups"]["skeleton"]["file"]).triangulate()
    sk = rotate_z180(sk)
    keep = np.abs(sk.cell_centers().points[:, 0]) < 0.45      # odsiew znacznikow grup
    if not keep.all():
        sk = sk.extract_cells(np.flatnonzero(keep)).extract_surface()
    sk = decimate(sk, args.skeleton_polys)
    sk.save(OUT / "skeleton.vtp", binary=True)

    doc = {
        "source": "Z-Anatomy (CC BY-SA 4.0), eksport i obrobka: tools/blender_export.py + prepare_meshes.py",
        "coordinate_system": "+X = prawa strona, +Y = przod, +Z = gora, metry",
        "skeleton": {"file": "skeleton.vtp", "polys": int(sk.n_cells)},
        "muscles": sorted(muscles, key=lambda m: m["name"]),
    }
    (ROOT / "data" / "muscles.json").write_text(
        json.dumps(doc, ensure_ascii=False, indent=1), encoding="utf-8")

    print(f"miesnie: {len(muscles)} plikow (pominieto {skipped} kaletek/powiezi)")
    print(f"trojkaty miesni: {total_in} -> {total_out}")
    print(f"szkielet: {sk.n_cells} trojkatow")
    print(f"katalog -> data/muscles.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
