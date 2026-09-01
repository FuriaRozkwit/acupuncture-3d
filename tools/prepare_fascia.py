"""Warstwa lacznotkankowa: powiezi, przegrody, troczki, rozciegna, pochewki, kaletki.

Te struktury byly swiadomie pomijane w prepare_meshes.py, bo zasmiecaly wybor
miesnia klikiem. Do mapy napiec sa jednak potrzebne: terapeuta palpuje wlasnie
te tkanki, a przegrody miedzymiesniowe wyznaczaja bruzdy na skorze.

Wejscie:  assets/mesh/       (surowy eksport z Blendera)
Wyjscie:  assets/ready/fascia/          pojedyncze siatki .vtp
          assets/ready/fascia_merged.vtp   jedna siatka z tablica `fid`
          data/fascia.json               katalog

    python tools/prepare_fascia.py [--max-polys 2500]
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

# klasa tkanki -> wzorce w nazwie; kolejnosc ma znaczenie (pierwsze trafienie wygrywa)
CLASSES = (
    ("septum",      ("septum",)),                  # przegroda miedzymiesniowa
    ("retinaculum", ("retinaculum",)),             # troczek
    ("aponeurosis", ("aponeurosis", "raphe", "tendinous arch")),
    ("sheath",      ("sheath",)),                  # pochewka sciegnista
    ("bursa",       ("bursa",)),                   # kaletka
    ("fascia",      ("fascia",)),                  # powiez wlasciwa
)
PL = {"septum": "przegroda miedzymiesniowa", "retinaculum": "troczek",
      "aponeurosis": "rozciegno", "sheath": "pochewka sciegnista",
      "bursa": "kaletka", "fascia": "powiez"}


def classify(name: str) -> str | None:
    low = name.lower()
    for cls, pats in CLASSES:
        if any(p in low for p in pats):
            return cls
    return None


def rotate_z180(mesh: pv.PolyData) -> pv.PolyData:
    mesh.points[:, 0] *= -1.0
    mesh.points[:, 1] *= -1.0
    return mesh


def sheet_axes(points: np.ndarray):
    """Powiez to plat, nie walec: interesuje nas jego NORMALNA i rozpietosc.

    Zwraca (srodek, normalna, grubosc, rozpietosc) - normalna to kierunek
    najmniejszej wariancji, czyli w poprzek platu.
    """
    c = points.mean(axis=0)
    u, s, vt = np.linalg.svd(points - c, full_matrices=False)
    n = vt[2]
    if n[0] < 0 or (abs(n[0]) < 1e-9 and n[1] < 0):
        n = -n                                    # znak umowny, dla powtarzalnosci
    proj = (points - c) @ n
    return c, n, float(proj.max() - proj.min()), float(s[0] / max(len(points), 1) ** 0.5)


def decimate(mesh: pv.PolyData, target: int) -> pv.PolyData:
    if mesh.n_cells <= target:
        return mesh
    return mesh.decimate_pro(1.0 - target / mesh.n_cells, feature_angle=45.0,
                             preserve_topology=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-polys", type=int, default=2500)
    args = ap.parse_args()

    cat_path = RAW / "catalog.json"
    if not cat_path.exists():
        print(f"brak {cat_path} - najpierw uruchom blender_export.py", file=sys.stderr)
        return 2
    cat = json.loads(cat_path.read_text(encoding="utf-8"))
    (OUT / "fascia").mkdir(parents=True, exist_ok=True)

    recs, parts = [], []
    total_in = total_out = 0
    for part in cat["parts"]:
        if part["raw"].endswith(".g") or abs(part["centroid"][0]) > 0.45:
            continue
        cls = classify(part["name"])
        if cls is None:
            continue
        mesh = pv.read(RAW / part["file"]).triangulate()
        if mesh.n_points < 4:
            continue
        total_in += mesh.n_cells
        mesh = decimate(rotate_z180(mesh), args.max_polys)
        total_out += mesh.n_cells
        c, n, thick, span = sheet_axes(np.asarray(mesh.points))
        fn = Path(part["file"]).with_suffix(".vtp").name
        mesh.save(OUT / "fascia" / fn, binary=True)
        recs.append({
            "id": Path(fn).stem,
            "name": part["name"],
            "class": cls,
            "side": part["side"],
            "file": f"fascia/{fn}",
            "polys": int(mesh.n_cells),
            "centroid": [round(float(x), 5) for x in c],
            "normal": [round(float(x), 5) for x in n],
            "thickness": round(thick, 5),
            "span": round(span, 5),
        })
        parts.append(mesh)

    order = sorted(range(len(recs)), key=lambda i: (recs[i]["class"], recs[i]["name"]))
    recs = [recs[i] for i in order]
    parts = [parts[i] for i in order]
    for i, (r, m) in enumerate(zip(recs, parts)):
        m.point_data["fid"] = np.full(m.n_points, i, dtype=np.int32)
        m.cell_data["fid"] = np.full(m.n_cells, i, dtype=np.int32)
        r["index"] = i
    merged = parts[0].merge(parts[1:]) if len(parts) > 1 else parts[0]
    merged.save(OUT / "fascia_merged.vtp", binary=True)

    counts = {}
    for r in recs:
        counts[r["class"]] = counts.get(r["class"], 0) + 1
    doc = {
        "source": "Z-Anatomy (CC BY-SA 4.0), eksport: tools/blender_export.py + prepare_fascia.py",
        "coordinate_system": "+X = prawa strona, +Y = przod, +Z = gora, metry",
        "note": "normal = kierunek najmniejszej wariancji, czyli w poprzek platu",
        "classes": {k: PL[k] for k in counts},
        "counts": counts,
        "merged": {"file": "fascia_merged.vtp", "polys": int(merged.n_cells)},
        "items": recs,
    }
    (ROOT / "data" / "fascia.json").write_text(
        json.dumps(doc, ensure_ascii=False, indent=1), encoding="utf-8")

    for k in sorted(counts):
        print(f"  {counts[k]:4d}  {PL[k]}")
    print(f"razem: {len(recs)} struktur, trojkaty {total_in} -> {total_out}")
    print(f"katalog -> data/fascia.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
