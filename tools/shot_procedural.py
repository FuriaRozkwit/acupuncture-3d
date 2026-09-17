"""Zrzuty modelu proceduralnego do ogladania (bez okna).

    python tools/shot_procedural.py -o docs/proc

Sluzy do sprawdzenia struktur z acu3d/structures_proc.py: czy punkty palcow,
stopy i twarzy siadaja tam, gdzie powinny. Anatomiczny model nie jest tu dobrym
wzorcem - ma dlon i stope w innej pozycji - wiec kontrola jest wzrokowa.
"""
from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
warnings.simplefilter("ignore")

import pyvista as pv  # noqa: E402

from acu3d.model import Atlas  # noqa: E402

VIEWS = {
    "body":  dict(pos=(0.0, -3.2, 0.95), focus=(0.0, 0.0, 0.95), up=(0, 0, 1)),
    "hand":  dict(pos=(0.75, -0.45, 0.74), focus=(0.31, 0.0, 0.74), up=(0, 0, 1)),
    "foot":  dict(pos=(0.55, 0.30, 0.16), focus=(0.07, 0.13, 0.05), up=(0, 0, 1)),
    "head":  dict(pos=(0.55, -0.75, 1.60), focus=(0.0, 0.0, 1.58), up=(0, 0, 1)),
}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--out", default="docs/proc")
    ap.add_argument("--height", type=float, default=1.72)
    ap.add_argument("--views", default="body,hand,foot,head")
    args = ap.parse_args()

    atlas = Atlas(args.height, procedural=True)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    struct_codes = set()
    for p in atlas.points:
        if "struct" in p["anchor"]:
            struct_codes.add(p["code"])

    for name in args.views.split(","):
        view = VIEWS[name]
        pl = pv.Plotter(off_screen=True, window_size=(900, 1100))
        pl.set_background("white")
        pl.add_mesh(atlas.layers["skin"], color="#e8d9c9", opacity=1.0,
                    smooth_shading=True, specular=0.2)
        plain = np.array([m["pos"] for m in atlas.markers
                          if m["code"] not in struct_codes])
        struc = np.array([m["pos"] for m in atlas.markers
                          if m["code"] in struct_codes])
        if len(plain):
            pl.add_points(plain, color="#2b6cb0", point_size=7,
                          render_points_as_spheres=True)
        if len(struc):
            pl.add_points(struc, color="#c53030", point_size=11,
                          render_points_as_spheres=True)
        pl.camera_position = [view["pos"], view["focus"], view["up"]]
        path = out.parent / f"{out.name}-{name}.png"
        pl.screenshot(str(path))
        pl.close()
        print(f"-> {path}")
    print("niebieskie = kotwice w cun, czerwone = kotwice strukturalne")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
