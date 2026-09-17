"""Eksport pakietu dla klienta mobilnego (Terrapka, renderer meridian-map).

    python tools/export_terrapka.py --height 1.72 --body anatomy -o dist/meridian-bundle.json

Pakiet jest samowystarczalny: punkty, przebiegi kanalow i siatka ciala w jednym
pliku JSON. Wspolrzedne sa w metrach, uklad +X = prawa strona, +Y = przod, +Z = gora,
0 = podloga - ten sam, w ktorym pracuje reszta projektu.

Siatka ciala idzie jako bufory base64 (Float32 / Uint16 / Uint32), bo tablice JSON
z 30 tys. liczb waza kilka razy tyle co te same dane binarnie, a three.js i tak
potrzebuje ich jako typed arrays.

UWAGA LICENCYJNA. Zawartosc pakietu zalezy od --body:

    procedural  geometria z wzorow (acu3d/geometry.py) - kod wlasny, MIT
    anatomy     geometria zmierzona na siatkach Z-Anatomy i BodyParts3D,
                a wiec CC BY-SA: wymaga atrybucji w aplikacji i share-alike
                na sam plik pakietu

Opisy z podrecznika (data/point_meta.json) NIE trafiaja do pakietu w zadnym trybie -
to material chroniony. Jedyne opisy w pakiecie to wlasne note_pl z data/points.json.
"""
from __future__ import annotations

import argparse
import base64
import json
import sys
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
warnings.simplefilter("ignore")

from acu3d.model import Atlas  # noqa: E402

LICENSES = {
    "procedural": {
        "geometry": "MIT",
        "note": "Geometria proceduralna i kotwice punktow - praca wlasna.",
        "attribution": [],
    },
    "anatomy": {
        "geometry": "CC BY-SA 4.0",
        "note": ("Geometria wyliczona z siatek anatomicznych na licencji CC BY-SA. "
                 "Pakiet musi byc rozpowszechniany na tej samej licencji, "
                 "z atrybucja widoczna w aplikacji."),
        "attribution": [
            "BodyParts3D, (c) The Database Center for Life Science, "
            "licensed under CC Attribution-Share Alike 2.1 Japan.",
            "Z-Anatomy, (c) Gauthier Kervyn, "
            "licensed under CC Attribution-Share Alike 4.0 International.",
        ],
    },
}


def b64(arr: np.ndarray) -> str:
    return base64.b64encode(np.ascontiguousarray(arr).tobytes()).decode("ascii")


def thin(line: np.ndarray, keep_every: int) -> np.ndarray:
    """Przerzedza probki krzywej, zawsze zostawiajac oba konce."""
    if keep_every <= 1 or len(line) <= 2:
        return line
    idx = list(range(0, len(line), keep_every))
    if idx[-1] != len(line) - 1:
        idx.append(len(line) - 1)
    return line[idx]


def export_body(atlas: Atlas, target_tris: int) -> dict:
    skin = atlas.layers["skin"].triangulate() if atlas.mode == "procedural" \
        else atlas.anatomy.skin.triangulate()
    tris_before = skin.n_cells
    if target_tris and tris_before > target_tris:
        skin = skin.decimate(1.0 - target_tris / tris_before).triangulate().clean()
    faces = skin.faces.reshape(-1, 4)[:, 1:]
    pos = np.asarray(skin.points, dtype=np.float32)
    idx_dtype = np.uint16 if len(pos) < 65536 else np.uint32
    return {
        "encoding": "base64",
        "positionType": "Float32Array",
        "indexType": "Uint16Array" if idx_dtype is np.uint16 else "Uint32Array",
        "vertexCount": int(len(pos)),
        "triangleCount": int(len(faces)),
        "trianglesBeforeDecimation": int(tris_before),
        "positions": b64(pos),
        "indices": b64(faces.astype(idx_dtype)),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--height", type=float, default=1.72, help="wzrost modelu w metrach")
    ap.add_argument("--body", choices=("anatomy", "procedural"), default="anatomy")
    ap.add_argument("--target-tris", type=int, default=8000,
                    help="docelowa liczba trojkatow skory (0 = bez decymacji)")
    ap.add_argument("--keep-every", type=int, default=3,
                    help="co ktora probke krzywej kanalu zostawic")
    ap.add_argument("-o", "--out", default="dist/meridian-bundle.json")
    args = ap.parse_args()

    atlas = Atlas(args.height, procedural=(args.body == "procedural"))
    if atlas.mode != args.body:
        print(f"! zadano --body {args.body}, a Atlas wystartowal w trybie {atlas.mode}",
              file=sys.stderr)

    points = [{
        "key": m["key"], "code": m["code"], "channel": m["channel"], "index": m["index"],
        "side": m["side"], "pinyin": m["pinyin"], "english": m["english"],
        "notePl": m["note_pl"], "accuracy": m["accuracy"],
        "p": [round(float(v), 4) for v in m["pos"]],
    } for m in atlas.markers]

    paths: dict[str, dict[str, list]] = {}
    samples = 0
    for ch in atlas.channels:
        for side in ("R", "L"):
            lines = [thin(l, args.keep_every) for l in atlas.meridian_lines(ch, side)]
            lines = [l for l in lines if len(l) >= 2]
            if not lines:
                continue
            samples += sum(len(l) for l in lines)
            paths.setdefault(ch, {})[side] = [
                [[round(float(v), 4) for v in p] for p in l] for l in lines
            ]

    bundle = {
        "version": 1,
        "generator": "acupuncture-3d tools/export_terrapka.py",
        "bodySource": atlas.mode,
        "heightM": round(float(atlas.rig.H), 3),
        "coordinateSystem": "+X = prawa strona, +Y = przod, +Z = gora; metry, 0 = podloga",
        "license": LICENSES[atlas.mode],
        "channels": atlas.channels,
        "points": points,
        "paths": paths,
        "body": export_body(atlas, args.target_tris),
    }

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(bundle, ensure_ascii=False, separators=(",", ":")),
                   encoding="utf-8")

    kb = out.stat().st_size / 1024
    print(f"tryb ciala:  {atlas.mode}  ({bundle['license']['geometry']})")
    print(f"punkty:      {len(points)}")
    print(f"kanaly:      {len(paths)}, probek krzywych: {samples}")
    print(f"skora:       {bundle['body']['vertexCount']} wierzcholkow, "
          f"{bundle['body']['triangleCount']} trojkatow "
          f"(z {bundle['body']['trianglesBeforeDecimation']})")
    print(f"-> {out}  {kb:.0f} kB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
