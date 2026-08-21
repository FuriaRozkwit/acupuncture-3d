"""Eksport rozwiazanych wspolrzednych 3D (dla innych narzedzi, np. fiszek).

    python tools/export_xyz.py --height 1.72 -o data/resolved_1.72.json
    python tools/export_xyz.py --csv punkty.csv

Kotwice sa zrodlem prawdy - ten plik to migawka dla konkretnego wzrostu.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
warnings.simplefilter("ignore")

from acu3d.model import Atlas  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--height", type=float, default=None, help="wzrost w metrach")
    ap.add_argument("-o", "--out", default=None, help="plik JSON")
    ap.add_argument("--csv", default=None, help="plik CSV")
    args = ap.parse_args()

    atlas = Atlas(args.height)
    rows = [{
        "key": m["key"], "code": m["code"], "channel": m["channel"], "index": m["index"],
        "side": m["side"], "pinyin": m["pinyin"], "english": m["english"],
        "accuracy": m["accuracy"],
        "x": round(float(m["pos"][0]), 4),
        "y": round(float(m["pos"][1]), 4),
        "z": round(float(m["pos"][2]), 4),
        "note_pl": m["note_pl"],
    } for m in atlas.markers]

    out = Path(args.out or f"data/resolved_{atlas.rig.H:.2f}.json")
    out.write_text(json.dumps({
        "height_m": atlas.rig.H,
        "coordinate_system": "+X = prawa strona, +Y = przod, +Z = gora, metry, 0 = podloga",
        "points": rows,
    }, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"{len(rows)} markerow -> {out}")

    if args.csv:
        with open(args.csv, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        print(f"{len(rows)} markerow -> {args.csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
