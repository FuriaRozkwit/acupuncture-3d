"""Kontrola spojnosci atlasu: czy wszystkie punkty sie rozwiazuja, leza na skorze
i czy sasiednie punkty kanalu nie stoja w nieprawdopodobnych odleglosciach.

    python tools/validate.py [--height 1.6]
"""
from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pyvista as pv

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
warnings.simplefilter("ignore")

from acu3d.model import Atlas  # noqa: E402

MAX_GAP_CM = 22.0        # sasiedzi w kanale dalej niz to = podejrzane
MIN_GAP_CM = 0.4         # ...i zbyt blisko siebie
MAX_SKIN_MM = 20.0       # marker ma lezec na skorze; 20 mm to zapas na rozstaw
#                          wierzcholkow siatki (ok. 6 mm) i odsuniecie markera


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--height", type=float, default=None)
    ap.add_argument("--gap", type=float, default=MAX_GAP_CM)
    args = ap.parse_args()

    atlas = Atlas(args.height)
    problems = 0
    print(f"wzrost modelu: {atlas.rig.H:.2f} m   punktow: {len(atlas.points)}   "
          f"markerow: {len(atlas.markers)}")

    pos = np.array([m["pos"] for m in atlas.markers])
    if np.isnan(pos).any():
        print("BLAD: NaN we wspolrzednych")
        problems += 1

    # siatka skory z BodyParts3D jest otwarta (oczodoly, usta), wiec test
    # "wewnatrz bryly" daje falszywe alarmy - mierzymy odleglosc od powierzchni
    skin = atlas.layers["skin"]
    dist = np.array([float(np.linalg.norm(np.asarray(skin.points)[skin.find_closest_point(p)] - p))
                     for p in pos])
    far = dist > MAX_SKIN_MM / 1000.0
    if far.any():
        problems += int(far.sum())
        order = np.argsort(-dist)
        print(f"BLAD: {int(far.sum())} markerow dalej niz {MAX_SKIN_MM:.0f} mm od skory: "
              + ", ".join(f"{atlas.markers[i]['key']} ({dist[i] * 1000:.0f} mm)"
                          for i in order[:12]))
    else:
        print(f"OK: wszystkie markery leza na skorze "
              f"(najdalszy {dist.max() * 1000:.1f} mm, sredni {dist.mean() * 1000:.1f} mm)")

    asym = [m["code"] for m in atlas.markers if m["side"] == "R"
            and not np.allclose(m["pos"] * [-1, 1, 1], atlas.by_key[m["code"] + "·L"]["pos"])]
    if asym:
        problems += len(asym)
        print(f"BLAD: niesymetryczne pary: {', '.join(asym)}")
    else:
        print("OK: strony L/R symetryczne")

    bad = []
    for ch in atlas.channels:
        codes = sorted((p for p in atlas.points if p["channel"] == ch), key=lambda r: r["index"])
        for a, b in zip(codes, codes[1:]):
            ka = atlas.by_key.get(a["code"], atlas.by_key.get(a["code"] + "·R"))
            kb = atlas.by_key.get(b["code"], atlas.by_key.get(b["code"] + "·R"))
            d = np.linalg.norm(ka["pos"] - kb["pos"]) * 100
            if d > args.gap or d < MIN_GAP_CM:
                bad.append(f"{a['code']}->{b['code']} {d:.1f} cm")
    if bad:
        print(f"UWAGA: {len(bad)} nietypowych odstepow miedzy sasiadami "
              f"(>{args.gap:.0f} cm lub <{MIN_GAP_CM} cm):")
        for b in bad:
            print("   " + b)
    else:
        print("OK: odstepy miedzy sasiednimi punktami w normie")

    acc = {}
    for p in atlas.points:
        acc[p["accuracy"]] = acc.get(p["accuracy"], 0) + 1
    print("dokladnosc kotwic: " + ", ".join(f"{k}={v}" for k, v in sorted(acc.items())))
    print("WYNIK:", "OK" if problems == 0 else f"{problems} bledow")
    return 0 if problems == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
