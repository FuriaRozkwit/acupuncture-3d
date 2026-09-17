"""Porownanie trybu proceduralnego z anatomicznym.

    python tools/check_procedural.py

Punkty strukturalne (kotwica `struct`) sa w trybie proceduralnym stawiane na
odpowiednikach z `acu3d/structures_proc.py`, a nie na siatkach Z-Anatomy.
Pytanie brzmi, czy trafiaja tam, gdzie trafialy na anatomii.

Sama odleglosc w milimetrach na to nie odpowiada, bo to dwa rozne ciala: stylizowany
manekin i skan. Roznica rosnie wzdluz konczyn niezaleznie od kotwic - na tulowiu jest
rzedu 25 mm, na udzie juz 75 mm. Dlatego mierzymy dwie rzeczy:

  swiat   odleglosc w mm, zestawiona z TLEM tego samego segmentu, czyli z punktami
          opisanymi w cun, ktore w obu trybach powstaja tak samo;
  ramka   polozenie w ukladzie wlasnego segmentu (os + poprzeczka), znormalizowane
          jego dlugoscia. To jest miara samej tabeli struktur: nie zalezy od tego,
          gdzie model trzyma reke ani jak gruba ma lydke.
"""
from __future__ import annotations

import collections
import json
import sys
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
warnings.simplefilter("ignore")

from acu3d.model import Atlas  # noqa: E402

DATA = Path(__file__).resolve().parent.parent / "data"
REGION_SEG = {"hand": "hand", "foot": "foot"}


def local_frame(seg):
    """Ortonormalna ramka segmentu, liczona tak samo dla obu modeli."""
    ez = np.asarray(seg.p1, dtype=float) - np.asarray(seg.p0, dtype=float)
    length = float(np.linalg.norm(ez))
    ez = ez / length
    up = np.array([0.0, 1.0, 0.0])
    if abs(float(np.dot(up, ez))) > 0.95:
        up = np.array([0.0, 0.0, 1.0])
    ex = np.cross(up, ez)
    ex = ex / np.linalg.norm(ex)
    ey = np.cross(ez, ex)
    return np.asarray(seg.p0, dtype=float), np.vstack([ex, ey, ez]), length


def local_coords(seg, p):
    p0, basis, length = local_frame(seg)
    return basis @ (np.asarray(p, dtype=float) - p0) / length


def seg_for(anchor: dict):
    if anchor.get("region") in REGION_SEG:
        return REGION_SEG[anchor["region"]]
    return anchor.get("seg")


def main() -> int:
    pts = json.loads((DATA / "points.json").read_text(encoding="utf-8"))["points"]
    anchor_of = {p["code"]: p["anchor"] for p in pts}

    anat = Atlas(1.72)
    if anat.mode != "anatomy":
        print("brak zasobow anatomicznych - nie ma z czym porownywac", file=sys.stderr)
        return 1
    proc = Atlas(1.72, procedural=True)

    bg = collections.defaultdict(list)
    world = collections.defaultdict(list)
    local = collections.defaultdict(list)
    worst = []

    for m in anat.markers:
        if m["side"] == "L":
            continue
        o = proc.by_key.get(m["key"])
        if o is None:
            continue
        a = anchor_of[m["code"]]
        d = float(np.linalg.norm(np.asarray(m["pos"]) - np.asarray(o["pos"]))) * 1000.0
        seg = seg_for(a)
        if "struct" not in a:
            if seg:
                bg[seg].append(d)
            continue
        seg = seg or "head"
        world[seg].append(d)
        try:
            la = local_coords(anat.body.seg(seg), m["pos"])
            lp = local_coords(proc.body.seg(seg), o["pos"])
            ps = proc.body.seg(seg)
            scale = float(np.linalg.norm(np.asarray(ps.p1) - np.asarray(ps.p0)))
        except Exception:
            continue
        dl = float(np.linalg.norm(la - lp)) * scale * 1000.0
        local[seg].append(dl)
        worst.append((dl, d, m["code"], a["struct"]))

    def med(xs):
        return float(np.median(xs)) if len(xs) else float("nan")

    print("segment    n      tlo    swiat   /tlo    ramka")
    for seg in sorted(world, key=lambda s: -med(local[s])):
        b = med(bg.get(seg, []))
        ratio = med(world[seg]) / b if b == b else float("nan")
        print(f"{seg:9s} {len(world[seg]):3d}  {b:6.1f}  {med(world[seg]):6.1f}"
              f"  {ratio:5.2f}  {med(local[seg]):7.1f}")
    allx = [v for xs in local.values() for v in xs]
    print()
    print(f"w ramce lokalnej razem: mediana {med(allx):.1f} mm, "
          f"p90 {np.percentile(allx, 90):.1f} mm  (n={len(allx)})")
    print()
    print("najgorsze w ramce lokalnej (ramka / swiat):")
    for dl, d, code, name in sorted(worst, reverse=True)[:12]:
        print(f"  {dl:6.1f} / {d:6.1f} mm  {code:8s} {name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
