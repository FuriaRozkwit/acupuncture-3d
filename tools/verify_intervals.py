"""Czy odleglosci MIEDZY kolejnymi punktami kanalu zgadzaja sie z kanonem.

Dotychczasowe testy sprawdzaly kazdy punkt osobno, wobec jego wlasnego opisu -
a to obejmowalo tylko te opisy, w ktorych jest liczba do wylowienia (123 z 361
punktow). Kanal jest jednak galezia o znanych interpunkcjach: LU-7 lezy 1.5 cun
nad zgieciem nadgarstka NIEZALEZNIE od tego, jak opisano LU-6 i LU-8, a odstep
LU-6 -> LU-7 jest przez to ustalony.

Odstepy kanoniczne to fakty tradycji, nie tekst z podrecznika - dlatego
`data/cun_intervals.json` moze byc w repo i mozna go pokazac komukolwiek,
w odrozieniu od `data/point_meta.json`.

Mierzymy skladowa WZDLUZ segmentu (tak liczy sie cun na konczynie) albo
roznice na odcinkowej linijce tulowia. Odstepy "anatomiczna" i "vertebra"
nie maja kanonicznej liczby i sa tylko liczone w pokryciu.

    python tools/verify_intervals.py [--tol 0.5]
"""
from __future__ import annotations

import argparse
import json
import sys
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
warnings.simplefilter("ignore")

from acu3d.model import Atlas  # noqa: E402

TRUNK = ("trunk", "neck", "head")


class Intervals:
    def __init__(self, atlas: Atlas):
        self.a = atlas
        self.an = atlas.anatomy
        self.anchor = {p["code"]: p["anchor"] for p in atlas.points}

    def pos(self, code: str):
        m = self.a.by_key.get(code) or self.a.by_key.get(f"{code}\u00b7R")
        return None if m is None else np.asarray(m["pos"], float)

    def measure(self, a: str, b: str):
        """Odleglosc w cun i sposob pomiaru, albo (None, powod)."""
        p, q = self.pos(a), self.pos(b)
        if p is None or q is None:
            return None, "brak punktu w modelu"
        sa = self.anchor.get(a, {}).get("seg")
        sb = self.anchor.get(b, {}).get("seg")
        if sa and sa == sb and sa not in TRUNK:
            seg = self.an.seg(sa)
            # SP, KI i LR na podudziu maja WLASNA linijke (przysrodkowa, 13 cun
            # zamiast 16 na tej samej dlugosci), wiec ich cun jest dluzszy.
            # Mierzenie ich standardowym cun segmentu dawalo +10..17% na kazdym
            # odstepie i wygladalo na blad atlasu, a bylo bledem miary.
            ra = self.anchor.get(a, {}).get("ruler")
            rb = self.anchor.get(b, {}).get("ruler")
            cun_m = seg.cun_m
            how = f"wzdluz segmentu {sa}"
            if ra and ra == rb:
                r = self.an.fit["rulers"][ra]
                length = abs(float(np.dot(np.asarray(r["p1"]) - np.asarray(r["p0"]),
                                          seg.e_ax)))
                cun_m = length / float(r["cun_length"])
                how = f"linijka {ra}"
            d = abs(float(np.dot(p - q, seg.e_ax))) / cun_m
            return d, how
        if (sa in TRUNK or sa is None) and (sb in TRUNK or sb is None):
            f = getattr(self.an, "cun_of_z", None)
            if f is not None:
                return abs(f(float(p[2])) - f(float(q[2]))), "linijka tulowia"
        return None, f"punkty w roznych segmentach ({sa} / {sb})"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tol", type=float, default=0.5, help="tolerancja w cun")
    args = ap.parse_args()

    doc = json.loads((ROOT / "data" / "cun_intervals.json").read_text(encoding="utf-8"))
    atlas = Atlas()
    iv = Intervals(atlas)

    rows, skipped = [], []
    pairs = sum(len(v) for v in doc["channels"].values())
    canon = 0
    for ch, entries in doc["channels"].items():
        for e in entries:
            if e["cun"] is None:
                continue
            canon += 1
            got, how = iv.measure(e["from"], e["to"])
            if got is None:
                skipped.append((e["from"], e["to"], how))
                continue
            rows.append((ch, e["from"], e["to"], float(e["cun"]), got,
                         got - float(e["cun"]), how))

    print(f"tabela: {pairs} par w 14 kanalach, {canon} z kanoniczna liczba "
          f"({canon * 100 // max(pairs, 1)}%)")
    print(f"zmierzonych: {len(rows)}, pominietych: {len(skipped)}")
    if not rows:
        return 1
    err = np.array([r[5] for r in rows])
    ok = np.abs(err) <= args.tol
    print(f"\nzgodnych w tolerancji {args.tol} cun: {ok.sum()}/{len(rows)} "
          f"({ok.sum() * 100 // len(rows)}%)")
    print(f"blad: mediana {np.median(np.abs(err)):.2f}, sredni ze znakiem {err.mean():+.2f}, "
          f"max {np.abs(err).max():.2f} cun")

    print("\nnajwieksze rozbieznosci:")
    for ch, a, b, want, got, e, how in sorted(rows, key=lambda r: -abs(r[5]))[:14]:
        print(f"   {a:7s} -> {b:7s} kanon {want:5.2f}  model {got:5.2f}  "
              f"({e:+5.2f} cun)  {how}")

    by = {}
    for ch, a, b, want, got, e, how in rows:
        by.setdefault(ch, []).append(e)
    print("\nwg kanalu:")
    for ch, es in sorted(by.items(), key=lambda x: -len(x[1])):
        es = np.array(es)
        n_ok = int((np.abs(es) <= args.tol).sum())
        print(f"   {ch:4s} {n_ok:2d}/{len(es):2d}  sredni blad {es.mean():+.2f} cun")

    if skipped:
        print("\npominiete (nie da sie zmierzyc jedna miara):")
        for a, b, why in skipped[:8]:
            print(f"   {a:7s} -> {b:7s}  {why}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
