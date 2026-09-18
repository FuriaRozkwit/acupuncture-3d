"""Poziom szumu testow weryfikatora - zmierzony na TYCH siatkach.

Prog ponizej szumu pomiarowego nie wykrywa bledow, tylko rozdzielczosc siatki.
Sprawdzamy, ile "zalamania" i ile rozjazdu Frecheta daja same poprawne krzywe:
geodezyjna jest z definicji gladka, wiec wszystko, co na niej widac, to szum.

Stad biora sie progi KINK_WARN_REL/KINK_FAIL_REL w verify_meridians.py oraz
liczby w polu "noise_floor" raportu - nie sa zgadniete.

    python tools/measure_noise.py
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))
warnings.simplefilter("ignore")

import numpy as np  # noqa: E402

import verify_meridians as V  # noqa: E402
from acu3d.model import Atlas  # noqa: E402

a = Atlas()
P = V.Paths(a, frechet=True)
MER = json.load(open(ROOT / "data" / "meridians.json", encoding="utf-8"))["paths"]

# --- 1. ile "zalamania" ma sama geodezyjna (krzywa z definicji gladka) --------
kk, rel, pairs = [], [], []
seen = set()
for ch, groups in MER.items():
    for g in groups:
        for x, y in zip(g[:-1], g[1:]):
            if a.by_key.get(x) or a.by_key.get(x + "·R"):
                if a.by_key.get(y) or a.by_key.get(y + "·R"):
                    pairs.append((x, y))
np.random.seed(0)
sel = [pairs[i] for i in np.random.choice(len(pairs), 60, replace=False)]
for x, y in sel:
    g = P.geodesic(P.pos(x), P.pos(y))
    if g is None or len(g) < 6:
        continue
    q = V.resample(g, 24)
    _, kink, krel = V.kinks(q)
    if kink.size:
        kk.append(float(kink.max()))
        rel.append(float(krel.max()))
kk, rel = np.array(kk), np.array(rel)
print(f"zalamanie samej geodezyjnej (mm):  mediana {np.median(kk):.2f}  "
      f"90. centyl {np.percentile(kk,90):.2f}  max {kk.max():.2f}   n={len(kk)}")
print(f"  prog w narzedziu: WARN {V.KINK_WARN_MM}  FAIL {V.KINK_FAIL_MM}")
print(f"zalamanie wzgledne (x odstep probek): mediana {np.median(rel):.3f}  "
      f"p90 {np.percentile(rel, 90):.3f}  p99 {np.percentile(rel, 99):.3f}  max {rel.max():.3f}")
print(f"  prog w narzedziu: WARN {V.KINK_WARN_REL}  FAIL {V.KINK_FAIL_REL}")

# --- 2. powtarzalnosc geodezyjnej: dwie liczone inaczej, ta sama para ---------
d = []
for x, y in sel[:30]:
    g1 = P.geodesic(P.pos(x), P.pos(y), n=48, rounds=40)
    g2 = P.geodesic(P.pos(x), P.pos(y), n=72, rounds=90)
    if g1 is None or g2 is None:
        continue
    fr = V.discrete_frechet(V.resample(g1, 48), V.resample(g2, 48))
    val, unit, _ = P._measure(fr, x)
    d.append((val, unit))
cun = np.array([v for v, u in d if u == "cun"])
mm = np.array([v for v, u in d if u == "mm"])
if cun.size:
    print(f"\npowtarzalnosc geodezyjnej (cun): mediana {np.median(cun):.2f}  "
          f"90. centyl {np.percentile(cun,90):.2f}  max {cun.max():.2f}   n={cun.size}")
    print(f"  prog w narzedziu: WARN {V.WARN_CUN}  FAIL {V.TOL_CUN}")
if mm.size:
    print(f"powtarzalnosc geodezyjnej (mm):  mediana {np.median(mm):.2f}  "
          f"90. centyl {np.percentile(mm,90):.2f}  max {mm.max():.2f}   n={mm.size}")
    print(f"  prog w narzedziu: WARN {V.WARN_MM}  FAIL {V.TOL_MM}")
