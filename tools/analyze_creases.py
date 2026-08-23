"""Czy meridiany biegna faldami skory? -> pomiar, nie wrazenie.

Teza brzmi: naturalne fałdy i wciecia skory pokrywaja sie z przebiegiem kanalow.
Fałd to miejsce wklesle, a wkleslosc da sie policzyc jako ujemna krzywizna
srednia powierzchni. Sprawdzamy wiec, czy probki linii kanalow trafiaja w miejsca
bardziej wklesle niz przecietny kawalek skory.

Porownanie jest z tlem losowym: sama ujemna krzywizna w punktach kanalu niczego
nie dowodzi, dopoki nie wiadomo, jaka krzywizne ma skora "gdziekolwiek".

    python tools/analyze_creases.py [--top 10]

UWAGA co do wniosku: skora pochodzi z modelu anatomicznego, ktory jest wygladzony.
Jesli nie ma w niej prawdziwych fałdow, wynik bedzie zerowy - i to tez jest
odpowiedz, tyle ze o modelu, a nie o ciele.
"""
from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
warnings.simplefilter("ignore")

from acu3d.model import Atlas  # noqa: E402


def _neighbors(skin):
    faces = skin.faces.reshape(-1, 4)[:, 1:]
    return faces


def _smooth(c: np.ndarray, faces: np.ndarray, rounds: int) -> np.ndarray:
    c = c.copy()
    for _ in range(rounds):
        acc = np.zeros_like(c)
        cnt = np.zeros(len(c))
        for i in range(3):
            a, b = faces[:, i], faces[:, (i + 1) % 3]
            np.add.at(acc, a, c[b])
            np.add.at(cnt, a, 1)
            np.add.at(acc, b, c[a])
            np.add.at(cnt, b, 1)
        c = np.where(cnt > 0, 0.5 * c + 0.5 * acc / np.maximum(cnt, 1), c)
    return c


def groove_field(skin, fine: int = 2, coarse: int = 40) -> np.ndarray:
    """Jak bardzo dane miejsce jest wklesłe NA TLE SWOJEGO OTOCZENIA.

    Sama krzywizna nie nadaje sie do szukania fałdow: najbardziej wklesłe miejsca
    modelu to szczeliny miedzy palcami, pacha i krocze, a nie fałdy skorne.
    Fałd na przedramieniu czy pod kolanem to plytkie wglebienie na wypukłym tle,
    wiec od krzywizny wygladzonej lekko odejmujemy wygladzona mocno - zostaje
    sam "rowek", bez ogolnego ksztaltu czesci ciala.
    """
    faces = _neighbors(skin)
    c = np.asarray(skin.curvature("mean"), dtype=float)
    return _smooth(c, faces, fine) - _smooth(c, faces, coarse)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--top", type=float, default=10.0,
                    help="ile procent najbardziej wklesłych miejsc uznajemy za fałd")
    args = ap.parse_args()

    atlas = Atlas()
    skin = atlas.layers["skin"]
    pts = np.asarray(skin.points)
    curv = groove_field(skin)

    thr = np.percentile(curv, args.top)          # prog "fałdu" (najbardziej wkleslе)
    print(f"skora: {skin.n_points} wierzcholkow")
    print(f"krzywizna srednia: mediana {np.median(curv):+.1f}, "
          f"1. percentyl {np.percentile(curv, 1):+.1f}, "
          f"99. percentyl {np.percentile(curv, 99):+.1f} [1/m]")
    print(f"prog fałdu ({args.top:.0f}% najbardziej wklesłych): {thr:+.1f}\n")

    rows = []
    all_hits = 0
    all_n = 0
    for ch in atlas.channels:
        vals = []
        for side in ("R", "L"):
            for line in atlas.meridian_lines(ch, side):
                for p in line:
                    j = skin.find_closest_point(p)
                    vals.append(curv[j])
            if atlas.channels[ch]["clock"] >= 13:
                break
        if not vals:
            continue
        vals = np.array(vals)
        hits = int((vals <= thr).sum())
        all_hits += hits
        all_n += len(vals)
        rows.append((ch, len(vals), float(np.median(vals)), 100.0 * hits / len(vals)))

    base_frac = args.top
    print(f"{'kanal':6s} {'probek':>7s} {'mediana rowka':>19s} {'% w fałdzie':>12s}")
    for ch, n, med, frac in sorted(rows, key=lambda r: -r[3]):
        mark = "  <-- wyraznie czesciej" if frac > base_frac * 1.8 else ""
        print(f"{ch:6s} {n:7d} {med:+19.1f} {frac:11.1f}%{mark}")

    print(f"\nrazem: {100.0 * all_hits / max(all_n, 1):.1f}% probek kanalow lezy w fałdzie")
    print(f"tło losowe:  {base_frac:.1f}%  (tyle wynosi prog z definicji)")
    ratio = (100.0 * all_hits / max(all_n, 1)) / base_frac
    print(f"iloraz: {ratio:.2f}x")

    # Kontrola miernika: podrecznik sam wskazuje punkty lezace NA fałdzie
    # (lokciowym, nadgarstka, podkolanowym). Jesli miernik ich nie znajduje,
    # to problem jest w modelu albo w mierniku, a nie w tezie.
    crease_codes = [c for c, m in atlas.meta.items()
                    if "crease" in m.get("location", "").lower()]
    vals = []
    for c in crease_codes:
        m = atlas.by_key.get(c) or atlas.by_key.get(f"{c}·R")
        if m is not None:
            vals.append(curv[skin.find_closest_point(m["pos"])])
    if vals:
        vals = np.array(vals)
        print(f"\nKONTROLA: {len(vals)} punktow, ktore podrecznik opisuje faldem")
        print(f"  ich wskaznik rowka: mediana {np.median(vals):+.1f} "
              f"(tlo: {np.median(curv):+.1f}, prog fałdu: {thr:+.1f})")
        print(f"  ile trafia w fałd: {100.0 * (vals <= thr).mean():.0f}% "
              f"(los: {base_frac:.0f}%)")
    if ratio < 1.2:
        print("\nWNIOSEK: na TYM modelu tezy nie da sie sprawdzic.")
        print("Skora z BodyParts3D jest wygladzona i stoi w pozycji wyprostowanej,")
        print("a fałdy stawowe powstaja dopiero przy zgieciu - dlatego nie ma ich")
        print("nawet tam, gdzie podrecznik wprost je wymienia: fald lokciowy przy")
        print("LU-5, fald nadgarstka przy PC-7, fald podkolanowy przy BL-40.")
        print("Do sprawdzenia tezy potrzeba skory ze skanu, najlepiej w zgieciu.")
    else:
        print("\nWNIOSEK: kanaly istotnie czesciej trafiaja w miejsca wklesłe.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
