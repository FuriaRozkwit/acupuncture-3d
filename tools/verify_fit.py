"""Czy kosci (Z-Anatomy) siedza w skorze (BodyParts3D) - audyt dopasowania.

Oba zbiory dopasowano na tulowiu, wiec na konczynach rozjezdzaja sie. Sprawdzamy
to na kosciach PODSKORNYCH: nad mostkiem, wyrostkami kolczystymi, obojczykiem czy
grzebieniem piszczeli jest kilka milimetrow tkanki. Kilkanascie oznacza, ze siatki
sie rozjechaly i pomiary liczone od tej kosci sa niepewne.

    python tools/verify_fit.py
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
warnings.simplefilter("ignore")

import pyvista as pv  # noqa: E402

from acu3d.model import Atlas  # noqa: E402

REF = ROOT / "assets" / "mesh" / "ref"


def bone_points(an, name: str):
    """Kosc z assets/ready/struct, a gdy jej tam nie ma - z surowej siatki ref."""
    try:
        return np.asarray(an.struct(name).mesh.points)
    except KeyError:
        pass
    import re
    slug = re.sub(r"[^a-z0-9]+", "-", name.strip().lower().replace("(", "").replace(")", "")).strip("-")
    f = REF / f"{slug}.ply"
    if not f.exists():
        return None
    P = np.asarray(pv.read(f).points).copy()
    P[:, 0] *= -1.0                      # uklad Z-Anatomy -> uklad acu3d
    P[:, 1] *= -1.0
    return P

LIMIT = 12.0     # mm - powyzej tego uznajemy, ze kosc nie siedzi w skorze

# struktury podskorne: nazwa -> (struktura, kierunek 'na zewnatrz', ile skory
# spodziewamy sie u zywego czlowieka)
PROBES = [
    ("mostek", "Manubrium of sternum", [0.0, 1.0, 0.0], "2-6 mm"),
    ("wyrostek mieczykowaty", "Xiphoid process", [0.0, 1.0, 0.0], "3-8 mm"),
    ("obojczyk", "Clavicle.r", [0.0, 1.0, 0.3], "3-8 mm"),
    ("lopatka - grzebien", "Scapula.r", [0.0, -1.0, 0.0], "5-15 mm"),
    ("zuchwa", "Mandible", [0.0, 0.3, -1.0], "3-8 mm"),
    ("kosc czolowa", "Frontal bone", [0.0, 1.0, 0.2], "4-8 mm"),
    ("rzepka", "Patella.r", [0.0, 1.0, 0.0], "3-8 mm"),
    ("kosc krzyzowa", "Sacrum", [0.0, -1.0, 0.0], "5-15 mm"),
    ("kosc biodrowa - talerz", "Hip bone.r", [0.6, 0.0, 1.0], "8-20 mm"),
    ("kretarz wiekszy", "Femur.r", [1.0, 0.0, 0.0], "8-20 mm"),
]


def skin_gap(atlas, pts: np.ndarray) -> np.ndarray:
    _, foot = atlas.anatomy.skin_cast.find_closest_cell(np.atleast_2d(pts),
                                                        return_closest_point=True)
    return np.linalg.norm(np.atleast_2d(foot) - np.atleast_2d(pts), axis=1) * 1000


# Sondy miednicy i uda skierowane w ROZNE strony. Przesuniecie kosci topi jedna
# krawedz i wynurza przeciwlegla; jednostajny nadmiar na wszystkich oznacza, ze
# kosciec jest po prostu wezszy od skory - i zadne przesuniecie tego nie naprawi.
PELVIS = [
    ("krzyzowa z tylu",    "Sacrum",     [0.0, -1.0,  0.0], 12.0),
    ("spojenie lonowe",    "Hip bone.r", [0.0,  1.0, -0.3], 12.0),
    ("talerz biodrowy",    "Hip bone.r", [0.6,  0.0,  1.0], 12.0),
    ("kolec biodr. przedni", "Hip bone.r", [0.4, 1.0,  0.3], 10.0),
    ("guz kulszowy",       "Hip bone.r", [0.2, -0.5, -1.0], 25.0),
    ("kretarz wiekszy",    "Femur.r",    [1.0,  0.0,  0.0], 15.0),
    ("rzepka",             "Patella.r",  [0.0,  1.0,  0.0],  5.0),
    ("nadklykiec przysr.", "Femur.r",    [-1.0, 0.0,  0.0],  6.0),
]


def shift_or_size(atlas) -> None:
    """Czy miednice da sie wpasowac przesunieciem (jak podudzie), czy nie."""
    an = atlas.anatomy
    rows = []
    for name, struct, d, want in PELVIS:
        P = bone_points(an, struct)
        if P is None:
            continue
        d = np.asarray(d, float); d /= np.linalg.norm(d)
        sel = P[P @ d > np.percentile(P @ d, 98)]
        rows.append((name, d, want, float(np.median(skin_gap(atlas, sel)))))
    if len(rows) < 4:
        return
    D = np.array([r[1] for r in rows])
    r = np.array([r[3] - r[2] for r in rows])          # nadmiar tkanki (mm)
    resid = r - (-D @ np.linalg.lstsq(-D, r, rcond=None)[0])
    print("\nmiednica i udo - przesuniecie czy rozmiar?")
    for (name, _, want, got), e in zip(rows, r):
        print(f"   {name:24s} {got:6.1f} mm   (oczekiwane {want:5.1f})   nadmiar {e:+6.1f}")
    print(f"   nadmiar na wszystkich sondach: srednio {np.mean(r):+.1f} mm "
          f"(odchylenie {np.std(r):.1f})")
    print(f"   po najlepszym przesunieciu zostaje {np.sqrt(np.mean(resid**2)):.1f} mm "
          f"z {np.sqrt(np.mean(r**2)):.1f} mm  ->  "
          + ("przesuniecie pomaga" if np.sqrt(np.mean(resid**2)) < 0.5 * np.sqrt(np.mean(r**2))
             else "przesuniecie NIE pomaga: kosciec jest wezszy od skory, nie przesuniety"))


def main() -> int:
    atlas = Atlas()
    an = atlas.anatomy
    print(f"{'struktura':26s} {'skory nad kescia':>18s}   oczekiwane   wynik")
    bad = 0
    for name, struct, direction, expect in PROBES:
        P = bone_points(an, struct)
        if P is None:
            print(f"{name:26s} {'brak siatki':>18s}")
            continue
        d = np.asarray(direction, dtype=float)
        d /= np.linalg.norm(d)
        # bierzemy 5% wierzcholkow najdalej wysunietych w zadanym kierunku:
        # to one leza tuz pod skora
        proj = P @ d
        sel = P[proj > np.percentile(proj, 95)]
        gap = float(np.median(skin_gap(atlas, sel)))
        ok = gap <= LIMIT
        bad += 0 if ok else 1
        print(f"{name:26s} {gap:15.1f} mm   {expect:>9s}   "
              f"{'ok' if ok else 'ZLE WPASOWANA'}")

    shift_or_size(atlas)

    # wyrostki kolczyste: zmierzone w rig_fitted, sa podskorne na calej dlugosci
    sp = np.array([v["spinous"] for k, v in an.fit["vertebra"].items()
                   if k.startswith(("T", "L"))])
    gap = skin_gap(atlas, sp)
    print(f"{'wyrostki kolczyste T/L':26s} {np.median(gap):15.1f} mm   {'5-15 mm':>9s}   "
          f"{'ok' if np.median(gap) <= 20 else 'ZLE WPASOWANE'}")
    print(f"\nstruktur zle wpasowanych: {bad}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
