"""Linie pracy miesni: os, wzdluz ktorej miesien ciagnie.

Miesien nie jest bryla, ktora "sie kurczy" - ma kierunek dzialania. Linia pracy
biegnie srodkiem brzuscia od przyczepu do przyczepu, a groty na koncach pokazuja,
w ktora strone idzie ciag: OBA konce ciagna do brzuscia, bo to on sie skraca.

Linie liczymy z samej geometrii miesnia, a nie z zalozen: wierzcholki rzutujemy
na os glowna, dzielimy na przedzialy i bierzemy srodek kazdego przedzialu.
Dzieki temu linia idzie srodkiem miesnia takze wtedy, gdy miesien jest wygiety
(dwuglowy uda, mostkowo-obojczykowo-sutkowy), zamiast przecinac go na skros.
"""
from __future__ import annotations

import numpy as np
import pyvista as pv

BINS = 9            # ile odcinkow linii srodkowej
END_FRAC = 0.12     # jaka czesc dlugosci uznajemy za przyczep


def centerline(points: np.ndarray, axis: np.ndarray, bins: int = BINS) -> np.ndarray:
    """Linia srodkowa miesnia: srodki ciezkosci kolejnych plastrow wzdluz osi."""
    if len(points) < bins * 3:
        return np.array([points.mean(axis=0)])
    u = points @ axis
    edges = np.linspace(u.min(), u.max(), bins + 1)
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = points[(u >= lo) & (u <= hi)]
        if len(sel) >= 3:
            out.append(sel.mean(axis=0))
    if len(out) < 2:
        return np.array([points.mean(axis=0)])
    line = np.array(out)
    # jedno wygladzenie, zeby plastry o nierownej gestosci nie robily zygzaka
    if len(line) >= 3:
        line[1:-1] = 0.5 * line[1:-1] + 0.25 * (line[:-2] + line[2:])
    return line


def attachments(points: np.ndarray, axis: np.ndarray, frac: float = END_FRAC):
    """Srodki obu przyczepow: skrajne plastry miesnia wzdluz jego osi."""
    u = points @ axis
    span = u.max() - u.min()
    lo = points[u <= u.min() + span * frac]
    hi = points[u >= u.max() - span * frac]
    a = lo.mean(axis=0) if len(lo) else points[np.argmin(u)]
    b = hi.mean(axis=0) if len(hi) else points[np.argmax(u)]
    return a, b


def build(points: np.ndarray, axis: np.ndarray, arrow_len: float = 0.022,
          bins: int = BINS) -> pv.PolyData:
    """Siatka jednej linii pracy: rurka wzdluz brzuscia + dwa groty ku srodkowi."""
    line = centerline(points, axis, bins)
    parts = []
    if len(line) >= 2:
        parts.append(pv.lines_from_points(line).tube(radius=0.0022, n_sides=8))
    a, b = attachments(points, axis)
    mid = line.mean(axis=0)
    for end in (a, b):
        d = mid - end
        n = np.linalg.norm(d)
        if n < 1e-6:
            continue
        d = d / n
        tip = min(arrow_len, n * 0.45)
        parts.append(pv.Cone(center=end + d * tip * 0.5, direction=d,
                             height=tip, radius=tip * 0.38, resolution=10))
    if not parts:
        return pv.PolyData()
    return parts[0].merge(parts[1:]) if len(parts) > 1 else parts[0]


def build_all(muscle_set, muscles=None, min_volume_l: float = 0.008) -> pv.PolyData:
    """Linie pracy dla zbioru miesni. Drobne pomijamy, bo tworza tylko szum."""
    out = []
    for m in (muscles if muscles is not None else muscle_set.muscles):
        if m.volume_l < min_volume_l:
            continue
        pts = muscle_set.rest[m.slice]
        if len(pts) < 40:
            continue
        mesh = build(pts, m.axis)
        if mesh.n_points:
            out.append(mesh)
    if not out:
        return pv.PolyData()
    return out[0].merge(out[1:]) if len(out) > 1 else out[0]
