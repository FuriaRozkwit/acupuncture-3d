"""Gladkie przebiegi kanalow: krzywe Beziera przechodzace przez punkty.

Kanal narysowany odcinkami prostymi wyglada jak lamana, a biegnie plynnie.
Uzywamy splajnu Catmulla-Roma zapisanego jako lancuch szesciennych krzywych
Beziera - w odroznieniu od zwyklego Beziera taki splajn PRZECHODZI przez wezly,
a wezly to punkty akupunkturowe, wiec krzywa nie moze ich omijac.

Parametryzacja dosrodkowa (alpha = 0.5) zamiast rownomiernej ma tu znaczenie
praktyczne: odstepy miedzy punktami wahaja sie od 1 cun do kilkunastu centymetrow,
a przy rownomiernej parametryzacji tak nierowne odcinki daja petle i ostrza.
"""
from __future__ import annotations

import numpy as np

ALPHA = 0.5          # 0 = rownomierna, 0.5 = dosrodkowa, 1 = cieciwowa
TENSION = 1.0        # 1 = pelne wybrzuszenie, 0 = odcinki proste


def _knots(points: np.ndarray, alpha: float) -> np.ndarray:
    d = np.linalg.norm(np.diff(points, axis=0), axis=1)
    d = np.maximum(d, 1e-9) ** alpha
    return np.concatenate([[0.0], np.cumsum(d)])


def bezier_segments(points: np.ndarray, alpha: float = ALPHA,
                    tension: float = TENSION) -> list:
    """Lancuch segmentow Beziera (P0, C1, C2, P3) przechodzacy przez `points`."""
    P = np.asarray(points, dtype=float)
    if len(P) < 2:
        return []
    if len(P) == 2:
        return [(P[0], P[0] + (P[1] - P[0]) / 3, P[1] - (P[1] - P[0]) / 3, P[1])]
    # konce przedluzamy odbiciem, zeby pierwszy i ostatni odcinek tez mial styczna
    ext = np.vstack([2 * P[0] - P[1], P, 2 * P[-1] - P[-2]])
    t = _knots(ext, alpha)
    turn = _turn_scale(ext)
    out = []
    for i in range(1, len(ext) - 2):
        p0, p1, p2, p3 = ext[i - 1], ext[i], ext[i + 1], ext[i + 2]
        t0, t1, t2, t3 = t[i - 1], t[i], t[i + 1], t[i + 2]
        dt = t2 - t1
        m1 = (p2 - p0) / max(t2 - t0, 1e-9)      # styczna w p1
        m2 = (p3 - p1) / max(t3 - t1, 1e-9)      # styczna w p2
        c1 = p1 + m1 * dt / 3.0 * tension * turn[i]
        c2 = p2 - m2 * dt / 3.0 * tension * turn[i + 1]
        out.append((p1, c1, c2, p2))
    return out


def _turn_scale(points: np.ndarray) -> np.ndarray:
    """O ile skrocic styczna w wezle, w ktorym kanal ostro zawraca.

    Przy zakosie (kanal pecherzyka zolciowego na skroni, biodro, zebra) pelna
    styczna Catmulla-Roma wyrzuca krzywa poza zakret - splajn wybrzusza sie
    kilka do kilkunastu milimetrow obok najkrotszej drogi po skorze. Im ostrzejszy
    zakret, tym krotsza styczna; przy linii prostej nic sie nie zmienia.
    """
    v = np.diff(points, axis=0)
    n = np.linalg.norm(v, axis=1, keepdims=True)
    u = v / np.maximum(n, 1e-12)
    scale = np.ones(len(points))
    cos = (u[:-1] * u[1:]).sum(axis=1)
    scale[1:-1] = np.clip(np.sqrt(np.clip((1.0 + cos) / 2.0, 0.0, 1.0)), 0.25, 1.0)
    return scale


def sample(points: np.ndarray, per_segment: int = 14, alpha: float = ALPHA,
           tension: float = TENSION) -> np.ndarray:
    """Gesto probkowana krzywa. Wezly wejsciowe leza dokladnie na wyniku."""
    segs = bezier_segments(points, alpha, tension)
    if not segs:
        return np.asarray(points, dtype=float)
    u = np.linspace(0.0, 1.0, per_segment, endpoint=False)[:, None]
    out = []
    for p0, c1, c2, p3 in segs:
        b = ((1 - u) ** 3 * p0 + 3 * (1 - u) ** 2 * u * c1
             + 3 * (1 - u) * u ** 2 * c2 + u ** 3 * p3)
        out.append(b)
    out.append(np.asarray(segs[-1][3])[None, :])
    return np.vstack(out)


def sample_varied(points: np.ndarray, counts, alpha: float = ALPHA,
                  tension: float = TENSION):
    """Krzywa, w ktorej kazdy odcinek ma wlasna liczbe probek (dluzszy - gestsza).

    Splajn liczony raz dla calej lamanej (styczne zaleza od sasiadow), probkowany
    osobno odcinek po odcinku. Zwraca (probki, indeksy wezlow).
    """
    segs = bezier_segments(points, alpha, tension)
    if not segs:
        P = np.asarray(points, dtype=float)
        return P, np.arange(len(P))
    out, nodes, start = [], [0], 0
    for (p0, c1, c2, p3), n in zip(segs, counts):
        u = np.linspace(0.0, 1.0, int(n), endpoint=False)[:, None]
        out.append((1 - u) ** 3 * p0 + 3 * (1 - u) ** 2 * u * c1
                   + 3 * (1 - u) * u ** 2 * c2 + u ** 3 * p3)
        start += int(n)
        nodes.append(start)
    out.append(np.asarray(segs[-1][3])[None, :])
    return np.vstack(out), np.array(nodes)


def sample_pair(points: np.ndarray, refs: np.ndarray, per_segment: int = 14,
                alpha: float = ALPHA, tension: float = TENSION):
    """Krzywa punktow i rownolegla krzywa punktow odniesienia (wnetrze ciala).

    Obie probkowane tak samo, wiec i-ta probka jednej odpowiada i-tej drugiej -
    dzieki temu kazda probke mozna potem wypchnac na skore wzdluz jej wlasnej
    normalnej, a nie wzdluz jednego, wspolnego kierunku.
    """
    return (sample(points, per_segment, alpha, tension),
            sample(refs, per_segment, alpha, tension))


def relax(curve: np.ndarray, pinned: np.ndarray, rounds: int = 2,
          weight: float = 0.5) -> np.ndarray:
    """Wygladza drgania powstale przy rzutowaniu na skore.

    Probki miedzywezlowe usredniamy z sasiadami, ale wezly (czyli punkty
    akupunkturowe) zostaja przypiete - krzywa ma je laczyc, a nie wygladzac.
    """
    c = np.array(curve, dtype=float, copy=True)
    if len(c) < 3:
        return c
    for _ in range(rounds):
        inner = 0.5 * (c[:-2] + c[2:])
        free = ~pinned[1:-1]
        c[1:-1][free] = (1 - weight) * c[1:-1][free] + weight * inner[free]
    return c
