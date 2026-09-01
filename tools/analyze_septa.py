"""Czy kanaly biegna bruzdami miedzy grupami miesniowymi?

Hipoteza (P-01): naturalne bruzdy na konczynach powstaja tam, gdzie przegroda
miedzymiesniowa dochodzi do powiezi glebokiej i skory. Jesli przebieg kanalow
odpowiada fizjonomii, punkty i linie kanalow powinny lezec blizej sladu przegrody
niz przecietne miejsce w tej samej okolicy.

Metoda
------
1. Slad przegrody na skorze: bierzemy najplytsze wierzcholki przegrody (dolny
   kwartyl odleglosci do skory) i rzutujemy je na skore.
2. Dla kazdego badanego miejsca liczymy odleglosc do najblizszego sladu.
3. TLO: dla tego samego miejsca losujemy punkty skory w promieniu `--halo`
   i liczymy te sama odleglosc. Bez tego porownania liczba w milimetrach
   nic nie znaczy - na cienkim przedramieniu wszystko jest blisko wszystkiego.
4. Wynik na miejsce to percentyl: jaki ulamek tla lezy blizej niz badane miejsce.
   Gdyby hipoteza byla falszywa, mediana percentyli wyszlaby 0.50.
5. Istotnosc: test permutacyjny na medianie percentyli.

    python tools/analyze_septa.py [--halo 0.09] [--null 400] [--json out.json]
"""
from __future__ import annotations

import argparse
import json
import sys
import warnings
from pathlib import Path

import numpy as np

warnings.simplefilter("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from acu3d.fascia import FasciaSet          # noqa: E402
from acu3d.model import Atlas               # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
# przegrody istnieja tylko dla tych odcinkow - reszty ciala ta analiza nie dotyczy
SEGMENTS = ("upper_arm", "thigh", "shank")


def trace_on_skin(sheet_pts, skin_pts, tree, quantile: float):
    """Rzut najplytszej czesci platu na skore."""
    d, idx = tree.query(sheet_pts)
    if len(d) == 0:
        return np.empty((0, 3))
    cut = np.quantile(d, quantile)
    keep = d <= cut
    return skin_pts[idx[keep]]


def percentiles(targets, traces_tree, skin_pts, skin_tree, halo, n_null, rng):
    """Dla kazdego miejsca: odleglosc do sladu i jej pozycja w lokalnym tle."""
    out = []
    for p in targets:
        d = float(traces_tree.query(p)[0])
        near = skin_tree.query_ball_point(p, halo)
        if len(near) < 30:
            continue
        pick = rng.choice(near, size=min(n_null, len(near)), replace=False)
        dn = traces_tree.query(skin_pts[pick])[0]
        out.append((d, float((dn < d).mean())))
    return out


def perm_test(pcts, rng, n=20000):
    """Czy mediana percentyli odbiega od 0.50? Losujemy z rozkladu jednostajnego."""
    obs = float(np.median(pcts))
    draws = np.median(rng.random((n, len(pcts))), axis=1)
    return obs, float((np.abs(draws - 0.5) >= abs(obs - 0.5)).mean())


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--halo", type=float, default=0.09, help="promien tla w metrach")
    ap.add_argument("--null", type=int, default=400, help="losowan tla na miejsce")
    ap.add_argument("--quantile", type=float, default=0.25, help="ktora czesc platu jest 'plytka'")
    ap.add_argument("--json", type=str, default="")
    args = ap.parse_args()

    from scipy.spatial import cKDTree

    rng = np.random.default_rng(20260901)
    atlas = Atlas()
    fascia = FasciaSet()
    skin_pts = np.asarray(atlas.layers["skin"].points)
    skin_tree = cKDTree(skin_pts)

    septa = fascia.of_class("septum")
    traces, per_septum = [], {}
    for s in septa:
        pts = np.asarray(fascia.sheet_mesh(s).points)
        tr = trace_on_skin(pts, skin_pts, skin_tree, args.quantile)
        per_septum[s.id] = {"name": s.name, "side": s.side, "trace_pts": len(tr),
                            "tree": cKDTree(tr) if len(tr) else None}
        if len(tr):
            traces.append(tr)
    trace_pts = np.vstack(traces)
    tree = cKDTree(trace_pts)
    print(f"przegrod: {len(septa)},  punktow sladu na skorze: {len(trace_pts)}")

    # --- punkty akupunkturowe na odcinkach z przegrodami ------------------
    seg_of = {p["code"]: p["anchor"].get("seg") for p in atlas.points}
    sel = [m for m in atlas.markers if seg_of.get(m["code"]) in SEGMENTS]
    res = percentiles([m["pos"] for m in sel], tree, skin_pts, skin_tree,
                      args.halo, args.null, rng)
    dists = np.array([r[0] for r in res]) * 1000.0
    pcts = np.array([r[1] for r in res])
    med, p = perm_test(pcts, rng)
    print(f"\nPUNKTY  (n={len(res)} z {len(sel)} kotwiczonych na ramieniu/udzie/podudziu)")
    print(f"  odleglosc do sladu przegrody: mediana {np.median(dists):.1f} mm, "
          f"kwartyle {np.percentile(dists,25):.1f}-{np.percentile(dists,75):.1f} mm")
    print(f"  percentyl wzgledem tla: mediana {med:.3f}  (0.50 = brak zwiazku)   p = {p:.4f}")
    print(f"  udzial punktow blizej niz polowa tla: {(pcts < 0.5).mean():.0%}")

    # --- linie kanalow na tych samych odcinkach ---------------------------
    seg_objs = [atlas.anatomy.seg(n) for n in SEGMENTS]

    def in_segment(q):
        for s in seg_objs:
            v = s.p1 - s.p0
            t = float(np.dot(q - s.p0, v) / max(np.dot(v, v), 1e-9))
            if 0.0 <= t <= 1.0:
                rad = np.linalg.norm(q - (s.p0 + t * v))
                if rad < 0.13:
                    return True
        return False

    line_pts = []
    for ch in atlas.channels:
        for side in ("R", "L"):
            for ln in atlas.meridian_lines(ch, side):
                for q in np.asarray(ln)[::3]:
                    if in_segment(q):
                        line_pts.append(q)
    line_pts = np.array(line_pts)
    lres = percentiles(line_pts, tree, skin_pts, skin_tree, args.halo, 120, rng)
    ldist = np.array([r[0] for r in lres]) * 1000.0
    lpct = np.array([r[1] for r in lres])
    lmed, lp = perm_test(lpct, rng)
    print(f"\nLINIE KANALOW  (n={len(lres)} probek na tych samych odcinkach)")
    print(f"  odleglosc do sladu przegrody: mediana {np.median(ldist):.1f} mm")
    print(f"  percentyl wzgledem tla: mediana {lmed:.3f}   p = {lp:.4f}")

    # --- rozbicie na kanaly ------------------------------------------------
    # Hipoteza w wersji mocnej ("wszystkie kanaly biegna bruzdami") jest za szeroka.
    # Sensowna wersja brzmi: NIEKTORE kanaly trzymaja sie konkretnej przegrody.
    # Kanal jest tu jednostka niezalezna - probki jednej linii sa ze soba zwiazane,
    # wiec liczenie p z 464 probek zawyzaloby istotnosc.
    print()
    print("WEDLUG KANALU  (punkty na ramieniu/udzie/podudziu)")
    by_ch = {}
    for m, (d, q) in zip(sel, res):
        by_ch.setdefault(m["channel"], []).append((d * 1000, q))
    rows = []
    for ch, vals in by_ch.items():
        if len(vals) < 3:
            continue
        rows.append((float(np.median([v[1] for v in vals])),
                     float(np.median([v[0] for v in vals])), ch, len(vals)))
    for q, d, ch, n in sorted(rows):
        bar = "#" * int(round((1.0 - q) * 24))
        print(f"  {ch:<3} n={n:<3} mediana {d:5.1f} mm   percentyl {q:.2f}  {bar}")
    ch_med = np.array([r[0] for r in rows])
    cm, cp = perm_test(ch_med, rng)
    print(f"  jednostka = kanal (n={len(rows)}): mediana percentyli {cm:.3f}   p = {cp:.4f}")

    # --- ktory punkt lezy na ktorej przegrodzie ---------------------------
    print("\nNAJBLIZSZE PARY  (punkt -> przegroda, do 15 mm)")
    pairs = []
    for m in sel:
        best = (1e9, None)
        for sid, rec in per_septum.items():
            if rec["tree"] is None:
                continue
            d = float(rec["tree"].query(m["pos"])[0])
            if d < best[0]:
                best = (d, rec)
        if best[1] and best[0] < 0.015:
            pairs.append((best[0] * 1000, m["code"], m["key"][-1], best[1]["name"],
                          best[1]["side"]))
    for d, code, side, name, sside in sorted(pairs):
        if side == sside or sside == "C":
            print(f"  {d:5.1f} mm   {code:<7} {name} ({sside})")

    if args.json:
        out = {"points": {"n": len(res), "median_mm": float(np.median(dists)),
                          "median_pct": med, "p": p},
               "lines": {"n": len(lres), "median_mm": float(np.median(ldist)),
                         "median_pct": lmed, "p": lp},
               "channels": [{"channel": ch, "n": n, "median_mm": round(d, 1),
                             "median_pct": round(q, 3)} for q, d, ch, n in sorted(rows)],
               "pairs": [{"dist_mm": round(d, 1), "code": c, "septum": n, "side": s}
                         for d, c, _, n, s in sorted(pairs)]}
        Path(args.json).write_text(json.dumps(out, ensure_ascii=False, indent=1),
                                   encoding="utf-8")
        print(f"\nzapisano {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
