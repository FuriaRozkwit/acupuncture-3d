"""Badanie struktury kanalow: co w niej jest naprawde, a co tylko sie wydaje.

Model ma teraz dwie niezalezne warstwy opisujace to samo cialo: TOPOLOGIE
(kto z kim sie laczy, wedlug tradycji) i GEOMETRIE (gdzie to lezy, zmierzone
na siatkach). Dopiero zestawienie ich ze soba daje cos, czego zadna z osobna
nie mowi - bo tradycja nigdzie nie podaje odleglosci w milimetrach, a siatka
nic nie wie o kanalach.

Kazdy test ma tlo. Liczba bez porownania z rozkladem losowym nic nie znaczy,
a "wyglada na zgodne" nie jest wynikiem.

    python tools/study_topology.py [--json out.json]
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

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"

# Obieg dobowy: kolejnosc, w ktorej wedlug tradycji qi przechodzi z kanalu
# do kanalu. To jest twierdzenie o CIAGLOSCI - da sie je zmierzyc.
CLOCK = ["LU", "LI", "ST", "SP", "HT", "SI", "BL", "KI", "PC", "SJ", "GB", "LR"]
FIVE_SHU = ["jing_well", "ying_spring", "shu_stream", "jing_river", "he_sea"]


def _safe(text: str) -> None:
    enc = getattr(sys.stdout, "encoding", None) or "utf-8"
    sys.stdout.write(text.encode(enc, errors="replace").decode(enc, errors="replace") + "\n")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", default="")
    args = ap.parse_args()

    from acu3d.channels import ChannelMap
    from acu3d.model import Atlas

    atlas = Atlas()
    cm = ChannelMap(atlas)
    rng = np.random.default_rng(20260902)
    out: dict = {}

    def pos(code, side="R"):
        m = atlas.by_key.get(f"{code}·{side}") or atlas.by_key.get(code)
        return None if m is None else np.asarray(m["pos"])

    # =================================================================== A
    _safe("A. GRAF KANALOW - kto z kim graniczy i ktore kanaly sa wezlami")
    _safe("=" * 74)
    deg: dict = {}
    for code, chans in cm.meets.items():
        for c in chans:
            deg.setdefault(c, set()).update(x for x in chans if x != c)
    shared: dict = {}
    for code, chans in cm.meets.items():
        for i, a in enumerate(chans):
            for b in chans[i + 1:]:
                shared[tuple(sorted((a, b)))] = shared.get(tuple(sorted((a, b))), 0) + 1
    hubs = sorted(((len(v), k) for k, v in deg.items()), reverse=True)
    for n, c in hubs[:8]:
        _safe(f"   {c:4s} styka sie z {n:2d} kanalami: {', '.join(sorted(deg[c]))}")
    top = sorted(shared.items(), key=lambda kv: -kv[1])[:6]
    _safe("   najczestsze pary (ile wspolnych punktow skrzyzowania):")
    for (a, b), n in top:
        _safe(f"      {a}-{b}: {n}")
    out["hubs"] = [{"channel": c, "degree": n} for n, c in hubs]

    # =================================================================== B
    _safe("")
    _safe("B. OBIEG DOBOWY - czy koniec kanalu lezy przy poczatku nastepnego")
    _safe("=" * 74)
    _safe("   Tradycja mowi, ze qi przechodzi z kanalu do kanalu na koncach.")
    _safe("   Jesli to opis ciala, a nie metafora, styki musza byc BLISKIE.")
    gaps, rand = [], []
    ends = {}
    for ch in CLOCK:
        seq = cm.order.get(ch, [])
        if seq:
            ends[ch] = (seq[0], seq[-1])
    for i, ch in enumerate(CLOCK):
        nx = CLOCK[(i + 1) % len(CLOCK)]
        if ch not in ends or nx not in ends:
            continue
        a, b = pos(ends[ch][1]), pos(ends[nx][0])
        if a is None or b is None:
            continue
        d = float(np.linalg.norm(a - b))
        gaps.append((d, ch, nx, ends[ch][1], ends[nx][0]))
    # Tlo musi byc rozkladem TEJ SAMEJ statystyki. Porownywanie mediany
    # dwunastu stykow z pula pojedynczych odleglosci losowych zanizalo istotnosc:
    # losujemy wiec cale sparowania i bierzemy ich mediany.
    tails = [pos(ends[c][1]) for c in ends if pos(ends[c][1]) is not None]
    heads = [pos(ends[c][0]) for c in ends if pos(ends[c][0]) is not None]
    obs = float(np.median([g[0] for g in gaps]))
    for _ in range(20000):
        q = rng.permutation(len(heads))
        d = [float(np.linalg.norm(tails[i] - heads[q[i]]))
             for i in range(len(tails)) if q[i] != i]
        if d:
            rand.append(np.median(d))
    med = obs * 1000
    medr = float(np.median(rand)) * 1000
    p = float(np.mean(np.array(rand) <= obs))
    for d, ch, nx, ca, cb in sorted(gaps):
        _safe(f"      {ch:3s} -> {nx:3s}   {ca:7s} -> {cb:7s}   {d * 1000:6.0f} mm")
    _safe(f"   mediana styku: {med:.0f} mm,  tlo (losowe sparowania konca z poczatkiem):"
          f" {medr:.0f} mm,  p = {p:.4f}")
    # Ktore styki sa dalekie - i czy to te, ktore tradycja opisuje jako przejscia
    # WEWNETRZNE, niewidoczne na skorze. Model nic o nich nie wie.
    INTERNAL = {("SP", "HT"), ("KI", "PC"), ("LR", "LU")}
    order = sorted(gaps, reverse=True)
    _safe("   przejscia wewnetrzne wg tradycji sa na miejscach: " + ", ".join(
        f"{i + 1}. {c}->{n} ({d * 1000:.0f} mm)"
        for i, (d, c, n, _, _) in enumerate(order) if (c, n) in INTERNAL))
    ext = [d for d, c, n, _, _ in gaps if (c, n) not in INTERNAL]
    intr = [d for d, c, n, _, _ in gaps if (c, n) in INTERNAL]
    if ext and intr:
        _safe(f"   mediana: przejscia zewnetrzne {np.median(ext) * 1000:.0f} mm, "
              f"wewnetrzne {np.median(intr) * 1000:.0f} mm")
    out["clock"] = {"median_mm": round(med, 1), "random_mm": round(medr, 1), "p": round(p, 4),
                    "pairs": [{"from": c, "to": n, "mm": round(d * 1000, 1)}
                              for d, c, n, _, _ in sorted(gaps)]}

    # =================================================================== C
    _safe("")
    _safe("C. PIEC PUNKTOW SHU - czy to samopodobienstwo, czy tylko kolejnosc")
    _safe("=" * 74)
    _safe("   Piec punktow transportu odwzorowuje caly obieg na jednej konczynie:")
    _safe("   od studni w koncu palca po morze w lokciu albo kolanie. Jesli to")
    _safe("   odwzorowanie skali, ich POLOZENIA WZGLEDNE musza byc podobne")
    _safe("   we wszystkich kanalach - nie tylko ich kolejnosc.")
    role_pos: dict = {}
    for code, rs in cm.roles.items():
        for r in rs:
            if r["role"] in FIVE_SHU:
                role_pos.setdefault(r["of"], {})[r["role"]] = code
    rows = []
    for ch, m in sorted(role_pos.items()):
        if len(m) < 5 or ch not in CLOCK:
            continue
        seq = cm.order.get(ch, [])
        if not seq:
            continue
        pts = [pos(m[r]) for r in FIVE_SHU]
        if any(p_ is None for p_ in pts):
            continue
        # odleglosc od studni, znormalizowana dlugoscia calego ukladu piatki
        d = np.array([float(np.linalg.norm(p_ - pts[0])) for p_ in pts])
        if d[-1] < 1e-6:
            continue
        rows.append((ch, d / d[-1], d[-1]))
    if rows:
        M = np.array([r[1] for r in rows])
        _safe(f"   kanalow z pelna piatka: {len(rows)}")
        _safe("        studnia  zrodlo  strumien  rzeka   morze")
        for ch, v, tot in rows:
            _safe(f"   {ch:4s}  " + "  ".join(f"{x:6.2f}" for x in v)
                  + f"   (calosc {tot * 1000:.0f} mm)")
        sd = M.std(axis=0)
        _safe(f"   odchylenie miedzy kanalami: "
              + "  ".join(f"{x:6.3f}" for x in sd))
        # tlo: pieciu losowych punktow tego samego kanalu
        null = []
        for _ in range(2000):
            ch = rows[rng.integers(0, len(rows))][0]
            seq = [c for c in cm.order.get(ch, []) if pos(c) is not None]
            if len(seq) < 5:
                continue
            pick = sorted(rng.choice(len(seq), 5, replace=False))
            pp = [pos(seq[k]) for k in pick]
            dd = np.array([float(np.linalg.norm(x - pp[0])) for x in pp])
            if dd[-1] > 1e-6:
                null.append(dd / dd[-1])
        if null:
            nsd = np.array(null).std(axis=0)
            _safe(f"   to samo dla piatek losowych:  "
                  + "  ".join(f"{x:6.3f}" for x in nsd))
            tighter = float(np.mean(sd[1:4] < nsd[1:4]))
            _safe(f"   piatka shu jest ciasniejsza niz losowa w "
                  f"{tighter:.0%} pozycji srodkowych")
        out["five_shu"] = {"channels": [r[0] for r in rows],
                           "mean": M.mean(axis=0).round(3).tolist(),
                           "sd": sd.round(3).tolist()}

    # =================================================================== D
    _safe("")
    _safe("D. GDZIE LEZA SKRZYZOWANIA - geografia spotkan")
    _safe("=" * 74)
    seg_of = {p["code"]: (p["anchor"].get("seg") or "struct")
              for p in atlas.points}
    cross = set(cm.crossings)
    tally: dict = {}
    for code, sg in seg_of.items():
        t = tally.setdefault(sg, [0, 0])
        t[1] += 1
        if code in cross:
            t[0] += 1
    for sg, (c, n) in sorted(tally.items(), key=lambda kv: -kv[1][0] / max(kv[1][1], 1)):
        _safe(f"   {sg:10s} {c:3d} / {n:3d} punktow to skrzyzowania  ({c / max(n,1):5.0%})")
    zc = [pos(c)[2] for c in cross if pos(c) is not None]
    za = [pos(p["code"])[2] for p in atlas.points if pos(p["code"]) is not None]
    _safe(f"   mediana wysokosci: skrzyzowania {np.median(zc):.2f} m, "
          f"wszystkie punkty {np.median(za):.2f} m")
    out["crossing_geography"] = {k: {"cross": v[0], "all": v[1]} for k, v in tally.items()}

    # =================================================================== E
    _safe("")
    _safe("E. TOPOLOGIA A GEOMETRIA - czy bliskosc w grafie to bliskosc w ciele")
    _safe("=" * 74)
    _safe("   Gdyby kanaly laczyly sie po prostu z sasiadami, topologia bylaby")
    _safe("   tylko innym zapisem geometrii i nic by nie wnosila.")
    cen = {}
    for ch in CLOCK:
        ps = [pos(c) for c in cm.order.get(ch, []) if pos(c) is not None]
        if ps:
            cen[ch] = np.mean(ps, axis=0)
    linked, unlinked = [], []
    for i, a in enumerate(CLOCK):
        for b in CLOCK[i + 1:]:
            if a not in cen or b not in cen:
                continue
            d = float(np.linalg.norm(cen[a] - cen[b]))
            key = tuple(sorted((a, b)))
            (linked if key in shared else unlinked).append(d)
    if linked and unlinked:
        ml, mu = np.median(linked) * 1000, np.median(unlinked) * 1000
        allv = np.array(linked + unlinked)
        lab = np.array([1] * len(linked) + [0] * len(unlinked))
        obs = np.median(allv[lab == 1]) - np.median(allv[lab == 0])
        perm = []
        for _ in range(20000):
            q = rng.permutation(lab)
            perm.append(np.median(allv[q == 1]) - np.median(allv[q == 0]))
        pv = float(np.mean(np.abs(perm) >= abs(obs)))
        _safe(f"   pary kanalow ze wspolnym skrzyzowaniem ({len(linked)}): "
              f"srodki oddalone o {ml:.0f} mm")
        _safe(f"   pary bez wspolnego skrzyzowania ({len(unlinked)}):     {mu:.0f} mm")
        _safe(f"   roznica {obs * 1000:+.0f} mm,  p = {pv:.4f}")
        out["topology_vs_geometry"] = {"linked_mm": round(ml, 1), "unlinked_mm": round(mu, 1),
                                       "p": round(pv, 4)}

    if args.json:
        Path(args.json).write_text(json.dumps(out, ensure_ascii=False, indent=1),
                                   encoding="utf-8")
        _safe(f"\nzapisano {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
