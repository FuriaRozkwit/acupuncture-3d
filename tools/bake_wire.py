"""Wypieka wektorowa siatke cunowa - cale cialo jako kilka tysiecy wierzcholkow.

Po co
-----
Warstwy (skora 90 tys. trojkatow, miesnie 1,08 mln, powiezi 237 tys.) sluzyly do
tego, zeby MODEL WIEDZIAL, co gdzie lezy. Do pracy z napieciami ta wiedza jest
potrzebna, ale same siatki juz nie - wystarczy jej wynik. Ten skrypt przenosi
wiedze z siatek do skalarow na rzadkiej siatce cunowej:

    dla kazdego wierzcholka:  grubosc tkanki miekkiej, glebokosc miesnia,
                              odleglosc do powiezi, najblizszy miesien

Potem przegladarka w trybie `--wire` nie wczytuje ani jednej siatki anatomicznej,
a mimo to potrafi powiedziec, co lezy pod oznaczeniem napiecia.

Budowa siatki
-------------
Kazdy odcinek ciala dostaje pierscienie co `--step` cun, a kazdy pierscien
`--rays` promieni rozlozonych po azymucie. Wierzcholki lacza sie w dwie rodziny
linii: obwodowe (pierscien) i podluzne (ten sam azymut przez kolejne pierscienie).
To jest ten sam uklad, w ktorym opisane sa punkty - cun wzdluz, azymut wokol -
wiec siatka nie jest ozdoba, tylko widocznym ukladem wspolrzednych atlasu.

    python tools/bake_wire.py [--step 0.5] [--rays 24]
"""
from __future__ import annotations

import argparse
import json
import sys
import warnings
from pathlib import Path

import math

import numpy as np

warnings.simplefilter("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"

# odcinki parzyste (maja lewy odpowiednik) i posrodkowe (pelny obrot)
LIMB = ("upper_arm", "forearm", "hand", "thigh", "shank", "foot")
CORE = ("trunk", "neck", "head")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--step", type=float, default=0.5, help="odstep pierscieni w cun")
    ap.add_argument("--rays", type=int, default=24, help="promieni na pierscien")
    ap.add_argument("--max-edge", type=float, default=0.085,
                    help="metry - dluzsza krawedz to przeskok, nie sasiedztwo")
    ap.add_argument("--out", default=str(DATA / "wire.npz"))
    args = ap.parse_args()

    import pyvista as pv
    from scipy.spatial import cKDTree

    from acu3d.anatomy import Anatomy
    from acu3d.fascia import FasciaSet
    from acu3d.muscles import MuscleSet

    an = Anatomy()
    ms = MuscleSet()
    try:
        fa = FasciaSet()
        fascia_pts = np.asarray(fa.mesh.points)
        fascia_tree = cKDTree(fascia_pts)
    except (FileNotFoundError, OSError):
        fascia_tree = None
        print("brak warstwy powieziowej - pomijam odleglosc do powiezi")

    bone_tree = cKDTree(np.asarray(an.skeleton.points))
    mus_pts = np.asarray(ms.mesh.points)
    mus_tree = cKDTree(mus_pts)
    mus_mid = np.asarray(ms.mesh.point_data["mid"])

    az = np.arange(args.rays) * (360.0 / args.rays)
    verts, rings, seg_id, cun_val, az_val, hit = [], [], [], [], [], []
    seg_names: list[str] = []

    def add_segment(name: str, seg, mirror: bool) -> None:
        """Pierscienie jednego odcinka; `mirror` odbija gotowe wierzcholki na lewa strone."""
        n_cun = seg.cun_length
        levels = np.arange(0.0, n_cun + 1e-6, args.step)
        sid = len(seg_names)
        seg_names.append(("L·" if mirror else "R·" if seg.name in LIMB else "") + name)
        for c in levels:
            t = seg.t_from_cun(float(c))
            if not (-0.02 <= t <= 1.02):
                continue
            ring = []
            tc = float(np.clip(t, 0.0, 1.0))
            base = np.asarray(seg.interior(tc), dtype=float)
            for a in az:
                p = np.asarray(seg.surface_by_azimuth(tc, float(a)), dtype=float)
                # cast() przy nietrafieniu siada na najblizszej skorze, wiec punkt
                # schodzi z wlasnego promienia - po tym poznajemy ucieczke
                rad = math.radians(float(a))
                d = np.cos(rad) * seg.e_ant + np.sin(rad) * seg.e_lat
                d = d / max(np.linalg.norm(d), 1e-9)
                w = p - base
                off = float(np.linalg.norm(w - d * float(np.dot(w, d))))
                if mirror:
                    p = p * np.array([-1.0, 1.0, 1.0])
                ring.append(len(verts))
                verts.append(p)
                seg_id.append(sid)
                cun_val.append(float(c))
                az_val.append(float(a))
                hit.append(off < 0.002)
            rings.append(ring)

    for name in CORE:
        add_segment(name, an.seg(name), mirror=False)
    for name in LIMB:
        seg = an.seg(name)
        add_segment(name, seg, mirror=False)
    for name in LIMB:                      # lewa strona to odbicie, nie drugi pomiar
        add_segment(name, an.seg(name), mirror=True)

    # Promien z tulowia w bok trafia w przylegajace ramie, nie w tulow. Poznajemy
    # to po tym, ze trafienie lezy blizej osi konczyny niz osi wlasnego odcinka.
    P = np.array(verts, dtype=np.float32)
    hit = np.array(hit, dtype=bool)
    core_ids = {i for i, n in enumerate(seg_names) if n in CORE}
    limb_axes = []
    for nm in ("upper_arm", "forearm"):
        sg = an.seg(nm)
        for sx in (1.0, -1.0):
            limb_axes.append((np.array(sg.p0) * [sx, 1, 1], np.array(sg.p1) * [sx, 1, 1]))
    sid_arr = np.array(seg_id)
    core_mask = np.isin(sid_arr, list(core_ids)) & hit

    def dist_to_axis(pts, p0, p1):
        v = p1 - p0
        t = np.clip((pts - p0) @ v / max(float(v @ v), 1e-9), 0.0, 1.0)
        return np.linalg.norm(pts - (p0 + t[:, None] * v), axis=1)

    idx = np.flatnonzero(core_mask)
    if len(idx):
        pts = P[idx].astype(float)
        own = np.full(len(idx), np.inf)
        for cid in core_ids:
            sg = an.seg(seg_names[cid])
            sel = sid_arr[idx] == cid
            if sel.any():
                own[sel] = dist_to_axis(pts[sel], np.array(sg.p0), np.array(sg.p1))
        near_limb = np.min([dist_to_axis(pts, a, b) for a, b in limb_axes], axis=0)
        stolen = near_limb < own
        hit[idx[stolen]] = False
        print(f"promieni tulowia przechwyconych przez konczyne: {int(stolen.sum())}")

    print(f"wierzcholkow: {len(P)}  pierscieni: {len(rings)}  odcinkow: {len(seg_names)}")

    # --- skalary: cala wiedza z siatek sprowadzona do liczb -----------------
    flesh_mm = bone_tree.query(P)[0] * 1000.0                 # grubosc tkanki miekkiej
    d_mus, i_mus = mus_tree.query(P)
    muscle_mm = d_mus * 1000.0                                # jak glaboko lezy brzusiec
    muscle_id = mus_mid[i_mus].astype(np.int16)
    if fascia_tree is not None:
        fascia_mm = fascia_tree.query(P)[0] * 1000.0
    else:
        fascia_mm = np.full(len(P), -1.0)

    # --- linie: obwodowe i podluzne ----------------------------------------
    edges = []
    for ring in rings:
        for k in range(len(ring)):
            edges.append((ring[k], ring[(k + 1) % len(ring)]))
    for a, b in zip(rings, rings[1:]):
        if len(a) == len(b) and seg_id[a[0]] == seg_id[b[0]]:
            edges += list(zip(a, b))
    E = np.array(edges, dtype=np.int32)

    # Dwa odsiewy. Promien, ktory uciekl przez otwor w siatce (oczodol, szpara
    # miedzy palcami), nie opisuje powierzchni - jego wierzcholek zostaje, bo ma
    # wlasciwe cun i azymut, ale nie laczymy go z sasiadami. Krawedz dluzsza niz
    # --max-edge to przeskok (promien z tulowia, ktory trafil w ramie), a nie
    # sasiedztwo - takie wlasnie robily "gwiazdy" w srodku sylwetki.
    H = hit
    keep = H[E[:, 0]] & H[E[:, 1]]
    seg_a = np.array(seg_id, dtype=np.int32)
    keep &= seg_a[E[:, 0]] == seg_a[E[:, 1]]
    P_all = np.array(verts, dtype=float)
    elen = np.linalg.norm(P_all[E[:, 0]] - P_all[E[:, 1]], axis=1)
    keep &= elen <= args.max_edge
    print(f"krawedzie: {len(E)} -> {int(keep.sum())} "
          f"(odsiane: {int((~H).sum())} ucieczek promienia, "
          f"{int((elen > args.max_edge).sum())} przeskokow)")
    E = E[keep]

    np.savez_compressed(
        args.out, points=P, edges=E,
        seg_id=np.array(seg_id, dtype=np.int16),
        cun=np.array(cun_val, dtype=np.float32),
        az=np.array(az_val, dtype=np.float32),
        flesh_mm=flesh_mm.astype(np.float32),
        muscle_mm=muscle_mm.astype(np.float32),
        muscle_id=muscle_id,
        hit=H,
        fascia_mm=fascia_mm.astype(np.float32))

    meta = {
        "source": "wyliczone z siatek Z-Anatomy i BodyParts3D (CC BY-SA), tools/bake_wire.py",
        "step_cun": args.step, "rays": args.rays,
        "segments": seg_names,
        "vertices": int(len(P)), "edges": int(len(E)),
        "scalars": {
            "flesh_mm": "odleglosc od skory do najblizszej kosci - grubosc tkanki miekkiej",
            "muscle_mm": "odleglosc do najblizszej powierzchni miesnia",
            "muscle_id": "indeks miesnia w data/muscles.json",
            "fascia_mm": "odleglosc do najblizszej struktury lacznotkankowej",
        },
        "muscles": [{"id": m["id"], "name": m["name"], "side": m["side"],
                     "axis": m["axis"], "centroid": m["centroid"]}
                    for m in json.loads((DATA / "muscles.json").read_text(
                        encoding="utf-8"))["muscles"]],
    }
    (DATA / "wire.json").write_text(json.dumps(meta, ensure_ascii=False, indent=1),
                                    encoding="utf-8")

    size = Path(args.out).stat().st_size / 1024
    print(f"krawedzi: {len(E)}")
    print(f"tkanka miekka: {flesh_mm.min():.0f}-{flesh_mm.max():.0f} mm, "
          f"mediana {np.median(flesh_mm):.0f} mm")
    print(f"zapisano {args.out} ({size:.0f} kB) i data/wire.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
