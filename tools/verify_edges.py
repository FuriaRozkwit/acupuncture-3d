"""Walidacja punktow konczyn wzgledem KRAWEDZI KOSCI, po skorze.

Opisy punktow mowia "jeden palec bocznie od grzebienia piszczeli", "na przednim
brzegu strzalki", "przy tylnym brzegu kosci lokciowej" - to odleglosci PO CIELE
od wyczuwalnej krawedzi kosci, a nie katy. Tu mierzymy dokladnie to.

Dwa kroki, bo nie kazda krawedz da sie odczytac z siatki:
  1. WIARYGODNOSC - krawedz liczona na 60 poziomach; jesli surowy pomiar odbiega
     od wygladzonej linii o wiecej niz kilka stopni, siatka jest za rzadka i taka
     krawedz nie nadaje sie na odniesienie (piszczel: grzebien tak, brzeg
     przysrodkowy nie).
  2. POMIAR - odleglosc punktu od krawedzi po skorze, w cun tego segmentu.

    python tools/verify_edges.py [--all]
"""
from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
warnings.simplefilter("ignore")

from acu3d.model import Atlas  # noqa: E402

# krawedzie do sprawdzenia: nazwa -> (segment, struktura, strona)
EDGES = {
    "tibial_crest": ("shank", "Tibia", "anterior"),
    "tibia_medial": ("shank", "Tibia", "medial"),
    "fibula_anterior": ("shank", "Fibula", "anterior"),
    "fibula_posterior": ("shank", "Fibula", "posterior"),
    "ulna_posterior": ("forearm", "Ulna", "posterior"),
    "ulna_anterior": ("forearm", "Ulna", "anterior"),
    "ulna_medial": ("forearm", "Ulna", "medial"),
    "radius_lateral": ("forearm", "Radius", "lateral"),
    "fibula_lateral": ("shank", "Fibula", "lateral"),
}

# Czy KOSC jest dobrze wpasowana w skore, sprawdzamy na jej fragmencie PODSKORNYM
# (grzebien piszczeli, boczna strona strzalki, tyl kosci lokciowej, bok promieniowej).
# Przedni brzeg strzalki lezy pod prostownikami - z niego nic o dopasowaniu nie wynika.
SUBCUT = {"Tibia": "tibial_crest", "Fibula": "fibula_lateral",
          "Ulna": "ulna_posterior", "Radius": "radius_lateral"}

# Czego spodziewamy sie wg OPISU ZRODLOWEGO tego punktu (nie wg pamieci - opisy
# roznia sie miedzy szkolami: GB-39 wg WHO lezy przed strzalka, a wg podrecznika,
# na ktorym stoja te dane, za jej tylnym brzegiem). Dodatnie = w strone rosnacego
# azymutu (bocznie/do przodu), ujemne = przysrodkowo/do tylu.
EXPECT = {
    "ST-36": ("tibial_crest", 1.0, "jeden palec bocznie od grzebienia"),
    "ST-37": ("tibial_crest", 1.0, "jeden palec bocznie"),
    "ST-38": ("tibial_crest", 1.0, "jeden palec bocznie"),
    "ST-39": ("tibial_crest", 1.0, "jeden palec bocznie"),
    "ST-40": ("tibial_crest", 1.5, "dwa palce bocznie"),
    "GB-36": ("fibula_anterior", 0.2, "przy przednim brzegu strzalki"),
    "GB-37": ("fibula_anterior", 0.2, "przy przednim brzegu strzalki"),
    "GB-38": ("fibula_anterior", 0.2, "przy przednim brzegu strzalki"),
    "GB-35": ("fibula_posterior", -0.2, "w zaglebieniu przy tylnym brzegu strzalki"),
    "GB-39": ("fibula_posterior", -0.3, "miedzy tylnym brzegiem strzalki a sciegnami strzalkowymi"),
    "SP-6": ("tibia_medial", 0.0, "przy brzegu przysrodkowym piszczeli"),
    "SP-8": ("tibia_medial", -0.3, "tuz za brzegiem przysrodkowym"),
    "LR-5": ("tibia_medial", -0.3, "tuz za brzegiem przysrodkowym"),
    "LR-6": ("tibia_medial", -0.3, "tuz za brzegiem przysrodkowym"),
    "SI-7": ("ulna_anterior", -0.3, "bruzda miedzy przednim brzegiem kosci lokciowej a FCU"),
    "SI-6": ("ulna_posterior", 0.0, "szczelina na grzbiecie glowy kosci lokciowej"),
}
TOL = 0.5           # cun - tyle wolno sie roznic od opisu
NOISE_OK = 6.0      # stopni - powyzej tego krawedz uznajemy za nieczytelna z siatki
SKIN_OK = 12.0      # mm - krawedz podskorna z grubsza tyle ma nad soba; wiecej znaczy,
#                     ze kosc jest zle wpasowana w skore i nie wolno od niej mierzyc


def edge_noise(an, name: str) -> tuple:
    """(szum w stopniach, zakres skretu) - jak bardzo surowy pomiar skacze."""
    seg_name, struct, side = EDGES[name]
    seg = an.seg(seg_name)
    mesh = an.struct(struct).mesh
    pick = {"anterior": seg.e_ant, "posterior": -seg.e_ant,
            "medial": -seg.e_lat, "lateral": seg.e_lat}[side]
    shaft = (an.EDGE_LEVELS > 0.15) & (an.EDGE_LEVELS < 0.85)
    raw = np.full(len(an.EDGE_LEVELS), np.nan)
    for i, t in enumerate(an.EDGE_LEVELS):
        c = seg.interior(float(t))
        sl = mesh.slice(normal=seg.e_ax, origin=c)
        if sl.n_points < 8:
            continue
        v = np.asarray(sl.points) - c
        v -= np.outer(v @ seg.e_ax, seg.e_ax)
        k = int(np.argmax(v @ pick))
        raw[i] = np.degrees(np.arctan2(float(v[k] @ seg.e_lat), float(v[k] @ seg.e_ant)))
    sm = np.array([an.bone_edge_azimuth(name, float(t)) for t in an.EDGE_LEVELS])
    ok = (~np.isnan(raw)) & shaft       # same konce kosci to klykcie i kostki
    noise = float(np.std(raw[ok] - sm[ok])) if ok.sum() > 5 else float("inf")
    return noise, float(np.ptp(sm[shaft])), int(ok.sum())


def edge_depth(an, name: str) -> float:
    """Ile 'skory' jest nad krawedzia. Krawedzie z tabeli sa podskorne, wiec
    kilkanascie milimetrow oznacza rozjazd dopasowania kosci do skory, nie tkanke."""
    seg_name, struct, side = EDGES[name]
    seg = an.seg(seg_name)
    mesh = an.struct(struct).mesh
    pick = {"anterior": seg.e_ant, "posterior": -seg.e_ant,
            "medial": -seg.e_lat, "lateral": seg.e_lat}[side]
    pts = []
    for t in np.linspace(0.2, 0.85, 14):
        c = seg.interior(float(t))
        sl = mesh.slice(normal=seg.e_ax, origin=c)
        if sl.n_points < 8:
            continue
        v = np.asarray(sl.points) - c
        v -= np.outer(v @ seg.e_ax, seg.e_ax)
        pts.append(np.asarray(sl.points)[int(np.argmax(v @ pick))])
    if not pts:
        return float("inf")
    _, foot = an.skin_cast.find_closest_cell(np.array(pts), return_closest_point=True)
    return float(np.median(np.linalg.norm(np.atleast_2d(foot) - np.array(pts), axis=1)) * 1000)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--all", action="store_true", help="wypisz takze punkty zgodne z opisem")
    args = ap.parse_args()
    atlas = Atlas()
    an = atlas.anatomy
    if not hasattr(an, "bone_edge_azimuth"):
        print("ten rig nie zna krawedzi kosci", file=sys.stderr)
        return 2

    an.BONE_EDGES.update(EDGES)
    print("KRAWEDZIE - czy da sie je odczytac z siatki")
    usable = {}
    fitted = {}
    for bone, probe in SUBCUT.items():
        d = edge_depth(an, probe)
        fitted[bone] = d <= SKIN_OK
        print(f"  kosc {bone:8s} podskornie ({probe}): {d:5.1f} mm skory"
              f"  -> {'wpasowana' if fitted[bone] else 'ZLE WPASOWANA (rozjazd siatek)'}")
    print()
    for name in EDGES:
        noise, span, n = edge_noise(an, name)
        bone = EDGES[name][1]
        ok = noise <= NOISE_OK and fitted.get(bone, True)
        usable[name] = ok
        why = ("uzyteczna" if ok else
               "ZA RZADKA SIATKA" if noise > NOISE_OK else f"kosc {bone} zle wpasowana")
        print(f"  {name:18s} szum {noise:5.1f} st., skret {span:5.1f} st.  -> {why}")

    print("\nPUNKTY - odleglosc po skorze od krawedzi (cun)")
    bad = 0
    for code, (edge, want, says) in EXPECT.items():
        key = code + "·R"
        if key not in atlas.by_key:
            continue
        seg_name = EDGES[edge][0]
        seg = an.seg(seg_name)
        p = np.asarray(atlas.by_key[key]["pos"])
        t = float(np.clip(np.dot(p - seg.p0, seg.e_ax) / seg.length, 0.0, 1.0))
        c = seg.interior(t)
        v = p - c
        v -= seg.e_ax * float(np.dot(v, seg.e_ax))
        az = float(np.degrees(np.arctan2(np.dot(v, seg.e_lat), np.dot(v, seg.e_ant))))
        az0 = an.bone_edge_azimuth(edge, t)
        got = an._arc_len(seg, t, az0, az) / seg.cun_m * (1.0 if az > az0 else -1.0)
        err = got - want
        mark = "" if abs(err) <= TOL else "  <-- ROZNI SIE OD OPISU"
        if not usable[edge]:
            mark = "  (odniesienie niepewne - wynik orientacyjny)"
        elif abs(err) > TOL:
            bad += 1
        if args.all or mark:
            print(f"  {code:6s} wg opisu {want:+5.2f}, zmierzone {got:+5.2f} cun "
                  f"({err:+.2f})  [{says}]{mark}")
    print(f"\npunktow niezgodnych z opisem na czytelnych krawedziach: {bad}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
