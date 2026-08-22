"""Wyznaczenie rigu z prawdziwej anatomii -> data/rig_fitted.json.

Zamiast proporcji z tabeli bierzemy punkty kostne wprost z siatek Z-Anatomy:
wyrostki kolczyste (poziomy dla punktow BL i DU), wciecie szyjne mostka, wyrostek
mieczykowaty, gorny brzeg spojenia lonowego, kretarz wiekszy, szpara stawu
kolanowego, kostki, wyrostek barkowy. To sa dokladnie te punkty, o ktore opiera
sie klasyczny system proporcjonalny (gu du fen cun).

    python tools/calibrate_rig.py
"""
from __future__ import annotations

import json
import warnings
from pathlib import Path

import numpy as np
import pyvista as pv

warnings.simplefilter("ignore")
ROOT = Path(__file__).resolve().parent.parent
REF = ROOT / "assets" / "mesh" / "ref"
READY = ROOT / "assets" / "ready"

VERTEBRAE = (["Atlas (C1)", "Axis (C2)"] + [f"Vertebra C{i}" for i in range(3, 8)]
             + [f"Vertebra T{i}" for i in range(1, 13)]
             + [f"Vertebra L{i}" for i in range(1, 6)])
VERT_KEY = (["C1", "C2"] + [f"C{i}" for i in range(3, 8)]
            + [f"T{i}" for i in range(1, 13)] + [f"L{i}" for i in range(1, 6)])


def slug(name: str) -> str:
    import re
    s = name.strip().lower().replace("(", "").replace(")", "")
    return re.sub(r"[^a-z0-9]+", "-", s).strip("-")


def load(name: str) -> np.ndarray:
    """Punkty kosci w ukladzie acu3d (+X prawa, +Y przod, +Z gora)."""
    p = np.asarray(pv.read(REF / f"{slug(name)}.ply").points).copy()
    p[:, 0] *= -1.0
    p[:, 1] *= -1.0
    return p


def top_region(pts: np.ndarray, frac: float = 0.06) -> np.ndarray:
    z = pts[:, 2]
    return pts[z > z.max() - frac * (z.max() - z.min() + 1e-9) - 0.004]


def bottom_region(pts: np.ndarray, frac: float = 0.06) -> np.ndarray:
    z = pts[:, 2]
    return pts[z < z.min() + frac * (z.max() - z.min() + 1e-9) + 0.004]


def find_navel(skin_pts: np.ndarray, z_pub: float, z_xi: float):
    """Pepek zmierzony na skorze, a nie zalozony z proporcji.

    Pepek to DOLEK: na linii posrodkowej skora jest cofnieta wzgledem sasiedztwa
    kilka centymetrow w bok. Mierzymy wlasnie te roznice, bo jest odporna na to,
    jak wybierzemy zakres (odchylka od wygladzonego trendu potrafila wskazac
    zamiast pepka miejsce, gdzie brzuch sie zweza - pomylka o 2 cm).
    """
    fallback = z_pub + (z_xi - z_pub) * 5.0 / 13.0

    def profile(z, x0, x1, tol=0.008):
        b = skin_pts[(np.abs(skin_pts[:, 2] - z) < tol) & (np.abs(skin_pts[:, 0]) >= x0)
                     & (np.abs(skin_pts[:, 0]) <= x1) & (skin_pts[:, 1] > 0.05)]
        return float(np.percentile(b[:, 1], 90)) if len(b) > 4 else np.nan

    zs = np.linspace(z_pub + 0.08, z_xi - 0.04, 120)
    dent = np.array([(profile(z, 0.028, 0.048) - profile(z, 0.0, 0.010)) for z in zs])
    if np.all(np.isnan(dent)) or np.nanmax(dent) < 0.002:
        return fallback, "proporcja 5/13 (nie znaleziono dolka)"
    i = int(np.nanargmax(dent))
    return float(zs[i]), f"zmierzony na skorze, dolek {dent[i] * 1000:.1f} mm"


def main() -> int:
    out = {"units": "m",
           "coordinate_system": "+X = prawa strona, +Y = przod, +Z = gora",
           "source": "Z-Anatomy (CC BY-SA 4.0) - punkty kostne zmierzone z siatek"}

    skin = pv.read(READY / "skin.vtp")
    sp = np.asarray(skin.points)
    out["height"] = round(float(sp[:, 2].max() - sp[:, 2].min()), 4)
    out["floor_z"] = round(float(sp[:, 2].min()), 4)

    # --- kregi: wierzcholek wyrostka kolczystego = punkt najbardziej ku tylowi
    verts = {}
    for name, key in zip(VERTEBRAE, VERT_KEY):
        f = REF / f"{slug(name)}.ply"
        if not f.exists():
            continue
        p = load(name)
        spinous = p[np.argmin(p[:, 1])]              # najdalej ku tylowi (-Y)
        body = p[p[:, 1] > np.percentile(p[:, 1], 85)].mean(axis=0)   # trzon, ku przodowi
        verts[key] = {"spinous": [round(float(x), 4) for x in spinous],
                      "body": [round(float(x), 4) for x in body],
                      "z": round(float(spinous[2]), 4)}
    sac = load("Sacrum")
    cox = load("Coccyx")
    zt, zb = sac[:, 2].max(), sac[:, 2].min()
    for i in range(1, 5):                            # S1..S4 rownomiernie po kosci krzyzowej
        z = zt - (zt - zb) * (i - 0.5) / 4.0
        band = sac[np.abs(sac[:, 2] - z) < 0.006]
        pt = band[np.argmin(band[:, 1])] if len(band) else np.array([0, 0, z])
        verts[f"S{i}"] = {"spinous": [round(float(x), 4) for x in pt], "z": round(float(z), 4)}
    verts["coccyx"] = {"spinous": [round(float(x), 4) for x in cox[np.argmin(cox[:, 2])]],
                       "z": round(float(cox[:, 2].min()), 4)}
    out["vertebra"] = verts

    # --- mostek, spojenie lonowe
    man = load("Manubrium of sternum")
    mid = man[np.abs(man[:, 0]) < 0.007]
    notch = mid[np.argmax(mid[:, 2])]                # wciecie szyjne (REN-22)
    xi = load("Xiphoid process")
    xiphoid_top = xi[np.argmax(xi[:, 2])]            # kat mostkowo-zebrowy (REN-15/16)
    hip = load("Hip bone.r")
    pub = hip[(hip[:, 0] < 0.015) & (hip[:, 1] > 0.0)]
    pubis = pub[np.argmax(pub[:, 2])] if len(pub) else hip[np.argmin(hip[:, 2])]

    L = {"sternal_notch": notch, "xiphoid": xiphoid_top, "pubis": pubis}

    # --- konczyna dolna (prawa)
    fem = load("Femur.r")
    fem_top = top_region(fem, 0.12)
    trochanter = fem_top[np.argmax(fem_top[:, 0])]   # krętarz wiekszy: najbardziej bocznie
    head = fem_top[fem_top[:, 0] < np.percentile(fem_top[:, 0], 25)].mean(axis=0)
    tib = load("Tibia.r")
    knee = top_region(tib, 0.03).mean(axis=0)        # plateau piszczeli = szpara stawu
    med_mall = tib[np.argmin(tib[:, 2])]
    fib = load("Fibula.r")
    lat_mall = fib[np.argmin(fib[:, 2])]
    tal = load("Talus.r")
    ankle = tal.mean(axis=0)
    cal = load("Calcaneus.r")

    # --- konczyna gorna (prawa)
    hum = load("Humerus.r")
    hum_head = top_region(hum, 0.10).mean(axis=0)
    elbow = bottom_region(hum, 0.05).mean(axis=0)
    rad = load("Radius.r")
    uln = load("Ulna.r")
    wrist = (bottom_region(rad, 0.04).mean(axis=0) + bottom_region(uln, 0.04).mean(axis=0)) / 2
    sca = load("Scapula.r")
    sca_top = top_region(sca, 0.12)
    acromion = sca_top[np.argmax(sca_top[:, 0])]
    cla = load("Clavicle.r")

    # --- konce dloni i stopy z powierzchni skory
    hand_zone = sp[(np.abs(sp[:, 0] - wrist[0]) < 0.10) & (sp[:, 2] < wrist[2])
                   & (sp[:, 2] > wrist[2] - 0.25) & (np.abs(sp[:, 1] - wrist[1]) < 0.12)]
    fingertip = hand_zone[np.argmin(hand_zone[:, 2])] if len(hand_zone) else wrist
    foot_zone = sp[(sp[:, 2] < 0.12) & (sp[:, 0] > 0.0)]
    toetip = foot_zone[np.argmax(foot_zone[:, 1])] if len(foot_zone) else ankle

    L.update({"trochanter": trochanter, "femoral_head": head, "knee": knee,
              "medial_malleolus": med_mall, "lateral_malleolus": lat_mall, "ankle": ankle,
              "acromion": acromion, "humeral_head": hum_head, "elbow": elbow, "wrist": wrist,
              "fingertip": fingertip, "toetip": toetip,
              "heel": cal[np.argmin(cal[:, 1])], "clavicle_lat": cla[np.argmax(cla[:, 0])],
              "vertex": sp[np.argmax(sp[:, 2])], "chin": None})
    # broda = dolny brzeg zuchwy (z kosci, nie z obrysu skory - obrys mylil sie o 6 cm)
    man_pts = load("Mandible")
    chin_bone = man_pts[np.argmin(man_pts[:, 2])]
    L["chin"] = chin_bone - np.array([0.0, 0.0, 0.004])   # cienka warstwa tkanek miekkich
    out["landmarks"] = {k: [round(float(x), 4) for x in v] for k, v in L.items() if v is not None}

    # --- linijka cun tulowia: pubis 0, xiphoid 13, wciecie szyjne 22
    z_pub, z_xi, z_notch = float(pubis[2]), float(xiphoid_top[2]), float(notch[2])
    z_umb, umb_src = find_navel(sp, z_pub, z_xi)
    out["trunk_ruler"] = {"knots": [[round(z_pub, 4), 0.0], [round(z_umb, 4), 5.0],
                                    [round(z_xi, 4), 13.0], [round(z_notch, 4), 22.0]],
                          "umbilicus_source": umb_src}
    out["landmarks"]["umbilicus"] = [0.0, 0.0, round(z_umb, 4)]

    # --- krzywa kregoslupa: wysokosc -> polozenie trzonu (punkt na pewno wewnatrz tulowia)
    curve = sorted(([v["body"][2], v["body"][1]] for v in verts.values() if "body" in v))
    out["spine_curve"] = [[round(z, 4), round(y, 4)] for z, y in curve]

    # --- segmenty: os + dlugosc w cun
    def seg(p0, p1, cun, **kw):
        d = {"p0": [round(float(x), 4) for x in p0], "p1": [round(float(x), 4) for x in p1],
             "cun_length": cun, "length_m": round(float(np.linalg.norm(np.array(p1) - np.array(p0))), 4)}
        d.update(kw)
        return d

    out["segments"] = {
        "upper_arm": seg(hum_head, elbow, 9, cun_t0=0.18),
        "forearm": seg(elbow, wrist, 12),
        "hand": seg(wrist, fingertip, 7.5),
        "thigh": seg(head, knee, 19),   # os anatomiczna idzie przez glowe kosci udowej,
        #                              linijka 19 cun liczona jest od krętarza
        "shank": seg(knee, lat_mall, 16),
        "foot": seg(ankle, toetip, 10),
    }
    out["cun"] = {
        "trunk_vertical_m": round((z_notch - z_pub) / 22.0, 5),
        "trunk_lateral_m": None,     # uzupelniane nizej z rozstawu brodawek
    }

    # linia srodkowoobojczykowa (= linia sutkowa) lezy 4 cun od posrodkowej;
    # bierzemy srodek obojczyka miedzy koncem mostkowym a barkowym
    cl_med = cla[np.argmin(cla[:, 0])]
    cl_lat = cla[np.argmax(cla[:, 0])]
    midclav_x = float((cl_med[0] + cl_lat[0]) / 2.0)
    out["cun"]["trunk_lateral_m"] = round(midclav_x / 4.0, 5)
    out["landmarks"]["midclavicular"] = [round(midclav_x, 4),
                                         round(float(cl_lat[1]), 4), round(float(cl_lat[2]), 4)]
    out["cun"]["note"] = ("pion tulowia: 22 cun od spojenia do wciecia szyjnego; "
                          "poziom: 4 cun od linii posrodkowej do brodawki")

    (ROOT / "data" / "rig_fitted.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")

    print(f"wzrost modelu: {out['height']:.3f} m")
    print(f"wciecie szyjne z={z_notch:.3f}  mieczyk z={z_xi:.3f}  spojenie z={z_pub:.3f}")
    print(f"1 cun tulowia w pionie: {out['cun']['trunk_vertical_m'] * 100:.2f} cm, "
          f"w poziomie: {out['cun']['trunk_lateral_m'] * 100:.2f} cm")
    for k, s in out["segments"].items():
        print(f"  {k:10s} {s['length_m']:.3f} m / {s['cun_length']} cun "
              f"= {s['length_m'] / s['cun_length'] * 100:.2f} cm/cun")
    print(f"kregow zmierzonych: {len(verts)}  ->  data/rig_fitted.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
