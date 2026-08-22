"""Rozwiazywanie kotwic punktow do wspolrzednych 3D.

Kotwica opisuje punkt tak, jak robi to podrecznik: *ktory segment ciala*,
*ile cun* wzdluz niego i *po ktorej stronie*. Dopiero rig zamienia to na metry,
wiec ten sam plik danych dziala dla kazdego wzrostu i kazdych proporcji modelu.

Pola kotwicy (`anchor`):
    seg        nazwa segmentu: trunk | neck | head | upper_arm | forearm | hand |
               thigh | shank | foot
  polozenie wzdluz osi (jedno z):
    t          ulamek dlugosci segmentu 0..1
    cun        odleglosc w cun od konca `from` ("p0" = proksymalny, "p1" = dystalny)
    z_frac     wysokosc jako ulamek wzrostu (segmenty pionowe)
    z_from     nazwa punktu odniesienia tulowia ("pubis", "umbilicus", "xiphoid",
               "sternal_notch", "nipple") + `z_cun` = ile cun powyzej (ujemne = ponizej)
    vertebra   "C7", "T1".."T12", "L1".."L5", "S1".."S4", "coccyx" (+ `dz_cun`)
    arc_cun    dla skalpu: odleglosc w cun po luku strzalkowym od przedniej linii
               wlosow (0) przez czubek glowy do tylnej linii wlosow (12)
  polozenie na obwodzie (jedno z):
    az         azymut w stopniach: 0 = przod, +90 = bok, 180 = tyl, -90 = przysrodkowo
    lat        odleglosc w cun od linii posrodkowej + `face`: anterior|posterior|lateral
  dodatkowo:
    out        wypchniecie ponad powierzchnie w cun (domyslnie 0)
"""
from __future__ import annotations

import numpy as np

from .rig import Rig

SCALP_ARC_CUN = 12.0          # przednia linia wlosow -> tylna, po luku przez czubek


def _t_for_z(rig: Rig, seg_name: str, z: float) -> float:
    seg = rig.seg(seg_name)
    dz = seg.p1[2] - seg.p0[2]
    if abs(dz) < 1e-9:
        return 0.5
    return float(np.clip((z - seg.p0[2]) / dz, 0.0, 1.0))


def _axial_t(rig: Rig, a: dict) -> float:
    seg = rig.seg(a["seg"])
    if "t" in a:
        return float(a["t"])
    if "cun" in a:
        return seg.t_from_cun(float(a["cun"]), a.get("from", "p0"))
    if "z_frac" in a:
        return _t_for_z(rig, a["seg"], float(a["z_frac"]) * rig.H)
    if "z_from" in a:
        cun = rig.landmark_cun[a["z_from"]] + float(a.get("z_cun", 0.0))
        return _t_for_z(rig, a["seg"], rig.z_of_cun(cun))
    if "vertebra" in a:
        z = rig.vertebra_z[a["vertebra"]]
        z += float(a.get("dz_cun", 0.0)) * rig.cun_trunk_lat
        return _t_for_z(rig, a["seg"], z)
    raise ValueError(f"kotwica bez polozenia osiowego: {a}")


def scalp_arc(rig: Rig, s_cun: float, lat_cun: float = 0.0, out_cun: float = 0.0) -> np.ndarray:
    """Punkt na sklepieniu czaszki wg luku strzalkowego (0 = przednia linia wlosow)."""
    head = rig.seg("head")
    z_front = rig.levels["hairline_ant"]
    t_front = _t_for_z(rig, "head", z_front)
    # luk: od przedniej linii wlosow (az 0) przez czubek (t=1) do tylu (az 180)
    ts_f = np.linspace(t_front, 1.0, 40)
    front = np.array([head.surface_by_azimuth(float(t), 0.0) for t in ts_f])
    ts_b = np.linspace(1.0, t_front - 0.10, 40)
    back = np.array([head.surface_by_azimuth(float(t), 180.0) for t in ts_b])
    arc = np.vstack([front, back[1:]])
    d = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(arc, axis=0), axis=1))])
    s = float(np.clip(s_cun / SCALP_ARC_CUN, 0.0, 1.0)) * d[-1]
    p = np.array([np.interp(s, d, arc[:, i]) for i in range(3)])
    if lat_cun:
        t_here = _t_for_z(rig, "head", p[2])
        face = "anterior" if p[1] >= head.axis_point(t_here)[1] else "posterior"
        # na samym czubku przekroj jest maly - trzymamy sie promienia pierscienia
        p = head.surface_by_lateral(t_here, lat_cun, face)
        if p[2] < 1.0:
            p[2] = max(p[2], 0.0)
    if out_cun:
        c = head.axis_point(_t_for_z(rig, "head", p[2]))
        n = p - c
        p = p + n / (np.linalg.norm(n) + 1e-9) * out_cun * head.cun_lat
    return p


def struct_base(body, a: dict) -> np.ndarray:
    """Punkt startowy kotwicy strukturalnej - lezy wewnatrz ciala, przy kosci."""
    region = a.get("region")
    ref = body.struct(a["struct"])
    along = float(a.get("along", 0.5))
    long_dir = body.frame(region)["long"] if region else None
    base = ref.at(along, long_dir)
    if a.get("toward"):
        other = body.struct(a["toward"])
        base = base + (other.at(along, long_dir) - base) * float(a.get("frac", 0.4))
    for name, mm in (a.get("shift") or {}).items():
        base = base + body.direction(name, region) * (float(mm) / 1000.0)
    return base


def resolve_struct(body, a: dict) -> np.ndarray:
    """Kotwica przy konkretnej strukturze: paliczek, kosc srodrecza, skrawek ucha.

    Zamiast zgadywac kat na uproszczonym stozku, startujemy z realnej kosci
    i szukamy skory w anatomicznie nazwanym kierunku (grzbietowo, promieniowo...).
    """
    base = struct_base(body, a)
    d = body.direction(a.get("dir", "dorsal"), a.get("region"))
    return body.cast_side(base, d)


def resolve(rig: Rig, anchor: dict) -> np.ndarray:
    """Zwraca wspolrzedne punktu po PRAWEJ stronie ciala (lewa = odbicie X)."""
    a = anchor
    if "struct" in a:
        return resolve_struct(rig, a)
    if "arc_cun" in a:
        return scalp_arc(rig, float(a["arc_cun"]), float(a.get("lat", 0.0)),
                         float(a.get("out", 0.0)))
    seg = rig.seg(a["seg"])
    t = _axial_t(rig, a)
    out = float(a.get("out", 0.0))
    if "az" in a:
        return seg.surface_by_azimuth(t, float(a["az"]), out)
    if "lat" in a:
        return seg.surface_by_lateral(t, float(a["lat"]), a.get("face", "anterior"), out)
    raise ValueError(f"kotwica bez polozenia obwodowego: {a}")


def resolve_sides(rig: Rig, point: dict) -> list:
    """Lista (side, xyz): dla punktow parzystych zwraca prawa i lewa strone."""
    p = resolve(rig, point["anchor"])
    if point.get("side") == "midline":
        p = p.copy()
        p[0] = 0.0
        return [("C", p)]
    left = p.copy()
    left[0] *= -1.0
    return [("R", p), ("L", left)]


def normal_at(rig: Rig, anchor: dict, pos: np.ndarray) -> np.ndarray:
    """Przyblizona normalna powierzchni w punkcie - do odsuwania etykiet."""
    seg = rig.seg(anchor.get("seg", "head"))
    t = _axial_t(rig, anchor) if "arc_cun" not in anchor else 1.0
    axis = seg.axis_point(t)
    n = pos - axis
    n[2] *= 0.35
    norm = np.linalg.norm(n)
    return n / norm if norm > 1e-9 else np.array([0.0, 1.0, 0.0])
