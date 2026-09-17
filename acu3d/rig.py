"""Szkielet proporcjonalny ciala.

Geometria budowana jest z ulamkow wzrostu (`data/rig.json`), a *mierzona* w cun -
klasycznych jednostkach proporcjonalnych, w ktorych opisane sa lokalizacje punktow.
Kazdy segment zna swoja dlugosc w cun, wiec kotwica typu "3 cun ponizej ST-35"
przelicza sie jednoznacznie na metry niezaleznie od wzrostu modelu.

Uklad wspolrzednych: +X = prawa strona osoby, +Y = przod, +Z = gora. Metry.
Wszystko liczone dla prawej strony; lewa powstaje przez odbicie X.
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

DATA = Path(__file__).resolve().parent.parent / "data"


def _unit(v: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(v)
    return v / n if n > 1e-12 else np.array([0.0, 0.0, 1.0])


@dataclass
class Segment:
    """Odcinek ciala o eliptycznym (lub superliptycznym) przekroju wzdluz osi."""

    name: str
    p0: np.ndarray                      # koniec blizszy (proksymalny)
    p1: np.ndarray                      # koniec dalszy (dystalny)
    profile: list                       # (t, rx, ry, dy) w metrach
    cun_length: float                   # dlugosc segmentu w cun
    exponent: float = 2.0               # 2 = elipsa, >2 = przekroj bardziej kanciasty
    cun_lat_m: float | None = None      # metry na 1 cun boczny (domyslnie z dlugosci)
    cun_t0: float = 0.0                 # t, w ktorym linijka cun ma 0 (np. fald pachowy)
    cun_t1: float = 1.0                 # t, w ktorym linijka cun ma cun_length
    e_ax: np.ndarray = field(init=False)
    e_lat: np.ndarray = field(init=False)
    e_ant: np.ndarray = field(init=False)

    def __post_init__(self) -> None:
        self.e_ax = _unit(self.p1 - self.p0)
        ant = np.array([0.0, 1.0, 0.0])
        ant = ant - self.e_ax * float(np.dot(ant, self.e_ax))
        if np.linalg.norm(ant) < 1e-6:            # os prawie pozioma (stopa)
            ant = np.array([0.0, 0.0, 1.0])
            ant = ant - self.e_ax * float(np.dot(ant, self.e_ax))
        self.e_ant = _unit(ant)
        self.e_lat = _unit(np.cross(self.e_ant, self.e_ax))
        if self.e_lat[0] < 0:                     # e_lat ma wskazywac na zewnatrz (w prawo)
            self.e_lat = -self.e_lat

    # --- metryka ---------------------------------------------------------
    @property
    def length(self) -> float:
        return float(np.linalg.norm(self.p1 - self.p0))

    @property
    def cun_m(self) -> float:
        """Metry na 1 cun mierzony wzdluz osi segmentu."""
        return self.length * (self.cun_t1 - self.cun_t0) / self.cun_length

    def t_from_cun(self, cun: float, frm: str = "p0") -> float:
        """t dla podanej odleglosci w cun od konca proksymalnego (`p0`) lub dystalnego (`p1`)."""
        span = self.cun_t1 - self.cun_t0
        if frm == "p1":
            return self.cun_t1 - span * (cun / self.cun_length)
        return self.cun_t0 + span * (cun / self.cun_length)

    @property
    def cun_lat(self) -> float:
        """Metry na 1 cun mierzony w poprzek (bocznie)."""
        return self.cun_lat_m if self.cun_lat_m else self.cun_m

    # --- geometria -------------------------------------------------------
    def axis_point(self, t: float) -> np.ndarray:
        t = float(np.clip(t, 0.0, 1.0))
        return self.p0 + (self.p1 - self.p0) * t

    def radii(self, t: float):
        """(rx, ry, dy) - polosie przekroju i przesuniecie srodka do przodu."""
        t = float(np.clip(t, 0.0, 1.0))
        ts = [p[0] for p in self.profile]
        rx = float(np.interp(t, ts, [p[1] for p in self.profile]))
        ry = float(np.interp(t, ts, [p[2] for p in self.profile]))
        dy = float(np.interp(t, ts, [p[3] for p in self.profile]))
        return rx, ry, dy

    def ring(self, t: float, n: int = 48, scale: float = 1.0) -> np.ndarray:
        """Zamkniety przekroj poprzeczny jako chmura punktow 3D."""
        rx, ry, dy = self.radii(t)
        rx, ry = rx * scale, ry * scale
        a = np.linspace(0.0, 2.0 * math.pi, n, endpoint=False)
        ex = 2.0 / self.exponent
        x = rx * np.sign(np.sin(a)) * np.abs(np.sin(a)) ** ex
        y = ry * np.sign(np.cos(a)) * np.abs(np.cos(a)) ** ex + dy
        c = self.axis_point(t)
        return c + np.outer(x, self.e_lat) + np.outer(y, self.e_ant)

    def surface_by_azimuth(self, t: float, az_deg: float, out_cun: float = 0.0) -> np.ndarray:
        """Punkt na powierzchni; azymut 0 = przod, +90 = bok, 180 = tyl, -90 = przysrodkowo."""
        rx, ry, dy = self.radii(t)
        a = math.radians(az_deg)
        ex = 2.0 / self.exponent
        sa, ca = math.sin(a), math.cos(a)
        x = rx * math.copysign(abs(sa) ** ex, sa)
        y = ry * math.copysign(abs(ca) ** ex, ca)
        base = self.axis_point(t) + dy * self.e_ant
        p = base + x * self.e_lat + y * self.e_ant
        if out_cun:
            p = p + _unit(p - base) * out_cun * self.cun_lat
        return p

    def surface_by_lateral(self, t: float, lat_cun: float, face: str = "anterior",
                           out_cun: float = 0.0) -> np.ndarray:
        """Punkt oddalony o `lat_cun` od linii posrodkowej segmentu, na wskazanej scianie."""
        rx, ry, dy = self.radii(t)
        x = lat_cun * self.cun_lat
        base = self.axis_point(t) + dy * self.e_ant
        if face == "lateral":
            sgn = 1.0 if x >= 0 else -1.0
            p = base + sgn * rx * self.e_lat
        else:
            xr = min(abs(x) / rx, 0.985)
            y = ry * (1.0 - xr ** self.exponent) ** (1.0 / self.exponent)
            if face == "posterior":
                y = -y
            p = base + x * self.e_lat + y * self.e_ant
        if out_cun:
            p = p + _unit(p - base) * out_cun * self.cun_lat
        return p


class Rig:
    """Zestaw segmentow + linijki cun (pionowa tulowia i poziomy kregow)."""

    def __init__(self, height_m: float | None = None, spec_path: Path | None = None):
        self.spec = json.loads((spec_path or DATA / "rig.json").read_text(encoding="utf-8"))
        self.H = float(height_m or self.spec["default_height_m"])
        self.levels = {k: v * self.H for k, v in self.spec["levels"].items()}
        self.W = {k: v * self.H for k, v in self.spec["widths"].items()}
        self.D = {k: v * self.H for k, v in self.spec["depths"].items()}
        self.R = {k: v * self.H for k, v in self.spec["limb_radii"].items()}
        self.cun_trunk_lat = self.spec["cun_trunk_lateral_ref"]["cun_per_H"] * self.H
        self.segments = {}
        self._structs = {}
        self._frames = {}
        self._skin_mesh = None
        self._build()

    # --- budowa ----------------------------------------------------------
    def _build(self) -> None:
        H, L, W, R = self.H, self.levels, self.W, self.R
        cl = self.spec["cun_lengths"]
        pose = self.spec["pose"]

        # --- tulow: pionowy loft o przekroju superliptycznym
        trunk_secs = [
            (0.440, 0.074, 0.050, -0.004), (0.470, 0.092, 0.058, -0.003),
            (0.500, 0.098, 0.061, -0.002),
            (0.530, 0.099, 0.060, -0.002), (0.600, 0.086, 0.056, 0.000),
            (0.625, 0.086, 0.055, 0.000), (0.680, 0.096, 0.062, 0.002),
            (0.720, 0.104, 0.065, 0.002), (0.780, 0.110, 0.062, 0.000),
            (0.815, 0.118, 0.056, -0.002), (0.838, 0.098, 0.050, -0.004),
        ]
        z0, z1 = 0.440 * H, 0.838 * H
        prof = [((z * H - z0) / (z1 - z0), rx * H, ry * H, dy * H) for z, rx, ry, dy in trunk_secs]
        self.segments["trunk"] = Segment(
            "trunk", np.array([0.0, 0.0, z0]), np.array([0.0, 0.0, z1]),
            prof, cun_length=cl["trunk"], exponent=2.6, cun_lat_m=self.cun_trunk_lat)

        # --- szyja i glowa
        self.segments["neck"] = Segment(
            "neck", np.array([0.0, 0.0, L["sternal_notch"]]),
            np.array([0.0, 0.0, L["chin"] + 0.004 * H]),
            [(0.0, W["neck"] * 1.2, self.D["neck"] * 1.15, 0.0),
             (1.0, W["neck"], self.D["neck"], 0.004 * H)],
            cun_length=cl["neck"], cun_lat_m=self.cun_trunk_lat)

        hz0, hz1 = L["chin"], L["vertex"]
        head_secs = [(0.00, 0.026, 0.030, 0.010), (0.18, 0.040, 0.052, 0.012),
                     (0.42, 0.045, 0.062, 0.008), (0.62, 0.046, 0.062, 0.004),
                     (0.82, 0.040, 0.052, 0.000), (1.00, 0.016, 0.020, -0.004)]
        self.segments["head"] = Segment(
            "head", np.array([0.0, 0.0, hz0]), np.array([0.0, 0.0, hz1]),
            [(t, rx * H, ry * H, dy * H) for t, rx, ry, dy in head_secs],
            cun_length=cl["head"], exponent=2.2, cun_lat_m=0.0125 * H)

        # --- konczyna gorna (prawa; lewa = odbicie)
        a1 = math.radians(pose["arm_abduction_deg"])
        a2 = math.radians(pose["forearm_abduction_deg"])
        sh = np.array([W["shoulder_joint_x"], 0.0, L["acromion"]])
        d1 = np.array([math.sin(a1), 0.0, -math.cos(a1)])
        el = sh + d1 * ((L["acromion"] - L["elbow"]) / math.cos(a1))
        d2 = _unit(np.array([math.sin(a2), 0.02, -math.cos(a2)]))
        wr = el + d2 * ((L["elbow"] - L["wrist"]) / abs(d2[2]))
        hand_tip = wr + d2 * (pose["hand_length"] * H)

        self.segments["upper_arm"] = Segment(
            "upper_arm", sh, el,
            [(0.0, R["upper_arm_prox"], R["upper_arm_prox"] * 1.05, 0.0),
             (0.25, R["upper_arm_prox"] * 0.92, R["upper_arm_prox"] * 0.95, 0.0),
             (1.0, R["upper_arm_dist"], R["upper_arm_dist"] * 1.05, 0.0)],
            cun_length=cl["upper_arm"], cun_t0=0.19)   # 0 cun = przedni fald pachowy
        self.segments["forearm"] = Segment(
            "forearm", el, wr,
            [(0.0, R["forearm_prox"], R["forearm_prox"], 0.0),
             (0.35, R["forearm_prox"] * 0.93, R["forearm_prox"] * 0.93, 0.0),
             (1.0, R["forearm_dist"], R["forearm_dist"] * 0.85, 0.0)],
            cun_length=cl["forearm"])
        self.segments["hand"] = Segment(
            "hand", wr, hand_tip,
            [(0.0, R["hand_prox"], R["hand_prox"] * 0.6, 0.0),
             (0.45, R["hand_prox"] * 1.15, R["hand_prox"] * 0.55, 0.0),
             (1.0, R["hand_dist"], R["hand_dist"] * 0.4, 0.0)],
            cun_length=cl["hand"])   # umowna skala dloni

        # --- konczyna dolna
        hip = np.array([W["hip_joint_x"], 0.0, L["trochanter"]])
        knee = np.array([W["knee_x"], 0.008 * H, L["knee"]])
        ankle = np.array([W["ankle_x"], 0.0, L["ankle"]])
        toe = ankle + np.array([0.0, pose["foot_length"] * H, -L["ankle"] * 0.55])
        self.segments["thigh"] = Segment(
            "thigh", hip, knee,
            [(0.0, R["thigh_prox"], R["thigh_prox"] * 1.02, 0.0),
             (0.3, R["thigh_prox"] * 0.95, R["thigh_prox"], 0.0),
             (1.0, R["thigh_dist"], R["thigh_dist"] * 1.05, 0.0)],
            cun_length=cl["thigh"])
        self.segments["shank"] = Segment(
            "shank", knee, ankle,
            [(0.0, R["shank_prox"], R["shank_prox"] * 1.05, 0.0),
             (0.22, R["shank_prox"] * 1.02, R["shank_prox"] * 1.1, -0.004 * H),
             (1.0, R["shank_dist"], R["shank_dist"], 0.0)],
            cun_length=cl["shank"])
        self.segments["foot"] = Segment(
            "foot", ankle, toe,
            [(0.0, R["foot_prox"] * 0.8, R["foot_prox"] * 0.8, 0.0),   # nie zachodzi na kostke
             (0.5, R["foot_prox"] * 1.05, R["foot_prox"] * 0.75, 0.0),
             (1.0, R["foot_dist"], R["foot_dist"] * 0.45, 0.0)],
            cun_length=cl["foot"])

        self._build_rulers()

    def _build_rulers(self) -> None:
        """Pionowa linijka cun tulowia + poziomy wyrostkow kolczystych."""
        L = self.levels
        knots = [("pubis", L["pubis"], 0.0), ("umbilicus", L["umbilicus"], 5.0),
                 ("xiphoid", L["xiphoid"], 13.0), ("sternal_notch", L["sternal_notch"], 22.0),
                 ("chin", L["chin"], 26.0)]
        self._ruler_z = [k[1] for k in knots]
        self._ruler_c = [k[2] for k in knots]
        self.landmark_cun = {k[0]: k[2] for k in knots}
        self.landmark_cun["nipple"] = self.cun_of_z(L["nipple"])
        self.landmark_cun["iliac_crest"] = self.cun_of_z(L["iliac_crest"])

        v = {"C7": 0.830}
        for i in range(1, 13):                       # T1..T12
            v[f"T{i}"] = 0.815 - (i - 1) * 0.013
        for i in range(1, 6):                        # L1..L5
            v[f"L{i}"] = 0.656 - (i - 1) * 0.016
        for i in range(1, 5):                        # S1..S4 (otwory krzyzowe)
            v[f"S{i}"] = 0.575 - (i - 1) * 0.016
        v["coccyx"] = 0.495
        self.vertebra_z = {k: z * self.H for k, z in v.items()}

    # --- linijki ---------------------------------------------------------
    def z_of_cun(self, cun: float) -> float:
        return float(np.interp(cun, self._ruler_c, self._ruler_z))

    def cun_of_z(self, z: float) -> float:
        return float(np.interp(z, self._ruler_z, self._ruler_c))

    def trunk_t_for_z(self, z: float) -> float:
        s = self.segments["trunk"]
        return float((z - s.p0[2]) / (s.p1[2] - s.p0[2]))

    def seg(self, name: str) -> Segment:
        return self.segments[name]

    # --- drobne struktury i kierunki -------------------------------------
    # anchors.resolve_struct wola struct(), frame(), direction() i cast_side()
    # nie wiedzac, czy pod spodem jest siatka anatomiczna, czy wzor. Ponizsze
    # cztery metody odtwarzaja ten sam interfejs co Anatomy, tylko z rigu.

    MAX_RAY = 0.6

    def use_skin(self, mesh) -> None:
        """Podstawia gotowa siatke skory, zeby nie budowac jej drugi raz."""
        self._skin_mesh = mesh

    @property
    def skin(self):
        if self._skin_mesh is None:
            from .geometry import soft_layer          # import lokalny: geometry uzywa rigu
            self._skin_mesh = soft_layer(self, 1.0).triangulate()
        return self._skin_mesh

    def struct(self, name: str):
        from . import structures_proc
        if name not in self._structs:
            self._structs[name] = structures_proc.build(self, name)
        return self._structs[name]

    def frame(self, region: str) -> dict:
        """Ramka kierunkow dloni albo stopy, wprost z wektorow segmentu.

        Anatomia wylicza ja z kosci, bo tam zwrot osi nie jest z gory znany.
        Tutaj segment ma juz swoje e_ax / e_lat / e_ant, wiec ramka jest ich
        przemianowaniem: w postawie anatomicznej kciuk jest bocznie (promieniowo),
        dlon zwrocona do przodu, a grzbiet stopy do gory.
        """
        if region in self._frames:
            return self._frames[region]
        seg = self.seg("hand" if region == "hand" else "foot")
        long, lat, ant = seg.e_ax, seg.e_lat, seg.e_ant
        if region == "hand":
            f = {"long": long, "radial": lat, "ulnar": -lat,
                 "palmar": ant, "dorsal": -ant}
        else:
            f = {"long": long, "lateral": lat, "medial": -lat,
                 "dorsal": ant, "plantar": -ant}
        f["proximal"], f["distal"] = -long, long
        self._frames[region] = f
        return f

    def direction(self, name: str, region: str = None) -> np.ndarray:
        if region and name in self.frame(region):
            return self.frame(region)[name]
        base = {"up": [0, 0, 1], "down": [0, 0, -1], "anterior": [0, 1, 0],
                "posterior": [0, -1, 0], "lateral": [1, 0, 0], "medial": [-1, 0, 0],
                "dorsal": [0, -1, 0], "palmar": [0, 1, 0]}
        if name not in base:
            raise KeyError("nie znam kierunku: " + name)
        return np.array(base[name], dtype=float)

    def cast(self, origin: np.ndarray, direction: np.ndarray,
             limit: float = None) -> np.ndarray:
        """Pierwsze trafienie w skore, idac z wnetrza na zewnatrz."""
        end = origin + direction * (self.MAX_RAY if limit is None else limit)
        pts, _ = self.skin.ray_trace(origin, end, first_point=True)
        pts = np.asarray(pts)
        if pts.size >= 3:
            return pts.reshape(-1)[:3]
        probe = origin + direction * 0.03           # promien uciekl - siadamy najblizej
        return np.asarray(self.skin.points)[self.skin.find_closest_point(probe)]

    def cast_side(self, base: np.ndarray, d: np.ndarray, region: str = None) -> np.ndarray:
        """Promien z punktu bocznego - nie moze przejsc przez linie posrodkowa."""
        limit = self.MAX_RAY
        if abs(base[0]) > 0.02 and abs(d[0]) > 1e-6 and d[0] * base[0] < 0:
            limit = float(np.clip(-base[0] / d[0], 0.01, self.MAX_RAY))
        return self.cast(base, d, limit)
