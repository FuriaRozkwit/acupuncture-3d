"""Cialo zbudowane z prawdziwych siatek anatomicznych (Z-Anatomy + BodyParts3D).

Zastepuje model proceduralny, zachowujac ten sam interfejs, ktorego uzywa
`anchors.resolve` - dzieki temu wszystkie 361 kotwic dziala bez zmiany w danych.

Roznica jest w tym, jak znajdujemy powierzchnie:
    model proceduralny  - wzor na elipsie przekroju segmentu,
    anatomia            - promien puszczony z wnetrza ciala na siatke skory.

Punkt odniesienia wewnatrz ciala bierze sie z osi segmentu (konczyny) albo
z krzywej rdzenia tulowia liczonej z trzonow kregow (tulow, szyja, glowa).
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pyvista as pv

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
READY = ROOT / "assets" / "ready"

MAX_RAY = 0.45          # metry - dalej niz to nie ma juz ciala


def _unit(v: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(v)
    return v / n if n > 1e-12 else np.array([0.0, 0.0, 1.0])


class FittedSegment:
    """Odcinek ciala oparty na kosciach, z powierzchnia szukana promieniem."""

    def __init__(self, name, p0, p1, cun_length, body, cun_t0=0.0, cun_t1=1.0,
                 cun_lat_m=None, core=False):
        self.name = name
        self.p0 = np.asarray(p0, dtype=float)
        self.p1 = np.asarray(p1, dtype=float)
        self.cun_length = float(cun_length)
        self.cun_t0 = float(cun_t0)
        self.cun_t1 = float(cun_t1)
        self.cun_lat_m = cun_lat_m
        self.body = body
        self.core = core                    # czy wnetrze brac z krzywej rdzenia
        self.e_ax = _unit(self.p1 - self.p0)
        ant = np.array([0.0, 1.0, 0.0])
        ant = ant - self.e_ax * float(np.dot(ant, self.e_ax))
        if np.linalg.norm(ant) < 1e-6:
            ant = np.array([0.0, 0.0, 1.0])
            ant = ant - self.e_ax * float(np.dot(ant, self.e_ax))
        self.e_ant = _unit(ant)
        self.e_lat = _unit(np.cross(self.e_ant, self.e_ax))
        if self.e_lat[0] < 0:
            self.e_lat = -self.e_lat

    # --- metryka ---------------------------------------------------------
    @property
    def length(self) -> float:
        return float(np.linalg.norm(self.p1 - self.p0))

    @property
    def cun_m(self) -> float:
        return self.length * (self.cun_t1 - self.cun_t0) / self.cun_length

    @property
    def cun_lat(self) -> float:
        return self.cun_lat_m if self.cun_lat_m else self.cun_m

    def t_from_cun(self, cun: float, frm: str = "p0") -> float:
        span = self.cun_t1 - self.cun_t0
        if frm == "p1":
            return self.cun_t1 - span * (cun / self.cun_length)
        return self.cun_t0 + span * (cun / self.cun_length)

    def axis_point(self, t: float) -> np.ndarray:
        t = float(np.clip(t, 0.0, 1.0))
        return self.p0 + (self.p1 - self.p0) * t

    def interior(self, t: float) -> np.ndarray:
        """Punkt na pewno wewnatrz ciala, z ktorego strzelamy promieniem."""
        p = self.axis_point(t)
        if self.core:
            p = p.copy()
            p[1] = self.body.core_y(float(p[2]))
        return p

    # --- powierzchnia ----------------------------------------------------
    def _limit(self, base: np.ndarray, d: np.ndarray) -> float:
        """Ile promien moze przelecieć, zanim przetnie linie posrodkowa.

        Bez tego punkt na przysrodkowej stronie uda potrafil przelecieć przez
        krocze i wyladowac na przeciwnej nodze - promien startuje wtedy z osi,
        ktora u gory uda biegnie jeszcze wewnatrz miednicy.
        """
        if abs(self.p0[0]) < 0.02 or abs(d[0]) < 1e-6:
            return MAX_RAY                       # segment posrodkowy albo promien pionowy
        if d[0] * base[0] >= 0:
            return MAX_RAY                       # promien idzie na zewnatrz ciala
        return float(np.clip(-base[0] / d[0], 0.01, MAX_RAY))

    def surface_by_azimuth(self, t: float, az_deg: float, out_cun: float = 0.0) -> np.ndarray:
        a = math.radians(az_deg)
        d = _unit(math.cos(a) * self.e_ant + math.sin(a) * self.e_lat)
        base = self.interior(t)
        p = self.body.cast(base, d, self._limit(base, d))
        return p + d * out_cun * self.cun_lat

    def surface_by_lateral(self, t: float, lat_cun: float, face: str = "anterior",
                           out_cun: float = 0.0) -> np.ndarray:
        base = self.interior(t) + self.e_lat * (lat_cun * self.cun_lat)
        if face == "posterior":
            d = -self.e_ant
        elif face == "lateral":
            d = self.e_lat if lat_cun >= 0 else -self.e_lat
            base = self.interior(t)
        else:
            d = self.e_ant
        p = self.body.cast(base, d, self._limit(base, d))
        return p + d * out_cun * self.cun_lat


class Anatomy:
    """Cialo anatomiczne: skora, szkielet, miesnie + linijka cun z kosci."""

    def __init__(self, fit_path: Path | None = None):
        self.fit = json.loads((fit_path or DATA / "rig_fitted.json").read_text(encoding="utf-8"))
        self.H = float(self.fit["height"])
        self.skin = pv.read(READY / "skin.vtp").triangulate()
        # BodyParts3D zapisuje skore z normalnymi do wewnatrz - bez tego jest czarna
        self.skin = self.skin.compute_normals(auto_orient_normals=True,
                                              consistent_normals=True, flip_normals=True)
        self.skeleton = pv.read(READY / "skeleton.vtp")
        self.cun_trunk_lat = float(self.fit["cun"]["trunk_lateral_m"])
        self.cun_trunk_vert = float(self.fit["cun"]["trunk_vertical_m"])
        self._core = self._build_core_curve()
        self._ruler_z = [k[0] for k in self.fit["trunk_ruler"]["knots"]]
        self._ruler_c = [k[1] for k in self.fit["trunk_ruler"]["knots"]]
        self.vertebra_z = {k: v["z"] for k, v in self.fit["vertebra"].items()}
        self.landmarks = {k: np.array(v) for k, v in self.fit["landmarks"].items()}
        self.levels = self._levels()
        self.segments = self._build_segments()
        self.landmark_cun = {
            "pubis": 0.0, "umbilicus": 5.0, "xiphoid": 13.0, "sternal_notch": 22.0,
            "nipple": self.cun_of_z(self.levels["nipple"]),
            "iliac_crest": self.cun_of_z(self.levels["iliac_crest"]),
        }

    # --- geometria pomocnicza -------------------------------------------
    def _build_core_curve(self):
        """Dla kazdej wysokosci: srodek ciala w osi przod-tyl (wnetrze tulowia)."""
        p = np.asarray(self.skin.points)
        band = p[np.abs(p[:, 0]) < 0.035]
        zs = np.linspace(band[:, 2].min(), band[:, 2].max(), 160)
        ys = []
        for z in zs:
            s = band[np.abs(band[:, 2] - z) < 0.012]
            ys.append(float((s[:, 1].min() + s[:, 1].max()) / 2.0) if len(s) > 6
                      else (ys[-1] if ys else 0.0))
        ys = np.convolve(np.array(ys), np.ones(5) / 5.0, mode="same")
        ys[:2], ys[-2:] = ys[2], ys[-3]
        return zs, ys

    def core_y(self, z: float) -> float:
        zs, ys = self._core
        return float(np.interp(z, zs, ys))

    def cast(self, origin: np.ndarray, direction: np.ndarray,
             limit: float = None) -> np.ndarray:
        """Pierwsze trafienie w skore, idac z wnetrza na zewnatrz."""
        end = origin + direction * (MAX_RAY if limit is None else limit)
        pts, _ = self.skin.ray_trace(origin, end, first_point=True)
        pts = np.asarray(pts)
        if pts.size >= 3:
            return pts.reshape(-1)[:3]
        # Promien uciekl - albo przez otwor w siatce (oczodol, usta), albo w szczeline
        # miedzy palcami. Zamiast zostawiac punkt w powietrzu, siadamy na najblizszym
        # kawalku skory w tym kierunku.
        probe = origin + direction * 0.03
        j = self.skin.find_closest_point(probe)
        return np.asarray(self.skin.points)[j]

    # --- linijki ---------------------------------------------------------
    def z_of_cun(self, cun: float) -> float:
        return float(np.interp(cun, self._ruler_c, self._ruler_z))

    def cun_of_z(self, z: float) -> float:
        return float(np.interp(z, self._ruler_z, self._ruler_c))

    def trunk_t_for_z(self, z: float) -> float:
        s = self.segments["trunk"]
        return float((z - s.p0[2]) / (s.p1[2] - s.p0[2]))

    def seg(self, name: str):
        return self.segments[name]

    def _levels(self) -> dict:
        lm = self.landmarks
        lv = {k: float(v[2]) for k, v in lm.items()}
        lv["vertex"] = float(lm["vertex"][2])
        lv["chin"] = float(lm["chin"][2]) if "chin" in lm else lv["vertex"] - 0.13 * self.H
        lv["nipple"] = self.z_of_cun(18.0) if False else float(
            lv["xiphoid"] + (lv["sternal_notch"] - lv["xiphoid"]) * 0.42)
        lv["iliac_crest"] = float(self.vertebra_z.get("L4", lv["pubis"] + 0.16))
        lv["hairline_ant"] = lv["chin"] + 0.66 * (lv["vertex"] - lv["chin"])
        return lv

    def _build_segments(self) -> dict:
        f = self.fit["segments"]
        segs = {}
        for name in ("upper_arm", "forearm", "hand", "thigh", "shank", "foot"):
            s = f[name]
            segs[name] = FittedSegment(name, s["p0"], s["p1"], s["cun_length"], self,
                                       cun_t0=s.get("cun_t0", 0.0))
        lv = self.levels
        z_lo = lv["pubis"] - 0.05
        z_hi = lv["sternal_notch"] + 0.03
        segs["trunk"] = FittedSegment(
            "trunk", [0.0, 0.0, z_lo], [0.0, 0.0, z_hi], 22.0, self,
            cun_lat_m=self.cun_trunk_lat, core=True)
        segs["neck"] = FittedSegment(
            "neck", [0.0, 0.0, lv["sternal_notch"]], [0.0, 0.0, lv["chin"] + 0.004],
            4.0, self, cun_lat_m=self.cun_trunk_lat, core=True)
        segs["head"] = FittedSegment(
            "head", [0.0, 0.0, lv["chin"]], [0.0, 0.0, lv["vertex"]],
            10.0, self, cun_lat_m=0.0125 * self.H, core=True)
        return segs

    # --- miesnie ---------------------------------------------------------
    def muscle_catalog(self) -> list:
        return json.loads((DATA / "muscles.json").read_text(encoding="utf-8"))["muscles"]
