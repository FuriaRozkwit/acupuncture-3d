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
import re
from pathlib import Path

import numpy as np
import pyvista as pv

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
READY = ROOT / "assets" / "ready"

MAX_RAY = 0.45          # metry - dalej niz to nie ma juz ciala
DOUBLE_SKIN = 0.010     # podwojny plat skory (faldy dloni) nie jest grubszy niz tyle


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
        shift = getattr(self.body, "section_shift", None)
        return p + shift(self.name, t) if shift else p

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


class StructRef:
    """Pojedyncza struktura anatomiczna (kosc, paliczek, skrawek malzowiny).

    Sluzy jako punkt odniesienia dla punktow, ktorych podrecznik nie opisuje w cun,
    tylko przez anatomie: "kacik paznokcia palca wskazujacego", "miedzy I a II koscia
    srodrecza", "przed skrawkiem ucha".
    """

    def __init__(self, name: str, mesh: pv.PolyData, meta: dict):
        self.name = name
        self.mesh = mesh
        self.points = np.asarray(mesh.points)
        self.centroid = np.array(meta["centroid"])
        self.axis = np.array(meta["axis"])
        self.extent = tuple(meta["extent"])

    def at(self, along: float, ref_dir: np.ndarray = None) -> np.ndarray:
        """Punkt na osi struktury: 0 = koniec blizszy, 1 = koniec dalszy.

        Zwrot osi z rozkladu PCA jest przypadkowy, wiec gdy podamy kierunek
        odniesienia (np. dlugosc dloni), os jest do niego dopasowywana - inaczej
        "koniec dalszy" paliczka raz oznaczalby czubek, a raz nasade.
        """
        axis, (lo, hi) = self.axis, self.extent
        if ref_dir is not None and float(np.dot(axis, ref_dir)) < 0:
            axis, lo, hi = -axis, -hi, -lo
        return self.centroid + axis * (lo + (hi - lo) * float(np.clip(along, 0.0, 1.0)))

    def extreme(self, direction: np.ndarray) -> np.ndarray:
        """Najdalej wysuniety punkt struktury we wskazanym kierunku."""
        d = _unit(np.asarray(direction, dtype=float))
        return self.points[int(np.argmax(self.points @ d))]


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
        self._structs_meta_early()
        # Malzowina sluzy tylko do RYSOWANIA. Promienie szukajace skory lecą na
        # powierzchnie bez uszu, bo zadnego z 361 punktow nie klada na malzowinie -
        # wszystkie leza na skorze glowy wokol niej, a cienki platek ucha
        # przechwytywalby promienie i przyciagal punkty na siebie.
        self.skin_cast = self.skin
        self.skin = self._with_auricles(self.skin)
        self.cun_trunk_lat = float(self.fit["cun"]["trunk_lateral_m"])
        self.cun_trunk_vert = float(self.fit["cun"]["trunk_vertical_m"])
        self._core = self._build_core_curve()
        self._ruler_z = [k[0] for k in self.fit["trunk_ruler"]["knots"]]
        self._ruler_c = [k[1] for k in self.fit["trunk_ruler"]["knots"]]
        self.vertebra_z = {k: v["z"] for k, v in self.fit["vertebra"].items()}
        self.landmarks = {k: np.array(v) for k, v in self.fit["landmarks"].items()}
        self.sagittal = self._build_sagittal()
        self.levels = self._levels()
        self.segments = self._build_segments()
        self._soft = self._soft_fit()
        self.landmark_cun = {
            "pubis": 0.0, "umbilicus": 5.0, "xiphoid": 13.0, "sternal_notch": 22.0,
            "nipple": self.cun_of_z(self.levels["nipple"]),
            "iliac_crest": self.cun_of_z(self.levels["iliac_crest"]),
        }

    def _structs_meta_early(self) -> None:
        sp = DATA / "structures.json"
        self._structs_meta = (json.loads(sp.read_text(encoding="utf-8"))["structures"]
                              if sp.exists() else {})
        self._struct_cache: dict = {}
        self._frames: dict = {}

    AURICLE = ("Helix", "Antihelix", "Concha of auricle", "Tragus", "Antitragus",
               "Lobule of auricle", "Apex of auricle", "Anterior notch of auricle")

    def _with_auricles(self, skin: pv.PolyData) -> pv.PolyData:
        """Dokleja malzowiny do powierzchni ciala.

        Skora z BodyParts3D ma w miejscu ucha tylko zawiazek (siega 7.8 cm od osi,
        gdy obraczka prawdziwego ucha siega 9.3 cm), wiec punkty uszne rzutowaly sie
        na plaski bok glowy. Malzowine bierzemy z Z-Anatomy i odbijamy na obie strony.
        """
        parts = [skin]
        for name in self.AURICLE:
            meta = self._structs_meta.get(self.slug(name + ".r"))
            if meta is None:
                continue
            m = pv.read(READY / meta["file"]).triangulate()
            parts.append(m)
            left = m.copy()
            left.points[:, 0] *= -1.0
            parts.append(left.flip_faces())
        if len(parts) == 1:
            return skin
        merged = parts[0].merge(parts[1:], merge_points=False)
        return merged.compute_normals(auto_orient_normals=False, consistent_normals=False)

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
             limit: float = None, probe: float = 0.03) -> np.ndarray:
        """Pierwsze trafienie w skore, idac z wnetrza na zewnatrz."""
        end = origin + direction * (MAX_RAY if limit is None else limit)
        pts, _ = self.skin_cast.ray_trace(origin, end, first_point=False)
        pts = np.asarray(pts).reshape(-1, 3)
        if len(pts):
            # Siatka skory w dloni jest miejscami podwojna (faldy palcow wchodza w dlon).
            # Pierwsze przeciecie bywa wtedy wewnetrznym platem, a prawdziwa skora lezy
            # kilka mm dalej - bierzemy pierwsze trafienie, za ktorym jest powietrze.
            # Tylko w obrebie kilku milimetrow: tam, gdzie reka przylega do tulowia,
            # siatki sie przenikaja i za skora ramienia tez jest "wnetrze" - skok do
            # dalszego trafienia przerzucalby punkt z reki na tulow.
            pts = pts[np.argsort(np.linalg.norm(pts - origin, axis=1))]
            near = pts[np.linalg.norm(pts - pts[0], axis=1) <= DOUBLE_SKIN]
            for hit in near:
                if not self.is_inside(hit + direction * 0.002):
                    return hit
            return pts[0]
        # Promien uciekl - albo przez otwor w siatce (oczodol, usta), albo w szczeline
        # miedzy palcami. Zamiast zostawiac punkt w powietrzu, siadamy na najblizszym
        # kawalku skory w tym kierunku.
        j = self.skin_cast.find_closest_point(origin + direction * probe)
        return np.asarray(self.skin_cast.points)[j]

    # --- linijki ---------------------------------------------------------
    # --- dopasowanie tkanek Z-Anatomy do skory BodyParts3D na konczynach -----
    LIMBS = ("upper_arm", "forearm", "thigh", "shank")
    # Przekroj przy biodrze lapie miednice i druga noge (promien 17-19 cm), a przy
    # kostce stope - takie przekroje odrzucamy po promieniu i po odlegloscii srodka
    # od osi. Tam, gdzie nie ma pomiaru: przy KONCU BLIZSZYM poprawka wygasa do zera
    # (nie zgadujemy przy pachwinie i pasze), przy DALSZYM trzyma ostatnia zmierzona.
    SECTION_MAX_R = {"upper_arm": 0.085, "forearm": 0.085, "thigh": 0.115, "shank": 0.080}
    SECTION_MAX_OFF = 0.035
    SECTION_T = np.linspace(0.0, 1.0, 41)

    def _soft_fit(self) -> dict:
        """Poprawki na to, ze kosci i miesnie (Z-Anatomy) nie siedza w skorze (BodyParts3D).

        Oba zbiory dopasowano na tulowiu (kosc krzyzowa, zuchwa), a konczyny maja
        w nich troche inna poze. Zmierzone: tkanki przedramienia i podudzia leza
        12-23 mm za srodkiem przekroju skory, kosci srodrecza wystaja 8-10 mm poza
        skore grzbietu dloni. Promien z osi kosci o azymucie 60-90 st. trafial
        przez to w grzbiet przedramienia (LU), a punkty liczone od kosci dloni
        ladowaly poza reka.

        Dwie poprawki: (1) na konczynach promienie wychodza ze srodka przekroju
        SKORY - azymut znaczy wtedy to, co w podreczniku; (2) punkty liczone od
        kosci dloni i stopy przesuwamy tak, zeby grzbiet kosci lezal 4 mm pod
        skora grzbietu. Wynik liczy sie raz i siedzi w assets/ready.
        """
        import hashlib
        h = hashlib.md5((DATA / "rig_fitted.json").read_bytes())
        h.update(str((READY / "skin.vtp").stat().st_size).encode())
        h.update(b"soft-fit-4")
        path = READY / f"soft_fit_{h.hexdigest()[:12]}.npz"
        if path.exists():
            with np.load(path) as z:
                return {k: z[k] for k in z.files}
        out = {}
        for name in self.LIMBS:
            out[name] = self._section_offsets(self.segments[name])
        out["hand"] = self._bone_shift("hand")
        try:
            for old in READY.glob("soft_fit_*.npz"):
                old.unlink()
            np.savez(path, **out)
        except Exception:
            pass
        return out

    def _skin_loop(self, center: np.ndarray, normal: np.ndarray):
        """Obrys skory w plaszczyznie przekroju - ten kawalek, ktory otacza `center`."""
        sl = self.skin_cast.slice(normal=normal, origin=center)
        if sl.n_points == 0:
            return None
        conn = sl.connectivity()
        best, bd = None, 1e9
        for rid in np.unique(conn["RegionId"]):
            P = np.asarray(conn.extract_points(conn["RegionId"] == rid).points)
            if len(P) < 10:
                continue
            d = np.linalg.norm(P.mean(axis=0) - center)
            if d < bd:
                best, bd = P, d
        return best

    def _section_offsets(self, seg) -> np.ndarray:
        """Wektor od osi kosci do srodka przekroju skory, dla kazdego t."""
        off = np.full((len(self.SECTION_T), 3), np.nan)
        for i, t in enumerate(self.SECTION_T):
            c = seg.axis_point(float(t))
            loop = self._skin_loop(c, seg.e_ax)
            if loop is None:
                continue
            sc = loop.mean(axis=0)
            radius = np.max(np.linalg.norm(loop - sc, axis=1))
            d = sc - c
            d -= seg.e_ax * float(np.dot(d, seg.e_ax))
            if radius < self.SECTION_MAX_R.get(seg.name, 0.085) \
                    and np.linalg.norm(d) < self.SECTION_MAX_OFF:
                off[i] = d
        ok = ~np.isnan(off[:, 0])
        if not ok.any():
            return np.zeros_like(off)
        idx = np.arange(len(off))
        first, last = int(idx[ok][0]), int(idx[ok][-1])
        for k in range(3):
            off[:, k] = np.interp(idx, idx[ok], off[ok, k])
        off[first:] = off[first:]
        off[:first] *= np.linspace(0.0, 1.0, first + 1)[:-1, None] if first else 1.0
        kernel = np.ones(5) / 5.0
        pad = np.pad(off, ((2, 2), (0, 0)), mode="edge")
        out = np.column_stack([np.convolve(pad[:, k], kernel, mode="valid") for k in range(3)])
        out[last:] = out[last]
        return out

    def section_shift(self, name: str, t: float) -> np.ndarray:
        soft = getattr(self, "_soft", None)
        if not soft or name not in soft:
            return np.zeros(3)
        off = soft[name]
        return np.array([np.interp(t, self.SECTION_T, off[:, k]) for k in range(3)])

    def _bone_shift(self, region: str) -> np.ndarray:
        """Przesuniecie kosci dloni/stopy: grzbiet kosci 4 mm pod skora grzbietu."""
        f = self.frame(region)
        bone = "metacarpal" if region == "hand" else "metatarsal"
        refs = [self.struct(f"{o} {bone} bone") for o in ("Second", "Third", "Fourth")]
        pts = np.vstack([r.points for r in refs])
        shifts = []
        for along in (0.3, 0.5, 0.7):
            base = np.mean([r.at(along, f["long"]) for r in refs], axis=0)
            loop = self._skin_loop(base, f["long"])
            if loop is None:
                continue
            s = pts[np.abs((pts - base) @ f["long"]) < 0.006]
            if len(s) < 5:
                continue
            skin_dorsal = float(np.max(loop @ f["dorsal"]))
            bone_dorsal = float(np.max(s @ f["dorsal"]))
            shifts.append(skin_dorsal - 0.004 - bone_dorsal)
        if not shifts:
            return np.zeros(3)
        # poprawiamy tylko kosci, ktore WYSTAJA poza skore; glebiej lezace zostaja
        return f["dorsal"] * min(0.0, float(np.median(shifts)))

    def soft_shift(self, region: str) -> np.ndarray:
        soft = getattr(self, "_soft", None)
        if not soft or region not in soft:
            return np.zeros(3)
        return np.asarray(soft[region])

    SAGITTAL_CUN = 18.0          # gladzizna -> C7: 3 + 12 + 3

    def _build_sagittal(self):
        """Linijka glowy po skorze w linii posrodkowej: od gladzizny do wyrostka C7.

        Model nie ma wlosow, wiec linii wlosow nie widac. Norma ma na to regule:
        gdy linie wlosow sa niewyrazne, od gladzizny do wyrostka kolczystego C7
        liczy sie 18 cun (3 do przedniej linii wlosow, 12 miedzy liniami, 3 od
        tylnej do C7). Oba konce to kosci, wiec linijka nie zgaduje niczego.
        Wczesniej przednia linia wlosow byla proporcja wysokosci glowy, a tylna
        wypadala na wysokosci guzowatosci potylicznej - o 9 cm za wysoko.
        """
        if "glabella" not in self.landmarks:
            return None
        g = self.landmarks["glabella"]
        c7 = np.array(self.fit["vertebra"]["C7"]["spinous"], dtype=float)
        center = np.array([0.0, (g[1] + c7[1]) / 2.0, g[2]])
        a0 = math.atan2(g[2] - center[2], g[1] - center[1])
        a1 = math.atan2(c7[2] - center[2], c7[1] - center[1])
        if a1 <= a0:
            a1 += 2.0 * math.pi                  # przez czubek glowy, nie przez twarz
        tree = self._obb(self.skin_cast)
        # koniec linijki to skora TUZ ZA wyrostkiem C7. Promien ze srodka glowy w strone
        # C7 idzie stromo w dol i trafia w skore dopiero na plecach, 5 cm nizej
        c7_skin = self._first_hit(tree, c7, c7 + np.array([0.0, -MAX_RAY, 0.0]))
        if c7_skin is None:
            c7_skin = c7
        pts = []
        for a in np.linspace(a0, a1 + 0.2, 480):
            d = np.array([0.0, math.cos(a), math.sin(a)])
            hit = self._first_hit(tree, center, center + d * MAX_RAY)
            if hit is None:
                continue
            if a > a0 + math.pi and hit[2] <= c7_skin[2]:
                break                                # doszlismy do poziomu C7 na karku
            pts.append(hit)
        pts = np.array(pts + [c7_skin])
        s = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))])
        # Cun po luku NIE jest rowny: jednolite 18 cun kladlo czubek (linia szczytow
        # malzowin, DU-20) na 6.2 cun zamiast 8. Linijka jest wiec odcinkowa jak na
        # tulowiu: gladzizna 0, linia szczytow uszu 8 (3 + 5), C7 18.
        knots_s, knots_c = [0.0, s[-1]], [0.0, self.SAGITTAL_CUN]
        apex = self._auricle_apex()
        if apex is not None:
            over = np.flatnonzero(pts[:, 2] > np.percentile(pts[:, 2], 70))
            k = int(over[np.argmin(np.abs(pts[over, 1] - apex[1]))])
            knots_s, knots_c = [0.0, s[k], s[-1]], [0.0, 8.0, self.SAGITTAL_CUN]
        cun = np.interp(s, knots_s, knots_c)
        return {"points": pts, "cun": cun, "knots_m": knots_s, "knots_cun": knots_c,
                "cun_m": s[-1] / self.SAGITTAL_CUN}

    def _auricle_apex(self):
        meta = self._structs_meta.get(self.slug("Apex of auricle.r"))
        if meta is None:
            return None
        pts = np.asarray(pv.read(READY / meta["file"]).points)
        return pts[np.argmax(pts[:, 2])]

    @staticmethod
    def _obb(mesh):
        import vtk
        tree = vtk.vtkOBBTree()
        tree.SetDataSet(mesh)
        tree.BuildLocator()
        return tree

    @staticmethod
    def _first_hit(tree, a, b):
        import vtk
        hits = vtk.vtkPoints()
        tree.IntersectWithLine(list(map(float, a)), list(map(float, b)), hits, None)
        if hits.GetNumberOfPoints() == 0:
            return None
        return np.array(hits.GetPoint(0))

    def _raw_normal(self, p: np.ndarray):
        """Normalna trojkata skory (zwrot dowolny) i najblizszy punkt powierzchni."""
        cell, foot = self.skin_cast.find_closest_cell(np.atleast_2d(p), return_closest_point=True)
        f = np.asarray(self.skin_cast.faces).reshape(-1, 4)[int(np.atleast_1d(cell)[0]), 1:]
        v = np.asarray(self.skin_cast.points)[f]
        n = np.cross(v[1] - v[0], v[2] - v[0])
        return n / max(float(np.linalg.norm(n)), 1e-12), np.asarray(foot).reshape(3)

    def surface_normal(self, p: np.ndarray) -> np.ndarray:
        """Normalna skory w najblizszym punkcie powierzchni, zwrocona na zewnatrz."""
        n, foot = self._raw_normal(p)
        return -n if self.is_inside(foot + n * 0.003) else n

    def walk_on_skin(self, p: np.ndarray, direction: np.ndarray, distance: float,
                     step: float = 0.002) -> np.ndarray:
        """Przejscie po powierzchni skory `distance` metrow w zadanym kierunku.

        Odstep boczny na sklepieniu czaszki (BL 1.5 cun, GB 2.25) norma mierzy PO
        GLOWIE od linii posrodkowej, a nie w poprzek przekroju - przy czubku glowy
        przekroj jest maly i przeliczenie przez niego rozrzucalo punkty.
        """
        if abs(distance) < 1e-9:
            return np.asarray(p, dtype=float)
        pos = np.asarray(p, dtype=float)
        d = np.asarray(direction, dtype=float) * np.sign(distance)
        left = abs(distance)
        while left > 1e-6:
            h = min(step, left)
            n, _ = self._raw_normal(pos)             # do rzutu stycznego zwrot nie ma znaczenia
            d = d - n * float(np.dot(d, n))
            d /= max(float(np.linalg.norm(d)), 1e-12)
            nxt = pos + d * h
            _, foot = self.skin_cast.find_closest_cell(np.atleast_2d(nxt),
                                                       return_closest_point=True)
            foot = np.asarray(foot).reshape(3)
            d = foot - pos
            d /= max(float(np.linalg.norm(d)), 1e-12)
            pos = foot
            left -= h
        return pos

    def scalp_side(self, p: np.ndarray, forward: np.ndarray) -> np.ndarray:
        """Kierunek w bok od linii posrodkowej czaszki (w strone prawej strony)."""
        n, _ = self._raw_normal(p)
        side = np.cross(forward, n)
        side -= n * float(np.dot(side, n))
        side /= max(float(np.linalg.norm(side)), 1e-12)
        return side if side[0] > 0 else -side

    # krawedzie kosci, do ktorych odwoluja sie opisy punktow: (segment, struktura, strona)
    BONE_EDGES = {
        "tibial_crest": ("shank", "Tibia", "anterior"),
        "tibia_medial": ("shank", "Tibia", "medial"),
        "fibula_anterior": ("shank", "Fibula", "anterior"),
        "fibula_posterior": ("shank", "Fibula", "posterior"),
        "fibula_posterior": ("shank", "Fibula", "posterior"),
        "ulna_posterior": ("forearm", "Ulna", "posterior"),
    }
    EDGE_LEVELS = np.linspace(0.02, 0.98, 60)

    # Ile tkanki ma nad soba dana krawedz u zywego czlowieka (mm) - na tym opiera
    # sie wpasowanie kosci w skore. Krawedzie podskorne maja kilka milimetrow.
    BONE_FIT_PROBES = {
        "shank": [("Tibia", "anterior", 4.0), ("Tibia", "medial", 4.0),
                  ("Fibula", "lateral", 7.0), ("Tibia", "posterior", 25.0)],
        "forearm": [("Ulna", "posterior", 5.0), ("Radius", "lateral", 4.0)],
    }

    def bone_shift(self, seg_name: str) -> np.ndarray:
        """Przesuniecie, ktore wsadza kosci tego segmentu z powrotem w skore.

        Skora jest z BodyParts3D, kosci z Z-Anatomy - to dwie rozne osoby, a
        dopasowanie liczono na tulowiu. Na podudziu nad grzebieniem piszczeli
        wychodzilo 24 mm "skory" zamiast 3-5, a nad powierzchnia przysrodkowa
        28 mm. Jedno przesuniecie (30 mm do przodu) godzi wszystkie cztery sondy
        z bledem ok. 3 mm - czyli kosci sa przesuniete, a nie inaczej zbudowane.
        """
        cache = getattr(self, "_bone_shift_cache", None)
        if cache is None:
            cache = self._bone_shift_cache = {}
        if seg_name in cache:
            return cache[seg_name]
        probes = self.BONE_FIT_PROBES.get(seg_name)
        if not probes:
            cache[seg_name] = np.zeros(3)
            return cache[seg_name]
        seg = self.seg(seg_name)
        pts, want = [], []
        for struct, side, expect in probes:
            d = {"anterior": seg.e_ant, "posterior": -seg.e_ant,
                 "medial": -seg.e_lat, "lateral": seg.e_lat}[side]
            P = self._edge_points(seg, self.struct(struct).mesh, d)
            if len(P):
                pts.append(P)
                want.append(expect)
        if not pts:
            cache[seg_name] = np.zeros(3)
            return cache[seg_name]

        def err(shift):
            total = 0.0
            for P, expect in zip(pts, want):
                _, foot = self.skin_cast.find_closest_cell(P + shift, return_closest_point=True)
                got = float(np.median(np.linalg.norm(np.atleast_2d(foot) - (P + shift), axis=1)))
                total += (got * 1000 - expect) ** 2
            return total

        grid = np.arange(-0.030, 0.0301, 0.003)
        zero = err(np.zeros(3))
        if (zero / max(len(pts), 1)) ** 0.5 <= 5.0:
            cache[seg_name] = np.zeros(3)      # kosci i tak siedza w skorze
            return cache[seg_name]
        best, best_err = np.zeros(3), zero
        for da in grid:
            for dl in grid:
                shift = seg.e_ant * da + seg.e_lat * dl
                e = err(shift)
                if e < best_err:
                    best, best_err = shift, e
        cache[seg_name] = best
        return best

    def _edge_points(self, seg, mesh, direction, ts=None) -> np.ndarray:
        out = []
        for t in (np.linspace(0.25, 0.8, 10) if ts is None else ts):
            c = seg.interior(float(t))
            sl = mesh.slice(normal=seg.e_ax, origin=c)
            if sl.n_points < 8:
                continue
            v = np.asarray(sl.points) - c
            v -= np.outer(v @ seg.e_ax, seg.e_ax)
            out.append(np.asarray(sl.points)[int(np.argmax(v @ direction))])
        return np.array(out)

    def bone_edge_azimuth(self, ref: str, t: float) -> float:
        """Azymut krawedzi kosci na danej wysokosci segmentu.

        Pojedynczy przekroj rzadkiej siatki daje kilkanascie stopni losu, wiec
        krawedz liczymy na 60 poziomach i wygladzamy: szum sie usrednia, a skret
        kosci (na goleni ok. 20 st. na dlugosci) zostaje.
        """
        table = getattr(self, "_edges", None)
        if table is None:
            table = self._edges = {}
        if ref not in table:
            table[ref] = self._edge_line(ref)
        return float(np.interp(float(np.clip(t, 0.0, 1.0)), self.EDGE_LEVELS, table[ref]))

    PALPABLE_MM = 8.0        # dopad tkanki, przy ktorym kosc przestaje byc wyczuwalna

    def _palpable_border(self, seg, mesh, start_az: float, direction: int,
                         t: float, shift: np.ndarray) -> float:
        """Azymut, na ktorym kosc chowa sie pod miesnie - tak samo szuka jej palec.

        Brzeg przysrodkowy piszczeli to nie "najbardziej przysrodkowy punkt kosci"
        (ten obraca sie razem z trojkatnym przekrojem), tylko miejsce, gdzie plaska
        powierzchnia podskorna sie konczy.
        """
        prev = start_az
        for step in range(1, 60):
            az = start_az + direction * step * 3.0
            p = seg.surface_by_azimuth(t, az)
            sl = mesh.slice(normal=seg.e_ax, origin=seg.interior(t) - shift)
            if sl.n_points < 4:
                return prev
            d = np.min(np.linalg.norm(np.asarray(sl.points) + shift - p, axis=1)) * 1000
            if d > self.PALPABLE_MM:
                return prev
            prev = az
        return prev

    def _edge_line(self, ref: str) -> np.ndarray:
        seg_name, struct, side = self.BONE_EDGES[ref]
        seg = self.seg(seg_name)
        mesh = self.struct(struct).mesh
        pick = {"anterior": seg.e_ant, "posterior": -seg.e_ant,
                "medial": -seg.e_lat, "lateral": seg.e_lat}[side]
        shift = self.bone_shift(seg_name)          # kosci wpasowane w skore
        out = np.full(len(self.EDGE_LEVELS), np.nan)
        for i, t in enumerate(self.EDGE_LEVELS):
            c = seg.interior(float(t))
            sl = mesh.slice(normal=seg.e_ax, origin=c - shift)
            if sl.n_points < 8:
                continue
            v = np.asarray(sl.points) + shift - c
            v -= np.outer(v @ seg.e_ax, seg.e_ax)
            k = int(np.argmax(v @ pick))
            az = np.degrees(np.arctan2(float(v[k] @ seg.e_lat), float(v[k] @ seg.e_ant)))
            if side in ("medial", "lateral"):
                # brzeg boczny/przysrodkowy szukamy palpacyjnie, idac od strony,
                # po ktorej kosc jest podskorna
                az = self._palpable_border(seg, mesh, az, -1 if side == "medial" else 1,
                                           float(t), shift)
            out[i] = az
        ok = ~np.isnan(out)
        if not ok.any():
            return np.zeros_like(out)
        out = np.interp(np.arange(len(out)), np.flatnonzero(ok), out[ok])
        pad = np.pad(out, 4, mode="edge")
        return np.convolve(pad, np.ones(9) / 9.0, mode="valid")

    def azimuth_at_surface_cun(self, seg_name: str, t: float, az0: float,
                               cun: float, tol: float = 0.02) -> float:
        """Azymut oddalony o `cun` PO SKORZE od azymutu `az0` na tej samej wysokosci.

        Tak mierzy podrecznik: "jeden palec bocznie od grzebienia piszczeli" to
        odleglosc po ciele, a nie kat. Przy owalnym przekroju ten sam kat znaczy
        raz 8, raz 14 mm - dlatego szukamy kata, ktory daje zadana odleglosc.
        """
        seg = self.seg(seg_name)
        target = abs(cun) * seg.cun_m
        step = 1.0 if cun >= 0 else -1.0
        lo, hi = 0.0, 5.0
        for _ in range(40):                       # zakres: ile stopni starcza
            if self._arc_len(seg, t, az0, az0 + step * hi) >= target or hi > 180.0:
                break
            hi *= 1.6
        for _ in range(30):
            mid = 0.5 * (lo + hi)
            if self._arc_len(seg, t, az0, az0 + step * mid) < target:
                lo = mid
            else:
                hi = mid
            if hi - lo < tol:
                break
        return az0 + step * 0.5 * (lo + hi)

    @staticmethod
    def _arc_len(seg, t: float, az_a: float, az_b: float, n: int = 16) -> float:
        pts = np.array([seg.surface_by_azimuth(t, float(x))
                        for x in np.linspace(az_a, az_b, n)])
        return float(np.linalg.norm(np.diff(pts, axis=0), axis=1).sum())

    def sagittal_point(self, cun_from_glabella: float) -> np.ndarray:
        sg = self.sagittal
        c = float(np.clip(cun_from_glabella, 0.0, self.SAGITTAL_CUN))
        return np.array([np.interp(c, sg["cun"], sg["points"][:, i]) for i in range(3)])

    def ruler_t(self, name: str, cun: float, frm: str = "p0") -> float:
        """Polozenie na osi segmentu wg osobnej linijki (np. przysrodkowa podudzia).

        Linijka ma wlasne konce kostne, ale korzysta z osi i azymutow segmentu -
        zmienia sie tylko to, gdzie wzdluz osi wypada dana liczba cun.
        """
        r = self.fit["rulers"][name]
        seg = self.seg(r["seg"])
        t0 = float(np.dot(np.asarray(r["p0"]) - seg.p0, seg.e_ax)) / seg.length
        t1 = float(np.dot(np.asarray(r["p1"]) - seg.p0, seg.e_ax)) / seg.length
        f = cun / float(r["cun_length"])
        return t1 + (t0 - t1) * f if frm == "p1" else t0 + (t1 - t0) * f

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
        if self.sagittal is not None:
            lv["hairline_ant"] = float(self.sagittal_point(3.0)[2])
            lv["hairline_post"] = float(self.sagittal_point(15.0)[2])
        else:                                   # stary rig bez gladzizny - proporcja
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

    # Na dloni i stopie skora lezy kilka milimetrow od kosci. Kosci (Z-Anatomy)
    # i skora (BodyParts3D) sa z roznych zrodel i miejscami kosc wystaje poza
    # skore - promien startuje wtedy juz NA ZEWNATRZ reki i pierwsza skora, jaka
    # trafia, jest udo (SI-2 i SI-3 lezaly 8 cm od malego palca). Krotki zasieg
    # i proba zastepcza tuz przy kosci trzymaja punkt na wlasnej czesci ciala.
    SHORT_RAY = {"hand": (0.045, 0.006), "foot": (0.065, 0.008)}

    def is_inside(self, p: np.ndarray) -> bool:
        """Czy punkt lezy w ciele: parzystosc przeciec ze skora, glosowanie 7 promieni
        (siatka ma otwory - oczy, usta - wiec jeden promien moze sie pomylic)."""
        import vtk
        if getattr(self, "_parity", None) is None:
            self._parity = self._obb(self.skin_cast)
            d = np.random.default_rng(7).normal(size=(7, 3))
            self._parity_rays = d / np.linalg.norm(d, axis=1, keepdims=True)
        hits = vtk.vtkPoints()
        v = 0
        for d in self._parity_rays:
            self._parity.IntersectWithLine(list(map(float, p)), list(map(float, p + d * 2.5)),
                                           hits, None)
            v += hits.GetNumberOfPoints() % 2
        return v >= 4

    def cast_side(self, base: np.ndarray, d: np.ndarray, region: str = None) -> np.ndarray:
        """Promien z punktu bocznego - nie moze przejsc przez linie posrodkowa."""
        limit = MAX_RAY
        if abs(base[0]) > 0.02 and abs(d[0]) > 1e-6 and d[0] * base[0] < 0:
            limit = float(np.clip(-base[0] / d[0], 0.01, MAX_RAY))
        if region in self.SHORT_RAY:
            reach, probe = self.SHORT_RAY[region]
            return self.cast(base, d, min(limit, reach), probe)
        return self.cast(base, d, limit)

    # --- drobne struktury -------------------------------------------------
    @staticmethod
    def slug(name: str) -> str:
        s = name.strip().lower().replace("(", "").replace(")", "")
        return re.sub(r"[^a-z0-9]+", "-", s).strip("-")

    def struct(self, name: str) -> StructRef:
        """Struktura po nazwie (domyslnie prawa strona). Wczytywana raz."""
        key = self.slug(name)
        if key not in self._structs_meta and self.slug(name + ".r") in self._structs_meta:
            key = self.slug(name + ".r")
        if key in self._struct_cache:
            return self._struct_cache[key]
        meta = self._structs_meta.get(key)
        if meta is None:
            raise KeyError(f"nie znam struktury: {name} ({key})")
        mesh = pv.read(READY / meta["file"])
        ref = StructRef(name, mesh, meta)
        self._struct_cache[key] = ref
        return ref

    def frame(self, region: str) -> dict:
        """Uklad kierunkow dloni albo stopy, wyliczony z kosci - nie zalozony.

        Kierunek promieniowy bierze sie z linii od V do I kosci srodrecza,
        a strone grzbietowa rozstrzyga kosc grochowata, ktora lezy po stronie
        dloniowej. Dla stopy strone grzbietowa rozstrzyga pion.
        """
        if region in self._frames:
            return self._frames[region]
        ordn = ("First", "Second", "Third", "Fourth", "Fifth")
        bone = "metacarpal" if region == "hand" else "metatarsal"
        refs = [self.struct(f"{o} {bone} bone") for o in ordn]
        base = np.mean([r.at(0.0) for r in refs], axis=0)
        tip = np.mean([r.at(1.0) for r in refs], axis=0)
        long = _unit(tip - base)
        first_to_fifth = refs[0].centroid - refs[4].centroid
        side = _unit(first_to_fifth - long * float(np.dot(first_to_fifth, long)))
        normal = _unit(np.cross(side, long))
        if region == "hand":
            pis = self.struct("Pisiform bone").centroid    # kosc grochowata = strona dloniowa
            if float(np.dot(pis - base, normal)) > 0:
                normal = -normal
            f = {"long": long, "radial": side, "ulnar": -side,
                 "dorsal": normal, "palmar": -normal}
        else:
            if normal[2] < 0:                              # grzbiet stopy patrzy w gore
                normal = -normal
            f = {"long": long, "medial": side, "lateral": -side,
                 "dorsal": normal, "plantar": -normal}
        f["proximal"], f["distal"] = -long, long
        self._frames[region] = f
        return f

    def direction(self, name: str, region: str = None) -> np.ndarray:
        """Kierunek po nazwie: lokalny (dloni/stopy) albo ogolny."""
        if region and name in self.frame(region):
            return self.frame(region)[name]
        base = {"up": [0, 0, 1], "down": [0, 0, -1], "anterior": [0, 1, 0],
                "posterior": [0, -1, 0], "lateral": [1, 0, 0], "medial": [-1, 0, 0],
                "dorsal": [0, -1, 0], "palmar": [0, 1, 0]}
        if name not in base:
            raise KeyError(f"nie znam kierunku: {name}")
        return np.array(base[name], dtype=float)

    # --- miesnie ---------------------------------------------------------
    def muscle_catalog(self) -> list:
        return json.loads((DATA / "muscles.json").read_text(encoding="utf-8"))["muscles"]
