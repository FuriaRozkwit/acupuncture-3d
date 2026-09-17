"""Weryfikacja linijki cun i przebiegow kanalow (wg VERIFY_MERIDIANS.md).

Trzy czesci, kazda z jawnym poziomem pewnosci:

  A. LINIJKA  - wektory kostne z data/cun_vectors.json (proporcje B-cun z normy).
               Kazdy wektor mierzymy na kosciach modelu i pytamy linijke, ile to cun.
               Rozstrzygaja wektory, na ktorych rig NIE byl kalibrowany.
  B. PRZEBIEGI - asercje bez zewnetrznej geometrii toru (bo normy toru nie ma):
               kolejnosc punktow na krzywej, kierunek przeplywu, przekazanie kanalu
               w cyklu, krzywa przy skorze (bez cieciw przez powietrze), nie w kosci,
               gladkosc stycznej i odchylenie Frecheta od geodezyjnej po skorze
               miedzy sasiednimi punktami (samozgodnosc - NIE podrecznik L1).
  C. PUNKTY   - skrot z tools/verify_points.py (opisy lokalne, L1). Norma L0
               (wspolrzedne WHO/GB) nie jest w repozytorium - raport to mowi.

    python tools/verify_meridians.py [--channels LU,LI] [--no-frechet]

Wynik: data/verify/verify_report.json i verify_report.csv (lokalne, nie w repo).
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import time
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pyvista as pv

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))
warnings.simplefilter("ignore")

from acu3d import anchors  # noqa: E402
from acu3d.model import Atlas, SKIN_OFFSET  # noqa: E402

DATA = ROOT / "data"
OUT = DATA / "verify"
REF = ROOT / "assets" / "mesh" / "ref"

# progi - konwencje robocze z notki (§4.5), nie wartosci znormalizowane
TOL_CUN, WARN_CUN = 0.5, 0.3
TOL_MM, WARN_MM = 5.0, 3.0                 # glowa i szyja
RULER_FAIL = 1.0                           # linijka: >1 cun FAIL, >0.5 cun WARN
AIR_MM = 10.0                              # krzywa ma lezec 4 mm nad skora
BONE_MM = 1.0                              # blizej kosci = krzywa w kosci
ANGLE_WARN, ANGLE_FAIL = 30.0, 60.0        # kat miedzy kolejnymi probkami
KINK_WARN_MM, KINK_FAIL_MM = 1.0, 2.0      # zalamanie liczy sie, gdy probka odstaje od cieciwy
HANDOFF_ADJ = 1.5                          # przekazanie "sasiednie": konce blizej niz tyle cun

ARM_YIN, ARM_YANG = ("LU", "HT", "PC"), ("LI", "SI", "SJ")
LEG_YANG, LEG_YIN = ("ST", "BL", "GB"), ("SP", "KI", "LR")

# cykl 12 kanalow (Lingshu 10). kind: adjacent = konce kanalow stykaja sie na
# powierzchni; branch = galaz od punktu `from` do poczatku nastepnego kanalu;
# internal = przejscie wewnetrzne, geometrycznie nieweryfikowalne
HANDOFF = [
    ("LU", "LI", "branch", "LU-7", "LI-1"),
    ("LI", "ST", "branch", "LI-20", "ST-1"),      # ST zaczyna sie przy LI-20, idzie do BL-1 i schodzi do ST-1
    ("ST", "SP", "branch", "ST-42", "SP-1"),
    ("SP", "HT", "internal", "SP-21", "HT-1"),
    ("HT", "SI", "adjacent", "HT-9", "SI-1"),
    ("SI", "BL", "branch", "SI-18", "BL-1"),
    ("BL", "KI", "branch", "BL-67", "KI-1"),
    ("KI", "PC", "internal", "KI-27", "PC-1"),
    ("PC", "SJ", "branch", "PC-8", "SJ-1"),
    ("SJ", "GB", "adjacent", "SJ-23", "GB-1"),
    ("GB", "LR", "branch", "GB-41", "LR-1"),
    ("LR", "LU", "internal", "LR-14", "LU-1"),
]


def status(value: float, fail: float, warn: float) -> str:
    if value is None or not np.isfinite(value):
        return "SKIP"
    return "FAIL" if value > fail else "WARN" if value > warn else "PASS"


def load_ref(name: str) -> np.ndarray:
    """Kosc z assets/mesh/ref w ukladzie acu3d (jak tools/calibrate_rig.py)."""
    p = np.asarray(pv.read(REF / f"{name}.ply").points).copy()
    p[:, 0] *= -1.0
    p[:, 1] *= -1.0
    return p


# --- A. linijka ---------------------------------------------------------------
class Ruler:
    def __init__(self, atlas: Atlas):
        self.a = atlas
        self.an = atlas.anatomy
        self.lm = {k: np.asarray(v, float) for k, v in self.an.landmarks.items()}
        self.head = self.an.seg("head")
        self.sag = self.an.sagittal
        self.cun_arc = self.sag["cun_m"]

    def _on_sagittal(self, p: np.ndarray):
        """Wspolrzedna cun i indeks najblizszej probki linijki strzalkowej glowy."""
        k = int(np.argmin(np.linalg.norm(self.sag["points"] - p, axis=1)))
        return float(self.sag["cun"][k]), k

    def _sagittal_len(self, p: np.ndarray, q: np.ndarray) -> float:
        (_, i), (_, j) = self._on_sagittal(p), self._on_sagittal(q)
        i, j = min(i, j), max(i, j)
        seg = self.sag["points"][i:j + 1]
        return float(np.linalg.norm(np.diff(seg, axis=0), axis=1).sum())

    def landmark(self, name: str) -> np.ndarray:
        if name in self.lm:
            return self.lm[name]
        if name == "midline":
            return np.zeros(3)
        if name == "scapula_medial_border":
            s = load_ref("scapula-r")
            z0, z1 = np.percentile(s[:, 2], [30, 70])
            band = s[(s[:, 2] > z0) & (s[:, 2] < z1)]
            return band[np.argmin(band[:, 0])]
        if name == "femur_medial_epicondyle":
            f = self.an.struct("Femur").mesh.points
            f = np.asarray(f)
            low = f[f[:, 2] < f[:, 2].min() + 0.12 * np.ptp(f[:, 2])]
            return low[np.argmin(low[:, 0])]
        if name == "tibia_medial_condyle_lower":
            t = np.asarray(self.an.struct("Tibia").mesh.points)
            top = t[t[:, 2] > t[:, 2].max() - 0.12 * np.ptp(t[:, 2])]
            rim = top[top[:, 0] < top[:, 0].min() + 0.003]
            return rim[np.argmin(rim[:, 2])]
        if name == "sole":
            return np.array([0.0, 0.0, float(self.an.fit["floor_z"])])
        if name == "glabella":
            fr = np.asarray(self.an.struct("Frontal bone").mesh.points)
            mid = fr[(np.abs(fr[:, 0]) < 0.006) & (fr[:, 2] < np.percentile(fr[:, 2], 25))]
            return mid[np.argmax(mid[:, 1])]
        if name == "hairline_ant":
            return self.an.sagittal_point(3.0)
        if name == "hairline_post":
            return self.an.sagittal_point(15.0)
        if name == "ear_apex_line":
            ap = np.asarray(self.an.struct("Apex of auricle").points)
            apex = ap[np.argmax(ap[:, 2])]
            P = self.sag["points"]
            over = np.flatnonzero(P[:, 2] > np.percentile(P[:, 2], 70))
            return P[over[np.argmin(np.abs(P[over, 1] - apex[1]))]]
        if name == "occipital_protuberance":
            sk = np.asarray(self.an.skeleton.points)
            occ = sk[(np.abs(sk[:, 0]) < 0.004) & (sk[:, 2] > 1.50) & (sk[:, 2] < 1.66)
                     & (sk[:, 1] < -0.03)]
            return occ[np.argmin(occ[:, 1])]
        if name == "C7":
            return np.array(self.an.fit["vertebra"]["C7"]["spinous"], dtype=float)
        if name in ("mastoid_R", "mastoid_L"):
            tb = np.asarray(self.an.struct("Temporal bone").mesh.points)
            back = tb[tb[:, 1] < np.median(tb[:, 1])]
            tip = back[np.argmin(back[:, 2])].copy()
            if name == "mastoid_L":
                tip[0] *= -1.0
            return tip
        raise KeyError(name)

    def measure(self, v: dict) -> dict:
        """Ile cun ma odcinek wg linijki modelu, oraz jak to zmierzono."""
        p, q = self.landmark(v["from"]), self.landmark(v["to"])
        reg = v["region"]
        if reg == "trunk_vertical":
            got = abs(self.an.cun_of_z(float(p[2])) - self.an.cun_of_z(float(q[2])))
            unit_m, how = self.an.cun_trunk_vert, "roznica wspolrzednych linijki odcinkowej tulowia"
            length = abs(p[2] - q[2])
        elif reg == "trunk_lateral":
            length = abs(abs(q[0]) - abs(p[0]))
            unit_m, how = self.an.cun_trunk_lat, "odleglosc od linii posrodkowej (os X)"
            got = length / unit_m
        elif reg in ("forearm", "thigh", "shank"):
            seg = self.an.seg(reg)
            axis = (seg.p1 - seg.p0) / np.linalg.norm(seg.p1 - seg.p0)
            length = abs(float(np.dot(q - p, axis)))
            unit_m, how = seg.cun_m, f"rzut na os segmentu {reg}"
            got = length / unit_m
        elif reg == "shank_medial":
            r = self.an.fit["rulers"]["shank_medial"]
            seg = self.an.seg(r["seg"])
            span = abs(float(np.dot(np.subtract(r["p1"], r["p0"]), seg.e_ax)))
            length = abs(float(np.dot(q - p, seg.e_ax)))
            unit_m, how = span / r["cun_length"], "rzut na os podudzia, linijka przysrodkowa"
            got = length / unit_m
        elif reg == "head_sagittal":
            got = abs(self._on_sagittal(p)[0] - self._on_sagittal(q)[0])
            length = self._sagittal_len(p, q)
            unit_m = self.cun_arc
            how = "po skorze w linii posrodkowej; linijka gladzizna -> C7 = 18 cun"
        elif reg == "head_lateral":
            length = self._behind_head(p, q)
            unit_m = self.head.cun_lat
            how = "geodezyjna po skorze za glowa, jednostka = cun boczny glowy (0.0125*H)"
            got = length / unit_m
        else:
            raise ValueError(reg)
        err = got - v["cun"]
        return {"id": v["id"], "region": reg, "from": v["from"], "to": v["to"],
                "norm_cun": v["cun"], "model_cun": round(got, 2), "err_cun": round(err, 2),
                "length_mm": round(length * 1000, 1), "cun_mm_model": round(unit_m * 1000, 2),
                "cun_mm_implied": round(length / v["cun"] * 1000, 2),
                "fitted": v["fitted"], "confirmed": v["confirmed"], "how": how,
                "status": status(abs(err), v.get("fail", RULER_FAIL), v.get("tol", TOL_CUN))}

    def _behind_head(self, p: np.ndarray, q: np.ndarray) -> float:
        """Geodezyjna po skorze na lokalnym wycinku (caly tulow liczylby sie minutami)."""
        skin = self.an.skin_cast
        lo = np.minimum(p, q) - np.array([0.03, 0.09, 0.05])
        hi = np.maximum(p, q) + np.array([0.03, 0.03, 0.05])
        patch = skin.clip_box([lo[0], hi[0], lo[1], hi[1], lo[2], hi[2]], invert=False)
        patch = patch.extract_surface().triangulate().clean()
        i, j = patch.find_closest_point(p), patch.find_closest_point(q)
        try:
            return float(patch.geodesic_distance(i, j))
        except Exception:
            return float(np.linalg.norm(p - q))

    def run(self) -> list:
        doc = json.loads((DATA / "cun_vectors.json").read_text(encoding="utf-8"))
        out = []
        for v in doc["vectors"]:
            try:
                out.append(self.measure(v))
            except Exception as e:                      # brak kosci = brak pomiaru, nie awaria
                out.append({"id": v["id"], "status": "SKIP", "note": f"{type(e).__name__}: {e}"})
        return out


# --- B. przebiegi ----------------------------------------------------------------
class Paths:
    def __init__(self, atlas: Atlas, frechet: bool = True):
        self.a = atlas
        self.an = atlas.anatomy
        self.frechet = frechet
        self.skin = self.an.skin_cast
        self.skin_loc = self._locator(self.skin)
        self.bone_loc = self._locator(atlas.layers["skeleton"])
        self.anchor = {p["code"]: p["anchor"] for p in atlas.points}

    @staticmethod
    def _locator(mesh):
        import vtk
        loc = vtk.vtkStaticCellLocator()
        loc.SetDataSet(mesh)
        loc.BuildLocator()
        return loc

    @staticmethod
    def _closest(loc, pts: np.ndarray):
        import vtk
        foot = np.empty_like(pts)
        dist = np.empty(len(pts))
        c = [0.0, 0.0, 0.0]
        cid, sub, d2 = vtk.reference(0), vtk.reference(0), vtk.reference(0.0)
        for k, p in enumerate(pts):
            loc.FindClosestPoint([float(x) for x in p], c, cid, sub, d2)
            foot[k] = c
            dist[k] = float(d2) ** 0.5
        return dist, foot

    def pos(self, code: str) -> np.ndarray:
        m = self.a.by_key.get(code) or self.a.by_key.get(f"{code}·R")
        return np.asarray(m["pos"], float)

    def region(self, code: str) -> str:
        a = self.anchor.get(code, {})
        if "arc_cun" in a:
            return "head"
        if "struct" in a:
            return a.get("region") or "head"
        return a.get("seg", "trunk")

    def unit(self, code: str):
        """(metry na 1 cun, czy liczyc w mm) dla rejonu punktu."""
        reg = self.region(code)
        if reg in ("head", "neck"):
            return None, True
        if reg == "trunk":
            return self.an.cun_trunk_lat, False
        return self.an.seg(reg).cun_m, False

    def _measure(self, meters: float, code: str):
        u, mm = self.unit(code)
        if mm:
            return meters * 1000, "mm", status(meters * 1000, TOL_MM, WARN_MM)
        return meters / u, "cun", status(meters / u, TOL_CUN, WARN_CUN)

    # --- asercje -------------------------------------------------------------
    def flow(self, ch: str, group: list) -> dict:
        """Kierunek przeplywu z numeracji: czy punkty ida tam, dokad plynie kanal."""
        P = np.array([self.pos(c) for c in group])
        notch = np.asarray(self.an.landmarks["sternal_notch"])
        if ch in ARM_YIN + ARM_YANG:
            key = np.linalg.norm(P - notch, axis=1)
            sign, rule = (1, "klatka -> dlon") if ch in ARM_YIN else (-1, "dlon -> glowa")
        elif ch in LEG_YANG + LEG_YIN:
            key = P[:, 2]
            sign, rule = (-1, "glowa -> stopa") if ch in LEG_YANG else (1, "stopa -> tulow")
        else:
            key = P[:, 2]
            sign, rule = 1, "w gore linii posrodkowej"
        steps = np.sign(np.diff(key)) * sign
        ends_ok = bool((key[-1] - key[0]) * sign > 0)
        return {"rule": rule, "ends_ok": ends_ok,
                "steps_forward": round(float((steps > 0).mean()), 2),
                "status": "PASS" if ends_ok else "FAIL"}

    def order(self, group: list, curve: np.ndarray, nodes: np.ndarray) -> list:
        """Kazdy punkt ma lezec przy SWOIM wezle krzywej - inaczej krzywa zawraca."""
        per = self.a.PER_SEGMENT
        bad = []
        for k, code in enumerate(group):
            node = int(nodes[min(k, len(nodes) - 1)])
            d = np.linalg.norm(curve - self.pos(code), axis=1)
            j = int(np.argmin(d))
            if abs(j - node) > per // 2 and d[j] < d[node] - 0.002:
                bad.append({"point": code, "node": node, "closest_sample": j,
                            "gap_mm": round(float(d[node] - d[j]) * 1000, 1)})
        return bad

    def segment(self, a: str, b: str, piece: np.ndarray) -> dict:
        # cieciwa przez powietrze chowa sie MIEDZY probkami (krok potrafi miec 3 cm),
        # wiec odleglosc od skory mierzymy tez w 1/4, 1/2 i 3/4 kazdego kroku
        t = np.array([0.25, 0.5, 0.75])[None, :, None]
        mids = (piece[:-1, None, :] + (piece[1:] - piece[:-1])[:, None, :] * t).reshape(-1, 3)
        mid_d, _ = self._closest(self.skin_loc, mids)
        skin_d, foot = self._closest(self.skin_loc, piece)
        skin_d = np.concatenate([skin_d, mid_d])
        bone_d, _ = self._closest(self.bone_loc, piece)
        ang, kink = kinks(piece)
        air = float(skin_d.max() * 1000)
        hole = False
        if air > AIR_MM:
            # siatka skory ma otwory (krocze, oczodoly, usta). Srodek kroku, ktory
            # wypada w otworze, lezy W CIELE - to dziura w modelu, nie cieciwa przez
            # powietrze. Prawdziwa cieciwa biegnie na zewnatrz.
            far = mids[int(np.argmax(mid_d))]
            hole = self.a.anatomy.is_inside(far)
        worst = int(np.argmax(kink)) if len(kink) else 0
        a_w = float(ang[worst]) if len(ang) else 0.0
        k_w = float(kink[worst]) if len(kink) else 0.0
        row = {"segment": f"{a}-{b}", "air_max_mm": round(air, 1), "hole": hole,
               "air_status": ("WARN" if hole else "FAIL") if air > AIR_MM else "PASS",
               # kosc i skora sa z roznych zrodel (Z-Anatomy, BodyParts3D) i miejscami
               # na siebie nachodza - bez odleglosci ze znakiem nie da sie rozstrzygnac,
               # czy krzywa jest W kosci, wiec to tylko ostrzezenie
               "bone_min_mm": round(float(bone_d.min() * 1000), 1),
               "bone_status": "WARN" if bone_d.min() * 1000 < BONE_MM else "PASS",
               "angle_max_deg": round(a_w, 0), "kink_mm": round(k_w, 1),
               "angle_status": ("FAIL" if a_w > ANGLE_FAIL and k_w > KINK_FAIL_MM else
                                "WARN" if a_w > ANGLE_WARN and k_w > KINK_WARN_MM else "PASS")}
        if self.frechet:
            row.update(self.against_geodesic(a, b, foot))
            if row["angle_status"] == "FAIL":
                # czy w tym samym miejscu zagina sie takze najkrotsza droga po skorze?
                # jesli tak, to fald ciala (pacha, pachwina), a nie usterka rysowania
                geo = self.geodesic(self.pos(a), self.pos(b))
                if geo is not None:
                    g_ang, g_kink = kinks(resample(geo, len(piece)))
                    if g_kink.size and g_kink.max() > 0.6 * row["kink_mm"]:
                        row["angle_status"] = "WARN"
                        row["fold"] = round(float(g_kink.max()), 1)
        return row

    def ruler_line(self, a: str, b: str, n: int = 40):
        """Odniesienie dla skalpu: linia wprost z linijki cun.

        Kanaly na sklepieniu czaszki biegna w STALEJ odleglosci od linii posrodkowej
        (BL 1.5 cun, GB 2.25, DU 0), a nie najkrotsza droga - geodezyjna scinalaby
        taka linie jak cieciwa. Odniesienie robimy wiec z samej linijki: kolejne
        kotwice miedzy jednym punktem a drugim.
        """
        A, B = self.anchor.get(a, {}), self.anchor.get(b, {})
        mid = self.a._anchor_line(A, B, "R", n)
        if mid is None:
            return None
        return np.vstack([self.pos(a), mid, self.pos(b)])

    def against_geodesic(self, a: str, b: str, on_skin: np.ndarray) -> dict:
        """Frechet i Hausdorff krzywej (zrzutowanej na skore) wzgledem odniesienia."""
        ruler = self.ruler_line(a, b)
        geo = ruler if ruler is not None else self.geodesic(self.pos(a), self.pos(b))
        if geo is None:
            return {"frechet": None, "frechet_status": "SKIP"}
        P, Q = resample(on_skin, 48), resample(geo, 48)
        fr = discrete_frechet(P, Q)
        hd = hausdorff(P, Q)
        length = float(np.linalg.norm(np.diff(Q, axis=0), axis=1).sum())
        val, unit, st = self._measure(fr, a)
        hval, _, _ = self._measure(hd, a)
        # Na dlugim odcinku bez punktow posrednich (kolano -> spojenie to 47 cm)
        # dwie rozne "najkrotsze drogi" po tej samej skorze rozchodza sie lokalnie
        # o kilkanascie milimetrow, choc biegna tak samo. Prog rosnie z dlugoscia.
        rel = self._measure(length * 0.05, a)[0]
        if st == "FAIL" and val <= rel:
            st = "WARN"
        return {"frechet": round(val, 2), "hausdorff": round(hval, 2), "unit": unit,
                "tol_len": round(rel, 2),
                "frechet_status": st,
                "ref": ("linia linijki cun (ta sama konstrukcja - sprawdza, czy rzut na skore "
                        "i omijanie usterek jej nie ruszyly)" if ruler is not None
                        else "geodezyjna po skorze (skracanie krzywej; samozgodnosc, nie L1)")}

    def geodesic(self, p: np.ndarray, q: np.ndarray, n: int = 48, rounds: int = 40):
        """Najkrotsza droga po skorze liczona skracaniem krzywej na trojkatach.

        Siatka skory ma wezly co 9-29 mm, wiec geodezyjna po krawedziach robi kroki
        po 3 cm i sama nie spelnia progu 0.5 cun. Tu zaczynamy od cieciwy zrzutowanej
        na skore i na kazdym kroku wygladzamy ja i rzutujemy z powrotem na najblizszy
        punkt powierzchni (nie wezel) - krzywa skraca sie do lokalnej geodezyjnej.
        """
        t = np.linspace(0.0, 1.0, n)[:, None]
        _, pts = self._closest(self.skin_loc, p + (q - p) * t)
        for _ in range(rounds):
            pts[1:-1] = 0.5 * pts[1:-1] + 0.25 * (pts[:-2] + pts[2:])
            _, pts = self._closest(self.skin_loc, pts)
        return pts

    def handoffs(self) -> list:
        out = []
        for src, dst, kind, a, b in HANDOFF:
            if kind == "internal":
                out.append({"from": src, "to": dst, "kind": kind, "status": "INFO",
                            "geom_verifiable": False})
                continue
            d = float(np.linalg.norm(self.pos(a) - self.pos(b)))
            u, mm = self.unit(a)
            if mm:
                u = self.an.cun_trunk_lat
            cun = d / u
            st = ("FAIL" if cun > HANDOFF_ADJ else "PASS") if kind == "adjacent" else "INFO"
            out.append({"from": src, "to": dst, "kind": kind, "points": f"{a}->{b}",
                        "dist_mm": round(d * 1000, 1), "dist_cun": round(cun, 2), "status": st})
        return out

    def channel(self, ch: str) -> dict:
        groups = self.a.paths.get(ch, [])
        curves = self.a.meridian_lines(ch, "R") or self.a.meridian_lines(ch, "C")
        node_sets = self.a.curve_nodes.get((ch, "R")) or self.a.curve_nodes.get((ch, "C")) or []
        rep = {"channel_id": ch, "path_checks": [], "order_errors": [], "flow": None}
        per = self.a.PER_SEGMENT
        ci = 0
        for group in groups:
            codes = [c for c in group if (c in self.a.by_key or f"{c}·R" in self.a.by_key)]
            if len(codes) < 2:
                continue
            curve = np.asarray(curves[ci])
            nodes = (np.asarray(node_sets[ci]) if ci < len(node_sets)
                     else np.minimum(np.arange(len(codes)) * per, len(curve) - 1))
            ci += 1
            if rep["flow"] is None or len(codes) > rep["flow"]["n"]:
                rep["flow"] = {**self.flow(ch, codes), "n": len(codes)}   # tor glowny = najdluzsza grupa
            rep["order_errors"] += self.order(codes, curve, nodes)
            for k in range(len(codes) - 1):
                lo, hi = int(nodes[k]), int(nodes[k + 1])
                rep["path_checks"].append(self.segment(codes[k], codes[k + 1], curve[lo:hi + 1]))
        return rep


def resample(poly: np.ndarray, n: int) -> np.ndarray:
    """Rowne odstepy po dlugosci luku - dyskretny Frechet zalezy od probkowania."""
    d = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(poly, axis=0), axis=1))])
    if d[-1] < 1e-9:
        return np.repeat(poly[:1], n, axis=0)
    s = np.linspace(0.0, d[-1], n)
    return np.column_stack([np.interp(s, d, poly[:, i]) for i in range(3)])


def kinks(poly: np.ndarray):
    """Kat skretu i wychylenie (mm) kazdej wewnetrznej probki od cieciwy sasiadow.

    Sam kat myli: haczyk 0.5 mm miedzy probkami co 1 mm daje 130 stopni i jest
    niewidoczny. Zalamanie to kat razem z wychyleniem.
    """
    if len(poly) < 5:
        return np.zeros(0), np.zeros(0)
    # pierwsza i ostatnia probka wewnetrzna leza przy punkcie - tam kanal ma prawo
    # skrecic (punkt jest zakretem toru); szukamy zalaman MIEDZY punktami
    poly = poly[1:-1] if len(poly) > 6 else poly
    a, m, b = poly[:-2], poly[1:-1], poly[2:]
    u, w = m - a, b - m
    nu, nw = np.linalg.norm(u, axis=1), np.linalg.norm(w, axis=1)
    cos = (u * w).sum(axis=1) / np.maximum(nu * nw, 1e-12)
    ang = np.degrees(np.arccos(np.clip(cos, -1, 1)))
    ab = b - a
    L2 = np.maximum((ab * ab).sum(axis=1), 1e-12)
    t = np.clip(((m - a) * ab).sum(axis=1) / L2, 0, 1)
    dev = np.linalg.norm(m - (a + t[:, None] * ab), axis=1) * 1000
    return ang, dev


def discrete_frechet(P: np.ndarray, Q: np.ndarray) -> float:
    """Dyskretna odleglosc Frecheta (Eiter-Mannila), iteracyjnie."""
    D = np.linalg.norm(P[:, None, :] - Q[None, :, :], axis=2)
    n, m = D.shape
    ca = np.empty((n, m))
    ca[0, 0] = D[0, 0]
    for i in range(1, n):
        ca[i, 0] = max(ca[i - 1, 0], D[i, 0])
    for j in range(1, m):
        ca[0, j] = max(ca[0, j - 1], D[0, j])
    for i in range(1, n):
        for j in range(1, m):
            ca[i, j] = max(min(ca[i - 1, j], ca[i - 1, j - 1], ca[i, j - 1]), D[i, j])
    return float(ca[-1, -1])


def hausdorff(P: np.ndarray, Q: np.ndarray) -> float:
    D = np.linalg.norm(P[:, None, :] - Q[None, :, :], axis=2)
    return float(max(D.min(axis=1).max(), D.min(axis=0).max()))


# --- C. punkty ------------------------------------------------------------------
def point_summary(atlas: Atlas) -> dict:
    if not atlas.meta:
        return {"status": "SKIP", "note": "brak lokalnego data/point_meta.json"}
    from verify_points import Verifier
    v = Verifier(atlas)
    v.run()
    rows = [{"point": c["code"], "kind": c["kind"], "claimed": c["claimed"],
             "measured": round(float(c["measured"]), 2), "err": round(float(c["err"]), 2),
             "status": status(abs(c["err"]), 1.5 if c["kind"] == "krag" else TOL_CUN, 1e9),
             "ref": "opis lokalny (L1)"} for c in v.claims]
    return {"L0": "brak - wspolrzedne WHO 2008 / GB/T 12346 nie sa w repozytorium",
            "checks": rows}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--channels", default="", help="np. LU,LI (domyslnie wszystkie)")
    ap.add_argument("--no-frechet", action="store_true", help="bez geodezyjnych (szybciej)")
    args = ap.parse_args()
    t0 = time.time()

    atlas = Atlas()
    if atlas.mode != "anatomy":
        print("weryfikacja wymaga trybu anatomii (assets/ready)", file=sys.stderr)
        return 2

    ruler = Ruler(atlas).run()
    print("A. LINIJKA CUN - wektory kostne")
    print(f"   {'wektor':44s} {'norma':>5s} {'model':>6s} {'blad':>6s}  {'cun mm: model/wynika':>21s}")
    for r in ruler:
        if r["status"] == "SKIP":
            print(f"   {r['id']:44s}  SKIP  {r.get('note', '')}")
            continue
        tag = " (kalibr.)" if r["fitted"] else ""
        print(f"   {r['id']:44s} {r['norm_cun']:5.1f} {r['model_cun']:6.2f} {r['err_cun']:+6.2f}  "
              f"{r['cun_mm_model']:9.1f} / {r['cun_mm_implied']:<9.1f} {r['status']}{tag}")

    paths = Paths(atlas, frechet=not args.no_frechet)
    wanted = [c for c in args.channels.split(",") if c] or list(atlas.paths)
    channels = []
    for ch in wanted:
        channels.append(paths.channel(ch))
    handoffs = paths.handoffs()
    points = point_summary(atlas)

    by_pt = {}
    for c in points.get("checks", []):
        by_pt.setdefault(c["point"].split("-")[0], []).append(c)
    total = {"PASS": 0, "WARN": 0, "FAIL": 0}
    print("\nB. PRZEBIEGI")
    for rep in channels:
        rep["point_checks"] = by_pt.get(rep["channel_id"], [])
        sts = []
        for s in rep["path_checks"]:
            sts += [s["air_status"], s["bone_status"], s["angle_status"], s.get("frechet_status", "SKIP")]
        if rep["flow"]:
            sts.append(rep["flow"]["status"])
        sts += ["FAIL"] * len(rep["order_errors"])
        worst = "FAIL" if "FAIL" in sts else "WARN" if "WARN" in sts else "PASS"
        rep["status"] = worst
        total[worst] += 1
        n = len(rep["path_checks"])
        cnt = lambda key: sum(s.get(key) in ("FAIL", "WARN") if key == "bone_status"  # noqa: E731
                              else s.get(key) == "FAIL" for s in rep["path_checks"])
        fr = [s["frechet"] for s in rep["path_checks"] if s.get("frechet") is not None]
        print(f"   {rep['channel_id']:4s} {worst:4s}  odcinkow {n:3d}  powietrze {cnt('air_status'):2d}  "
              f"kosc(W) {cnt('bone_status'):2d}  zalamania {cnt('angle_status'):2d}  "
              f"Frechet FAIL {cnt('frechet_status'):2d}"
              + (f" (mediana {np.median(fr):.2f})" if fr else "")
              + f"  kolejnosc {len(rep['order_errors'])}  przeplyw "
              + (rep["flow"]["status"] if rep["flow"] else "-"))
    print("\n   przekazania w cyklu:")
    for h in handoffs:
        extra = f"{h['points']} {h['dist_cun']:.1f} cun" if "points" in h else "wewnetrzne"
        print(f"     {h['from']:3s}->{h['to']:3s} {h['kind']:9s} {extra:28s} {h['status']}")

    OUT.mkdir(parents=True, exist_ok=True)
    report = {
        "run_id": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "confidence": {"points": "L1 (opis lokalny); L0 niedostepne",
                       "paths": "samozgodnosc geometryczna; toru nie definiuje zadna norma",
                       "ruler": "wektory B-cun z normy, do potwierdzenia z egzemplarzem"},
        "thresholds": {"cun": TOL_CUN, "mm_head": TOL_MM, "air_mm": AIR_MM, "bone_mm": BONE_MM,
                       "angle_warn": ANGLE_WARN, "angle_fail": ANGLE_FAIL,
                       "skin_offset_mm": SKIN_OFFSET * 1000},
        "summary": {"channels": len(channels), **{k.lower(): v for k, v in total.items()}},
        "ruler": ruler, "handoffs": handoffs, "channels": channels,
        "points_L0": points.get("L0"),
    }
    (OUT / "verify_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=1,
                                                       default=float), encoding="utf-8")
    with open(OUT / "verify_report.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["channel_id", "check_type", "segment_or_point", "metric", "value_cun",
                    "value_mm", "threshold", "status", "source", "note"])
        for r in ruler:
            if r["status"] != "SKIP":
                w.writerow(["-", "ruler", r["id"], "err_cun", r["err_cun"], r["length_mm"],
                            RULER_FAIL, r["status"], "B-cun (norma)", r["how"]])
        for rep in channels:
            ch = rep["channel_id"]
            for s in rep["path_checks"]:
                w.writerow([ch, "air", s["segment"], "max_mm", "", s["air_max_mm"], AIR_MM,
                            s["air_status"], "samozgodnosc", ""])
                w.writerow([ch, "bone", s["segment"], "min_mm", "", s["bone_min_mm"], BONE_MM,
                            s["bone_status"], "samozgodnosc", ""])
                w.writerow([ch, "kink", s["segment"], "kink_mm", "", s["kink_mm"], KINK_FAIL_MM,
                            s["angle_status"], "samozgodnosc", f"{s['angle_max_deg']} deg"])
                if s.get("frechet") is not None:
                    cun = s["frechet"] if s["unit"] == "cun" else ""
                    mm = s["frechet"] if s["unit"] == "mm" else ""
                    w.writerow([ch, "frechet", s["segment"], "frechet", cun, mm,
                                TOL_CUN if cun != "" else TOL_MM, s["frechet_status"],
                                "geodezyjna po skorze", ""])
            if rep["flow"]:
                w.writerow([ch, "flow", "-", "ends_ok", "", "", "", rep["flow"]["status"],
                            "numeracja + regula przeplywu", rep["flow"]["rule"]])
            for e in rep["order_errors"]:
                w.writerow([ch, "order", e["point"], "gap_mm", "", e["gap_mm"], "", "FAIL",
                            "samozgodnosc", ""])
            for c in rep["point_checks"]:
                w.writerow([ch, f"point:{c['kind']}", c["point"], "err", c["err"], "", TOL_CUN,
                            c["status"], c["ref"], ""])
    print(f"\nkanaly: PASS {total['PASS']}  WARN {total['WARN']}  FAIL {total['FAIL']}   "
          f"({time.time() - t0:.0f} s)  ->  data/verify/verify_report.json, .csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
