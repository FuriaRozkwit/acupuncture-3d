"""Spina cialo (anatomiczne albo proceduralne), dane punktow i meridiany.

Tryb `anatomy` uzywa siatek z Z-Anatomy i BodyParts3D oraz rigu wyliczonego
z kosci (data/rig_fitted.json). Tryb `procedural` to poprzedni manekin - zostaje
jako zapasowy, gdy nie ma pobranych zasobow.
"""
from __future__ import annotations

import hashlib
import json
import warnings
from pathlib import Path

import numpy as np
import pyvista as pv

from . import anchors, curves
from .anatomy import Anatomy
from .geometry import build_layers
from .rig import Rig

DATA = Path(__file__).resolve().parent.parent / "data"
READY = Path(__file__).resolve().parent.parent / "assets" / "ready"
SKIN_OFFSET = 0.004          # metry - o tyle marker odstaje od skory


class Atlas:
    def __init__(self, height_m: float | None = None, procedural: bool = False):
        self.anatomy = None
        if not procedural and (READY / "skin.vtp").exists():
            self.anatomy = Anatomy()
            self.body = self.anatomy
            self.rig = self.anatomy                      # zgodnosc ze starym kodem
            self.mode = "anatomy"
            self.layers = {"skin": self.anatomy.skin, "skeleton": self.anatomy.skeleton}
            self._skin_tri = self.anatomy.skin
        else:
            self.rig = Rig(height_m)
            self.body = self.rig
            self.mode = "procedural"
            self.layers = build_layers(self.rig)
            self._skin_tri = self.layers["skin"].triangulate()
            self.rig.use_skin(self._skin_tri)      # ta sama skora dla rzutowania kotwic

        doc = json.loads((DATA / "points.json").read_text(encoding="utf-8"))
        self.channels = doc["channels"]
        self.points = doc["points"]
        meta_path = DATA / "point_meta.json"        # opisy z podrecznika - lokalne,
        self.meta = (json.loads(meta_path.read_text(encoding="utf-8"))   # nie w repozytorium
                     if meta_path.exists() else {})
        self.paths = json.loads((DATA / "meridians.json").read_text(encoding="utf-8"))["paths"]
        self.curve_nodes: dict = {}      # (kanal, strona) -> indeksy punktow na krzywych
        self.markers = self._resolve_markers()
        self.by_key = {m["key"]: m for m in self.markers}
        self._curves = None

    # --- rozwiazywanie punktow ------------------------------------------
    def _ref_point(self, anchor: dict, pos: np.ndarray) -> np.ndarray:
        """Punkt wewnatrz ciala, z ktorego 'patrzymy' na punkt."""
        if "struct" in anchor:
            return anchors.struct_base(self.body, anchor)
        seg = self.body.seg(anchor.get("seg", "head"))
        t = 0.9 if "arc_cun" in anchor else anchors._axial_t(self.body, anchor)
        return seg.interior(t) if hasattr(seg, "interior") else seg.axis_point(t)

    # --- rzutowanie na powierzchnie (tylko tryb proceduralny) ------------
    STEPS = 30
    SPAN = 0.60

    def _march(self, refs: np.ndarray, pos: np.ndarray) -> np.ndarray:
        """Wypycha punkty na widoczna powierzchnie sumy bryl proceduralnych."""
        refs, pos = np.atleast_2d(refs), np.atleast_2d(pos)
        d = pos - refs
        n = np.linalg.norm(d, axis=1, keepdims=True)
        n = np.where(n < 1e-9, 1e-9, n)
        u = d / n
        fac = 1.0 + np.linspace(0.0, self.SPAN, self.STEPS)
        cand = refs[:, None, :] + u[:, None, :] * (n[:, :, None] * fac[None, None, :]).transpose(0, 2, 1)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sel = pv.PolyData(cand.reshape(-1, 3)).select_enclosed_points(
                self._skin_tri, tolerance=0.0, check_surface=False)
        inside = np.asarray(sel["SelectedPoints"]).astype(bool).reshape(len(pos), self.STEPS)
        idx = np.where(inside.all(axis=1), 0, inside.argmin(axis=1))
        return cand[np.arange(len(pos)), idx] + u * SKIN_OFFSET

    def place(self, refs: np.ndarray, pos: np.ndarray) -> np.ndarray:
        """Ostateczne polozenie markerow na skorze."""
        if self.mode == "procedural":
            return self._march(refs, pos)
        refs, pos = np.atleast_2d(refs), np.atleast_2d(pos)
        # punkt juz lezy na skorze - odsuwamy go wzdluz normalnej skierowanej NA ZEWNATRZ.
        # Kierunek "od kosci" zawodzil na dloni: kosc wystaje tam poza siatke skory,
        # wiec marker i krzywa lądowaly pod skora (znikajace PC-8, PC-9).
        cells, foot = self._on_skin(pos)
        n = self._tri_normals()[cells]
        n = n * self._outward_sign(foot, n, pos - refs)[:, None]
        return self._surface_out(foot + n * SKIN_OFFSET, n)

    def _surface_out(self, pts: np.ndarray, normals: np.ndarray, keep=None) -> np.ndarray:
        """Ostatnia kontrola: nic, co ma byc widac, nie moze lezec w ciele.

        W dloni siatka skory sama sie przenika (faldy palcow wchodza w dlon) - tam
        nawet normalna skierowana na zewnatrz konczy sie wewnatrz. Szukamy wtedy
        wzdluz normalnej, w obie strony, najblizszego miejsca poza cialem.
        """
        body = self.anatomy
        out = pts.copy()
        steps = np.arange(0.001, 0.0121, 0.001)
        for i, (p, n) in enumerate(zip(pts, normals)):
            # Linii z linijki nie ruszamy: ona z definicji lezy na skorze, a w szczelinie
            # (pacha, krocze) test parzystosci i tak mowi "wewnatrz", bo punkt jest
            # otoczony przez sasiednie czesci ciala.
            if (keep is not None and keep[i]) or not body.is_inside(p):
                continue
            for s in steps:
                if not body.is_inside(p + n * s):
                    out[i] = p + n * (s + 0.002)
                    break
                if not body.is_inside(p - n * s):
                    out[i] = p - n * (s + 0.002)
                    break
        return out

    def _resolve_markers(self) -> list:
        out, raws, refs = [], [], []
        for p in self.points:
            side, raw = anchors.resolve_sides(self.body, p)[0]      # "R" albo "C"
            r = self._ref_point(p["anchor"], raw).copy()
            if side == "C":
                r[0] = 0.0
            raws.append(raw)
            refs.append(r)
            out.append({
                "key": p["code"] if side == "C" else f"{p['code']}·R",
                "code": p["code"], "channel": p["channel"], "index": p["index"],
                "side": side, "ref": r, "accuracy": p["accuracy"],
                "pinyin": p.get("pinyin", ""), "english": p.get("english", ""),
                "note_pl": p.get("note_pl", ""),
            })
        pos = self.place(np.array(refs), np.array(raws))
        markers = []
        for rec, pt in zip(out, pos):
            rec["pos"] = pt
            markers.append(rec)
            if rec["side"] == "R":                                  # lewa strona = odbicie
                mir = dict(rec)
                mir["key"] = f"{rec['code']}·L"
                mir["side"] = "L"
                mir["pos"] = pt * np.array([-1.0, 1.0, 1.0])
                mir["ref"] = rec["ref"] * np.array([-1.0, 1.0, 1.0])
                markers.append(mir)
        return markers

    # --- meridiany -------------------------------------------------------
    PER_SEGMENT = 24        # najmniej probek miedzy sasiednimi punktami
    SAMPLE_MM = 5.0         # ...a przy dluzszych odcinkach mniej wiecej co tyle mm
    MAX_SAMPLES = 96
    SMOOTH_ROUNDS = 6       # wygladzanie + rzut: zdejmuje drgania, zostawia ksztalt splajnu
    SHORTEN_ROUNDS = 80     # skracanie po powierzchni: lokalna geodezyjna miedzy punktami
    AIR_MM = 6.0            # srodek kroku dalej od skory = cieciwa przez powietrze
    KINK_DEG, KINK_MM = 60.0, 2.0
    DETOUR = 1.25           # odcinek dluzszy od geodezyjnej o tyle (i o 1 cm) = objazd

    def _cache_path(self):
        h = hashlib.md5()
        for name in ("points.json", "meridians.json", "rig_fitted.json"):
            f = DATA / name
            if f.exists():
                h.update(f.read_bytes())
        # krzywe zaleza tez od kodu, ktory rozwiazuje kotwice - zmiana linijki
        # w samym kodzie musi uniewaznic stare krzywe
        for mod in ("anchors.py", "anatomy.py", "curves.py", "model.py"):
            h.update((Path(__file__).resolve().parent / mod).read_bytes())
        h.update(f"{self.mode}-{self.PER_SEGMENT}".encode())
        return READY / f"meridians_{h.hexdigest()[:12]}.npz"

    def _load_curves(self) -> dict:
        """Krzywe zaleza tylko od danych, wiec liczymy je raz i trzymamy na dysku."""
        path = self._cache_path()
        if path.exists():
            try:
                with np.load(path) as z:
                    out = {}
                    for k in sorted(z.files):
                        ch, side, i, kind = k.split("|")
                        if kind == "nodes":
                            self.curve_nodes.setdefault((ch, side), []).append(z[k])
                        else:
                            out.setdefault((ch, side), []).append(z[k])
                    return out
            except Exception:
                pass
        keys = [(ch, side) for ch in self.channels for side in ("R", "L")]
        if self.mode == "procedural":
            out = {k: self._build_meridian(*k) for k in keys}
        else:
            out = self._build_on_skin(keys)
        try:
            READY.mkdir(parents=True, exist_ok=True)
            for old in READY.glob("meridians_*.npz"):
                old.unlink()
            data = {f"{ch}|{side}|{i:02d}|line": arr
                    for (ch, side), lines in out.items() for i, arr in enumerate(lines)}
            data.update({f"{ch}|{side}|{i:02d}|nodes": arr
                         for (ch, side), ns in self.curve_nodes.items()
                         for i, arr in enumerate(ns)})
            np.savez_compressed(path, **data)
        except Exception:
            pass
        return out

    def meridian_lines(self, channel: str, side: str = "R") -> list:
        if self._curves is None:
            self._curves = self._load_curves()
        return self._curves.get((channel, side), [])

    def _build_meridian(self, channel: str, side: str = "R") -> list:
        """Przebiegi kanalu jako gladkie krzywe oblepione na powierzchni ciala.

        Kanal nie biegnie lamana od punktu do punktu, wiec zamiast odcinkow
        prowadzimy splajn Catmulla-Roma (lancuch krzywych Beziera). Splajn
        przechodzi przez wezly, a po zrzutowaniu na skore wezly sa jeszcze raz
        przypinane dokladnie do polozen punktow - krzywa ma je laczyc, a nie
        przesuwac.
        """
        lines = []
        for group in self.paths.get(channel, []):
            pts, refs = [], []
            for code in group:
                m = self.by_key.get(code) or self.by_key.get(f"{code}·{side}")
                if m is None:
                    continue
                pts.append(m["pos"])
                refs.append(m["ref"])
            if len(pts) < 2:
                continue
            P, R = np.array(pts), np.array(refs)
            dense_p, dense_r = curves.sample_pair(P, R, self.PER_SEGMENT)
            nodes = np.arange(len(P)) * self.PER_SEGMENT
            nodes[-1] = len(dense_p) - 1
            pinned = np.zeros(len(dense_p), dtype=bool)
            pinned[nodes] = True

            hug = self._hug(dense_r, dense_p)
            smoothed = curves.relax(hug, pinned)     # zdejmij drgania po rzutowaniu
            # Wygladzanie rusza wiekszosc probek o ulamek milimetra - te zostawiamy.
            # Rzutujemy ponownie tylko tam, gdzie ruch byl duzy, czyli gdzie promien
            # wczesniej uciekl i krzywa naprawde zeszla ze skory.
            moved = np.linalg.norm(smoothed - hug, axis=1) > 0.001
            if moved.any():
                idx = np.flatnonzero(moved)
                smoothed[idx] = self._hug(dense_r[idx], smoothed[idx])
            hug = smoothed
            hug[nodes] = P                           # wezly dokladnie w punktach
            lines.append(hug)
        return lines

    # --- przebiegi po powierzchni skory (tryb anatomii) ----------------------
    def _anchor_line(self, a: dict, b: dict, side: str, n: int):
        """Linia miedzy punktami prowadzona LINIJKA, a nie cieciwa splajnu.

        Kanal miedzy dwoma punktami nie idzie na skroty: na sklepieniu czaszki
        trzyma stala odleglosc od linii posrodkowej, na konczynie staly azymut,
        na tulowiu stala odleglosc od linii posrodkowej. Gdy obie kotwice sa tego
        samego rodzaju, probki miedzy nimi liczymy tak samo jak same punkty -
        splajn zostaje tam, gdzie kotwice sa rozne (np. przedramie -> dlon).
        """
        body = self.body
        f = np.linspace(0.0, 1.0, n, endpoint=False)[1:]
        mix = lambda k, d=0.0: float(a.get(k, d)) + (float(b.get(k, d)) - float(a.get(k, d))) * u  # noqa: E731
        out = []
        if "arc_cun" in a and "arc_cun" in b:
            for u in f:
                out.append(anchors.resolve(body, {"seg": "head", "arc_cun": mix("arc_cun"),
                                                  "lat": mix("lat"), "out": mix("out")}))
        elif (a.get("seg") == b.get("seg") and a.get("seg") is not None
              and "arc_cun" not in a and "struct" not in a and "struct" not in b):
            seg = body.seg(a["seg"])
            ta, tb = anchors._axial_t(body, a), anchors._axial_t(body, b)
            if "az" in a and "az" in b:
                az0 = float(a["az"])
                d_az = (float(b["az"]) - az0 + 180.0) % 360.0 - 180.0   # krotszym lukiem
                for u in f:
                    t = ta + (tb - ta) * u
                    out.append(seg.surface_by_azimuth(t, az0 + d_az * u, mix("out")))
            elif "lat" in a and "lat" in b and a.get("face") == b.get("face"):
                for u in f:
                    t = ta + (tb - ta) * u
                    out.append(seg.surface_by_lateral(t, mix("lat"), a.get("face", "anterior"),
                                                      mix("out")))
            else:
                # rozne style kotwic na tym samym segmencie (raz kat, raz odstep od
                # posrodkowej) - sprowadzamy oba konce do polozenia wzdluz segmentu
                # i kata wokol jego osi, a miedzy nimi interpolujemy
                pa, pb = anchors.resolve(body, a), anchors.resolve(body, b)
                az_a, az_b = self._seg_azimuth(seg, ta, pa), self._seg_azimuth(seg, tb, pb)
                d_az = (az_b - az_a + 180.0) % 360.0 - 180.0
                for u in f:
                    t = ta + (tb - ta) * u
                    out.append(seg.surface_by_azimuth(t, az_a + d_az * u))
        else:
            return None
        P = np.array(out)
        if side == "L":
            P = P * np.array([-1.0, 1.0, 1.0])
        return P

    @staticmethod
    def _seg_azimuth(seg, t: float, p: np.ndarray) -> float:
        v = np.asarray(p) - seg.interior(t)
        v = v - seg.e_ax * float(np.dot(v, seg.e_ax))
        return float(np.degrees(np.arctan2(np.dot(v, seg.e_lat), np.dot(v, seg.e_ant))))

    def _groups(self, channel: str, side: str) -> list:
        anchor_of = {p["code"]: p["anchor"] for p in self.points}
        out = []
        for group in self.paths.get(channel, []):
            pts, ancs = [], []
            for code in group:
                m = self.by_key.get(code) or self.by_key.get(f"{code}·{side}")
                if m is not None:
                    pts.append(m["pos"])
                    ancs.append(anchor_of.get(code, {}))
            if len(pts) >= 2:
                out.append((np.array(pts), ancs))
        return out

    def _on_skin(self, pts: np.ndarray):
        """Najblizszy punkt POWIERZCHNI skory (na trojkacie, nie w wezle siatki)."""
        cells, foot = self.anatomy.skin_cast.find_closest_cell(pts, return_closest_point=True)
        return np.atleast_1d(cells), np.atleast_2d(foot)

    def _outward_sign(self, foot: np.ndarray, normals: np.ndarray, hint: np.ndarray) -> np.ndarray:
        """+1 gdy normalna wskazuje na zewnatrz ciala, -1 gdy do srodka.

        Rozstrzyga parzystosc: z punktu na zewnatrz promien przecina skore parzysta
        liczbe razy. Liczymy to z obu stron skory w kilku kierunkach, bo siatka ma
        otwory (oczy, usta, szczatek ucha). Gdy glosowanie nie jest wyrazne, decyduje
        podpowiedz (kierunek od wnetrza ciala).
        """
        import vtk
        if getattr(self, "_obb", None) is None:
            self._obb = vtk.vtkOBBTree()
            self._obb.SetDataSet(self.anatomy.skin_cast)
            self._obb.BuildLocator()
            rng = np.random.default_rng(7)
            d = rng.normal(size=(7, 3))
            self._rays = d / np.linalg.norm(d, axis=1, keepdims=True)
        hits = vtk.vtkPoints()

        def inside(p) -> int:
            v = 0
            for d in self._rays:
                self._obb.IntersectWithLine(list(map(float, p)), list(map(float, p + d * 2.5)),
                                            hits, None)
                v += hits.GetNumberOfPoints() % 2
            return v

        sign = np.where((normals * hint).sum(axis=1) < 0, -1.0, 1.0)
        for i, (f, n) in enumerate(zip(foot, normals)):
            vote = inside(f - n * 0.003) - inside(f + n * 0.003)   # +7: n na zewnatrz
            if abs(vote) >= 5:
                sign[i] = 1.0 if vote > 0 else -1.0
        return sign

    def _tri_normals(self) -> np.ndarray:
        if getattr(self, "_normals", None) is None:
            sk = self.anatomy.skin_cast
            f = sk.faces.reshape(-1, 4)[:, 1:]
            v = np.asarray(sk.points)
            n = np.cross(v[f[:, 1]] - v[f[:, 0]], v[f[:, 2]] - v[f[:, 0]])
            self._normals = n / np.maximum(np.linalg.norm(n, axis=1, keepdims=True), 1e-12)
        return self._normals

    def _build_on_skin(self, keys: list) -> dict:
        """Wszystkie przebiegi naraz, prowadzone po powierzchni skory.

        Wczesniej kazda probke splajnu wypychalo sie promieniem z wnetrza ciala.
        Miedzy palcami, w pachwinie i przy uchu promien uciekal i probka siadala
        na sasiednim kawalku skory - stad cieciwy przez powietrze i zygzaki.
        Teraz probki rzutujemy na najblizszy punkt powierzchni, wygladzamy
        i rzutujemy znowu. Tam, gdzie mimo to odcinek odrywa sie, lamie albo
        robi objazd, zastepuje go krzywa skrocona po powierzchni (lokalna
        geodezyjna). Wszystkie kanaly ida w jednej tablicy, bo jedno zapytanie
        o najblizsze punkty kosztuje prawie tyle samo dla 100 co dla 20 000 probek.
        """
        per = self.PER_SEGMENT
        jobs, dense, pinned = [], [], []
        start = 0
        for key in keys:
            for P, ancs in self._groups(*key):
                # dluzszy odcinek dostaje wiecej probek: przy 24 na kazdy odcinek
                # krok siegal 24 mm i krzywa scinala pachwine zamiast ja obrysowac
                L = np.linalg.norm(np.diff(P, axis=0), axis=1)
                counts = np.clip(np.ceil(L * 1000 / self.SAMPLE_MM), per, self.MAX_SAMPLES)
                d, nodes = curves.sample_varied(P, counts)
                ruler = np.zeros(len(P) - 1, dtype=bool)
                for k in range(len(P) - 1):
                    line = self._anchor_line(ancs[k], ancs[k + 1], key[1], int(counts[k]))
                    if line is not None and len(line) == nodes[k + 1] - nodes[k] - 1:
                        d[nodes[k] + 1:nodes[k + 1]] = line
                        ruler[k] = True
                pin = np.zeros(len(d), dtype=bool)
                pin[nodes] = True
                jobs.append((key, start, len(d), nodes, P, ruler))
                dense.append(d)
                pinned.append(pin)
                start += len(d)
        out = {k: [] for k in keys}
        if not jobs:
            return out
        pin = np.concatenate(pinned)          # konce grup przypiete - grupy sie nie mieszaja
        _, foot = self._on_skin(np.vstack(dense))
        ruler_mask = np.zeros(len(foot), dtype=bool)
        for (ch, side), s0, n, nodes, P, ruler in jobs:
            for k, is_ruler in enumerate(ruler):
                if is_ruler:
                    ruler_mask[s0 + nodes[k]:s0 + nodes[k + 1] + 1] = True
        smooth = self._settle(foot.copy(), pin | ruler_mask, self.SMOOTH_ROUNDS)
        short = self._settle(foot.copy(), pin, self.SHORTEN_ROUNDS)
        mid = 0.5 * (smooth[:-1] + smooth[1:])
        _, mid_foot = self._on_skin(mid)
        air = np.linalg.norm(mid_foot - mid, axis=1)

        chosen = smooth.copy()
        self.shortcuts = []                  # do raportu: ktore odcinki zastapiono i czemu
        for (ch, side), s0, n, nodes, P, ruler in jobs:
            for k in range(len(nodes) - 1):
                lo, hi = s0 + nodes[k], s0 + nodes[k + 1]
                why = (self._defect(smooth[lo:hi + 1], short[lo:hi + 1], air[lo:hi], True)
                       if ruler[k] else "styk segmentow - najkrotsza droga po skorze")
                if why:
                    chosen[lo:hi + 1] = short[lo:hi + 1]
                    if side != "L" and ruler[k]:
                        self.shortcuts.append((ch, k, why))

        # krotka droga po powierzchni potrafi zostawic drobne haczyki na zagieciach
        # (pachwina, pacha) - kilka rund wygladzania z powrotnym rzutem je zdejmuje,
        # a linii z linijki nie rusza
        chosen = self._settle(chosen, pin | ruler_mask, 3)
        # rowne odstepy miedzy probkami: splajn dosrodkowy zageszcza je przy bliskich
        # punktach i rozrzedza na dlugich odcinkach (krok do 33 mm), a dlugi krok
        # przeskakuje wklesłosc - stad "powietrze" w kroczu i haczyki na udzie
        for key, s0, cnt, nodes, P, ruler in jobs:
            chosen[s0:s0 + cnt] = self._respace(chosen[s0:s0 + cnt], nodes)
        _, chosen = self._on_skin(chosen)
        cells, chosen = self._on_skin(chosen)
        normals = self._tri_normals()[cells]
        hint = np.empty_like(normals)             # podpowiedz zwrotu: odsuniecie markerow
        for key, s0, cnt, nodes, P, ruler in jobs:
            hint[s0:s0 + cnt] = self._blend_offsets(chosen[s0:s0 + cnt], P, nodes)
        normals = normals * self._outward_sign(chosen, normals, hint)[:, None]
        for key, s0, cnt, nodes, P, ruler in jobs:
            line = self._lift(chosen[s0:s0 + cnt], normals[s0:s0 + cnt], P, nodes,
                              self._surface_out, ruler_mask[s0:s0 + cnt])
            line = self._min_clearance(line, normals[s0:s0 + cnt])
            line[nodes] = P
            out[key].append(line)
            self.curve_nodes.setdefault(key, []).append(np.asarray(nodes))
        return out

    @staticmethod
    def _blend_offsets(foot: np.ndarray, P: np.ndarray, nodes: np.ndarray) -> np.ndarray:
        off = P - foot[nodes]
        idx = np.arange(len(foot))
        seg = np.clip(np.searchsorted(nodes, idx, side="right") - 1, 0, len(nodes) - 2)
        f = ((idx - nodes[seg]) / np.maximum(nodes[seg + 1] - nodes[seg], 1))[:, None]
        return off[seg] * (1 - f) + off[seg + 1] * f

    @staticmethod
    def _respace(line: np.ndarray, nodes: np.ndarray) -> np.ndarray:
        """Rowne odstepy po dlugosci luku w kazdym odcinku; wezly zostaja na miejscu."""
        out = line.copy()
        for k in range(len(nodes) - 1):
            lo, hi = nodes[k], nodes[k + 1]
            seg = line[lo:hi + 1]
            d = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(seg, axis=0), axis=1))])
            if d[-1] < 1e-9:
                continue
            u = np.linspace(0.0, d[-1], hi - lo + 1)
            out[lo:hi + 1] = np.column_stack([np.interp(u, d, seg[:, i]) for i in range(3)])
        return out

    def _settle(self, pts: np.ndarray, pinned: np.ndarray, rounds: int) -> np.ndarray:
        for _ in range(rounds):
            pts = curves.relax(pts, pinned, rounds=1)
            _, pts = self._on_skin(pts)
        return pts

    def _defect(self, smooth: np.ndarray, short: np.ndarray, air: np.ndarray,
                ruler: bool = False) -> str:
        """Czemu odcinek po wygladzeniu nie nadaje sie do pokazania ('' = nadaje sie)."""
        if air.size and air.max() * 1000 > self.AIR_MM:
            return f"powietrze {air.max() * 1000:.0f} mm"
        a, m, b = smooth[:-2], smooth[1:-1], smooth[2:]
        u, w = m - a, b - m
        cos = (u * w).sum(1) / np.maximum(np.linalg.norm(u, axis=1) * np.linalg.norm(w, axis=1),
                                          1e-12)
        ang = np.degrees(np.arccos(np.clip(cos, -1, 1)))
        ab = b - a
        t = np.clip(((m - a) * ab).sum(1) / np.maximum((ab * ab).sum(1), 1e-12), 0, 1)
        dev = np.linalg.norm(m - (a + t[:, None] * ab), axis=1) * 1000
        # linia z linijki biegnie tam, gdzie cialo ma fald (pacha, sciegno Achillesa,
        # krocze) i ma prawo sie zagiac - odrzucamy ja dopiero przy ostrym zalamaniu
        limit = self.KINK_MM * 3 if ruler else self.KINK_MM
        if np.any((ang > self.KINK_DEG) & (dev > limit)):
            return f"zalamanie {dev.max():.0f} mm"
        if ruler:
            return ""            # linia z linijki MA prawo byc dluzsza od najkrotszej drogi
        L_s = np.linalg.norm(np.diff(smooth, axis=0), axis=1).sum()
        L_g = np.linalg.norm(np.diff(short, axis=0), axis=1).sum()
        if L_s > self.DETOUR * L_g and L_s - L_g > 0.010:
            return f"objazd {L_s * 1000:.0f} zamiast {L_g * 1000:.0f} mm"
        return ""

    @staticmethod
    def _lift(foot: np.ndarray, normals: np.ndarray, P: np.ndarray, nodes: np.ndarray,
              surface_out=None, keep=None) -> np.ndarray:
        """Odsuwa krzywa od skory na zewnatrz i przypina wezly do markerow.

        `normals` maja juz zwrot na zewnatrz (patrz _outward_sign).
        """
        off = P - foot[nodes]
        nrm = normals.copy()
        # sasiednie trojkaty rzadkiej siatki potrafia miec normalne rozne o 40 stopni,
        # co przy 4 mm odsuniecia daje szpilki - usredniamy kierunek wzdluz krzywej
        pin = np.zeros(len(nrm), dtype=bool)
        pin[nodes] = True
        nrm[nodes] = off / np.maximum(np.linalg.norm(off, axis=1, keepdims=True), 1e-12)
        nrm = curves.relax(nrm, pin, rounds=8, weight=0.6)
        nrm /= np.maximum(np.linalg.norm(nrm, axis=1, keepdims=True), 1e-12)
        line = foot + nrm * SKIN_OFFSET
        line[nodes] = P
        return line if surface_out is None else surface_out(line, nrm, keep)

    def _min_clearance(self, line: np.ndarray, normals: np.ndarray) -> np.ndarray:
        """Krzywa ma wystawac nad skore. W bruzdach (nadgarstek, lopatka) usredniony
        kierunek odsuniecia kladzie sie prawie stycznie i linia ledwie wystaje -
        polowa rurki chowa sie wtedy w skorze."""
        _, foot = self._on_skin(line)
        d = np.linalg.norm(line - foot, axis=1)
        need = np.clip(SKIN_OFFSET * 0.75 - d, 0.0, SKIN_OFFSET)
        if not need.any():
            return line
        # podnosimy PLYNNIE i lagodnie: skokowa poprawka pojedynczych probek sama
        # robila szpilki, a zbyt mocna zaginala krzywa na zagieciach
        k = np.ones(11) / 11.0
        need = np.convolve(np.pad(need, 5, mode="edge"), k, mode="valid")
        return line + normals * need[:, None]

    def _hug(self, refs: np.ndarray, pts: np.ndarray) -> np.ndarray:
        if self.mode == "procedural":
            return self._march(refs, pts)
        out = []
        for r, p in zip(refs, pts):
            d = p - r
            n = np.linalg.norm(d)
            if n < 1e-9:
                out.append(p)
                continue
            u = d / n
            out.append(self.anatomy.cast(r, u) + u * SKIN_OFFSET)
        return np.array(out)

    def channel_codes(self, channel: str) -> list:
        return [p["code"] for p in self.points if p["channel"] == channel]

    def describe(self, key: str) -> dict:
        m = self.by_key[key]
        return {**m, "meta": self.meta.get(m["code"], {})}
