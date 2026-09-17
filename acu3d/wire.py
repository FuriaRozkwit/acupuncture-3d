"""Wektorowa siatka cunowa - cialo bez ani jednej siatki anatomicznej.

Wiedza o tkankach nie znika razem z siatkami: zostaje wypieczona w skalarach
na wierzcholkach (`tools/bake_wire.py`). Kazdy wierzcholek wie, ile ma nad soba
tkanki miekkiej, jak gleboko lezy brzusiec i ktory to miesien. Dzieki temu
9144 wierzcholki i 211 kB na dysku odpowiadaja na te same pytania, na ktore
wczesniej odpowiadalo 700 MB siatek - tylko od razu, bez szukania.

Siatka nie jest rysunkiem ciala. To WIDOCZNY UKLAD WSPOLRZEDNYCH atlasu:
pierscienie to poziomy w cun, linie podluzne to azymuty. Punkty akupunkturowe
sa opisane dokladnie w tym ukladzie, wiec siatka pokazuje siatke odniesienia,
w ktorej istnieja - a nie sylwetke, ktora je zaslania.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pyvista as pv

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"

# skala barw grubosci tkanki miekkiej: od kosci tuz pod skora po gruby brzusiec
FLESH_COLORS = ["#f2f6fb", "#8fb6d6", "#4f7fae", "#2b4a6b"]
FLESH_RANGE = (2.0, 55.0)          # mm


class Wire:
    """Siatka cunowa z wypieczona wiedza o tkankach."""

    def __init__(self, path: Path | None = None):
        npz = np.load(path or DATA / "wire.npz")
        self.points = npz["points"].astype(float)
        self.edges = npz["edges"]
        self.seg_id = npz["seg_id"]
        self.cun = npz["cun"]
        self.az = npz["az"]
        self.flesh_mm = npz["flesh_mm"]
        self.muscle_mm = npz["muscle_mm"]
        self.muscle_id = npz["muscle_id"]
        self.fascia_mm = npz["fascia_mm"]
        # wierzcholki z ucieczka promienia zostaja w danych (maja poprawne cun
        # i azymut), ale oznaczenie nie moze na nie usiasc - nie leza na skorze
        self.hit = npz["hit"] if "hit" in npz else np.ones(len(self.points), bool)

        meta = json.loads((DATA / "wire.json").read_text(encoding="utf-8"))
        self.segments = meta["segments"]
        self.muscles = meta["muscles"]
        self.muscle_axis = np.array([m["axis"] for m in self.muscles], dtype=float)
        self.step_cun = meta["step_cun"]
        self.rays = meta["rays"]

        self.mesh = self._build()
        self._valid = np.flatnonzero(self.hit)
        self._cloud = pv.PolyData(self.points[self._valid])
        # rzadsza kopia do testu "czy odcinek biegnie przy ciele" - przy progu
        # 5 cm co czwarty wierzcholek wystarcza, a liczy sie cztery razy szybciej
        self._coarse = self.points[self._valid][::4]
        self._graph = None          # budowany przy pierwszym trasowaniu
        self._routes: dict = {}     # trasa liczona raz na pare wierzcholkow
        self._ptree = None          # drzewo do szybkiego korytarza

    def _build(self) -> pv.PolyData:
        lines = np.hstack([np.full((len(self.edges), 1), 2, dtype=np.int32),
                           self.edges]).ravel()
        m = pv.PolyData(self.points)
        m.lines = lines
        m.point_data["tkanka_mm"] = self.flesh_mm
        m.point_data["miesien_mm"] = self.muscle_mm
        return m

    # --- odpytywanie ------------------------------------------------------
    def nearest(self, pos) -> int:
        j = self._cloud.find_closest_point(np.asarray(pos, dtype=float))
        return int(self._valid[j])

    def snap(self, pos) -> np.ndarray:
        """Przyciaga wskazanie do siatki - oznaczenie ma lezec w ukladzie cun."""
        return self.points[self.nearest(pos)]

    def describe(self, pos, index: int = -1) -> dict:
        """Co model wie o tym miejscu. Same liczby wypieczone z siatek.

        `index` pozwala pominac szukanie, gdy wywolujacy juz zna wierzcholek -
        przy podgladzie pod kursorem oszczedza to polowe pracy na kazdy ruch myszy.
        """
        i = index if index >= 0 else self.nearest(pos)
        mid = int(self.muscle_id[i])
        mus = self.muscles[mid] if 0 <= mid < len(self.muscles) else None
        return {
            "index": i,
            "pos": self.points[i],
            "segment": self.segments[int(self.seg_id[i])],
            "cun": float(self.cun[i]),
            "az": float(self.az[i]),
            "flesh_mm": float(self.flesh_mm[i]),
            "muscle": mus,
            "muscle_mm": float(self.muscle_mm[i]),
            "muscle_axis": self.muscle_axis[mid] if mus else None,
            "fascia_mm": float(self.fascia_mm[i]),
        }

    def label(self, pos) -> str:
        d = self.describe(pos)
        side = {"L": "lewy", "R": "prawy", "C": ""}.get(
            (d["muscle"] or {}).get("side", ""), "")
        mus = (d["muscle"] or {}).get("name", "-")
        return (f"{d['segment']}  {d['cun']:.1f} cun  az {d['az']:.0f}"
                f"   tkanka {d['flesh_mm']:.0f} mm"
                f"   {mus}{(' (' + side + ')') if side else ''} "
                f"{d['muscle_mm']:.0f} mm")

    def off_body(self, a, b, samples: int = 24, near: float = 0.05) -> float:
        """Jaka czesc odcinka biegnie POZA cialem (0 = caly przy ciele, 1 = caly w powietrzu).

        Prosta miedzy dwoma odleglymi punktami - powiedzmy dlonia i stopa - idzie
        w duzej czesci przez powietrze. Wszystko, co "lezy wzdluz" takiej prostej,
        to przypadkowi sasiedzi koncow, a nie sciezka w ciele. Bez tej liczby
        zestawienie zachecaloby do nadinterpretacji.
        """
        a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
        t = np.linspace(0.0, 1.0, samples)[:, None]
        line = a + (b - a) * t
        # jedno rozgloszenie na rzadszej kopii: bez petli i bez narzutu VTK,
        # bo to jest sciezka klikniecia i kazda milisekunda jest odczuwalna
        d2 = ((line[:, None, :] - self._coarse[None, :, :]) ** 2).sum(axis=2)
        return float((d2.min(axis=1) > near * near).mean())

    # --- korytarz wzdluz linii sily ---------------------------------------
    def along(self, a, b, radius: float = 0.03) -> np.ndarray:
        """Indeksy wierzcholkow lezacych w korytarzu miedzy dwoma miejscami."""
        a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
        v = b - a
        n2 = float(np.dot(v, v))
        if n2 < 1e-12:
            return np.array([self.nearest(a)])
        t = np.clip((self.points - a) @ v / n2, 0.0, 1.0)
        d = np.linalg.norm(self.points - (a + t[:, None] * v), axis=1)
        return np.flatnonzero(d <= radius)

    # --- trasa po ciele ----------------------------------------------------
    def _build_graph(self):
        """Graf siatki z dlugosciami krawedzi jako wagami.

        Krawedzie sa te same, ktore widac na ekranie - lacznie ze szwami miedzy
        odcinkami. Dlatego trasa nie moze przeskoczyc z dloni na stope przez
        powietrze: musi przejsc przez przedramie, ramie, tulow i noge.
        """
        from scipy.sparse import coo_matrix
        e = self.edges
        w = np.linalg.norm(self.points[e[:, 0]] - self.points[e[:, 1]], axis=1)
        n = len(self.points)
        self._graph = coo_matrix((np.r_[w, w],
                                  (np.r_[e[:, 0], e[:, 1]],
                                   np.r_[e[:, 1], e[:, 0]])), shape=(n, n)).tocsr()
        return self._graph

    def route(self, pos_a, pos_b):
        """Najkrotsza droga po powierzchni ciala miedzy dwoma miejscami.

        Zwraca lamana przez wierzcholki siatki. Gdy trasy nie ma (wierzcholek
        odciety od reszty), zwraca odcinek prosty - wywolujacy pozna to po tym,
        ze dostal dwa punkty.
        """
        from scipy.sparse.csgraph import dijkstra
        g = self._graph if self._graph is not None else self._build_graph()
        i, j = self.nearest(pos_a), self.nearest(pos_b)
        # ta sama para pyta o trase dwa razy: raz analiza, raz rysowanie
        cached = self._routes.get((i, j))
        if cached is not None:
            return cached
        dist, pred = dijkstra(g, indices=i, return_predecessors=True)
        if not np.isfinite(dist[j]):
            out = (np.array([self.points[i], self.points[j]]), float("inf"))
            self._routes[(i, j)] = out
            return out
        path, k = [j], j
        while k != i and k >= 0:
            k = pred[k]
            path.append(k)
        path = np.array(path[::-1])
        out = (self.points[path], float(dist[j]))
        self._routes[(i, j)] = out
        self._routes[(j, i)] = (out[0][::-1].copy(), out[1])
        if len(self._routes) > 400:            # dzienna sesja nie zrobi ich tylu
            self._routes.clear()
        return out

    def near_path(self, poly: np.ndarray, radius: float = 0.03):
        """Wierzcholki w korytarzu wzdluz lamanej, plus numer mijanego odcinka.

        Odleglosc liczymy do WEZLOW trasy, nie do jej odcinkow. Wezly sa gesto -
        to sa wierzcholki tej samej siatki, odlegle o 1-2 cm - wiec blad jest
        mniejszy niz polowa dlugosci krawedzi, a liczy sie to drzewem zamiast
        rozgloszeniem 9144 na 155 i schodzi z 40 ms do 3 ms.
        """
        from scipy.spatial import cKDTree
        tree = cKDTree(poly)
        d, k = tree.query(self.points, distance_upper_bound=radius)
        idx = np.flatnonzero(np.isfinite(d))
        return idx, np.clip(k[idx], 0, max(len(poly) - 2, 0))

    def along_path(self, poly: np.ndarray, radius: float = 0.03) -> np.ndarray:
        return self.near_path(poly, radius)[0]

    def tangents(self, poly: np.ndarray) -> np.ndarray:
        """Kierunek lamanej na kazdym jej odcinku - potrzebny, bo trasa skreca."""
        v = poly[1:] - poly[:-1]
        n = np.linalg.norm(v, axis=1, keepdims=True)
        return v / np.where(n < 1e-12, 1e-12, n)


def point_to_path(pts: np.ndarray, poly: np.ndarray):
    """Odleglosc kazdego punktu do lamanej i polozenie wzdluz niej (w metrach).

    Liczone naraz dla wszystkich punktow i wszystkich odcinkow lamanej -
    to jest sciezka klikniecia, wiec petla po odcinkach bylaby odczuwalna.
    """
    pts = np.atleast_2d(pts)
    a = poly[:-1]
    v = poly[1:] - a
    L = np.linalg.norm(v, axis=1)
    L = np.where(L < 1e-12, 1e-12, L)
    start = np.r_[0.0, np.cumsum(L)[:-1]]
    w = pts[:, None, :] - a[None, :, :]                    # (punkty, odcinki, 3)
    t = np.clip((w * v[None]).sum(axis=2) / (L ** 2)[None], 0.0, 1.0)
    foot = a[None] + t[..., None] * v[None]
    d = np.linalg.norm(pts[:, None, :] - foot, axis=2)
    k = np.argmin(d, axis=1)
    rows = np.arange(len(pts))
    return d[rows, k], start[k] + t[rows, k] * L[k], k
