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

    def describe(self, pos) -> dict:
        """Co model wie o tym miejscu. Same liczby wypieczone z siatek."""
        i = self.nearest(pos)
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
