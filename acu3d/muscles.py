"""Miesnie: wybor pojedynczego miesnia i modulacja jego parametrow.

Wszystkie miesnie trzymamy w JEDNEJ siatce z tablica `mid` (indeks miesnia)
na punktach i komorkach. Dzieki temu rysuje je jeden aktor (szybko), a mimo to
wiadomo, w ktory miesien trafil klikniety piksel, i mozna ruszac tylko jego
wierzcholkami.

Model napiecia
--------------
Miesien kurczy sie wzdluz swojej osi glownej, a przyczepy zostaja na kosciach:

    s   = s_max * (alfa * sztywnosc + tonus)      skrocenie
    w(t)= sin(pi * t)                             profil - 0 na koncach, 1 w brzuscu
    u'  = u * (1 - s*w)                           wzdluz osi
    r'  = r / sqrt(1 - s*w)                       w poprzek (objetosc bez zmian)

Przerost/zanik to osobna, jednorodna zmiana promienia - nie miesza sie
z napieciem, bo to inne zjawisko.
"""
from __future__ import annotations

import json
import math
import warnings
from pathlib import Path

import numpy as np
import pyvista as pv

warnings.simplefilter("ignore")
ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
READY = ROOT / "assets" / "ready"
MERGED = READY / "muscles_merged.vtp"

S_MAX = 0.16                  # maksymalne skrocenie brzusca przy pelnej aktywacji


class Muscle:
    """Jeden miesien: metryka, parametry i zakres wierzcholkow we wspolnej siatce."""

    def __init__(self, rec: dict, index: int, lo: int, hi: int):
        self.index = index
        self.id = rec["id"]
        self.name = rec["name"]
        self.side = rec["side"]
        self.latin = rec.get("latin", "")
        self.polys = rec["polys"]
        self.volume_l = rec.get("volume_l", 0.0)
        self.centroid = np.array(rec["centroid"])
        self.axis = np.array(rec["axis"])
        self.extent = tuple(rec["extent"])
        self.slice = slice(lo, hi)
        # parametry sterowane z panelu
        self.activation = 0.0
        self.tone = 0.0
        self.hypertrophy = 1.0
        self.stiffness = 1.0

    @property
    def label(self) -> str:
        s = {"L": "lewy", "R": "prawy", "C": ""}.get(self.side, "")
        return f"{self.name}{(' (' + s + ')') if s else ''}"

    @property
    def length_m(self) -> float:
        return float(self.extent[1] - self.extent[0])

    def shortening(self) -> float:
        """Wzgledne skrocenie brzusca dla obecnych parametrow."""
        drive = self.activation * self.stiffness + self.tone
        return float(np.clip(S_MAX * drive, 0.0, 0.45))

    def force_ratio(self) -> float:
        """Uproszczona krzywa sila-dlugosc Hilla wokol dlugosci optymalnej."""
        rel = 1.0 - self.shortening() * 0.5
        fl = math.exp(-((rel - 1.0) / 0.45) ** 2)
        return float(self.activation * fl)


class MuscleSet:
    """Zbior miesni w jednej siatce, z deformacja liczona na zywo."""

    def __init__(self, catalog_path: Path | None = None):
        doc = json.loads((catalog_path or DATA / "muscles.json").read_text(encoding="utf-8"))
        self.records = doc["muscles"]
        self.mesh, spans = self._build()
        self.rest = np.array(self.mesh.points, copy=True)
        self.muscles = [Muscle(r, i, *spans[i]) for i, r in enumerate(self.records)]
        self.by_id = {m.id: m for m in self.muscles}
        self._mid_points = np.asarray(self.mesh.point_data["mid"])
        self._on_cache: dict = {}          # miesien -> punkty nad nim (liczone raz)

    # --- budowa ----------------------------------------------------------
    def _build(self):
        if MERGED.exists():
            mesh = pv.read(MERGED)
            mid = np.asarray(mesh.point_data["mid"])
            spans = []
            for i in range(len(self.records)):
                idx = np.flatnonzero(mid == i)
                spans.append((int(idx[0]), int(idx[-1]) + 1) if len(idx) else (0, 0))
            return mesh, spans
        parts, spans, off = [], [], 0
        for i, rec in enumerate(self.records):
            m = pv.read(READY / rec["file"]).triangulate()
            m.clear_data()
            m.point_data["mid"] = np.full(m.n_points, i, dtype=np.int32)
            m.cell_data["mid"] = np.full(m.n_cells, i, dtype=np.int32)
            parts.append(m)
            spans.append((off, off + m.n_points))
            off += m.n_points
        mesh = parts[0].merge(parts[1:], merge_points=False)
        mesh.save(MERGED, binary=True)
        return mesh, spans

    # --- deformacja ------------------------------------------------------
    def deform(self, muscle: Muscle, rest: np.ndarray = None) -> np.ndarray:
        """Nowe polozenia wierzcholkow miesnia. Czysta funkcja - nic nie zapisuje.

        Liczy sie w ulamku milisekundy, wiec moze chodzic za suwakiem w czasie
        rzeczywistym; kosztowne bylo tylko to, co bylo wokol niej.
        """
        rest = self.rest[muscle.slice] if rest is None else rest
        if rest.size == 0:
            return rest
        a = muscle.axis
        c = muscle.centroid
        rel = rest - c
        u = rel @ a                                   # skladowa wzdluz osi
        radial = rel - np.outer(u, a)                 # skladowa w poprzek
        lo, hi = muscle.extent
        span = max(hi - lo, 1e-6)
        t = np.clip((u - lo) / span, 0.0, 1.0)
        w = np.sin(np.pi * t)                         # 0 na przyczepach, 1 w brzuscu
        s = muscle.shortening()
        k = 1.0 - s * w
        scale = muscle.hypertrophy / np.sqrt(np.maximum(k, 0.15))
        return c + np.outer(u * k, a) + radial * scale[:, None]

    def apply(self, muscle: Muscle) -> None:
        """Wpisuje deformacje do wspolnej siatki (uzywane poza trybem interaktywnym)."""
        if muscle.slice.stop > muscle.slice.start:
            self.mesh.points[muscle.slice] = self.deform(muscle)

    def reset(self, muscle: Muscle) -> None:
        muscle.activation = 0.0
        muscle.tone = 0.0
        muscle.hypertrophy = 1.0
        muscle.stiffness = 1.0
        self.mesh.points[muscle.slice] = self.rest[muscle.slice]

    def refresh(self) -> None:
        self.mesh.Modified()

    # --- wyszukiwanie ----------------------------------------------------
    def at_position(self, pos: np.ndarray) -> Muscle | None:
        """Miesien najblizszy wskazanemu punktowi w przestrzeni."""
        i = self.mesh.find_closest_point(np.asarray(pos, dtype=float))
        if i < 0:
            return None
        return self.muscles[int(self._mid_points[i])]

    def find(self, query: str) -> list:
        q = query.lower()
        return [m for m in self.muscles if q in m.name.lower()]

    def points_on(self, muscle: Muscle, markers: list, radius: float = 0.022) -> list:
        """Punkty akupunkturowe lezace nad tym miesniem.

        Wynik zalezy od pozycji spoczynkowej, a nie od parametrow, wiec liczymy go
        raz na miesien i trzymamy w pamieci - inaczej kazde drgniecie suwaka
        kosztowaloby 36 ms.
        """
        cached = self._on_cache.get(muscle.index)
        if cached is not None:
            return cached
        sub = self.rest[muscle.slice]
        if sub.size == 0:
            return []
        tree = pv.PolyData(sub)
        out = []
        for m in markers:
            j = tree.find_closest_point(m["pos"])
            d = float(np.linalg.norm(sub[j] - m["pos"]))
            if d <= radius:
                out.append((m, d))
        out = sorted(out, key=lambda x: x[1])
        self._on_cache[muscle.index] = out
        return out

    def submesh(self, muscle: Muscle) -> pv.PolyData:
        """Osobna siatka jednego miesnia - do rysowania i deformacji w locie."""
        rec = self.records[muscle.index]
        return pv.read(READY / rec["file"]).triangulate()

    def muscle_under(self, pos: np.ndarray, max_dist: float = 0.05):
        """Miesien lezacy pod wskazanym punktem skory (albo None)."""
        i = self.mesh.find_closest_point(np.asarray(pos, dtype=float))
        if i < 0:
            return None, None
        d = float(np.linalg.norm(self.rest[i] - pos))
        if d > max_dist:
            return None, None
        return self.muscles[int(self._mid_points[i])], d
