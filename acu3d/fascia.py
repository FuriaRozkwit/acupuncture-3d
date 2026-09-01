"""Warstwa lacznotkankowa: powiezi, przegrody, troczki, rozciegna, pochewki, kaletki.

Osobna warstwa, a nie domieszka do miesni - z dwoch powodow. Po pierwsze te
struktury otaczaja miesnie, wiec wrzucone do wspolnej siatki przechwytywalyby
kazde klikniecie. Po drugie palpuje sie je inaczej niz brzusiec: powiez jest
platem, nie walcem, i opisuje ja NORMALNA (kierunek w poprzek platu), a nie os.

Przegrody miedzymiesniowe maja tu status szczegolny: to one wyznaczaja bruzdy
miedzy grupami miesniowymi, wiec sluza do porownania z przebiegiem kanalow
(tools/analyze_septa.py).
"""
from __future__ import annotations

import json
import warnings
from pathlib import Path

import numpy as np
import pyvista as pv

warnings.simplefilter("ignore")
ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
READY = ROOT / "assets" / "ready"
MERGED = READY / "fascia_merged.vtp"

# kolejnosc = kolejnosc przelaczania klawiszem [h]
CLASS_ORDER = ("septum", "fascia", "retinaculum", "aponeurosis", "sheath", "bursa")
CLASS_PL = {"septum": "przegrody miedzymiesniowe", "fascia": "powiezi",
            "retinaculum": "troczki", "aponeurosis": "rozciegna",
            "sheath": "pochewki sciegniste", "bursa": "kaletki"}
CLASS_COLOR = {"septum": "#8fd4c2", "fascia": "#cfd8e3", "retinaculum": "#e0b46a",
               "aponeurosis": "#d9d2b0", "sheath": "#b6c9e0", "bursa": "#d59ac2"}


class Sheet:
    """Jedna struktura lacznotkankowa."""

    def __init__(self, rec: dict):
        self.index = rec["index"]
        self.id = rec["id"]
        self.name = rec["name"]
        self.cls = rec["class"]
        self.side = rec["side"]
        self.file = rec["file"]
        self.polys = rec["polys"]
        self.centroid = np.array(rec["centroid"])
        self.normal = np.array(rec["normal"])
        self.thickness = rec["thickness"]
        self.span = rec["span"]

    @property
    def label(self) -> str:
        s = {"L": "lewa", "R": "prawa", "C": ""}.get(self.side, "")
        return f"{self.name}{(' (' + s + ')') if s else ''}"


class FasciaSet:
    """Wszystkie platy w jednej siatce z tablica `fid`."""

    def __init__(self, catalog_path: Path | None = None):
        path = catalog_path or DATA / "fascia.json"
        doc = json.loads(path.read_text(encoding="utf-8"))
        self.sheets = [Sheet(r) for r in doc["items"]]
        self.by_id = {s.id: s for s in self.sheets}
        self.mesh = pv.read(MERGED)
        self._fid = np.asarray(self.mesh.point_data["fid"])
        self._fid_cells = np.asarray(self.mesh.cell_data["fid"])
        self.classes = [c for c in CLASS_ORDER if any(s.cls == c for s in self.sheets)]

    def __len__(self) -> int:
        return len(self.sheets)

    def of_class(self, cls: str) -> list[Sheet]:
        return [s for s in self.sheets if s.cls == cls]

    def submesh(self, classes) -> pv.PolyData:
        """Siatka zlozona z wybranych klas tkanki (pusta, gdy nic nie wybrano)."""
        want = {s.index for s in self.sheets if s.cls in set(classes)}
        if not want:
            return pv.PolyData()
        keep = np.isin(self._fid_cells, list(want))
        return self.mesh.extract_cells(np.flatnonzero(keep)).extract_surface()

    def sheet_mesh(self, sheet: Sheet) -> pv.PolyData:
        return pv.read(READY / sheet.file).triangulate()

    def at_position(self, pos, max_dist: float = 0.04, classes=None):
        """Najblizszy plat lacznotkankowy (albo None, gdy zaden nie jest w zasiegu)."""
        pos = np.asarray(pos, dtype=float)
        if classes is None:
            mesh, idx_map = self.mesh, self._fid
        else:
            want = {s.index for s in self.sheets if s.cls in set(classes)}
            sel = np.isin(self._fid, list(want))
            if not sel.any():
                return None, None
            pts = np.flatnonzero(sel)
            mesh, idx_map = pv.PolyData(self.mesh.points[pts]), self._fid[pts]
        i = mesh.find_closest_point(pos)
        if i < 0:
            return None, None
        d = float(np.linalg.norm(np.asarray(mesh.points)[i] - pos))
        if d > max_dist:
            return None, None
        return self.sheets[int(idx_map[i])], d

    def counts(self) -> dict:
        out = {}
        for s in self.sheets:
            out[s.cls] = out.get(s.cls, 0) + 1
        return out
