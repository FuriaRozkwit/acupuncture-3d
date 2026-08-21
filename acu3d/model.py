"""Spina cialo (anatomiczne albo proceduralne), dane punktow i meridiany.

Tryb `anatomy` uzywa siatek z Z-Anatomy i BodyParts3D oraz rigu wyliczonego
z kosci (data/rig_fitted.json). Tryb `procedural` to poprzedni manekin - zostaje
jako zapasowy, gdy nie ma pobranych zasobow.
"""
from __future__ import annotations

import json
import warnings
from pathlib import Path

import numpy as np
import pyvista as pv

from . import anchors
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

        doc = json.loads((DATA / "points.json").read_text(encoding="utf-8"))
        self.channels = doc["channels"]
        self.points = doc["points"]
        meta_path = DATA / "point_meta.json"        # opisy z podrecznika - lokalne,
        self.meta = (json.loads(meta_path.read_text(encoding="utf-8"))   # nie w repozytorium
                     if meta_path.exists() else {})
        self.paths = json.loads((DATA / "meridians.json").read_text(encoding="utf-8"))["paths"]
        self.markers = self._resolve_markers()
        self.by_key = {m["key"]: m for m in self.markers}

    # --- rozwiazywanie punktow ------------------------------------------
    def _ref_point(self, anchor: dict, pos: np.ndarray) -> np.ndarray:
        """Punkt wewnatrz ciala, z ktorego 'patrzymy' na punkt."""
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
        d = pos - refs
        n = np.linalg.norm(d, axis=1, keepdims=True)
        u = d / np.where(n < 1e-9, 1e-9, n)
        return pos + u * SKIN_OFFSET          # w trybie anatomii punkt juz lezy na skorze

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
    def meridian_lines(self, channel: str, side: str = "R") -> list:
        """Polilinie przebiegu kanalu, oblepione na powierzchni ciala."""
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
            dense_p, dense_r = [], []
            for i in range(len(pts) - 1):
                for f in np.linspace(0.0, 1.0, 7, endpoint=False):
                    dense_p.append(pts[i] * (1 - f) + pts[i + 1] * f)
                    dense_r.append(refs[i] * (1 - f) + refs[i + 1] * f)
            dense_p.append(pts[-1])
            dense_r.append(refs[-1])
            lines.append(self._hug(np.array(dense_r), np.array(dense_p)))
        return lines

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
