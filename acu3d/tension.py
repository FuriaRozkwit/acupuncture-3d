"""Mapa napiec: znaczniki na manekinie, ich zapis i analiza.

Zalozenia
---------
1. Znacznik jest KOTWICZONY ANATOMICZNIE, a nie we wspolrzednych. Zapisujemy
   miesien i polozenie wzgledem jego osi, wiec ten sam zapis da sie nalozyc na
   inny model i na inna sylwetke - dokladnie tak, jak punkty w tym projekcie.
2. Znacznik ma KIERUNEK. Napiecie ma zwrot (ciagnie, ogranicza ruch w osi),
   wiec sam punkt to za malo.
3. Analiza jest SPRAWDZALNA. Nie zgaduje przyczyn: podaje, co lezy pod znacznikiem
   i w jego zasiegu, z odlegloscia w milimetrach i z cytatem ze zrodla. Kazde
   zdanie da sie sprawdzic na modelu.

Analiza jest regulowa, liczona z danych - nie jest modelem jezykowym i niczego
nie wnioskuje poza tym, co da sie zmierzyc. Warstwa jezykowa moze to podsumowac
pozniej, ale wtedy dane opuszczaja komputer i to jest osobna decyzja.
"""
from __future__ import annotations

import json
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
SESSIONS = ROOT / "data" / "tension"

LEVELS = {0: "bez napiecia", 1: "lekkie", 2: "wyrazne", 3: "silne, bolesne"}
LEVEL_COLORS = {0: "#5fb98f", 1: "#d8c447", 2: "#e08a34", 3: "#cf3b3b"}
DEFAULT_RADIUS = 0.035          # metry - zasieg pojedynczego oznaczenia


class Mark:
    """Pojedyncze oznaczenie napiecia."""

    def __init__(self, pos, level=2, direction=None, radius=DEFAULT_RADIUS,
                 muscle=None, along=None, note="", ts=None, mid=0):
        self.id = mid
        self.pos = np.asarray(pos, dtype=float)
        self.level = int(level)
        self.direction = None if direction is None else np.asarray(direction, dtype=float)
        self.radius = float(radius)
        self.muscle = muscle
        self.along = along
        self.note = note
        self.ts = ts or datetime.now(timezone.utc).isoformat(timespec="seconds")

    def to_dict(self) -> dict:
        d = {"id": self.id, "pos": [round(float(x), 5) for x in self.pos],
             "level": self.level, "radius_m": round(self.radius, 4),
             "muscle": self.muscle, "along": self.along, "note": self.note, "ts": self.ts}
        if self.direction is not None:
            d["direction"] = [round(float(x), 4) for x in self.direction]
        return d

    @classmethod
    def from_dict(cls, d: dict) -> "Mark":
        return cls(d["pos"], d.get("level", 2), d.get("direction"),
                   d.get("radius_m", DEFAULT_RADIUS), d.get("muscle"), d.get("along"),
                   d.get("note", ""), d.get("ts"), d.get("id", 0))


class TensionMap:
    """Zbior oznaczen jednej sesji + pole napiecia rozlozone na skorze."""

    def __init__(self, atlas, muscles=None, name: str = None):
        self.atlas = atlas
        self.muscles = muscles
        self.name = name or datetime.now().strftime("sesja-%Y-%m-%d-%H%M")
        self.marks: list = []
        self._next_id = 1

    # --- zbieranie oznaczen ---------------------------------------------
    def add(self, pos, level: int = 2, direction=None, radius: float = DEFAULT_RADIUS,
            note: str = "") -> Mark:
        muscle_name, along = None, None
        if self.muscles is not None:
            m, dist = self.muscles.muscle_under(np.asarray(pos), max_dist=0.06)
            if m is not None:
                muscle_name = m.name + ("." + m.side.lower() if m.side in ("L", "R") else "")
                rel = np.asarray(pos) - m.centroid
                lo, hi = m.extent
                along = round(float((rel @ m.axis - lo) / max(hi - lo, 1e-6)), 3)
        mark = Mark(pos, level, direction, radius, muscle_name, along, note,
                    mid=self._next_id)
        self._next_id += 1
        self.marks.append(mark)
        return mark

    def remove_last(self):
        return self.marks.pop() if self.marks else None

    def clear(self) -> None:
        self.marks = []

    # --- pole napiecia ---------------------------------------------------
    def field(self, points: np.ndarray) -> np.ndarray:
        """Poziom napiecia w kazdym punkcie powierzchni.

        Kazde oznaczenie rozmywa sie gaussowsko na swoj promien, a pola sie sumuja -
        dwa sasiadujace oznaczenia daja jedno wieksze ognisko, a nie dwie kropki.
        """
        out = np.zeros(len(points), dtype=float)
        for m in self.marks:
            d2 = np.sum((points - m.pos) ** 2, axis=1)
            out += m.level * np.exp(-d2 / (2.0 * (m.radius / 2.0) ** 2))
        return out

    # --- zapis i odczyt ---------------------------------------------------
    def path(self) -> Path:
        return SESSIONS / f"{self.name}.json"

    def save(self) -> Path:
        SESSIONS.mkdir(parents=True, exist_ok=True)
        doc = {"session": self.name,
               "saved": datetime.now(timezone.utc).isoformat(timespec="seconds"),
               "model_height_m": round(float(self.atlas.rig.H), 3),
               "note": "znaczniki kotwiczone anatomicznie: miesien + polozenie na jego osi",
               "marks": [m.to_dict() for m in self.marks]}
        p = self.path()
        p.write_text(json.dumps(doc, ensure_ascii=False, indent=1), encoding="utf-8")
        return p

    def load(self, path: Path) -> int:
        doc = json.loads(Path(path).read_text(encoding="utf-8"))
        self.name = doc.get("session", self.name)
        self.marks = [Mark.from_dict(d) for d in doc.get("marks", [])]
        self._next_id = max((m.id for m in self.marks), default=0) + 1
        return len(self.marks)

    # --- analiza ---------------------------------------------------------
    def analyse(self, mark: Mark, max_points: int = 6) -> dict:
        """Co lezy pod oznaczeniem i w jego zasiegu. Same fakty z modelu i zrodla."""
        out = {"mark": mark, "muscle": None, "muscle_dist_mm": None, "points": [],
               "channels": {}}
        if self.muscles is not None:
            m, dist = self.muscles.muscle_under(mark.pos, max_dist=0.06)
            if m is not None:
                out["muscle"] = m
                out["muscle_dist_mm"] = round(dist * 1000, 1)
                if mark.direction is not None:
                    cos = float(np.dot(_unit(mark.direction), m.axis))
                    out["dir_vs_fiber_deg"] = round(float(np.degrees(np.arccos(
                        np.clip(abs(cos), -1, 1)))), 1)
        near = []
        for k in self.atlas.markers:
            d = float(np.linalg.norm(k["pos"] - mark.pos))
            if d <= mark.radius * 1.6:
                near.append((d, k))
        near.sort(key=lambda x: x[0])
        for d, k in near[:max_points]:
            meta = self.atlas.meta.get(k["code"], {})
            actions = meta.get("actions") or []
            if isinstance(actions, list):
                actions = " ".join(str(a) for a in actions if a)
            out["points"].append({
                "code": k["code"], "side": k["side"], "dist_mm": round(d * 1000, 1),
                "pinyin": k.get("pinyin", ""), "channel": k["channel"],
                "actions": (actions or "")[:220],
                "source": "Deadman" if actions else "",
            })
            out["channels"][k["channel"]] = out["channels"].get(k["channel"], 0) + 1
        return out

    def report(self) -> str:
        """Zestawienie calej mapy - kazda pozycja z odlegloscia i zrodlem."""
        if not self.marks:
            return "Mapa napiec jest pusta."
        lines = [f"MAPA NAPIEC: {self.name}   oznaczen: {len(self.marks)}",
                 "=" * 78]
        by_channel: dict = {}
        by_muscle: dict = {}
        for mark in self.marks:
            a = self.analyse(mark)
            head = f"[{mark.id}] poziom {mark.level} ({LEVELS.get(mark.level, '?')})"
            if a["muscle"] is not None:
                head += f"   miesien: {a['muscle'].label} ({a['muscle_dist_mm']} mm)"
                by_muscle[a["muscle"].label] = by_muscle.get(a["muscle"].label, 0) + mark.level
            if "dir_vs_fiber_deg" in a:
                head += f"   kierunek vs wlokna: {a['dir_vs_fiber_deg']} st."
            lines.append(head)
            if mark.note:
                lines.append(f"     notatka: {mark.note}")
            for p in a["points"]:
                lines.append(f"     {p['code']:8s} {p['pinyin']:<12s} {p['dist_mm']:5.0f} mm"
                             f"   [{p['channel']}]")
                if p["actions"]:
                    lines.append(f"        dzialanie wg zrodla: {p['actions']}")
            for ch, n in a["channels"].items():
                by_channel[ch] = by_channel.get(ch, 0) + n * mark.level
            lines.append("")
        if by_muscle:
            lines.append("SUMA WG MIESNI (poziom napiecia narastajaco):")
            for name, score in sorted(by_muscle.items(), key=lambda kv: -kv[1]):
                lines.append(f"   {score:4d}  {name}")
        if by_channel:
            lines.append("")
            lines.append("SUMA WG KANALOW (ile punktow kanalu lezy w polu napiecia):")
            for ch, score in sorted(by_channel.items(), key=lambda kv: -kv[1]):
                info = self.atlas.channels[ch]
                lines.append(f"   {score:4d}  {ch:4s} {info['pl']}")
        lines.append("")
        lines.append("Zestawienie jest wyliczone z geometrii modelu i opisow zrodlowych.")
        lines.append("Nie jest rozpoznaniem ani zaleceniem - to spis tego, co lezy pod")
        lines.append("oznaczeniami, z odlegloscia do sprawdzenia na modelu.")
        return "\n".join(lines)


def _unit(v: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(v)
    return v / n if n > 1e-12 else np.array([0.0, 0.0, 1.0])
