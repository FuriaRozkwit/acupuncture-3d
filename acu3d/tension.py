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
        self.anchor = None          # {seg, cun, az} - kotwica w ukladzie atlasu
        self.ts = ts or datetime.now(timezone.utc).isoformat(timespec="seconds")

    def to_dict(self) -> dict:
        d = {"id": self.id, "pos": [round(float(x), 5) for x in self.pos],
             "level": self.level, "radius_m": round(self.radius, 4),
             "muscle": self.muscle, "along": self.along, "note": self.note, "ts": self.ts}
        if getattr(self, "anchor", None):
            d["anchor"] = self.anchor
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

    def __init__(self, atlas, muscles=None, name: str = None, fascia=None, wire=None):
        self.atlas = atlas
        self.muscles = muscles
        # siatka cunowa zastepuje siatki anatomiczne tam, gdzie ich nie ma:
        # ma wypieczone to samo (miesien pod spodem, grubosc tkanki), tylko w liczbach
        self.wire = wire
        # warstwa lacznotkankowa jest opcjonalna, ale to wlasnie ona jest tym,
        # co palec wyczuwa jako "napiecie" - brzusiec miesnia lezy glebiej
        self.fascia = fascia
        self.name = name or datetime.now().strftime("sesja-%Y-%m-%d-%H%M")
        self.marks: list = []
        self._next_id = 1
        # linia sily to lancuch oznaczen w kolejnosci stawiania - tak, jak palec
        # prowadzi po ciele: od miejsca, ktore ciagnie, do miejsca, ktore boli
        self.chains: list = []

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
        anchor = None
        if muscle_name is None and self.wire is not None:
            d = self.wire.describe(pos)
            if d["muscle"] is not None:
                m = d["muscle"]
                muscle_name = m["name"] + ("." + m["side"].lower()
                                           if m["side"] in ("L", "R") else "")
            anchor = {"seg": d["segment"], "cun": round(d["cun"], 2),
                      "az": round(d["az"], 1), "flesh_mm": round(d["flesh_mm"], 1)}
        mark = Mark(pos, level, direction, radius, muscle_name, along, note,
                    mid=self._next_id)
        mark.anchor = anchor
        self._next_id += 1
        self.marks.append(mark)
        if self.chains:
            self.chains[-1].append(mark.id)
        else:
            self.chains.append([mark.id])
        return mark

    def break_chain(self) -> None:
        """Konczy biezaca linie sily; kolejne oznaczenia zaczna nowa."""
        if self.chains and self.chains[-1]:
            self.chains.append([])

    def chain_marks(self) -> list:
        """Lancuchy jako listy oznaczen, z pominieciem pustych i jednoelementowych."""
        by_id = {m.id: m for m in self.marks}
        out = []
        for ch in self.chains:
            got = [by_id[i] for i in ch if i in by_id]
            if len(got) > 1:
                out.append(got)
        return out

    def remove_last(self):
        if not self.marks:
            return None
        m = self.marks.pop()
        for ch in self.chains:
            if m.id in ch:
                ch.remove(m.id)
        while self.chains and not self.chains[-1]:
            self.chains.pop()
        return m

    def clear(self) -> None:
        self.marks = []
        self.chains = []

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
               "marks": [m.to_dict() for m in self.marks],
               "chains": [list(c) for c in self.chains if c]}
        p = self.path()
        p.write_text(json.dumps(doc, ensure_ascii=False, indent=1), encoding="utf-8")
        return p

    def load(self, path: Path) -> int:
        doc = json.loads(Path(path).read_text(encoding="utf-8"))
        self.name = doc.get("session", self.name)
        self.marks = [Mark.from_dict(d) for d in doc.get("marks", [])]
        self._next_id = max((m.id for m in self.marks), default=0) + 1
        self.chains = [list(c) for c in doc.get("chains", [])]
        if not self.chains and self.marks:            # starsze zapisy: jeden lancuch
            self.chains = [[m.id for m in self.marks]]
        return len(self.marks)

    # --- analiza ---------------------------------------------------------
    def analyse(self, mark: Mark, max_points: int = 6) -> dict:
        """Co lezy pod oznaczeniem i w jego zasiegu. Same fakty z modelu i zrodla."""
        out = {"mark": mark, "muscle": None, "muscle_dist_mm": None, "points": [],
               "channels": {}, "fascia": None, "fascia_dist_mm": None,
               "wire": None}
        if self.muscles is not None:
            m, dist = self.muscles.muscle_under(mark.pos, max_dist=0.06)
            if m is not None:
                out["muscle"] = m
                out["muscle_dist_mm"] = round(dist * 1000, 1)
                if mark.direction is not None:
                    cos = float(np.dot(_unit(mark.direction), m.axis))
                    out["dir_vs_fiber_deg"] = round(float(np.degrees(np.arccos(
                        np.clip(abs(cos), -1, 1)))), 1)
        if self.muscles is None and self.wire is not None:
            d = self.wire.describe(mark.pos)
            out["wire"] = d                       # segment, cun, azymut, grubosc tkanki
            if d["muscle"] is not None:
                out["muscle_dist_mm"] = round(d["muscle_mm"], 1)
                if mark.direction is not None:
                    cos = float(np.dot(_unit(mark.direction), d["muscle_axis"]))
                    out["dir_vs_fiber_deg"] = round(float(np.degrees(np.arccos(
                        np.clip(abs(cos), -1, 1)))), 1)
        if self.fascia is not None:
            sheet, fdist = self.fascia.at_position(mark.pos, max_dist=0.035)
            if sheet is not None:
                out["fascia"] = sheet
                out["fascia_dist_mm"] = round(fdist * 1000, 1)
                if mark.direction is not None:
                    # kat do NORMALNEJ platu: 90 st. = ciagniecie wzdluz powiezi,
                    # 0 st. = w poprzek, czyli napiecie ja rozciaga na przestrzal
                    cosn = float(np.dot(_unit(mark.direction), _unit(sheet.normal)))
                    out["dir_vs_sheet_deg"] = round(float(np.degrees(np.arccos(
                        np.clip(abs(cosn), -1, 1)))), 1)
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


    # --- linia sily -------------------------------------------------------
    def analyse_line(self, a: "Mark", b: "Mark", corridor: float = 0.03) -> dict:
        """Co lezy na odcinku miedzy dwoma oznaczeniami.

        Nie zgaduje mechanizmu. Podaje trzy mierzalne rzeczy: dlugosc i zwrot
        odcinka, punkty akupunkturowe w korytarzu wzdluz niego, oraz miesnie,
        ktorych OS PRACY jest z tym odcinkiem zgodna - bo to one moga ciagnac
        wzdluz tej linii, a nie te, ktore leza w poprzek.
        """
        v = b.pos - a.pos
        length = float(np.linalg.norm(v))
        u = _unit(v)
        out = {"a": a, "b": b, "length_mm": round(length * 1000, 1),
               "dir": u, "points": [], "muscles": []}

        for k in self.atlas.markers:                 # punkty w korytarzu
            w = k["pos"] - a.pos
            t = float(np.dot(w, u))
            if not (-0.01 <= t <= length + 0.01):
                continue
            d = float(np.linalg.norm(w - u * t))
            if d <= corridor:
                out["points"].append({"code": k["code"], "side": k["side"],
                                      "channel": k["channel"], "off_mm": round(d * 1000, 1),
                                      "at_mm": round(t * 1000, 1),
                                      "pinyin": k.get("pinyin", "")})
        out["points"].sort(key=lambda q: q["at_mm"])

        if self.wire is not None:                    # miesnie zgodne z osia odcinka
            idx = self.wire.along(a.pos, b.pos, radius=corridor)
            seen = {}
            for i in idx:
                mid = int(self.wire.muscle_id[i])
                if mid < 0 or mid >= len(self.wire.muscles):
                    continue
                cos = abs(float(np.dot(u, self.wire.muscle_axis[mid])))
                deg = float(np.degrees(np.arccos(np.clip(cos, -1, 1))))
                prev = seen.get(mid)
                if prev is None or deg < prev[0]:
                    seen[mid] = (deg, float(self.wire.muscle_mm[i]))
            for mid, (deg, mm) in sorted(seen.items(), key=lambda kv: kv[1][0]):
                m = self.wire.muscles[mid]
                out["muscles"].append({"name": m["name"], "side": m["side"],
                                       "angle_deg": round(deg, 1), "depth_mm": round(mm, 1)})
        return out

    def line_report(self, corridor: float = 0.03) -> str:
        """Zestawienie linii sil: co lezy wzdluz kazdej i co moze ja ciagnac."""
        chains = self.chain_marks()
        if not chains:
            return "Nie ma zadnej linii sily - potrzebne sa co najmniej dwa oznaczenia."
        lines = [f"LINIE SIL: {self.name}   lancuchow: {len(chains)}", "=" * 78]
        for n, ch in enumerate(chains, 1):
            total = sum(float(np.linalg.norm(b.pos - a.pos)) for a, b in zip(ch, ch[1:]))
            lines.append(f"[{n}] {len(ch)} oznaczen, dlugosc lacznie "
                         f"{total * 1000:.0f} mm, poziomy "
                         f"{'-'.join(str(m.level) for m in ch)}")
            for a, b in zip(ch, ch[1:]):
                r = self.analyse_line(a, b, corridor)
                lines.append(f"    [{a.id}] -> [{b.id}]  {r['length_mm']:.0f} mm")
                if r["points"]:
                    lines.append("       punkty na linii: " + ", ".join(
                        f"{q['code']} ({q['off_mm']:.0f} mm od osi)"
                        for q in r["points"][:8]))
                for m in r["muscles"][:3]:
                    side = {"L": "lewy", "R": "prawy", "C": ""}.get(m["side"], "")
                    lines.append(f"       wzdluz linii: {m['name']}"
                                 f"{(' (' + side + ')') if side else ''}"
                                 f"  os pod {m['angle_deg']} st., glebokosc {m['depth_mm']:.0f} mm")
            lines.append("")
        lines.append("Kat osi miesnia mowi, czy miesien MOZE ciagnac wzdluz tej linii:")
        lines.append("0 stopni = os rownolegla do linii, 90 = w poprzek. To geometria,")
        lines.append("nie rozpoznanie - mechanizm trzeba sprawdzic badaniem.")
        return "\n".join(lines)

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
            elif a["wire"] is not None:
                w = a["wire"]
                head += (f"   {w['segment']} {w['cun']:.1f} cun az {w['az']:.0f}"
                         f"   tkanka {w['flesh_mm']:.0f} mm")
                if mark.muscle:
                    head += f"   miesien: {mark.muscle} ({a['muscle_dist_mm']} mm)"
                    by_muscle[mark.muscle] = by_muscle.get(mark.muscle, 0) + mark.level
            if "dir_vs_fiber_deg" in a:
                head += f"   kierunek vs os miesnia: {a['dir_vs_fiber_deg']} st."
            lines.append(head)
            if a["fascia"] is not None:
                f = f"     tkanka laczna: {a['fascia'].label} ({a['fascia_dist_mm']} mm)"
                if "dir_vs_sheet_deg" in a:
                    f += f"   kierunek vs plat: {a['dir_vs_sheet_deg']} st."
                lines.append(f)
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
