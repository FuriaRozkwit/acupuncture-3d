"""Edycja kotwicy punktu - logika bez okna (gizmo i panel wlasciwosci korzystaja z niej).

Punkt nie ma wspolrzednych XYZ, tylko kotwice proporcjonalna. Uchwyt gizma nie moze
wiec przesuwac punktu "w X": przesuwa go wzdluz STOPNIA SWOBODY kotwicy - `cun`
wzdluz konczyny, `az` wokol niej, `lat` od linii posrodkowej - i zapisuje zmiane
w jednostkach tej kotwicy. Wymiana modelu ciala nadal niczego nie uniewaznia.

Kierunek kazdego uchwytu liczymy roznicowo: dokad przesunie sie punkt, gdy pole
wzrosnie o krok. Jeden mechanizm dziala dzieki temu dla kazdego rodzaju kotwicy,
a strzalka pokazuje dokladnie to, co zrobi przeciagniecie.
"""
from __future__ import annotations

import copy
import json
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from . import anchors, profiles

ROOT = Path(__file__).resolve().parent.parent
INTERVALS = ROOT / "data" / "cun_intervals.json"

# prefiks zrodla, ktory chroni recznie ustawione pole przed tools/fit_azimuth.py
EDIT_TAG = "edycja:"
# pola, ktore nie moga zejsc ponizej zera (odleglosc od konca odcinka, ulamek)
CLAMPED = ("cun", "along", "frac")
UNIT = ("along", "frac")                   # ...a te dodatkowo nie powyzej 1


@dataclass(frozen=True)
class Dof:
    """Jeden stopien swobody kotwicy: sciezka pola, krok i opis dla czlowieka."""
    path: tuple            # np. ("cun",) albo ("shift", "up")
    step: float            # jeden klik [+]/[-]
    unit: str
    label: str

    @property
    def name(self) -> str:
        return ".".join(self.path)


def dofs(anchor: dict) -> list:
    """Stopnie swobody, ktore ma sens ruszac w tej kotwicy (w kolejnosci uchwytow)."""
    a = anchor
    out = []
    if "struct" in a:
        if "along" in a:
            out.append(Dof(("along",), 0.02, "", "wzdluz kosci"))
        if "frac" in a:
            out.append(Dof(("frac",), 0.02, "", "miedzy kosciami"))
        for k, lab in (("up", "w gore"), ("posterior", "do tylu"), ("anterior", "do przodu"),
                       ("medial", "przysrodkowo"), ("lateral", "bocznie")):
            if k in a.get("shift", {}):
                out.append(Dof(("shift", k), 1.0, "mm", lab))
        return out
    if "arc_cun" in a:
        out.append(Dof(("arc_cun",), 0.1, "cun", "wzdluz luku czaszki"))
        out.append(Dof(("lat",), 0.1, "cun", "od linii posrodkowej"))
        return out
    # poziom: jedno z pol wysokosci, zaleznie od stylu kotwicy
    for key, step, lab in (("cun", 0.1, "wzdluz segmentu"), ("z_cun", 0.1, "w pionie"),
                           ("dz_cun", 0.1, "w pionie od kregu"), ("z_frac", 0.002, "w pionie")):
        if key in a:
            out.append(Dof((key,), step, "cun" if key != "z_frac" else "", lab))
            break
    # obwod: krawedz kosci, azymut albo odstep boczny
    if "edge" in a:
        out.append(Dof(("edge_cun",), 0.1, "cun", "od krawedzi kosci"))
    elif "az" in a:
        out.append(Dof(("az",), 2.0, "st.", "wokol segmentu"))
    elif "lat" in a:
        out.append(Dof(("lat",), 0.1, "cun", "od linii posrodkowej"))
    return out


def get(anchor: dict, dof: Dof) -> float:
    v = anchor
    for k in dof.path:
        v = v.get(k, 0.0) if isinstance(v, dict) else 0.0
    return float(v)


def nudge(anchor: dict, dof: Dof, delta: float) -> dict:
    """Nowa kotwica z polem zmienionym o `delta` (w jednostkach pola). Stara zostaje."""
    a = copy.deepcopy(anchor)
    node = a
    for k in dof.path[:-1]:
        node = node.setdefault(k, {})
    last = dof.path[-1]
    val = float(node.get(last, 0.0)) + delta
    if last in CLAMPED:
        val = max(val, 0.0)                      # nie wychodzimy poza koniec odcinka
    if last in UNIT:
        val = min(val, 1.0)                      # ulamek dlugosci kosci
    node[last] = round(val, 4)
    return a


class AnchorEditor:
    """Stan edycji jednego punktu: kotwica przed i po, pozycja, uchwyty, odstepy."""

    def __init__(self, atlas, code: str, side: str = "R"):
        self.atlas = atlas
        self.code = code
        self.side = side if side in ("R", "L", "C") else "R"
        self.point = next(p for p in atlas.points if p["code"] == code)
        self.original = copy.deepcopy(self.point["anchor"])     # to, co jest w uzyciu
        self.anchor = copy.deepcopy(self.original)
        fp = profiles.factory_point(code)                        # nigdy nie nadpisywana
        self.factory = fp["anchor"]
        self.factory_source = fp.get("az_source", "")
        self.dofs = dofs(self.anchor)
        if self.point.get("side") == "midline":
            # punkt linii posrodkowej ma X wymuszone na 0 - odstep boczny nic nie robi
            self.dofs = [d for d in self.dofs if d.path[-1] != "lat"]
        self._canon = None

    # --- geometria ----------------------------------------------------------
    def position(self, anchor: dict | None = None) -> np.ndarray:
        """Pozycja markera dla kotwicy - ta sama droga co w Atlas._resolve_markers."""
        a = self.anchor if anchor is None else anchor
        raw = np.asarray(anchors.resolve(self.atlas.body, a), float)
        ref = self.atlas._ref_point(a, raw).copy()
        midline = self.point.get("side") == "midline"
        if midline:
            ref[0] = 0.0
        pos = np.asarray(self.atlas.place(np.array([ref]), np.array([raw]))[0], float)
        if midline:
            pos[0] = 0.0
        elif self.side == "L":
            pos = pos * np.array([-1.0, 1.0, 1.0])
        return pos

    def handles(self) -> list:
        """Dla kazdego stopnia swobody: (dof, kierunek jednostkowy, mm na jeden krok).

        Roznica centralna: pozycja przy polu +krok minus pozycja przy -krok. Gdy pole
        jest na granicy (np. cun = 0), liczymy roznice jednostronna.
        """
        out = []
        here = self.position()
        for d in self.dofs:
            try:
                v = get(self.anchor, d)
                low = d.path[-1] in CLAMPED and v - d.step < 0
                high = d.path[-1] in UNIT and v + d.step > 1.0
                if high:                                 # na gornej granicy: tylko w dol
                    vec = here - self.position(nudge(self.anchor, d, -d.step))
                elif low:                                # na dolnej: tylko w gore
                    vec = self.position(nudge(self.anchor, d, d.step)) - here
                else:
                    plus = self.position(nudge(self.anchor, d, d.step))
                    minus = self.position(nudge(self.anchor, d, -d.step))
                    vec = (plus - minus) / 2.0
            except Exception:
                continue
            mm = float(np.linalg.norm(vec)) * 1000
            if mm < 1e-3:
                continue
            out.append((d, vec / np.linalg.norm(vec), mm))
        return out

    # --- zmiana -------------------------------------------------------------
    def set_value(self, dof: Dof, value: float) -> None:
        self.anchor = nudge(self.anchor, dof, value - get(self.anchor, dof))

    def step(self, dof: Dof, n: float) -> None:
        self.anchor = nudge(self.anchor, dof, dof.step * n)

    def reset(self) -> None:
        """Cofnij: wraca do kotwicy, ktora byla w uzyciu przed ta edycja."""
        self.anchor = copy.deepcopy(self.original)

    def to_factory(self) -> None:
        """Ustawia kotwice fabryczna. Zapis usunie wtedy nadpisanie z profilu."""
        self.anchor = copy.deepcopy(self.factory)

    @property
    def profile(self) -> str:
        return getattr(self.atlas, "profile", {}).get("name", profiles.FACTORY)

    @property
    def is_factory(self) -> bool:
        """Czy biezaca kotwica (z niezapisanymi zmianami) to fabryczna."""
        return self.anchor == self.factory

    @property
    def dirty(self) -> bool:
        return self.anchor != self.original

    # --- odstepy do sasiadow w kanale --------------------------------------
    def neighbour_intervals(self) -> list:
        """Odstep do poprzedniego i nastepnego punktu kanalu vs kanon z cun_intervals.

        Pozycja edytowanego punktu jest brana z BIEZACEJ kotwicy, wiec przy
        przeciaganiu widac na zywo, czy odstep zbliza sie do kanonu, czy oddala.
        """
        if self._canon is None:
            self._canon = {}
            if INTERVALS.exists():
                doc = json.loads(INTERVALS.read_text(encoding="utf-8"))
                for entries in doc["channels"].values():
                    for e in entries:
                        self._canon[(e["from"], e["to"])] = e
        tools = str(ROOT / "tools")
        if tools not in sys.path:
            sys.path.insert(0, tools)
        from verify_intervals import Intervals  # noqa: E402

        me, here = self.code, self.position()
        if self.side == "L":
            here = here * np.array([-1.0, 1.0, 1.0])     # odstepy mierzymy po prawej
        editor = self

        class _Live(Intervals):
            def pos(self, code):
                return here if code == me else super().pos(code)

        live = _Live(self.atlas)
        live.anchor[me] = editor.anchor                 # linijka i segment z nowej kotwicy
        out = []
        for (a, b), e in self._canon.items():
            if me not in (a, b):
                continue
            got, how = live.measure(a, b)
            out.append({"from": a, "to": b, "canon": e["cun"], "basis": e["basis"],
                        "model": got, "how": how})
        return sorted(out, key=lambda r: (r["to"] != me, r["from"]))

    # --- zapis ---------------------------------------------------------------
    def save(self, note: str = "") -> dict:
        """Zapisuje kotwice do PROFILU. Fabryczny data/points.json zostaje nietkniety.

        Gdy w uzyciu jest fabryczna, profiles.save_point zaklada profil "moja"
        i robi go domyslnym - edytowana wersja wlacza sie przy nastepnym starcie.
        Kotwica rowna fabrycznej usuwa nadpisanie z profilu (powrot do fabryki).
        Pole zrodla dostaje prefiks EDIT_TAG, ktory chroni je przed fit_azimuth.
        """
        changes = []
        for d in self.dofs:
            before, after = get(self.original, d), get(self.anchor, d)
            if abs(before - after) > 1e-9:
                changes.append(f"{d.name} {before:g} -> {after:g}")
        if self.is_factory:
            src = self.factory_source
        else:
            src = f"{EDIT_TAG} {', '.join(changes) or 'bez zmian'}"
            if note:
                src += f" - {note}"
            prev = self.point.get("az_source", "")
            if prev and changes:
                src += f" (wczesniej: {prev})"
            elif prev:
                src = prev
        name = profiles.save_point(self.profile, self.code, self.anchor, src)
        prof = getattr(self.atlas, "profile", None)
        if isinstance(prof, dict):
            prof["name"] = name                   # z fabrycznej mogl powstac "moja"
            changed = set(prof.get("changed", []))
            changed.discard(self.code) if self.is_factory else changed.add(self.code)
            prof["changed"] = sorted(changed)
        self.point["anchor"] = copy.deepcopy(self.anchor)
        if src:
            self.point["az_source"] = src
        self.original = copy.deepcopy(self.anchor)
        return {"code": self.code, "profile": name, "source": src,
                "factory": self.is_factory}
