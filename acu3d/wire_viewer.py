"""Silnik napiec: siatka cunowa, punkty i linie sil. Bez warstw i bez siatek.

Czym to sie rozni od `viewer.py`
--------------------------------
Tamten pokazuje CIALO. Ten pokazuje UKLAD ODNIESIENIA, w ktorym cialo jest
opisane, i sluzy do jednej rzeczy: precyzyjnego stawiania oznaczen napiecia
i laczenia ich w linie sil.

Zadna siatka anatomiczna nie jest tu wczytywana. Wiedza o tkankach jest
wypieczona w skalarach siatki cunowej (`tools/bake_wire.py`), wiec scena ma
9 tys. wierzcholkow zamiast 1,4 mln trojkatow, a mimo to po kazdym klikniecu
wiadomo, ile jest w tym miejscu tkanki miekkiej i ktory miesien lezy pod spodem.

Klikniecie stawia oznaczenie od razu - to jest funkcja glowna, a nie tryb,
do ktorego trzeba wejsc. Oznaczenia stawiane po kolei tworza LANCUCH, rysowany
jako wektor ze strzalka: tak wyglada kierunek ciagu, ktory sie bada.
"""
from __future__ import annotations

import sys

import numpy as np
import pyvista as pv

from .model import Atlas
from .tension import LEVEL_COLORS, LEVELS, TensionMap
from .wire import FLESH_COLORS, FLESH_RANGE, Wire

THEME = {"bg": "#0b0e13", "text": "#dfe6f0", "dim": "#7c879a",
         "point": "#5f6b7d", "sel": "#ffe36a", "force": "#ff8a5c"}


def _safe_print(text: str) -> None:
    enc = getattr(sys.stdout, "encoding", None) or "utf-8"
    sys.stdout.write(text.encode(enc, errors="replace").decode(enc, errors="replace") + "\n")
    sys.stdout.flush()


class WireViewer:
    def __init__(self, atlas: Atlas, off_screen: bool = False, window_size=(1500, 960)):
        self.atlas = atlas
        self.wire = Wire()
        self.tension = TensionMap(atlas, muscles=None, wire=self.wire)
        self.pl = pv.Plotter(off_screen=off_screen, window_size=window_size,
                             title="Silnik napiec - siatka cunowa")
        self.pl.set_background(THEME["bg"])
        self.level = 2
        self.show_points = True
        self.show_labels = False
        self.filtered = False        # filtr wiedzy: tylko punkty lezace na liniach sil
        self.point_actors: dict = {}
        self.label_actors: dict = {}
        self._build()

    # ------------------------------------------------------------------ scena
    def _build(self) -> None:
        self.pl.add_mesh(self.wire.mesh, scalars="tkanka_mm", cmap=FLESH_COLORS,
                         clim=list(FLESH_RANGE), show_scalar_bar=False,
                         line_width=1, opacity=0.55, name="wire")
        # punkty: male i szare, dopoki nie sa potrzebne - siatka ma byc czytelna
        for ch, info in self.atlas.channels.items():
            ms = [m for m in self.atlas.markers if m["channel"] == ch]
            pts = np.array([m["pos"] for m in ms])
            self.point_actors[ch] = self.pl.add_mesh(
                pv.PolyData(pts).glyph(geom=pv.Sphere(radius=0.0045, theta_resolution=10,
                                                      phi_resolution=10),
                                       scale=False, orient=False),
                color=info["color"], opacity=0.85, name=f"pts_{ch}")
            self.label_actors[ch] = self.pl.add_point_labels(
                pts, [m["code"] for m in ms], font_size=10, text_color=THEME["text"],
                shape=None, show_points=False, always_visible=True, name=f"lab_{ch}")
            self.label_actors[ch].SetVisibility(False)

        self._add_help()
        self.info("Klikaj w cialo - kazde klikniecie stawia oznaczenie napiecia.")
        self.pl.enable_point_picking(callback=self._on_pick, show_message=False,
                                     left_clicking=True, show_point=False, tolerance=0.03)
        self._keys()
        self._view("xz", 180)

    def _add_help(self) -> None:
        self.pl.add_text(
            "SILNIK NAPIEC\n"
            "[klik] oznaczenie   [,]/[.] poziom 0-3   [z] cofnij\n"
            "[spacja] koniec linii sily, nastepna zaczyna sie od nowa\n"
            "[q] punkty   [l] etykiety   [f] filtr: tylko punkty na liniach\n"
            "[p] zestawienie linii   [s] zapis\n"
            "[1-4] widoki   [0] domyslny",
            position="upper_right", font_size=8, color=THEME["dim"], name="help")

    def _keys(self) -> None:
        keys = {
            "period": lambda: self.set_level(self.level + 1),
            "comma": lambda: self.set_level(self.level - 1),
            "z": self.undo,
            "space": self.break_chain,
            "q": self.toggle_points,
            "f": self.toggle_filter,
            "l": self.toggle_labels,
            "p": self.report,
            "s": self.save,
            "0": lambda: self._view("xz", 180),
            "1": lambda: self._view("xz", 180),
            "2": lambda: self._view("xz", 0),
            "3": lambda: self._view("yz", 0),
            "4": lambda: self._view("xy", 0),
        }
        for k, fn in keys.items():
            self.pl.add_key_event(k, fn)

    # ---------------------------------------------------------------- sterowanie
    def info(self, text: str) -> None:
        self.pl.add_text(text, position="lower_left", font_size=9,
                         color=THEME["text"], name="info")

    def _view(self, cam: str, azimuth: float) -> None:
        self.pl.camera_position = cam
        if azimuth:
            self.pl.camera.azimuth = azimuth
        self.pl.reset_camera()
        self.pl.camera.zoom(1.35)
        self.pl.render()

    def set_level(self, level: int) -> None:
        self.level = int(np.clip(level, 0, 3))
        self.info(f"Poziom kolejnych oznaczen: {self.level} ({LEVELS[self.level]})")
        self.pl.render()

    def toggle_points(self) -> None:
        self.show_points = not self.show_points
        for a in self.point_actors.values():
            a.SetVisibility(self.show_points)
        self.pl.render()

    def toggle_filter(self) -> None:
        """Filtr wiedzy: z 722 markerow zostawia te, ktore linie sil faktycznie mijaja.

        To jest ten sam korytarz, ktorego uzywa zestawienie - wiec to, co widac,
        i to, co wypisane, jest jednym i tym samym zbiorem.
        """
        self.filtered = not self.filtered
        if not self.filtered:
            for ch, a in self.point_actors.items():
                a.SetVisibility(self.show_points)
            self.pl.remove_actor("filtered", render=False)
            self.info("Filtr wylaczony - wszystkie punkty z powrotem.")
            self.pl.render()
            return
        codes, keys = set(), []
        for chain in self.tension.chain_marks():
            for a, b in zip(chain, chain[1:]):
                for q in self.tension.analyse_line(a, b)["points"]:
                    codes.add(q["code"])
                    keys.append((q["code"], q["side"]))
        if not codes:
            self.filtered = False
            self.info("Nie ma jeszcze zadnej linii sily - filtr nie ma czego zawezic.")
            self.pl.render()
            return
        for a in self.point_actors.values():
            a.SetVisibility(False)
        sel = [m for m in self.atlas.markers if (m["code"], m["side"]) in set(keys)]
        cloud = pv.PolyData(np.array([m["pos"] for m in sel]))
        self.pl.add_mesh(
            cloud.glyph(geom=pv.Sphere(radius=0.008), scale=False, orient=False),
            color=THEME["sel"], name="filtered")
        by_ch: dict = {}
        for m in sel:
            by_ch[m["channel"]] = by_ch.get(m["channel"], 0) + 1
        top = ", ".join(f"{c} ({n})" for c, n in
                        sorted(by_ch.items(), key=lambda kv: -kv[1]))
        self.info(f"FILTR: {len(sel)} punktow na liniach sil.  Kanaly: {top}")
        _safe_print("\npunkty na liniach sil: " + ", ".join(sorted(codes)))
        self.pl.render()

    def toggle_labels(self) -> None:
        self.show_labels = not self.show_labels
        for a in self.label_actors.values():
            a.SetVisibility(self.show_labels)
        self.pl.render()

    def break_chain(self) -> None:
        self.tension.break_chain()
        self.info("Linia sily zamknieta. Nastepne klikniecie zaczyna nowa.")
        self.pl.render()

    def undo(self) -> None:
        m = self.tension.remove_last()
        if m is None:
            return
        self._refresh()
        self.info(f"Cofnieto oznaczenie [{m.id}]. Zostalo: {len(self.tension.marks)}")
        self.pl.render()

    def save(self) -> None:
        if not self.tension.marks:
            return
        path = self.tension.save()
        self.info(f"Zapisano {len(self.tension.marks)} oznaczen: {path.name}")
        _safe_print(f"\nzapisano: {path}")
        self.pl.render()

    def report(self) -> None:
        _safe_print("\n" + self.tension.line_report())
        self.info("Zestawienie linii sil wypisane na konsoli.")
        self.pl.render()

    # ------------------------------------------------------------------ klik
    def _on_pick(self, picked, *args) -> None:
        pos = np.asarray(getattr(picked, "points", picked)).reshape(-1, 3)
        if not len(pos):
            return
        # przyciagniecie do siatki: oznaczenie ma lezec w ukladzie cun, a nie
        # tam, gdzie akurat trafil piksel
        pos = self.wire.snap(pos[0])
        mark = self.tension.add(pos, level=self.level)
        self._refresh()
        a = self.tension.analyse(mark)
        w = a["wire"]
        txt = [f"[{mark.id}] poziom {mark.level} ({LEVELS[mark.level]})   "
               f"{w['segment']}  {w['cun']:.1f} cun  az {w['az']:.0f}",
               f"tkanka miekka {w['flesh_mm']:.0f} mm"
               + (f"   miesien: {mark.muscle} ({w['muscle_mm']:.0f} mm)"
                  if mark.muscle else "")]
        if a["points"]:
            txt.append("punkty w zasiegu: " + ", ".join(
                f"{q['code']} ({q['dist_mm']:.0f} mm)" for q in a["points"][:4]))
        chain = [c for c in self.tension.chain_marks() if mark in c]
        if chain:
            prev = chain[-1][-2]
            r = self.tension.analyse_line(prev, mark)
            head = f"linia [{prev.id}]->[{mark.id}] {r['length_mm']:.0f} mm"
            if r["points"]:
                head += "   na linii: " + ", ".join(q["code"] for q in r["points"][:5])
            if r["muscles"]:
                m0 = r["muscles"][0]
                head += f"   wzdluz: {m0['name']} ({m0['angle_deg']} st.)"
            txt.append(head)
        self.info("\n".join(txt))
        self.pl.render()

    # --------------------------------------------------------------- rysowanie
    def _refresh(self) -> None:
        for name in ("marks", "forces"):
            self.pl.remove_actor(name, render=False)
        marks = self.tension.marks
        if not marks:
            return
        cloud = pv.PolyData(np.array([m.pos for m in marks]))
        cloud["poziom"] = np.array([float(m.level) for m in marks])
        self.pl.add_mesh(
            cloud.glyph(geom=pv.Sphere(radius=0.009), scale=False, orient=False),
            scalars="poziom", clim=[0, 3], show_scalar_bar=False,
            cmap=[LEVEL_COLORS[i] for i in range(4)], name="marks")

        # linia sily: rurka przez lancuch + strzalka na kazdym odcinku
        parts = []
        for chain in self.tension.chain_marks():
            pts = np.array([m.pos for m in chain])
            parts.append(pv.lines_from_points(pts).tube(radius=0.0026))
            for a, b in zip(chain, chain[1:]):
                v = b.pos - a.pos
                n = float(np.linalg.norm(v))
                if n > 1e-6:
                    parts.append(pv.Arrow(start=a.pos + v * 0.5, direction=v / n,
                                          scale=min(0.05, n * 0.5)))
        if parts:
            merged = parts[0].merge(parts[1:]) if len(parts) > 1 else parts[0]
            self.pl.add_mesh(merged, color=THEME["force"], name="forces")

    # ------------------------------------------------------------------- start
    def show(self) -> None:
        self.pl.show()

    def screenshot(self, path: str) -> None:
        self.pl.screenshot(path)
