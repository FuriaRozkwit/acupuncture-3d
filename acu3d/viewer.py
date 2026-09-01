"""Interaktywna przegladarka 3D: warstwy ciala, kanaly, punkty, miesnie.

Wydajnosc suwakow
-----------------
Wybrany miesien rysowany jest JAKO OSOBNY AKTOR, a w wielkiej siatce wszystkich
miesni jego komorki stają się przezroczyste (wpis w tablicy kolorow z alfa = 0).
Dzieki temu ruch suwaka przelicza i wysyla na karte kilka tysiecy wierzcholkow,
a nie pol miliona. Mapowanie punkt-miesien liczone jest raz na miesien i
zapamietywane, bo nie zalezy od parametrow.
"""
from __future__ import annotations

import sys
import textwrap

import numpy as np
import pyvista as pv

from . import action_lines
from .fascia import CLASS_COLOR, CLASS_PL, FasciaSet
from .model import Atlas
from .tension import LEVEL_COLORS, LEVELS, TensionMap
from .muscles import MuscleSet

THEME = {
    "bg": "#0e1116",
    "skin": "#e8cbaa",
    "muscle": "#a8382c",
    "muscle_sel": "#ffcf40",
    "bone": "#efe9dc",
    "text": "#dfe6f0",
    "dim": "#8792a3",
    "accent": "#7fd1c8",
    "highlight": "#ffe36a",
}
LAYER_LABELS = {"skin": "skora", "muscle": "miesnie", "skeleton": "kosci",
                "fascia": "powiezi"}
SIDE_PL = {"R": "prawa", "L": "lewa", "C": "linia posrodkowa"}
PARAMS = (("napiecie", "activation", 0.0, 1.0, 0.0),
          ("tonus spoczynkowy", "tone", 0.0, 0.6, 0.0),
          ("przerost", "hypertrophy", 0.75, 1.35, 1.0),
          ("sztywnosc", "stiffness", 0.2, 1.6, 1.0))
# przezroczystosci w trybie fokusu: reszta ciala schodzi na drugi plan
FOCUS_OPACITY = {"skin": 0.07, "muscle": 0.10, "skeleton": 0.18, "fascia": 0.22}


def _safe_print(text: str) -> None:
    """Wypis odporny na konsole, ktora nie zna znakow spoza swojej strony kodowej."""
    enc = getattr(sys.stdout, "encoding", None) or "utf-8"
    sys.stdout.write(text.encode(enc, errors="replace").decode(enc, errors="replace") + "\n")
    sys.stdout.flush()


def _rgba(hex_color: str, alpha: int = 255) -> list:
    h = hex_color.lstrip("#")
    return [int(h[i:i + 2], 16) for i in (0, 2, 4)] + [alpha]


class Viewer:
    def __init__(self, atlas: Atlas, off_screen: bool = False, window_size=(1560, 980)):
        self.atlas = atlas
        self.muscles = MuscleSet() if atlas.mode == "anatomy" else None
        self.pl = pv.Plotter(off_screen=off_screen, window_size=window_size,
                             title="Atlas punktow akupunktury 3D")
        self.pl.set_background(THEME["bg"])
        self.layer_actors: dict = {}
        self.point_actors: dict = {}
        self.line_actors: dict = {}
        self.label_actors: dict = {}
        self.visible_channels = set(atlas.channels)
        self.show_labels = False
        self.show_meridians = True
        self.solo: str | None = None
        self.selected = None                 # wybrany miesien
        self.focused_key: str | None = None  # wybrany punkt
        self.focus_mode = False
        self._focus_saved: dict | None = None
        self._sel_mesh = None
        self._sel_rest = None
        self._sel_actor = None
        self._sliders: list = []
        self._buttons: list = []
        self.tension = None               # mapa napiec (tworzona przy pierwszym uzyciu)
        self.tension_mode = False
        self.tension_level = 2            # poziom nadawany kolejnym oznaczeniom
        self._await_direction = False
        self.show_action = False          # linie pracy miesni
        self._action_actor = None
        self.fascia = None
        self.fascia_step = 0            # 0 = wylaczona, potem kolejne klasy
        self._muscle_opacity = 1.0        # krycie warstwy miesni sprzed wlaczenia linii
        self._build()

    # ------------------------------------------------------------------ scena
    def _build(self) -> None:
        a = self.atlas
        self.layer_actors["skin"] = self.pl.add_mesh(
            a.layers["skin"], color=THEME["skin"], opacity=1.0, smooth_shading=True,
            specular=0.12, name="skin")
        if self.muscles is not None:
            m = self.muscles.mesh
            m.point_data["sel"] = np.zeros(m.n_points, dtype=np.float32)
            # 0 = miesien widoczny, 1 = zastapiony osobnym aktorem, wiec przezroczysty
            lut = pv.LookupTable(values=np.array([_rgba(THEME["muscle"]), [0, 0, 0, 0]],
                                                 dtype=np.uint8))
            lut.scalar_range = (0.0, 1.0)
            self.layer_actors["muscle"] = self.pl.add_mesh(
                m, scalars="sel", cmap=lut, clim=[0, 1], show_scalar_bar=False,
                smooth_shading=True, specular=0.2, name="muscle")
        else:
            self.layer_actors["muscle"] = self.pl.add_mesh(
                a.layers["muscle"], color=THEME["muscle"], smooth_shading=True, name="muscle")
        self.layer_actors["skeleton"] = self.pl.add_mesh(
            a.layers["skeleton"], color=THEME["bone"], smooth_shading=True, name="skeleton")
        self.layer_actors["fascia"] = self.pl.add_mesh(
            pv.PolyData(), color=CLASS_COLOR["septum"], smooth_shading=True,
            opacity=0.55, name="fascia")
        self.layer_actors["muscle"].SetVisibility(False)
        self.layer_actors["skeleton"].SetVisibility(False)
        self.layer_actors["fascia"].SetVisibility(False)

        r = 0.0058
        for ch, info in a.channels.items():
            ms = [m for m in a.markers if m["channel"] == ch]
            pts = np.array([m["pos"] for m in ms])
            glyph = pv.PolyData(pts).glyph(
                geom=pv.Sphere(radius=r, theta_resolution=14, phi_resolution=14),
                scale=False, orient=False)
            self.point_actors[ch] = self.pl.add_mesh(
                glyph, color=info["color"], smooth_shading=True, name=f"pts_{ch}")
            lines = []
            for side in ("R", "L"):
                for poly in a.meridian_lines(ch, side):
                    if len(poly) > 3:
                        # krzywa jest juz gladka i lezy na skorze - drugi splajn
                        # tylko odciagnalby ja od powierzchni
                        lines.append(pv.lines_from_points(poly).tube(radius=r * 0.33))
                if info["clock"] >= 13:              # REN / DU sa pojedyncze
                    break
            if lines:
                merged = lines[0].merge(lines[1:]) if len(lines) > 1 else lines[0]
                self.line_actors[ch] = self.pl.add_mesh(
                    merged, color=info["color"], opacity=0.65, name=f"mer_{ch}")
            self.label_actors[ch] = self.pl.add_point_labels(
                pts, [m["code"] for m in ms], font_size=11, text_color=THEME["text"],
                shape=None, show_points=False, always_visible=True, name=f"lab_{ch}")
            self.label_actors[ch].SetVisibility(False)

        self._add_widgets()
        self._add_help()
        self.info("Kliknij punkt albo miesien.  [j] warstwa miesni   [f] fokus na wybranym")
        self.pl.enable_point_picking(callback=self._on_pick, show_message=False,
                                     left_clicking=True, show_point=False, tolerance=0.02)
        self._view("xz", 180)

    # -------------------------------------------------------------- widgety
    def _button(self, label: str, fn, position, color=None) -> None:
        """Przycisk chwilowy - wykonuje akcje i sam sie odklika."""
        holder = {}

        def cb(state):
            fn()
            w = holder.get("w")
            if w is not None:
                w.GetRepresentation().SetState(False)
            self.pl.render()

        holder["w"] = self.pl.add_checkbox_button_widget(
            cb, value=False, position=position, size=16, border_size=1,
            color_on=color or THEME["accent"], color_off="#38414f")
        self.pl.add_text(label, position=(position[0] + 21, position[1] + 1),
                         font_size=7, color=THEME["text"], name=f"btn_{label}")
        self._buttons.append(holder["w"])

    def _add_widgets(self) -> None:
        y = 12
        for name in ("skin", "muscle", "skeleton"):
            self.pl.add_checkbox_button_widget(
                self._toggle_layer(name), value=(name == "skin"), position=(12, y),
                size=20, border_size=1, color_on=THEME["skin"], color_off="#38414f")
            self.pl.add_text(LAYER_LABELS[name], position=(38, y + 2), font_size=8,
                             color=THEME["text"], name=f"lt_{name}")
            y += 26

        # kanaly, od dolu do gory
        self._channel_widgets = {}
        x, y0 = 12, 200
        for i, (ch, info) in enumerate(self.atlas.channels.items()):
            y = y0 + 23 * (len(self.atlas.channels) - 1 - i)
            self._channel_widgets[ch] = self.pl.add_checkbox_button_widget(
                self._toggle_channel(ch), value=True, position=(x, y),
                size=17, border_size=1, color_on=info["color"], color_off="#38414f")
            self.pl.add_text(f"{ch}  {info['pl']}", position=(x + 22, y + 1), font_size=7,
                             color=THEME["text"], name=f"ct_{ch}")
        top = y0 + 23 * len(self.atlas.channels)
        self.pl.add_text("KANALY", position=(12, top + 84), font_size=8, color=THEME["dim"],
                         name="hdr_ch")
        self._button("wszystkie", lambda: self.set_channels(set(self.atlas.channels)),
                     (12, top + 60))
        self._button("zaden", lambda: self.set_channels(set()), (12, top + 38))
        self._button("odwroc", self.invert_channels, (12, top + 16))

        skin_slider = self.pl.add_slider_widget(
            self._set_skin_opacity, [0.0, 1.0], value=1.0, title="krycie skory",
            pointa=(0.74, 0.055), pointb=(0.96, 0.055), style="modern",
            title_height=0.015, fmt="%.2f")
        skin_slider.GetRepresentation().SetLabelHeight(0.013)

        if self.muscles is not None:
            for i, (title, attr, lo, hi, val) in enumerate(PARAMS):
                yy = 0.44 - i * 0.098
                sl = self.pl.add_slider_widget(
                    self._set_param(attr), [lo, hi], value=val, title=title,
                    pointa=(0.74, yy), pointb=(0.96, yy), style="modern",
                    title_height=0.015, fmt="%.2f")
                rep = sl.GetRepresentation()      # wartosc nie ma zachodzic na tytul
                rep.SetLabelHeight(0.013)
                self._sliders.append(sl)
            w, h = self.pl.window_size
            self._button("fokus [f]", self.toggle_focus,
                         (int(w * 0.74), int(h * 0.60)), color=THEME["highlight"])
            self._button("zeruj miesien [r]", self.reset_muscle,
                         (int(w * 0.74), int(h * 0.565)))
            self._button("linie pracy [w]", self.toggle_action_lines,
                         (int(w * 0.74), int(h * 0.53)))

    def _add_help(self) -> None:
        txt = ("[f] fokus na wybranym    [a]/[n] wszystkie kanaly / zaden   [i] odwroc\n"
               "[c]/[v] solo kanalu   [x] bez solo   [m] meridiany   [l] etykiety\n"
               "[k] skora   [j] miesnie   [b] kosci   [h] powiezi   [r] wyzeruj miesien\n"
               "[1] przod  [2] tyl  [3] bok  [4] gora   [0] widok domyslny   [klik] opis")
        self.pl.add_text(txt, position="upper_right", font_size=7, color=THEME["dim"],
                         name="help")
        keys = {
            "a": lambda: self.set_channels(set(self.atlas.channels)),
            "n": lambda: self.set_channels(set()),
            "i": self.invert_channels,
            "c": lambda: self.cycle_solo(-1),
            "v": lambda: self.cycle_solo(1),
            "x": lambda: self.set_solo(None),
            "m": self.toggle_meridians,
            "l": self.toggle_labels,
            "k": lambda: self._flip_layer("skin"),
            "j": lambda: self._flip_layer("muscle"),
            "b": lambda: self._flip_layer("skeleton"),
            "h": self.cycle_fascia,
            "r": self.reset_muscle,
            "w": self.toggle_action_lines,
            "t": self.toggle_tension_mode,
            "z": self.tension_undo,
            "g": self.tension_direction_next,
            "period": lambda: self.set_tension_level(self.tension_level + 1),
            "comma": lambda: self.set_tension_level(self.tension_level - 1),
            "p": self.tension_report,
            "s": self.tension_save,
            "f": self.toggle_focus,
            "0": lambda: self._view("xz", 180),
            "1": lambda: self._view("xz", 180),
            "2": lambda: self._view("xz", 0),
            "3": lambda: self._view("yz", 0),
            "4": lambda: self._view("xy", 0),
        }
        for k, fn in keys.items():
            self.pl.add_key_event(k, fn)

    # -------------------------------------------------------- widok i warstwy
    def _view(self, cam: str, azimuth: float) -> None:
        if self.focus_mode:
            self._focus_leave(keep_camera=False)
        self.pl.camera_position = cam
        if azimuth:
            self.pl.camera.azimuth = azimuth
        self.pl.reset_camera()
        self.pl.camera.zoom(1.35)
        self.pl.render()

    # --------------------------------------------------- warstwa lacznotkankowa
    def cycle_fascia(self) -> None:
        """[h] przelacza kolejne klasy tkanki lacznej.

        Ladujemy je dopiero przy pierwszym uzyciu - to 237 tysiecy trojkatow,
        ktore w atlasie punktow nie sa potrzebne, dopoki ktos o nie nie poprosi.
        Klasy dochodza kumulacyjnie: przegrody, potem powiezi, potem reszta,
        bo tak ubywa czytelnosci i tak rosnie szczegol palpacyjny.
        """
        if self.fascia is None:
            try:
                self.fascia = FasciaSet()
            except (FileNotFoundError, OSError):
                self.info("Brak warstwy powieziowej - uruchom tools/prepare_fascia.py")
                self.pl.render()
                return
        order = self.fascia.classes
        self.fascia_step = (self.fascia_step + 1) % (len(order) + 1)
        act = self.layer_actors["fascia"]
        if self.fascia_step == 0:
            act.SetVisibility(False)
            self.info("Warstwa lacznotkankowa wylaczona.")
            self.pl.render()
            return
        shown = order[:self.fascia_step]
        mesh = self.fascia.submesh(shown)
        self.pl.add_mesh(mesh, color=CLASS_COLOR[shown[-1]], smooth_shading=True,
                         opacity=0.55, name="fascia")
        self.layer_actors["fascia"] = self.pl.renderer.actors["fascia"]
        self.layer_actors["fascia"].SetVisibility(True)
        counts = self.fascia.counts()
        txt = ", ".join(f"{CLASS_PL[c]} ({counts[c]})" for c in shown)
        self.info(f"Warstwa lacznotkankowa: {txt}"
                  f"\n{mesh.n_cells} trojkatow.  [h] dalej")
        self.pl.render()

    def _toggle_layer(self, name: str):
        def cb(state: bool) -> None:
            self.layer_actors[name].SetVisibility(bool(state))
            self.pl.render()
        return cb

    def _flip_layer(self, name: str) -> None:
        act = self.layer_actors[name]
        act.SetVisibility(not act.GetVisibility())
        self.pl.render()

    def _set_skin_opacity(self, value: float) -> None:
        self.layer_actors["skin"].GetProperty().SetOpacity(float(value))
        self.pl.render()

    # ---------------------------------------------------------------- kanaly
    def _toggle_channel(self, ch: str):
        def cb(state: bool) -> None:
            (self.visible_channels.add if state else self.visible_channels.discard)(ch)
            self._apply_visibility()
        return cb

    def set_channels(self, chans: set) -> None:
        self.visible_channels = set(chans)
        self._sync_channel_buttons()
        self._apply_visibility()

    def invert_channels(self) -> None:
        self.set_channels(set(self.atlas.channels) - self.visible_channels)

    def _sync_channel_buttons(self) -> None:
        """Przyciski kanalow maja pokazywac stan, nawet gdy zmienil go klawisz."""
        for ch, w in getattr(self, "_channel_widgets", {}).items():
            try:
                w.GetRepresentation().SetState(ch in self.visible_channels)
            except AttributeError:
                pass

    def set_solo(self, ch: str | None) -> None:
        self.solo = ch
        self.visible_channels = set(self.atlas.channels) if ch is None else {ch}
        self.show_labels = ch is not None
        self._sync_channel_buttons()
        self._apply_visibility()
        info = self.atlas.channels.get(ch)
        self.info(f"Kanal {ch} - {info['pl']} ({info['en']}), zywiol: {info['element']}"
                  if info else "Wszystkie kanaly")

    def cycle_solo(self, step: int) -> None:
        order = list(self.atlas.channels)
        if self.solo is None:
            self.set_solo(order[0 if step > 0 else -1])
        else:
            self.set_solo(order[(order.index(self.solo) + step) % len(order)])

    def toggle_meridians(self) -> None:
        self.show_meridians = not self.show_meridians
        self._apply_visibility()

    def toggle_labels(self) -> None:
        self.show_labels = not self.show_labels
        self._apply_visibility()

    def _apply_visibility(self) -> None:
        for ch in self.atlas.channels:
            on = ch in self.visible_channels
            self.point_actors[ch].SetVisibility(on)
            if ch in self.line_actors:
                self.line_actors[ch].SetVisibility(on and self.show_meridians)
            self.label_actors[ch].SetVisibility(on and self.show_labels)
        self.pl.render()

    # --------------------------------------------------------------- miesnie
    def _set_param(self, attr: str):
        def cb(value: float) -> None:
            if self.selected is None or self.muscles is None:
                return
            setattr(self.selected, attr, float(value))
            self._deform_selected()
            self._muscle_info(self.selected, to_console=False)
            self.pl.render()
        return cb

    def _deform_selected(self) -> None:
        """Przelicza tylko wybrany miesien - kilka tysiecy wierzcholkow."""
        if self._sel_mesh is None:
            return
        pts = self.muscles.deform(self.selected, self._sel_rest)
        self._sel_mesh.points[:] = pts
        self._sel_mesh.compute_normals(inplace=True, point_normals=True,
                                       cell_normals=False, split_vertices=False)
        self._sel_mesh.Modified()

    def select_muscle(self, muscle, to_console: bool = True) -> None:
        if self.muscles is None or muscle is None:
            return
        sel = self.muscles.mesh.point_data["sel"]
        sel[:] = 0.0
        sel[muscle.slice] = 1.0
        self.muscles.mesh.Modified()
        self.selected = muscle
        self._sel_mesh = self.muscles.submesh(muscle)
        self._sel_rest = np.array(self._sel_mesh.points, copy=True)
        self._deform_selected()
        self._sel_actor = self.pl.add_mesh(
            self._sel_mesh, color=THEME["muscle_sel"], smooth_shading=True,
            specular=0.25, name="sel_muscle")
        for slider, (_, attr, _, _, _) in zip(self._sliders, PARAMS):
            slider.GetRepresentation().SetValue(getattr(muscle, attr))
        self._muscle_info(muscle, to_console=to_console)
        self._selected_action_line()
        if self.focus_mode:
            self._focus_camera()
        self.pl.render()

    def reset_muscle(self) -> None:
        if self.selected is None or self.muscles is None:
            return
        self.muscles.reset(self.selected)
        self._deform_selected()
        for slider, (_, attr, _, _, _) in zip(self._sliders, PARAMS):
            slider.GetRepresentation().SetValue(getattr(self.selected, attr))
        self._muscle_info(self.selected, to_console=False)
        self.pl.render()

    def toggle_action_lines(self) -> None:
        """Pokazuje kierunki ciagu miesni - os pracy z grotami ku brzuscowi."""
        if self.muscles is None:
            return
        self.show_action = not self.show_action
        # Linie biegna SRODKIEM miesni, wiec przy pelnym kryciu byly niewidoczne -
        # na czas ich pokazywania warstwa miesni schodzi do polprzezroczystej.
        mus = self.layer_actors["muscle"]
        if self.show_action:
            self._muscle_opacity = mus.GetProperty().GetOpacity()
            mus.SetVisibility(True)
            mus.GetProperty().SetOpacity(0.28)
        else:
            mus.GetProperty().SetOpacity(self._muscle_opacity)
        if self.show_action and self._action_actor is None:
            mesh = action_lines.build_all(self.muscles)
            if mesh.n_points:
                self._action_actor = self.pl.add_mesh(
                    mesh, color=THEME["accent"], name="action_lines")
        if self._action_actor is not None:
            self._action_actor.SetVisibility(self.show_action)
        self._selected_action_line()
        self.info("Linie pracy miesni: wlaczone" if self.show_action
                  else "Linie pracy miesni: wylaczone")
        self.pl.render()

    def _selected_action_line(self) -> None:
        """Linia pracy wybranego miesnia - grubsza i w kolorze zaznaczenia."""
        if self.muscles is None:
            return
        if self.selected is None or not self.show_action:
            self.pl.remove_actor("action_sel", render=False)
            return
        pts = self.muscles.rest[self.selected.slice]
        if len(pts) < 40:
            return
        mesh = action_lines.build(pts, self.selected.axis, arrow_len=0.030)
        if mesh.n_points:
            self.pl.add_mesh(mesh, color=THEME["muscle_sel"], name="action_sel")
        # wybrany miesien tez musi przepuscic swoja linie
        if self._sel_actor is not None:
            self._sel_actor.GetProperty().SetOpacity(0.35 if self.show_action else 1.0)

    def _muscle_info(self, muscle, to_console: bool = False) -> None:
        pts = self.muscles.points_on(
            muscle, [m for m in self.atlas.markers if m["side"] in (muscle.side, "C")])
        shown = [(p, d) for p, d in pts if p["channel"] in self.visible_channels][:8]
        lines = [f"MIESIEN: {muscle.label}",
                 f"dlugosc {muscle.length_m * 100:.1f} cm   objetosc {muscle.volume_l:.3f} l",
                 f"napiecie {muscle.activation:.2f}   tonus {muscle.tone:.2f}   "
                 f"przerost {muscle.hypertrophy:.2f}   sztywnosc {muscle.stiffness:.2f}",
                 f"skrocenie brzusca {muscle.shortening() * 100:.0f}%   "
                 f"sila wzgledna {muscle.force_ratio():.2f}"]
        if shown:
            lines.append("punkty nad miesniem: " + ", ".join(
                f"{p['code']} ({d * 1000:.0f} mm)" for p, d in shown))
        self.info("\n".join(lines))
        if to_console:
            _safe_print("\n" + "=" * 78 + f"\n{muscle.label}\n" + "\n".join(lines[1:]))


    # ------------------------------------------------------------ mapa napiec
    def _ensure_tension(self):
        if self.tension is None:
            if self.fascia is None:          # do palpacji powiez jest istotniejsza
                try:                         # niz brzusiec, wiec ladujemy ja tutaj
                    self.fascia = FasciaSet()
                except (FileNotFoundError, OSError):
                    self.fascia = None
            self.tension = TensionMap(self.atlas, self.muscles, fascia=self.fascia)
        return self.tension

    def toggle_tension_mode(self) -> None:
        """Tryb mapowania napiec: klikniecie w ciało stawia oznaczenie."""
        self._ensure_tension()
        self.tension_mode = not self.tension_mode
        self._await_direction = False
        self._refresh_tension()
        if self.tension_mode:
            self.info(f"MAPA NAPIEC — klikaj w cialo, poziom {self.tension_level} "
                      f"({LEVELS[self.tension_level]}).  [,]/[.] poziom   [g] kierunek\n"
                      f"[z] cofnij   [s] zapisz   [p] raport   [t] wyjscie")
        else:
            self.info("Mapa napiec wylaczona (oznaczenia zostaja).")
        self.pl.render()

    def set_tension_level(self, level: int) -> None:
        self.tension_level = int(np.clip(level, 0, 3))
        if self.tension_mode:
            self.info(f"Poziom oznaczenia: {self.tension_level} "
                      f"({LEVELS[self.tension_level]})")
            self.pl.render()

    def tension_direction_next(self) -> None:
        """Po wcisnieciu [g] kolejne klikniecie nadaje kierunek ostatniemu oznaczeniu."""
        if not self.tension_mode or not (self.tension and self.tension.marks):
            return
        self._await_direction = True
        self.info("Kliknij tam, DOKAD ciagnie napiecie - to nada kierunek "
                  "ostatniemu oznaczeniu.")
        self.pl.render()

    def tension_undo(self) -> None:
        if self.tension and self.tension.marks:
            m = self.tension.remove_last()
            self._refresh_tension()
            self.info(f"Cofnieto oznaczenie [{m.id}]. Zostalo: {len(self.tension.marks)}")
            self.pl.render()

    def tension_save(self) -> None:
        if not (self.tension and self.tension.marks):
            return
        path = self.tension.save()
        self.info(f"Zapisano {len(self.tension.marks)} oznaczen: {path.name}")
        _safe_print(f"\nmapa napiec zapisana: {path}")
        self.pl.render()

    def tension_report(self) -> None:
        if not (self.tension and self.tension.marks):
            self.info("Mapa napiec jest pusta - nie ma czego podsumowac.")
            self.pl.render()
            return
        _safe_print("\n" + self.tension.report())
        self.info(f"Raport z {len(self.tension.marks)} oznaczen wypisany na konsoli.")
        self.pl.render()

    def _add_tension_mark(self, pos) -> None:
        t = self._ensure_tension()
        if self._await_direction and t.marks:
            last = t.marks[-1]
            d = np.asarray(pos) - last.pos
            if np.linalg.norm(d) > 1e-4:
                last.direction = d / np.linalg.norm(d)
            self._await_direction = False
            self._refresh_tension()
            self.info(f"Nadano kierunek oznaczeniu [{last.id}].")
            self.pl.render()
            return
        mark = t.add(pos, level=self.tension_level)
        self._refresh_tension()
        a = t.analyse(mark)
        txt = [f"Oznaczenie [{mark.id}] poziom {mark.level} ({LEVELS[mark.level]})"]
        if a["muscle"] is not None:
            txt.append(f"miesien pod spodem: {a['muscle'].label} ({a['muscle_dist_mm']} mm)")
        if a["fascia"] is not None:
            txt.append(f"tkanka laczna: {a['fascia'].label} ({a['fascia_dist_mm']} mm)")
        if a["points"]:
            txt.append("punkty w zasiegu: " + ", ".join(
                f"{p['code']} ({p['dist_mm']:.0f} mm)" for p in a["points"][:5]))
        self.info("\n".join(txt))
        self.pl.render()

    def _refresh_tension(self) -> None:
        """Rysuje oznaczenia, ich kierunki i rozmyte pole napiecia na skorze."""
        for name in ("tension_marks", "tension_dirs"):
            self.pl.remove_actor(name, render=False)
        t = self.tension
        skin = self.atlas.layers["skin"]
        if not (self.tension_mode and t and t.marks):
            if "napiecie" in skin.point_data:
                del skin.point_data["napiecie"]
                self.layer_actors["skin"] = self.pl.add_mesh(
                    skin, color=THEME["skin"], smooth_shading=True, specular=0.12,
                    name="skin", opacity=self.layer_actors["skin"].GetProperty().GetOpacity())
            return

        # pole napiecia na skorze: oznaczenia rozmywaja sie i sumuja
        field = t.field(np.asarray(skin.points))
        skin.point_data["napiecie"] = field
        op = self.layer_actors["skin"].GetProperty().GetOpacity()
        self.layer_actors["skin"] = self.pl.add_mesh(
            skin, scalars="napiecie", clim=[0.0, 3.0], show_scalar_bar=False,
            cmap=[THEME["skin"], LEVEL_COLORS[1], LEVEL_COLORS[2], LEVEL_COLORS[3]],
            smooth_shading=True, specular=0.12, name="skin", opacity=op)

        pts = np.array([m.pos for m in t.marks])
        lev = np.array([float(m.level) for m in t.marks])
        cloud = pv.PolyData(pts)
        cloud["poziom"] = lev
        glyph = cloud.glyph(geom=pv.Sphere(radius=0.008), scale=False, orient=False)
        self.pl.add_mesh(glyph, scalars="poziom", clim=[0, 3], show_scalar_bar=False,
                         cmap=[LEVEL_COLORS[0], LEVEL_COLORS[1], LEVEL_COLORS[2],
                               LEVEL_COLORS[3]], name="tension_marks")
        arrows = [pv.Arrow(start=m.pos, direction=m.direction, scale=0.06)
                  for m in t.marks if m.direction is not None]
        if arrows:
            merged = arrows[0].merge(arrows[1:]) if len(arrows) > 1 else arrows[0]
            self.pl.add_mesh(merged, color=THEME["text"], name="tension_dirs")

    # ----------------------------------------------------------------- fokus
    def _focus_target(self):
        if self.focused_key:
            m = self.atlas.by_key[self.focused_key]
            return np.asarray(m["pos"]), f"{m['code']} ({SIDE_PL[m['side']]})"
        if self.selected is not None:
            return np.asarray(self.selected.centroid), self.selected.label
        return None, None

    def toggle_focus(self) -> None:
        if self.focus_mode:
            self._focus_leave()
        else:
            self._focus_enter()

    def _focus_enter(self) -> None:
        target, label = self._focus_target()
        if target is None:
            self.info("Fokus: najpierw kliknij punkt albo miesien.")
            self.pl.render()
            return
        self._focus_saved = {
            "opacity": {k: a.GetProperty().GetOpacity() for k, a in self.layer_actors.items()},
            "visible": {k: bool(a.GetVisibility()) for k, a in self.layer_actors.items()},
            "channels": set(self.visible_channels),
            "labels": self.show_labels,
            "camera": self.pl.camera_position,
        }
        for name, op in FOCUS_OPACITY.items():
            act = self.layer_actors[name]
            # powiezi wchodza w fokus tylko wtedy, gdy uzytkownik je wlaczyl -
            # inaczej [f] wskrzeszalby pusta warstwe
            if name == "fascia" and self.fascia_step == 0:
                continue
            act.SetVisibility(True)
            act.GetProperty().SetOpacity(op)
        if self.focused_key:                       # zostaw tylko kanal wybranego punktu
            self.visible_channels = {self.atlas.by_key[self.focused_key]["channel"]}
            self.show_labels = True
        self._apply_visibility()
        self.focus_mode = True
        self._focus_camera()
        self.info(f"FOKUS: {label}   [f] wyjscie   "
                  f"{'miesien pod punktem podswietlony' if self.selected else ''}")
        self.pl.render()

    def _focus_camera(self, distance: float = None) -> None:
        """Ustawia wybrany obiekt w srodku obrotu i przybliza tak, by go objac."""
        target, _ = self._focus_target()
        if target is None:
            return
        if distance is None:
            size = self.selected.length_m if self.selected is not None else 0.12
            distance = float(np.clip(2.6 * size, 0.30, 1.10))
        cam = self.pl.camera
        d = np.array(cam.position) - np.array(cam.focal_point)
        n = np.linalg.norm(d)
        u = d / n if n > 1e-9 else np.array([0.0, -1.0, 0.0])
        cam.focal_point = tuple(target)
        cam.position = tuple(target + u * distance)
        self.pl.render()

    def _focus_leave(self, keep_camera: bool = True) -> None:
        if not self.focus_mode or self._focus_saved is None:
            self.focus_mode = False
            return
        s = self._focus_saved
        for k, a in self.layer_actors.items():
            a.GetProperty().SetOpacity(s["opacity"][k])
            a.SetVisibility(s["visible"][k])
        self.visible_channels = s["channels"]
        self.show_labels = s["labels"]
        self._sync_channel_buttons()
        self._apply_visibility()
        if keep_camera:
            self.pl.camera_position = s["camera"]
        self.focus_mode = False
        self._focus_saved = None
        self.info("Fokus wylaczony.")
        self.pl.render()

    # ---------------------------------------------------------- wybor klikiem
    def _on_pick(self, picked, *args) -> None:
        pos = np.asarray(getattr(picked, "points", picked)).reshape(-1, 3)
        if not len(pos):
            return
        pos = pos[0]
        if self.tension_mode:
            self._add_tension_mark(pos)
            return
        cands = [m for m in self.atlas.markers if m["channel"] in self.visible_channels
                 and self.point_actors[m["channel"]].GetVisibility()]
        if cands:
            d = np.linalg.norm(np.array([m["pos"] for m in cands]) - pos, axis=1)
            i = int(np.argmin(d))
            if d[i] < 0.025:
                self.focus(cands[i]["key"])
                return
        if self.muscles is not None and self.layer_actors["muscle"].GetVisibility():
            m = self.muscles.at_position(pos)
            if m is not None:
                self.focused_key = None
                self.select_muscle(m)

    # ------------------------------------------------------------ opis punktu
    def focus(self, key: str) -> None:
        d = self.atlas.describe(key)
        meta = d["meta"]
        ch = self.atlas.channels[d["channel"]]
        self.focused_key = key
        lines = [f"{d['code']}  {d['pinyin']}  ({meta.get('english', '')})",
                 f"kanal: {ch['pl']} / {ch['en']}   strona: {SIDE_PL[d['side']]}"
                 f"   dokladnosc: {d['accuracy']}"]
        if meta.get("categories"):
            lines.append("typ: " + ", ".join(meta["categories"])[:110])
        lines += textwrap.wrap("lokalizacja: " + (d["note_pl"] or meta.get("location", "")), 96)[:3]
        if self.muscles is not None:
            mus, dist = self.muscles.muscle_under(d["pos"])
            if mus is not None:
                lines.append(f"miesien pod punktem: {mus.label}  ({dist * 1000:.0f} mm od skory)")
                self.select_muscle(mus, to_console=False)
                self.focused_key = key           # wybor miesnia nie kasuje punktu
        if meta.get("needling"):
            lines += textwrap.wrap("naklucie: " + meta["needling"], 96)[:2]
        lines.append("[f] fokus na tym punkcie")
        self.info("\n".join(lines))
        self.pl.add_mesh(pv.Sphere(radius=0.011, center=d["pos"], theta_resolution=20,
                                   phi_resolution=20),
                         color=THEME["highlight"], name="highlight", opacity=0.9)
        if self.focus_mode:
            self._focus_camera()
        self._print_full(d)
        self.pl.render()

    def _print_full(self, d: dict) -> None:
        """Pelny opis na konsole; nie moze wywrocic obslugi klikniecia."""
        try:
            meta = d["meta"]
            out = ["\n" + "=" * 78,
                   f"{d['code']}  {d['pinyin']}  {meta.get('english', '')}   [{d['accuracy']}]",
                   "-" * 78]
            for label, key in (("LOKALIZACJA (Deadman)", "location"), ("NAKLUCIE", "needling"),
                               ("DZIALANIE", "actions"), ("WSKAZANIA", "indications")):
                v = meta.get(key)
                if not v:
                    continue
                if isinstance(v, list):
                    v = " ".join(str(x) for x in v if x)
                out.append(f"{label}: " + textwrap.fill(str(v), 76, subsequent_indent="  ")[:900])
            if d["note_pl"]:
                out.append("NOTATKA: " + d["note_pl"])
            _safe_print("\n".join(out))
        except Exception as exc:                      # opis to dodatek, nie powod do awarii
            _safe_print(f"[nie udalo sie wypisac opisu {d.get('code')}: {exc}]")

    def info(self, text: str) -> None:
        self.pl.add_text(text, position=(250, 12), font_size=8, color=THEME["text"],
                         name="info")

    def show(self) -> None:
        self.pl.show()

    def screenshot(self, path: str) -> None:
        self.pl.screenshot(path)
