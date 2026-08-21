"""Interaktywna przegladarka 3D: warstwy ciala, kanaly, punkty, miesnie."""
from __future__ import annotations

import sys
import textwrap

import numpy as np
import pyvista as pv

from .model import Atlas
from .muscles import MuscleSet

THEME = {
    "bg": "#0e1116",
    "skin": "#e8cbaa",
    "muscle": "#a8382c",
    "muscle_sel": "#ffcf40",
    "bone": "#efe9dc",
    "text": "#dfe6f0",
    "dim": "#8792a3",
    "highlight": "#ffe36a",
}
LAYER_LABELS = {"skin": "skora", "muscle": "miesnie", "skeleton": "kosci"}
SIDE_PL = {"R": "prawa", "L": "lewa", "C": "linia posrodkowa"}


def _safe_print(text: str) -> None:
    """Wypis odporny na konsole, ktora nie zna znakow spoza swojej strony kodowej."""
    enc = getattr(sys.stdout, "encoding", None) or "utf-8"
    sys.stdout.write(text.encode(enc, errors="replace").decode(enc, errors="replace") + "\n")
    sys.stdout.flush()


class Viewer:
    def __init__(self, atlas: Atlas, off_screen: bool = False, window_size=(1560, 980)):
        self.atlas = atlas
        self.muscles = MuscleSet() if atlas.mode == "anatomy" else None
        self.selected = None                       # wybrany miesien
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
        self._sliders: list = []
        self._build()

    # --- scena -----------------------------------------------------------
    def _build(self) -> None:
        a = self.atlas
        self.layer_actors["skin"] = self.pl.add_mesh(
            a.layers["skin"], color=THEME["skin"], opacity=1.0, smooth_shading=True,
            specular=0.12, name="skin")
        if self.muscles is not None:
            m = self.muscles.mesh
            m.point_data["sel"] = np.zeros(m.n_points, dtype=np.float32)
            self.layer_actors["muscle"] = self.pl.add_mesh(
                m, scalars="sel", clim=[0, 1], show_scalar_bar=False,
                cmap=[THEME["muscle"], THEME["muscle_sel"]], smooth_shading=True,
                specular=0.2, name="muscle")
        else:
            self.layer_actors["muscle"] = self.pl.add_mesh(
                a.layers["muscle"], color=THEME["muscle"], smooth_shading=True, name="muscle")
        self.layer_actors["skeleton"] = self.pl.add_mesh(
            a.layers["skeleton"], color=THEME["bone"], smooth_shading=True, name="skeleton")
        self.layer_actors["muscle"].SetVisibility(False)
        self.layer_actors["skeleton"].SetVisibility(False)

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
                        lines.append(pv.Spline(poly, len(poly) * 3).tube(radius=r * 0.33))
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
        self.info("Kliknij punkt albo miesien (wlacz warstwe miesni klawiszem [j]).")
        self.pl.enable_point_picking(callback=self._on_pick, show_message=False,
                                     left_clicking=True, show_point=False, tolerance=0.02)
        self._view("xz", 180)

    def _add_widgets(self) -> None:
        y = 12
        for name in ("skin", "muscle", "skeleton"):
            self.pl.add_checkbox_button_widget(
                self._toggle_layer(name), value=(name == "skin"), position=(12, y),
                size=20, border_size=1, color_on=THEME["skin"], color_off="#3a4250")
            self.pl.add_text(LAYER_LABELS[name], position=(38, y + 2), font_size=8,
                             color=THEME["text"], name=f"lt_{name}")
            y += 26
        self.pl.add_slider_widget(
            self._set_skin_opacity, [0.0, 1.0], value=1.0, title="krycie skory",
            pointa=(0.74, 0.055), pointb=(0.96, 0.055), style="modern",
            title_height=0.018, fmt="%.2f")

        x, y = 12, 210 + 23 * (len(self.atlas.channels) - 1)
        for ch, info in self.atlas.channels.items():
            self.pl.add_checkbox_button_widget(
                self._toggle_channel(ch), value=True, position=(x, y),
                size=17, border_size=1, color_on=info["color"], color_off="#3a4250")
            self.pl.add_text(f"{ch}  {info['pl']}", position=(x + 22, y + 1), font_size=7,
                             color=THEME["text"], name=f"ct_{ch}")
            y -= 23

        if self.muscles is not None:
            specs = [("napiecie", "activation", 0.0, 1.0, 0.0),
                     ("tonus spoczynkowy", "tone", 0.0, 0.6, 0.0),
                     ("przerost", "hypertrophy", 0.75, 1.35, 1.0),
                     ("sztywnosc", "stiffness", 0.2, 1.6, 1.0)]
            for i, (title, attr, lo, hi, val) in enumerate(specs):
                yy = 0.40 - i * 0.075
                self._sliders.append(self.pl.add_slider_widget(
                    self._set_param(attr), [lo, hi], value=val, title=title,
                    pointa=(0.74, yy), pointb=(0.96, yy), style="modern",
                    title_height=0.016, fmt="%.2f"))

    def _add_help(self) -> None:
        txt = ("[a] wszystkie kanaly   [n] zaden   [c]/[v] solo kanalu   [x] bez solo\n"
               "[m] meridiany   [l] etykiety   [k] skora   [j] miesnie   [b] kosci\n"
               "[r] wyzeruj miesien   [1] przod  [2] tyl  [3] bok  [4] gora   [klik] opis")
        self.pl.add_text(txt, position="upper_right", font_size=7, color=THEME["dim"],
                         name="help")
        keys = {
            "a": lambda: self.set_channels(set(self.atlas.channels)),
            "n": lambda: self.set_channels(set()),
            "c": lambda: self.cycle_solo(-1),
            "v": lambda: self.cycle_solo(1),
            "x": lambda: self.set_solo(None),
            "m": self.toggle_meridians,
            "l": self.toggle_labels,
            "k": lambda: self._flip_layer("skin"),
            "j": lambda: self._flip_layer("muscle"),
            "b": lambda: self._flip_layer("skeleton"),
            "r": self.reset_muscle,
            "1": lambda: self._view("xz", 180),
            "2": lambda: self._view("xz", 0),
            "3": lambda: self._view("yz", 0),
            "4": lambda: self._view("xy", 0),
        }
        for k, fn in keys.items():
            self.pl.add_key_event(k, fn)

    # --- widok i warstwy -------------------------------------------------
    def _view(self, cam: str, azimuth: float) -> None:
        self.pl.camera_position = cam
        if azimuth:
            self.pl.camera.azimuth = azimuth
        self.pl.reset_camera()
        self.pl.camera.zoom(1.35)
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

    # --- kanaly ----------------------------------------------------------
    def _toggle_channel(self, ch: str):
        def cb(state: bool) -> None:
            (self.visible_channels.add if state else self.visible_channels.discard)(ch)
            self._apply_visibility()
        return cb

    def set_channels(self, chans: set) -> None:
        self.visible_channels = set(chans)
        self._apply_visibility()

    def set_solo(self, ch: str | None) -> None:
        self.solo = ch
        self.visible_channels = set(self.atlas.channels) if ch is None else {ch}
        self.show_labels = ch is not None
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

    # --- miesnie ---------------------------------------------------------
    def _set_param(self, attr: str):
        def cb(value: float) -> None:
            if self.selected is None or self.muscles is None:
                return
            setattr(self.selected, attr, float(value))
            self.muscles.apply(self.selected)
            self.muscles.refresh()
            self._muscle_info(self.selected)
            self.pl.render()
        return cb

    def select_muscle(self, muscle) -> None:
        if self.muscles is None or muscle is None:
            return
        sel = self.muscles.mesh.point_data["sel"]
        sel[:] = 0.0
        sel[muscle.slice] = 1.0
        self.muscles.mesh.Modified()
        self.selected = muscle
        for slider, attr in zip(self._sliders,
                                ("activation", "tone", "hypertrophy", "stiffness")):
            slider.GetRepresentation().SetValue(getattr(muscle, attr))
        self._muscle_info(muscle)
        self.pl.render()

    def reset_muscle(self) -> None:
        if self.selected is None or self.muscles is None:
            return
        self.muscles.reset(self.selected)
        self.muscles.refresh()
        self.select_muscle(self.selected)

    def _muscle_info(self, muscle) -> None:
        pts = self.muscles.points_on(
            muscle, [m for m in self.atlas.markers
                     if m["side"] in (muscle.side, "C") and m["channel"] in self.visible_channels])
        lines = [f"MIESIEN: {muscle.label}",
                 f"dlugosc {muscle.length_m * 100:.1f} cm   objetosc {muscle.volume_l:.3f} l   "
                 f"trojkatow {muscle.polys}",
                 f"napiecie {muscle.activation:.2f}   tonus {muscle.tone:.2f}   "
                 f"przerost {muscle.hypertrophy:.2f}   sztywnosc {muscle.stiffness:.2f}",
                 f"skrocenie brzusca {muscle.shortening() * 100:.0f}%   "
                 f"sila wzgledna {muscle.force_ratio():.2f}"]
        if pts:
            lines.append("punkty nad miesniem: " + ", ".join(
                f"{p['code']} ({d * 1000:.0f} mm)" for p, d in pts[:8]))
        self.info("\n".join(lines))
        _safe_print("\n" + "=" * 78 + f"\n{muscle.label}\n" + "\n".join(lines[1:]))

    # --- wybor klikniciem ------------------------------------------------
    def _on_pick(self, picked, *args) -> None:
        pos = np.asarray(getattr(picked, "points", picked)).reshape(-1, 3)
        if not len(pos):
            return
        pos = pos[0]
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
                self.select_muscle(m)

    # --- opis punktu -----------------------------------------------------
    def focus(self, key: str) -> None:
        d = self.atlas.describe(key)
        meta = d["meta"]
        ch = self.atlas.channels[d["channel"]]
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
        if meta.get("needling"):
            lines += textwrap.wrap("naklucie: " + meta["needling"], 96)[:2]
        self.info("\n".join(lines))
        self.pl.add_mesh(pv.Sphere(radius=0.010, center=d["pos"], theta_resolution=20,
                                   phi_resolution=20),
                         color=THEME["highlight"], name="highlight", opacity=0.85)
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
