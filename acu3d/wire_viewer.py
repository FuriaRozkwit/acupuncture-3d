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
import time

import numpy as np
import pyvista as pv

from . import combos
from .channels import ChannelMap
from .model import Atlas
from .tension import DEFAULT_RADIUS, LEVEL_COLORS, LEVELS, TensionMap
from .wire import FLESH_COLORS, FLESH_RANGE, Wire

# aktory oznaczen i linii startuja puste i zapelniaja sie przy pierwszym kliknieciu
pv.global_theme.allow_empty_mesh = True

THEME = {"bg": "#0b0e13", "text": "#dfe6f0", "dim": "#7c879a",
         "point": "#5f6b7d", "sel": "#ffe36a", "force": "#ff8a5c",
         "accent": "#7fd1c8"}

# Klasyfikacja punktow wedlug definicji kanalowych. Kolejnosc jest merytoryczna,
# nie alfabetyczna: najpierw piatka transportu (bo tworzy powtarzalny wzorzec na
# kazdej konczynie), potem punkty wiazace kanal z innymi, na koncu zbiorcze.
ROLE_ORDER = (
    ("jing_well", "studnia", "#ffe36a"),
    ("ying_spring", "zrodlo", "#ffc25c"),
    ("shu_stream", "strumien", "#ff9a5c"),
    ("jing_river", "rzeka", "#c47ad6"),
    ("he_sea", "morze", "#6ba8ff"),
    ("lower_he_sea", "dolne morze", "#4d7fd6"),
    ("yuan", "yuan (zrodlowy)", "#7fd1c8"),
    ("luo", "luo (laczacy)", "#5fd18a"),
    ("great_luo", "wielki luo", "#3fa86a"),
    ("xi", "xi (szczelina)", "#e0655c"),
    ("front_mu", "mu (przednie)", "#d98cc4"),
    ("back_shu", "shu grzbietu", "#a882d6"),
    ("confluent", "otwierajacy", "#ffab5c"),
    ("hui", "hui (zbiorczy)", "#c9d46a"),
    ("command", "rozkazujacy", "#8fd4c2"),
    ("window", "okno nieba", "#9fb6d6"),
    ("sea", "morze qi/krwi", "#6ad4c9"),
    ("crossing", "skrzyzowanie", "#ff7a9c"),
)


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
        # Aktory oznaczen i linii powstaja RAZ, a potem tylko podmieniamy w nich
        # geometrie. add_mesh kosztuje 5 ms, bo przebudowuje caly potok VTK -
        # przy dwoch wywolaniach na klik to byla polowa odczuwanej zwloki.
        self._marks_poly = pv.PolyData()
        self._force_poly = pv.PolyData()
        self._arrow = pv.Arrow(start=(0, 0, 0), direction=(1, 0, 0), scale=1.0)
        self._sphere = pv.Sphere(radius=0.009, theta_resolution=12, phi_resolution=12)
        # tryb kombinacji: klik wybiera PUNKT zamiast stawiac oznaczenie
        self.combo_mode = False
        self.combo: list = []            # klucze wybranych markerow, max 3
        self.book = combos.CombinationBook()
        self.chan = ChannelMap(atlas)     # meridianologia: czym punkty sie lacza
        self._browse = -1                # ktora zapisana kombinacja jest pokazana
        self._combo_poly = pv.PolyData()
        self._sel_sphere = pv.Sphere(radius=0.011, theta_resolution=14,
                                     phi_resolution=14)
        self._mpos = np.array([m["pos"] for m in atlas.markers], dtype=float)
        self.role_on: set = set()          # ktore klasy sa wlaczone w pasku
        self._role_actors: dict = {}
        self._role_widgets: dict = {}
        self._inspected: str | None = None
        self.hover_on = True             # podglad obszaru pod kursorem
        self._hover_idx = -1             # ostatni wierzcholek pod mysza
        self._hover_t = 0.0              # ostatnie odswiezenie (dlawienie)
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
            self.label_actors[ch].SetPickable(False)

        self._marks_poly.point_data["poziom"] = np.zeros(0)
        self.pl.add_mesh(self._marks_poly, scalars="poziom", clim=[0, 3],
                         show_scalar_bar=False,
                         cmap=[LEVEL_COLORS[i] for i in range(4)], name="marks")
        self.pl.add_mesh(self._force_poly, color=THEME["force"], name="forces")
        self.pl.add_mesh(self._combo_poly, color=THEME["sel"], name="combo")

        # Podglad obszaru pod kursorem. Kula ma DOKLADNIE promien oznaczenia
        # napiecia, wiec przed klikiem widac, co obejmie - a nie tylko gdzie
        # sie trafi. Aktory powstaja raz i tylko sie przesuwaja: ruch myszy
        # nie moze przebudowywac geometrii.
        # Siatka, nie bryla: przezroczysta kula nad siatka cunowa robila sie
        # nieczytelna, a rownolezniki czytaja sie na ciemnym tle od razu i nie
        # zaslaniaja tego, co pod nimi.
        self._hover_actor = self.pl.add_mesh(
            pv.Sphere(radius=DEFAULT_RADIUS, theta_resolution=20, phi_resolution=14),
            color=THEME["accent"], style="wireframe", line_width=1,
            opacity=0.55, name="hover")
        self._hover_dot = self.pl.add_mesh(
            pv.Sphere(radius=0.004, theta_resolution=10, phi_resolution=10),
            color=THEME["accent"], name="hover_dot")
        self._hover_actor.SetVisibility(False)
        self._hover_dot.SetVisibility(False)
        self._hover_actor.SetPickable(False)
        self._hover_dot.SetPickable(False)

        self._build_role_index()
        self._add_role_panel()
        self._add_help()
        self.info("Klikaj w cialo - kazde klikniecie stawia oznaczenie napiecia.")
        self.pl.enable_point_picking(callback=self._on_pick, show_message=False,
                                     left_clicking=True, show_point=False, tolerance=0.03)
        self._keys()
        self.pl.track_mouse_position()
        try:
            self.pl.iren.add_observer("MouseMoveEvent", self._on_hover)
        except Exception:                 # bez interaktora (zrzuty) podglad milczy
            pass
        self._view("xz", 180)
        self.wire._build_graph()      # 0,3 s tu zamiast przy pierwszym kliknieciu

    def _add_help(self) -> None:
        self.pl.add_text(
            "SILNIK NAPIEC\n"
            "[klik] oznaczenie   [,]/[.] poziom 0-3   [z] cofnij\n"
            "[spacja] koniec linii sily, nastepna zaczyna sie od nowa\n"
            "[q] punkty   [l] etykiety   [f] filtr: tylko punkty na liniach\n"
            "[c] KOMBINACJE: klik wybiera punkt (2-3)   [x] wyczysc\n"
            "[Enter] zapisz kombinacje   [n] przegladaj zapisane\n"
            "[p] zestawienie   [s] zapis sesji   [u] zdejmij klasyfikacje\n"
            "[1-4] widoki   [0] domyslny   [h] podglad obszaru pod kursorem\n"
            "pasek po lewej: klasy punktow wg definicji kanalowych",
            position="upper_right", font_size=8, color=THEME["dim"], name="help")

    def _keys(self) -> None:
        keys = {
            "period": lambda: self.set_level(self.level + 1),
            "comma": lambda: self.set_level(self.level - 1),
            "z": self.undo,
            "space": self.break_chain,
            "q": self.toggle_points,
            "f": self.toggle_filter,
            "c": self.toggle_combo,
            "x": self.combo_clear,
            "u": lambda: self.set_roles(set()),
            "h": self.toggle_hover,
            "Return": self.combo_save,
            "n": self.combo_browse,
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

    # -------------------------------------------------- podglad pod kursorem
    def toggle_hover(self) -> None:
        self.hover_on = not self.hover_on
        if not self.hover_on:
            self._hover_actor.SetVisibility(False)
            self._hover_dot.SetVisibility(False)
            self._hover_idx = -1
        self.info(f"Podglad obszaru pod kursorem: "
                  f"{'wlaczony' if self.hover_on else 'wylaczony'}")
        self.pl.render()

    def _on_hover(self, *args) -> None:
        """Wolane przy kazdym drgnieciu myszy - wiec musi byc tanie i dlawione."""
        if not self.hover_on:
            return
        now = time.perf_counter()
        if now - self._hover_t < 0.03:            # najwyzej ~30 razy na sekunde
            return
        self._hover_t = now
        try:
            pos = self.pl.pick_mouse_position()
        except Exception:
            return
        if pos is None:
            return
        self.hover_at(np.asarray(pos, dtype=float))

    def hover_at(self, pos: np.ndarray) -> None:
        """Przesuwa podglad na wierzcholek siatki najblizszy wskazaniu.

        Gdy kursor jest poza cialem, picker zwraca punkt na plaszczyznie kamery -
        dlatego odrzucamy trafienia dalsze niz `MAX` od siatki, zamiast rysowac
        kule w powietrzu.
        """
        MAX = 0.06
        i = self.wire.nearest(pos)
        q = self.wire.points[i]
        if float(np.linalg.norm(q - pos)) > MAX:
            if self._hover_idx != -1:
                self._hover_actor.SetVisibility(False)
                self._hover_dot.SetVisibility(False)
                self._hover_idx = -1
                self.pl.render()
            return
        if i == self._hover_idx:
            return
        self._hover_idx = i
        self._hover_actor.SetPosition(*q)
        self._hover_dot.SetPosition(*q)
        self._hover_actor.SetVisibility(True)
        self._hover_dot.SetVisibility(True)
        d = self.wire.describe(q, index=i)
        txt = [f"{d['segment']}  {d['cun']:.1f} cun  az {d['az']:.0f}"
               f"   tkanka {d['flesh_mm']:.0f} mm"]
        if d["muscle"] is not None:
            side = {"L": "lewy", "R": "prawy", "C": ""}.get(d["muscle"]["side"], "")
            txt.append(f"pod spodem: {d['muscle']['name']}"
                       f"{(' (' + side + ')') if side else ''} "
                       f"{d['muscle_mm']:.0f} mm")
        near = np.linalg.norm(self._mpos - q, axis=1)
        j = np.flatnonzero(near <= DEFAULT_RADIUS)
        if len(j):
            j = j[np.argsort(near[j])][:5]
            txt.append("w zasiegu: " + ", ".join(
                f"{self.atlas.markers[k]['code']} ({near[k] * 1000:.0f} mm)" for k in j))
        self.pl.add_text("\n".join(txt), position=(320, 896), font_size=8,
                         color=THEME["accent"], name="hoverinfo")
        self.pl.render()

    # ------------------------------------------------------ pasek wlasciwosci
    def _build_role_index(self) -> None:
        """Kod punktu -> zbior klas, do ktorych nalezy wedlug definicji kanalowych."""
        self.role_of_code: dict = {}
        for code, rs in self.chan.roles.items():
            for r in rs:
                self.role_of_code.setdefault(code, set()).add(r["role"])
        for code in self.chan.crossings:
            self.role_of_code.setdefault(code, set()).add("crossing")
        self.codes_of_role: dict = {}
        for code, rr in self.role_of_code.items():
            for r in rr:
                self.codes_of_role.setdefault(r, set()).add(code)

    def _add_role_panel(self) -> None:
        self.pl.add_text("KLASYFIKACJA PUNKTOW", position=(12, 940), font_size=8,
                         color=THEME["dim"], name="hdr_roles")
        y = 918
        for role, label, color in ROLE_ORDER:
            n = len(self.codes_of_role.get(role, ()))
            if not n:
                continue
            self._role_widgets[role] = self.pl.add_checkbox_button_widget(
                self._toggle_role(role), value=False, position=(12, y),
                size=15, border_size=1, color_on=color, color_off="#2b3340")
            self.pl.add_text(f"{label}  ({n})", position=(33, y + 1), font_size=7,
                             color=THEME["text"], name=f"rl_{role}")
            y -= 19
        self._button_role("wszystkie", lambda: self.set_roles(
            {r for r, _, _ in ROLE_ORDER if self.codes_of_role.get(r)}), (12, y - 6))
        self._button_role("zaden", lambda: self.set_roles(set()), (12, y - 28))

    def _button_role(self, label: str, fn, position) -> None:
        holder = {}

        def cb(state):
            fn()
            w = holder.get("w")
            if w is not None:
                w.GetRepresentation().SetState(False)
            self.pl.render()

        holder["w"] = self.pl.add_checkbox_button_widget(
            cb, value=False, position=position, size=14, border_size=1,
            color_on=THEME["dim"], color_off="#2b3340")
        self.pl.add_text(label, position=(position[0] + 20, position[1] + 1),
                         font_size=7, color=THEME["dim"], name=f"rb_{label}")

    def _toggle_role(self, role: str):
        def cb(state: bool) -> None:
            (self.role_on.add if state else self.role_on.discard)(role)
            self._apply_roles()
        return cb

    def set_roles(self, roles: set) -> None:
        self.role_on = set(roles)
        for r, w in self._role_widgets.items():
            w.GetRepresentation().SetState(r in self.role_on)
        self._apply_roles()

    def _apply_roles(self) -> None:
        """Zaznaczone klasy dostaja wlasny kolor; reszta punktow schodzi na drugi plan.

        Punkt moze nalezec do kilku klas naraz (LU-7 jest i luo, i otwierajacy) -
        dostaje wtedy kolor pierwszej z listy, ale w opisie widac wszystkie.
        """
        for name in list(self._role_actors):
            self.pl.remove_actor(name, render=False)
        self._role_actors.clear()
        dim = bool(self.role_on)
        for a in self.point_actors.values():
            a.GetProperty().SetOpacity(0.12 if dim else 0.85)
            a.SetVisibility(self.show_points)
        for role, label, color in ROLE_ORDER:
            if role not in self.role_on:
                continue
            codes = self.codes_of_role.get(role, set())
            pts = [m["pos"] for m in self.atlas.markers if m["code"] in codes]
            if not pts:
                continue
            name = f"role_{role}"
            cloud = pv.PolyData(np.array(pts))
            self.pl.add_mesh(cloud.glyph(geom=self._sel_sphere, scale=False,
                                         orient=False),
                             color=color, name=name)
            self._role_actors[name] = True
        if self.role_on:
            total = len({c for r in self.role_on for c in self.codes_of_role.get(r, ())})
            names = ", ".join(l for r, l, _ in ROLE_ORDER if r in self.role_on)
            self.info(f"KLASYFIKACJA: {names}\n{total} punktow. "
                      f"Kliknij punkt, zeby zobaczyc jego opis.")
        else:
            self.info("Klasyfikacja wylaczona.")
        self.pl.render()

    def _properties(self, key: str) -> str:
        """Opis punktu: czym jest w kanale, gdzie lezy i co ma pod soba."""
        m = self.atlas.by_key.get(key)
        if m is None:
            return ""
        code = m["code"]
        info = self.atlas.channels[m["channel"]]
        w = self.wire.describe(m["pos"])
        rows = [f"{code}   {m.get('pinyin', '')}",
                f"kanal: {info['pl']} ({m['channel']})   strona: "
                f"{ {'R': 'prawa', 'L': 'lewa', 'C': 'posrodkowa'}.get(m['side'], '') }",
                f"dokladnosc kotwicy: {m.get('accuracy', '?')}"]
        mine = self.role_of_code.get(code, set())
        if mine:
            names = [l for r, l, _ in ROLE_ORDER if r in mine]
            rows.append("klasy: " + ", ".join(names))
        for r in self.chan.roles.get(code, []):
            if r["role"] in ("luo", "yuan", "xi", "front_mu", "back_shu", "he_sea"):
                rows.append(f"   {r['role']} kanalu {r['of']}")
        if code in self.chan.crossings:
            rows.append("   spotykaja sie tu: "
                        + ", ".join(self.chan.meets.get(code, [])))
        if code in self.chan.opening:
            rows.append("   otwiera: " + ", ".join(self.chan.opening[code]))
        rows.append(f"polozenie: {w['segment']}  {w['cun']:.1f} cun  az {w['az']:.0f}")
        rows.append(f"tkanka miekka: {w['flesh_mm']:.0f} mm")
        if w["muscle"] is not None:
            rows.append(f"miesien: {w['muscle']['name']} ({w['muscle_mm']:.0f} mm)")
        if m.get("note_pl"):
            rows.append(m["note_pl"][:96])
        return "\n".join(rows)

    def inspect(self, key: str) -> None:
        self._inspected = key
        self.pl.add_text(self._properties(key), position=(1130, 470), font_size=8,
                         color=THEME["text"], name="props")
        self.pl.add_mesh(pv.Sphere(radius=0.013, center=self.atlas.by_key[key]["pos"]),
                         color=THEME["sel"], name="inspected")
        self.pl.render()

    # ------------------------------------------------------------ kombinacje
    def toggle_combo(self) -> None:
        """[c] przelacza to, CO robi klikniecie: oznaczenie albo wybor punktu."""
        self.combo_mode = not self.combo_mode
        if not self.combo_mode:
            self.combo_clear(quiet=True)
            self.info("Tryb oznaczen napiecia. Klik stawia oznaczenie.")
        else:
            n = len(self.book)
            self.info(f"TRYB KOMBINACJI - klikaj punkty (2 albo 3).  "
                      f"[Enter] zapis   [n] przegladaj ({n} w dzienniczku)   [x] wyczysc")
        self.pl.render()

    def combo_clear(self, quiet: bool = False) -> None:
        self.combo = []
        self._browse = -1
        self._refresh_combo()
        if not quiet:
            self.info("Wybor wyczyszczony.")
            self.pl.render()

    def _combo_pick(self, pos) -> None:
        """Wybiera marker najblizszy klikniecia; czwarty klik zaczyna od nowa."""
        d = np.linalg.norm(self._mpos - pos, axis=1)
        j = int(np.argmin(d))
        if d[j] > 0.05:
            self.info("Za daleko od punktu - kliknij blizej markera.")
            self.pl.render()
            return
        key = self.atlas.markers[j]["key"]
        if len(self.combo) >= 3 or key in self.combo:
            self.combo = []
        self.combo.append(key)
        self._refresh_combo()
        if len(self.combo) < 2:
            self.info(f"Wybrano {key}.  Kliknij drugi punkt kombinacji.")
            self.pl.render()
            return
        d = combos.describe(self.tension, self.atlas.by_key, self.combo,
                            atlas=self.atlas, chan=self.chan)
        txt = [" + ".join(k.split("\u00b7")[0] for k in self.combo)]
        saved = self.book.find(self.combo)
        if saved is not None:
            txt[0] += f"   [w dzienniczku: {saved.name}]"
            if saved.note:
                txt[0] += f" - {saved.note}"
        for r in d["pairs"]:
            line = f"{r['from']['code']}->{r['to']['code']} {r['length_mm']:.0f} mm"
            if r.get("plan") and r["plan"].get("via"):
                line += f"  [{r['plan']['why']}]"
            if r["between"]:
                line += "   po drodze: " + ", ".join(q["code"] for q in r["between"][:6])
            if r["muscles"]:
                m = r["muscles"][0]
                line += f"   wzdluz: {m['name']} ({m['angle_deg']} st.)"
            txt.append(line)
        self.info("\n".join(txt))
        self.pl.render()

    def _refresh_combo(self) -> None:
        if len(self.combo) < 1:
            self._combo_poly.copy_from(pv.PolyData())
            return
        pts = np.array([self.atlas.by_key[k]["pos"] for k in self.combo])
        parts = [pv.PolyData(pts).glyph(geom=self._sel_sphere, scale=False, orient=False)]
        for i in range(len(pts)):                    # kazda para, nie lancuch:
            for j in range(i + 1, len(pts)):         # trojka ma pokazac trojkat
                poly, _plan = combos.channel_route(self.tension, self.atlas,
                                                   self.chan, self.combo[i], self.combo[j])
                if poly is None:
                    poly, _ = self.wire.route(pts[i], pts[j])
                if len(poly) > 2:
                    parts.append(pv.lines_from_points(poly).tube(radius=0.002, n_sides=8))
                else:
                    parts.append(pv.Tube(pointa=pts[i], pointb=pts[j], radius=0.002,
                                         n_sides=8))
        merged = parts[0].merge(parts[1:]) if len(parts) > 1 else parts[0]
        self._combo_poly.copy_from(merged)

    def combo_save(self) -> None:
        if len(self.combo) < 2:
            self.info("Kombinacja to co najmniej dwa punkty.")
            self.pl.render()
            return
        if self.book.find(self.combo) is not None:
            self.info("Ta kombinacja jest juz w dzienniczku.")
            self.pl.render()
            return
        c = self.book.add(list(self.combo))
        path = self.book.save()
        self.info(f"Zapisano kombinacje {c.name}  ({len(self.book)} w dzienniczku).\n"
                  f"Opis wpisz w {path.name} - pole 'note'.")
        _safe_print(f"\n{combos.report(self.tension, self.atlas.by_key, self.combo)}")
        self.pl.render()

    def combo_browse(self) -> None:
        """[n] pokazuje kolejna zapisana kombinacje na modelu."""
        if not len(self.book):
            self.info("Dzienniczek jest pusty - wybierz punkty i zapisz [Enter].")
            self.pl.render()
            return
        if not self.combo_mode:
            self.toggle_combo()
        self._browse = (self._browse + 1) % len(self.book)
        c = self.book.items[self._browse]
        self.combo = [k for k in c.keys if k in self.atlas.by_key]
        self._refresh_combo()
        txt = f"[{self._browse + 1}/{len(self.book)}]  {c.name}"
        if c.note:
            txt += f"\n{c.note}"
        self.info(txt)
        self.pl.render()

    # ------------------------------------------------------------------ klik
    def _on_pick(self, picked, *args) -> None:
        pos = np.asarray(getattr(picked, "points", picked)).reshape(-1, 3)
        if not len(pos):
            return
        if self.combo_mode:
            self._combo_pick(pos[0])
            return
        # Wlaczona klasyfikacja zmienia sens klikniecia: blisko podswietlonego
        # punktu ogladamy jego opis zamiast stawiac oznaczenie. Inaczej nie da sie
        # przejrzec klasy, nie zasmiecajac mapy napiec.
        if self.role_on:
            codes = {c for r in self.role_on for c in self.codes_of_role.get(r, ())}
            cand = [m for m in self.atlas.markers if m["code"] in codes]
            if cand:
                d = np.linalg.norm(np.array([m["pos"] for m in cand]) - pos[0], axis=1)
                j = int(np.argmin(d))
                if d[j] < 0.03:
                    self.inspect(cand[j]["key"])
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
        """Podmienia geometrie w istniejacych aktorach zamiast tworzyc nowe.

        Strzalki sa jednym glifem: gotowy grot powielamy wzdluz wektorow, zamiast
        budowac `pv.Arrow` dla kazdego odcinka - jeden taki obiekt kosztuje 2 ms,
        wiec przy dziesieciu oznaczeniach same groty zjadaly 18 ms na klik.
        """
        marks = self.tension.marks
        if not marks:
            self._marks_poly.copy_from(pv.PolyData())
            self._force_poly.copy_from(pv.PolyData())
            return
        pts = np.array([m.pos for m in marks])
        lev = np.array([float(m.level) for m in marks])
        cloud = pv.PolyData(pts)
        cloud["poziom"] = lev
        glyph = cloud.glyph(geom=self._sphere, scale=False, orient=False)
        self._marks_poly.copy_from(glyph)

        segs, mids, vecs = [], [], []
        for chain in self.tension.chain_marks():
            cp = np.array([m.pos for m in chain])
            route = [cp[0]]
            for x, y in zip(cp, cp[1:]):
                poly, _ = self.wire.route(x, y)
                route.extend(poly[1:] if len(poly) > 1 else [y])
            segs.append(pv.lines_from_points(np.array(route)).tube(radius=0.0026,
                                                                   n_sides=8))
            for a, b in zip(chain, chain[1:]):
                v = b.pos - a.pos
                n = float(np.linalg.norm(v))
                if n > 1e-6:
                    mids.append(a.pos + v * 0.5)
                    vecs.append(v / n * min(0.05, n * 0.5))
        if mids:
            heads = pv.PolyData(np.array(mids))
            heads["v"] = np.array(vecs)
            segs.append(heads.glyph(geom=self._arrow, orient="v", scale="v",
                                    factor=1.0))
        if segs:
            merged = segs[0].merge(segs[1:]) if len(segs) > 1 else segs[0]
            self._force_poly.copy_from(merged)
        else:
            self._force_poly.copy_from(pv.PolyData())

    # ------------------------------------------------------------------- start
    def show(self) -> None:
        """Otwiera okno, wymuszajac jawne polozenie i rozmiar.

        WSLg potrafi zaparkowac nowe okno na wspolrzednych -25600,-25600 i
        zostawic z niego kikut wielkosci paska tytulu: proces dziala, VTK
        raportuje pelny rozmiar, a na ekranie nie ma nic. Jawna pozycja jest
        podpowiedzia dla X11, ktora menedzer okien respektuje.
        """
        try:
            rw = self.pl.ren_win
            rw.SetPosition(60, 60)
            rw.SetSize(*self.pl.window_size)
        except Exception:
            pass
        self.pl.show()

    def screenshot(self, path: str) -> None:
        self.pl.screenshot(path)
