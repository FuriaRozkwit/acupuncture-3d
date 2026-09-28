"""Tryb edycji punktu: gizmo transformacji + kontekstowy panel wlasciwosci.

Wpina sie w SurfaceViewer Z ZEWNATRZ, bez zmieniania jego kodu: dokleja wlasne
obserwatory myszy do stylu interakcji z wyzszym priorytetem i przejmuje zdarzenie
tylko wtedy, gdy kursor lapie uchwyt albo przycisk panelu. Kazdy inny klik idzie
do przegladarki jak dotad (obrot kamery, wybor punktu, prawy panel).

Uchwyty to stopnie swobody KOTWICY, nie osie X/Y/Z (patrz anchor_edit.py):
    bursztynowy  - poziom (cun wzdluz konczyny, wysokosc na tulowiu, luk czaszki)
    turkusowy    - obwod (azymut wokol konczyny, odstep od linii posrodkowej)
    fioletowy    - trzeci stopien (przesuniecie w mm przy kotwicach kostnych)

Przeciagniecie przyciaga sie do kroku pola (0.1 cun, 2 st., 1 mm), bo takimi
jednostkami punkty opisuje sie klinicznie. Przebiegi kanalow NIE przeliczaja
sie na zywo (to kilka minut) - zamiast nich cienkie linie do sasiadow w kanale.

    [t]  wlacz / wylacz tryb edycji; potem kliknij punkt
    [v]  konfiguracje: fabryczna i profile - wybor domyslnej

Zapis trafia do PROFILU (acu3d/profiles.py), nigdy do fabrycznego points.json.
"""
from __future__ import annotations

import json
import textwrap
from pathlib import Path

import numpy as np
import pyvista as pv

from . import profiles
from .anchor_edit import AnchorEditor, get

ROOT = Path(__file__).resolve().parent.parent

HANDLE_M = 0.030                  # dlugosc uchwytu w kazda strone (m)
GRAB_PX = 12                      # jak blisko uchwytu trzeba kliknac (piksele)
CENTER_PX = 9                     # srodek gizma (tam uchwyty sie krzyzuja) nie lapie
# Gizmo wyzej niz linie kanalow (te leza ok. 4 mm nad skora): na tej samej
# wysokosci bursztynowy uchwyt "wzdluz" znikal w pomaranczowej linii ST na goleni.
LIFT_M = 0.009
OUTLINE = "#0b0f14"               # ciemny obrys - uchwyt odcina sie od kazdego koloru kanalu
PRIORITY = 10.0                   # wyzej niz obserwatory przegladarki (0.0)
COLORS = ("#ffb020", "#22d3ee", "#e879f9")
C_POINT, C_GHOST, C_LINK = "#ffffff", "#8b95a5", "#64748b"


def _hex(h: str) -> tuple:
    h = h.lstrip("#")
    return tuple(int(h[i:i + 2], 16) / 255.0 for i in (0, 2, 4))


class EditMode:
    def __init__(self, viewer):
        self.v = viewer
        self.pl = viewer.pl
        self.on = False
        self.ed: AnchorEditor | None = None
        self.handles = []                  # (dof, kierunek, mm/krok, kolor)
        self.pos = None                    # biezaca pozycja punktu
        self.drag = None                   # stan przeciagania uchwytu
        self.buttons = []                  # (aktor, akcja)
        self.intervals = []
        self._names = set()
        self._silhouettes = {}             # sylwetki nie maja nazw w PyVista - trzymamy aktory
        self.prof_on = False               # panel profili [v]
        self.pbuttons = []
        self._pnames = set()
        self._panel_bottom = None
        self._neighbours = self._load_neighbours()

    # --- wpiecie -----------------------------------------------------------
    def attach(self) -> "EditMode":
        v = self.v
        orig_select, orig_clear = v.select_point, v.clear_selection

        def select_point(key_or_code, *a, **k):
            ok = orig_select(key_or_code, *a, **k)
            if ok and self.on and getattr(v, "selected_point", None):
                self.edit(v.selected_point)
            return ok

        def clear_selection(*a, **k):
            self.close()
            return orig_clear(*a, **k)

        v.select_point, v.clear_selection = select_point, clear_selection
        self.pl.add_key_event("t", self.toggle)
        self.pl.add_key_event("v", self.toggle_profiles)
        try:
            self.pl.ren_win.AddObserver("WindowResizeEvent", self._on_resize)
        except Exception:
            pass

        style = getattr(v, "_style", None)
        if style is not None:
            self._tags = {}
            for ev, fn in (("LeftButtonPressEvent", self._press),
                           ("MouseMoveEvent", self._move),
                           ("LeftButtonReleaseEvent", self._release)):
                self._tags[ev] = style.AddObserver(ev, fn, PRIORITY)
            self._style = style
        return self

    def _consume(self, event: str) -> None:
        """Zatrzymuje zdarzenie - obserwatory przegladarki go nie dostana."""
        self._style.GetCommand(self._tags[event]).SetAbortFlag(1)

    def _take_press(self) -> None:
        """Przejmuje wcisniecie RAZEM z puszczeniem, ktore po nim nastapi.

        Przegladarka rozpoznaje klik po parze wcisniecie/puszczenie. Gdy wcisniecie
        zabierzemy, a puszczenie oddamy, przegladarka liczy klik od STAREJ pozycji
        wcisniecia - i wybierala przypadkowy punkt obok (w tescie: SJ-5 zamiast
        zostac przy edytowanym LU-7, z utrata niezapisanej zmiany).
        """
        self._consume("LeftButtonPressEvent")
        self._swallow = True
        if hasattr(self.v, "_press_pos"):
            self.v._press_pos = None

    # --- wlacz / wylacz ---------------------------------------------------
    def toggle(self) -> None:
        self.on = not self.on
        if self.on:
            key = getattr(self.v, "selected_point", None)
            if key:
                self.edit(key)
            else:
                self._say("TRYB EDYCJI: kliknij punkt, zeby go przesuwac  [t] wylacza")
        else:
            self.close()
            self._say("tryb edycji wylaczony")

    def edit(self, key: str) -> None:
        code, side = self._split(key)
        if self.ed is not None and self.ed.dirty and self.ed.code != code:
            self._say(f"niezapisana zmiana {self.ed.code} odrzucona")
        try:
            self.ed = AnchorEditor(self.v.atlas, code, side)
        except StopIteration:
            return
        self._rebuild()

    def close(self) -> None:
        self.ed, self.drag, self.handles = None, None, []
        for act in self._silhouettes.values():
            self.pl.remove_actor(act)
        self._silhouettes.clear()
        for n in list(self._names):
            self.pl.remove_actor(n)
        self._names.clear()
        self.buttons.clear()
        self.pl.render()

    # --- geometria gizma ---------------------------------------------------
    def _rebuild(self) -> None:
        """Pelne przeliczenie: pozycja, kierunki uchwytow, rysunek, panel."""
        ed = self.ed
        self.pos = ed.position()
        self.handles = [(d, vec, mm, COLORS[i % len(COLORS)])
                        for i, (d, vec, mm) in enumerate(ed.handles())]
        try:
            self._normal = np.asarray(self.v.atlas.anatomy.surface_normal(self.pos), float)
        except Exception:
            self._normal = np.zeros(3)
        self.intervals = ed.neighbour_intervals()
        self._draw()
        self._panel()
        self.pl.render()

    def _add(self, name: str, mesh, **kw) -> None:
        self.pl.add_mesh(mesh, name=name, pickable=False, reset_camera=False, **kw)
        self._names.add(name)

    def _draw(self) -> None:
        p = self.pos + self._normal * LIFT_M
        for i, (d, vec, mm, color) in enumerate(self.handles):
            a, b = p - vec * HANDLE_M, p + vec * HANDLE_M
            shaft = pv.Line(a, b).tube(radius=0.0014, n_sides=16)
            tip = pv.Cone(center=b, direction=vec, height=0.009, radius=0.0042, resolution=18)
            self._add(f"edit_shaft_{i}", shaft, color=color)
            self._add(f"edit_tip_{i}", tip, color=color)
            # Obrys jako SYLWETKA liczona z kamery, nie druga, grubsza bryla: grubsza
            # rurka otaczala kolorowa i uchwyty wychodzily czarne.
            for part, mesh in (("s", shaft), ("t", tip)):
                self._outline(f"edit_line_{part}{i}", mesh)
            self._add(f"edit_tail_{i}", pv.Sphere(radius=0.0022, center=a), color=color)
        self._add("edit_point", pv.Sphere(radius=0.0048, center=p), color=C_POINT)
        # gdzie punkt byl przed edycja - zeby bylo widac, o ile sie ruszyl
        orig = self.ed.position(self.ed.original)
        if np.linalg.norm(orig - self.pos) > 1e-4:
            self._add("edit_ghost", pv.Sphere(radius=0.0035, center=orig + self._normal * LIFT_M),
                      color=C_GHOST, opacity=0.6)
        else:
            self._drop("edit_ghost")
        # proste linie do sasiadow w kanale: przebieg przeliczy sie dopiero po zapisie
        for j, other in enumerate(self._neighbours.get(self.ed.code, [])):
            q = self._marker_pos(other)
            if q is None:
                self._drop(f"edit_link_{j}")
                continue
            self._add(f"edit_link_{j}", pv.Line(p, q), color=C_LINK, line_width=2,
                      opacity=0.8)

    def _outline(self, name: str, mesh) -> None:
        """Ciemny kontur uchwytu - odcina go od kanalu w tym samym kolorze."""
        self._drop(name)
        try:
            act = self.pl.add_silhouette(mesh, color=OUTLINE, line_width=3)
        except Exception:
            return
        self._silhouettes[name] = act
        self._names.add(name)

    def _drop(self, name: str) -> None:
        act = self._silhouettes.pop(name, None)
        if act is not None:
            self.pl.remove_actor(act)
            self._names.discard(name)
            return
        if name in self._names:
            self.pl.remove_actor(name)
            self._names.discard(name)

    def _marker_pos(self, code: str):
        side = self.ed.side
        for key in (code, f"{code}·{'L' if side == 'L' else 'R'}"):
            m = self.v.atlas.by_key.get(key)
            if m is not None:
                return np.asarray(m["pos"], float)
        return None

    # --- panel wlasciwosci --------------------------------------------------
    # Panel stoi u gory posrodku: lewy gorny rog zajmuje legenda, lewy dolny karta
    # opisu punktu, prawa kolumna panel kliniczny. Na jasnej skorze tekst bez tla
    # byl nieczytelny, stad ciemne tlo i czcionka o stalej szerokosci - wtedy
    # wiersze maja rowne ramki, a przyciski staja tam, gdzie konczy sie tekst.
    ROW_H = 21
    BG = "#10151b"

    def _char_w(self) -> float:
        cw = getattr(self, "_cw", None)
        if cw:
            return cw
        act = self.pl.add_text("M" * 20, position=(0, 0), font_size=9, name="edit_probe")
        act.GetTextProperty().SetFontFamilyToCourier()
        bb = [0.0] * 4
        try:
            act.GetBoundingBox(self.pl.renderer, bb)
            cw = (bb[1] - bb[0]) / 20.0
        except Exception:
            cw = 0.0
        self.pl.remove_actor("edit_probe")
        self._cw = cw if cw > 3 else 7.5
        return self._cw

    def _text(self, name: str, txt: str, x: int, y: int, color="#e8edf3", bg=BG, frame=None):
        act = self.pl.add_text(txt, position=(x, y), font_size=9, color=color, name=name)
        self._names.add(name)
        tp = act.GetTextProperty()
        tp.SetFontFamilyToCourier()
        if bg is not None:
            tp.SetBackgroundColor(*_hex(bg))
            tp.SetBackgroundOpacity(0.9)
        if frame:
            tp.SetFrame(1)
            tp.SetFrameColor(*_hex(frame))
            tp.SetFrameWidth(1)
        return act

    def _button(self, name: str, label: str, x: int, y: int, action, accent="#5b6778"):
        act = self._text(name, f" {label} ", x, y, bg="#1f2630", frame=accent)
        self.buttons.append((act, action))
        return act

    def _origin(self):
        try:
            W, H = self.pl.ren_win.GetSize()
        except Exception:
            W, H = self.pl.window_size
        width = int(self.WIDTH_CH * self._char_w())
        x = int(np.clip(W / 2 - width / 2, 400, max(400, W - 280 - width)))
        return x, H - 34

    WIDTH_CH = 56

    def _fmt(self, d, val: float) -> str:
        if d.unit == "st.":
            return f"{val:.1f}"
        return f"{val:.3f}" if d.path[-1] == "z_frac" else f"{val:.2f}"

    def _dof_row(self, d) -> str:
        return f" {d.label:<24s}{self._fmt(d, get(self.ed.anchor, d)):>8s} {d.unit:<4s}"

    def _panel(self) -> None:
        for n in [n for n in self._names if n.startswith("edit_ui_")]:
            self.pl.remove_actor(n)
            self._names.discard(n)
        self.buttons.clear()

        ed, W = self.ed, self.WIDTH_CH
        x, y = self._origin()
        cw = self._char_w()
        pad = lambda t: t.ljust(W)[:W]          # noqa: E731 - rowne ramki wierszy

        star = "  *niezapisane" if ed.dirty else ""
        side = {"R": "prawa", "L": "lewa", "C": "srodek"}[ed.side]
        self._text("edit_ui_head", pad(f" EDYCJA KOTWICY  {ed.code} ({side}){star}"), x, y,
                   color="#ffe36a")
        y -= self.ROW_H
        prof = getattr(self.v.atlas, "profile", {}) or {}
        pname = prof.get("name", profiles.FACTORY)
        n_changed = len(prof.get("changed", []))
        if ed.code in prof.get("stale", []):
            state, scol = "fabryka zmienila sie od edycji!", "#fdba74"
        elif ed.is_factory:
            state, scol = "punkt fabryczny", "#8fa3b8"
        else:
            state, scol = "punkt zmieniony", "#c4b5fd"
        head = f" profil: {pname}" + (f" ({n_changed} zm.)" if pname != profiles.FACTORY else "")
        # Stan punktu po prawej stronie; nazwa profilu po lewej z lagodnym skracaniem,
        # dzieki czemu stan ("punkt fabryczny", "punkt zmieniony") NIGDY nie jest uciety.
        avail = W - len(state) - 2
        if len(head) > avail:
            suffix = f" ({n_changed} zm.)" if pname != profiles.FACTORY else ""
            pname_max = max(6, avail - len(" profil: ") - len(suffix) - 1)
            pname_short = pname[:pname_max] + "…"
            head = f" profil: {pname_short}{suffix}"
            if len(head) > avail:
                head = head[:avail]
        gap = max(1, W - len(head) - len(state))
        self._text("edit_ui_prof", pad(f"{head}{' ' * gap}{state}"), x, y, color=scol)
        y -= self.ROW_H
        for k, w in enumerate(textwrap.wrap(self._anchor_line(ed.anchor), width=W - 2)[:2]):
            self._text(f"edit_ui_anchor_{k}", pad(" " + w), x, y, color="#aab6c4")
            y -= self.ROW_H
        y -= 6

        # jeden wiersz na stopien swobody, przyciski tuz za tekstem
        for i, d in enumerate(ed.dofs):
            color = next((h[3] for h in self.handles if h[0] == d), "#6b7280")
            row = self._dof_row(d)
            self._text(f"edit_ui_dof_{i}", row, x, y, color=color)
            bx = int(x + len(row) * cw + 8)
            self._button(f"edit_ui_minus_{i}", "-", bx, y, lambda d=d: self._step(d, -1), color)
            self._button(f"edit_ui_plus_{i}", "+", bx + int(3 * cw + 10), y,
                         lambda d=d: self._step(d, +1), color)
            y -= self.ROW_H + 2
        y -= 6

        # odstepy do sasiadow w kanale, liczone na zywo z biezacej kotwicy
        if self.intervals:
            self._text("edit_ui_int_h", pad(" odstep w kanale        kanon  teraz"), x, y,
                       color="#8fa3b8")
            y -= self.ROW_H
            for k, r in enumerate(self.intervals):
                got = "  -" if r["model"] is None else f"{r['model']:5.2f}"
                if r["canon"] is None:
                    canon, verdict, color = "   - ", r["basis"][:10], "#8b95a5"
                else:
                    canon = f"{r['canon']:5.2f}"
                    diff = None if r["model"] is None else r["model"] - r["canon"]
                    verdict = "?" if diff is None else ("ok" if abs(diff) <= 0.5 else f"{diff:+.2f}")
                    color = "#86efac" if verdict == "ok" else "#fdba74"
                self._text(f"edit_ui_int_{k}",
                           pad(f"   {r['from']:>6s} -> {r['to']:<6s}   {canon}  {got}  {verdict}"),
                           x, y, color=color)
                y -= self.ROW_H
            y -= 6

        # skad sie wziela obecna kotwica
        src = self.ed.point.get("az_source") or "(brak opisu pochodzenia)"
        for k, w in enumerate(textwrap.wrap(f"zrodlo: {src}", width=W - 2)[:2]):
            self._text(f"edit_ui_src_{k}", pad(" " + w), x, y, color="#7c8796")
            y -= self.ROW_H
        y -= 8
        self._button("edit_ui_save", "Zapisz", x, y, self._save,
                     "#3fb950" if ed.dirty else "#5b6778")
        self._button("edit_ui_undo", "Cofnij", x + int(10 * cw), y, self._undo)
        self._button("edit_ui_factory", "Fabryczny", x + int(19 * cw), y, self._factory,
                     "#5b6778" if ed.is_factory else "#c4b5fd")
        self._button("edit_ui_close", "Zamknij [t]", x + int(31 * cw), y, self.toggle)
        self._panel_bottom = y
        if self.prof_on:
            self._profiles_panel()

    @staticmethod
    def _anchor_line(a: dict) -> str:
        skip = ("region",)
        parts = [f"{k} {v}" for k, v in a.items() if k not in skip and not isinstance(v, dict)]
        parts += [f"{k}.{kk} {vv}" for k, v in a.items() if isinstance(v, dict) for kk, vv in v.items()]
        return "kotwica: " + ", ".join(parts)

    def toggle_profiles(self) -> None:
        self.prof_on = not self.prof_on
        if self.prof_on:
            self._profiles_panel()
        else:
            self._clear_profiles()
        self.pl.render()

    def _clear_profiles(self) -> None:
        for n in list(self._pnames):
            self.pl.remove_actor(n)
        self._pnames.clear()
        self.pbuttons.clear()

    def _profiles_panel(self) -> None:
        """Lista konfiguracji: klik ustawia domyslna (wchodzi przy nastepnym starcie).

        Przelaczanie na zywo nie jest mozliwe bez przebudowy calej sceny - kanaly,
        etykiety i warstwy buduje przegladarka przy starcie - wiec uczciwiej jest
        ustawic domyslny profil i powiedziec, ze zadziala po ponownym uruchomieniu.
        """
        self._clear_profiles()
        W = self.WIDTH_CH
        cw = self._char_w()
        x, top = self._origin()
        y = (self._panel_bottom - 34) if (self.ed is not None and self._panel_bottom) else top
        loaded = (getattr(self.v.atlas, "profile", {}) or {}).get("name", profiles.FACTORY)
        default = profiles.active_name()

        def pad(t):
            return t.ljust(W)[:W]

        def row(name, txt, color, px=None, bg=None, frame=None, action=None):
            act = self.pl.add_text(txt, position=(px if px is not None else x, y),
                                   font_size=9, color=color, name=name)
            self._pnames.add(name)
            tp = act.GetTextProperty()
            tp.SetFontFamilyToCourier()
            tp.SetBackgroundColor(*_hex(bg or self.BG))
            tp.SetBackgroundOpacity(0.9)
            if frame:
                tp.SetFrame(1)
                tp.SetFrameColor(*_hex(frame))
                tp.SetFrameWidth(1)
            if action:
                self.pbuttons.append((act, action))
            return act

        def_head = f"domyslna: {default}"
        gap = max(2, W - len(" KONFIGURACJE") - len(def_head))
        row("prof_ui_head", pad(f" KONFIGURACJE{' ' * gap}{def_head}"), "#ffe36a")
        y -= self.ROW_H + 4
        for k, name in enumerate([profiles.FACTORY] + profiles.names()):
            n = 0 if name == profiles.FACTORY else len(profiles.load(name)["points"])
            mark = "*" if name == default else " "
            note = "wczytana" if name == loaded else ""
            count = "" if name == profiles.FACTORY else f"{n:3d} pkt"
            row(f"prof_ui_p{k}", pad(f" [{mark}] {name:<22s}{count:>8s}  {note}"),
                "#e8edf3" if name == default else "#aab6c4", bg="#1f2630",
                frame="#3fb950" if name == default else "#5b6778",
                action=lambda name=name: self._use_profile(name))
            y -= self.ROW_H + 3
        y -= 4
        row("prof_ui_new", " Nowy z wczytanego ", "#e8edf3", bg="#1f2630", frame="#5b6778",
            action=self._new_profile)
        row("prof_ui_close", " Zamknij [v] ", "#e8edf3", px=x + int(24 * cw), bg="#1f2630",
            frame="#5b6778", action=self.toggle_profiles)
        y -= self.ROW_H + 2
        hint = " klik = domyslna; wchodzi przy nastepnym uruchomieniu"
        hint_lines = [hint] if len(hint) <= W else textwrap.wrap(hint, width=W)
        for k, hline in enumerate(hint_lines):
            row(f"prof_ui_hint_{k}" if len(hint_lines) > 1 else "prof_ui_hint",
                pad(hline), "#7c8796")
            y -= self.ROW_H

    def _use_profile(self, name: str) -> None:
        profiles.set_active(name)
        self._profiles_panel()
        self.pl.render()
        self._say(f"domyslna konfiguracja: '{name}' - wczyta sie przy nastepnym uruchomieniu")

    def _new_profile(self) -> None:
        import time
        loaded = (getattr(self.v.atlas, "profile", {}) or {}).get("name", profiles.FACTORY)
        name = profiles.create(f"profil-{time.strftime('%m%d-%H%M%S')}", copy_from=loaded)
        profiles.set_active(name)
        # dalsze zapisy w tej sesji ida do nowego profilu
        if isinstance(getattr(self.v.atlas, "profile", None), dict):
            self.v.atlas.profile["name"] = name
        self._profiles_panel()
        if self.ed is not None:
            self._panel()
        self.pl.render()
        self._say(f"utworzono profil '{name}' (kopia '{loaded}') - jest teraz domyslny")

    def _on_resize(self, *_a) -> None:
        if self.ed is not None:
            self._panel()
        elif self.prof_on:
            self._profiles_panel()

    # --- akcje --------------------------------------------------------------
    def _step(self, dof, n: int) -> None:
        if self.ed is None:
            return
        self.ed.step(dof, n)
        self._rebuild()

    def _undo(self) -> None:
        if self.ed is None:
            return
        self.ed.reset()
        self._rebuild()
        self._say(f"{self.ed.code}: przywrocona kotwica sprzed edycji")

    def _save(self) -> None:
        if self.ed is None:
            return
        if not self.ed.dirty:
            self._say("nic do zapisania")
            return
        rec = self.ed.save("przesuniete w przegladarce")
        self._rebuild()
        what = "wrocil do fabrycznego" if rec["factory"] else "zapisany"
        self._say(f"{rec['code']} {what} w profilu '{rec['profile']}' - fabryczna konfiguracja "
                  f"nietknieta; przebieg kanalu przeliczy sie przy nastepnym uruchomieniu")
        print(f"[edycja] {rec['code']} -> profil '{rec['profile']}': {rec['source']}")

    def _factory(self) -> None:
        """Ustawia kotwice fabryczna punktu (do zapisania przyciskiem Zapisz)."""
        if self.ed is None:
            return
        self.ed.to_factory()
        self._rebuild()
        self._say(f"{self.ed.code}: kotwica fabryczna - [Zapisz] usunie zmiane z profilu")

    def _say(self, text: str) -> None:
        try:
            self.v.info(text)
        except Exception:
            print(text)

    # --- mysz ---------------------------------------------------------------
    def _mouse(self):
        return self.pl.iren.interactor.GetEventPosition()

    def _to_screen(self, p) -> np.ndarray:
        r = self.pl.renderer
        r.SetWorldPoint(float(p[0]), float(p[1]), float(p[2]), 1.0)
        r.WorldToDisplay()
        return np.array(r.GetDisplayPoint()[:2])

    def _hit_button(self, x, y):
        # GetBoundingBox zwraca ramke WZGLEDEM polozenia aktora, nie w pikselach
        # ekranu - bez dodania pozycji kazdy przycisk "lezal" w lewym dolnym rogu
        bbox = [0.0] * 4
        for act, action in self.buttons + self.pbuttons:
            try:
                act.GetBoundingBox(self.pl.renderer, bbox)
                px, py = act.GetPosition()
            except Exception:
                continue
            if px + bbox[0] - 2 <= x <= px + bbox[1] + 2 and py + bbox[2] - 2 <= y <= py + bbox[3] + 2:
                return action
        return None

    def _hit_handle(self, x, y):
        """Uchwyt, ktorego odcinek na ekranie przechodzi najblizej kursora."""
        m = np.array([x, y], float)
        p = self.pos + self._normal * LIFT_M
        # W srodku wszystkie uchwyty sie krzyzuja - tam wybor bylby losowy (decydowal
        # szum liczb), wiec srodek nie lapie zadnego uchwytu.
        if np.linalg.norm(m - self._to_screen(p)) < CENTER_PX:
            return None
        best, best_d = None, GRAB_PX
        for i, (d, vec, mm, color) in enumerate(self.handles):
            a, b = self._to_screen(p - vec * HANDLE_M), self._to_screen(p + vec * HANDLE_M)
            ab = b - a
            t = np.clip(np.dot(m - a, ab) / max(np.dot(ab, ab), 1e-9), 0.0, 1.0)
            dist = np.linalg.norm(m - (a + t * ab))
            if dist < best_d:
                best, best_d = i, dist
        return best

    def _press(self, obj, _ev) -> None:
        x, y = self._mouse()
        action = self._hit_button(x, y)
        if action is not None:
            self._take_press()
            action()
            return
        if not self.on or self.ed is None:
            return
        p = self.pos + self._normal * LIFT_M
        if np.linalg.norm(np.array([x, y], float) - self._to_screen(p)) < CENTER_PX:
            self._take_press()          # klik w sam edytowany punkt: nic nie rob
            return
        i = self._hit_handle(x, y)
        if i is None:
            return                      # nie nasz klik - przegladarka robi swoje
        self._take_press()
        d, vec, mm, color = self.handles[i]
        p = self.pos
        # ile pikseli ekranu to jeden krok pola, i w ktora strone
        s0 = self._to_screen(p)
        s1 = self._to_screen(p + vec * (mm / 1000.0))
        step_px = s1 - s0
        self.drag = {"i": i, "dof": d, "start": np.array([x, y], float),
                     "value": get(self.ed.anchor, d), "step_px": step_px,
                     "steps": 0, "vec": vec, "mm": mm, "p0": p.copy()}

    def _move(self, obj, _ev) -> None:
        if self.drag is None:
            return
        self._consume("MouseMoveEvent")
        g = self.drag
        m = np.array(self._mouse(), float)
        L = float(np.dot(g["step_px"], g["step_px"]))
        if L < 1e-6:
            return
        steps = int(round(float(np.dot(m - g["start"], g["step_px"])) / L))
        if steps == g["steps"]:
            return
        g["steps"] = steps
        self.ed.set_value(g["dof"], g["value"] + steps * g["dof"].step)
        # w trakcie ciagniecia tylko przesuwamy punkt; kierunki uchwytow
        # przeliczamy dopiero po puszczeniu, bo to kilka rzutow promieni
        self.pos = self.ed.position()
        self._draw()
        self._panel_value_only()
        self.pl.render()

    def _release(self, obj, _ev) -> None:
        if getattr(self, "_swallow", False):
            self._swallow = False
            self._consume("LeftButtonReleaseEvent")
        if self.drag is None:
            return
        self.drag = None
        self._rebuild()

    def _panel_value_only(self) -> None:
        """Podczas przeciagania odswiezamy tylko liczbe - pelny panel po puszczeniu."""
        d = self.drag["dof"]
        name = f"edit_ui_dof_{self.ed.dofs.index(d)}"
        act = self.pl.renderer.actors.get(name) if hasattr(self.pl.renderer, "actors") else None
        if act is not None:
            act.SetInput(self._dof_row(d))

    # --- pomocnicze ---------------------------------------------------------
    @staticmethod
    def _split(key: str):
        for sep in ("·", ":"):
            if sep in key:
                code, side = key.rsplit(sep, 1)
                return code, side
        return key, "C"

    @staticmethod
    def _load_neighbours() -> dict:
        """Poprzedni i nastepny punkt w kazdej grupie kanalu (z meridians.json)."""
        f = ROOT / "data" / "meridians.json"
        out = {}
        if not f.exists():
            return out
        for groups in json.loads(f.read_text(encoding="utf-8"))["paths"].values():
            for g in groups:
                for a, b in zip(g[:-1], g[1:]):
                    out.setdefault(a, []).append(b)
                    out.setdefault(b, []).append(a)
        return {k: sorted(set(v)) for k, v in out.items()}


def attach(viewer) -> EditMode:
    """Dokleja tryb edycji do dzialajacej przegladarki. Zwraca obiekt trybu."""
    return EditMode(viewer).attach()
