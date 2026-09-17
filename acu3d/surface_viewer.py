"""Podglad powierzchni: sama powloka ciala, do ogladania jak statuetki.

Bez warstw, bez siatki cunowej i bez stawiania oznaczen. Jest tylko skora
w kolorze bialej porcelany, kursor lezacy na powierzchni i punkty, ktore mozna
wlaczyc do sprawdzenia lokalizacji.

Kamera ma naśladowac badanie wlasnymi oczami:

  * kolko myszy przybliza DO MIEJSCA POD KURSOREM, a nie do srodka ekranu -
    tak jak nachyla sie glowe ku temu, co sie oglada;
  * [f] ustawia punkt pod kursorem jako os obrotu, wiec obraca sie wokol
    ogladanej struktury, a nie wokol srodka sylwetki;
  * lewy przycisk obraca, Shift + lewy przesuwa.

Kursor to pierscien ulozony stycznie do skory (wzdluz normalnej w miejscu
trafienia). Przesuwa sie go przez macierz aktora, bez przebudowy geometrii,
wiec nadaza za mysza, zamiast przeskakiwac po wezlach siatki.
"""
from __future__ import annotations

import sys
import time

import numpy as np
import pyvista as pv
import vtk

from .model import Atlas

THEME = {
    "skin": "#ece6dc",        # biala porcelana z cieplym odcieniem
    "bg_top": "#3a3f47",
    "bg_bottom": "#15181d",
    # kursor w kolorze spoza palety przemian, zeby nie mylil sie z zadnym kanalem
    "ring": "#36c5d6",
    "text": "#e6e9ee",
    "dim": "#9aa3b1",
    "off": "#3a404a",
}

# Wu Xing: przemiana -> (nazwa, odcien yin, odcien yang).
#
# Kolory klasyczne to drewno zielone, ogien czerwony, ziemia zolta, metal bialy,
# woda czarna. Dwa trzeba bylo przesunac, bo postac jest biala: metal jako srebro
# (biel na bieli bylaby niewidoczna), woda jako granat bliski czerni (czysta czern
# zlewala sie z tlem legendy). Kanal yin ma odcien ciemniejszy, yang jasniejszy -
# para wewnetrzno-zewnetrzna dzieli barwe, a rozni sie tonem.
ELEMENTS = {
    "Wood":  ("DREWNO", "#1f7a45", "#4cb574"),
    "Fire":  ("OGIEN",  "#b3261e", "#e25b4f"),
    "Earth": ("ZIEMIA", "#a8741a", "#dcab3f"),
    "Metal": ("METAL",  "#5f6b77", "#9aa5b1"),
    # woda: granat zamiast czerni - przy #16202b zaznaczone pole w legendzie
    # bylo ciemniejsze od odznaczonego i zaznaczenia nie bylo widac
    "Water": ("WODA",   "#2c4a6e", "#5f86b3"),
    "-":     ("POZA WU XING", "#6a4c9c", "#a58bd0"),
}
ELEMENT_ORDER = ("Wood", "Fire", "Earth", "Metal", "Water", "-")
YIN = frozenset({"LU", "SP", "HT", "KI", "PC", "LR", "REN"})
LINE_RADIUS = 0.0011          # m - rurka przebiegu kanalu
AIR_GAP = 0.010               # m - odcinek, ktorego wnetrze jest dalej od skory, idzie powietrzem


def split_through_air(poly: np.ndarray, distance_to_skin, gap: float = AIR_GAP) -> list:
    """Tnie przebieg tam, gdzie odcinek miedzy probkami przechodzi przez powietrze.

    Kryterium NIE jest dlugosc odcinka, tylko to, czy jego WNETRZE lezy na skorze.
    Konce kazdego odcinka sa dociagniete do skory, wiec sprawdzanie samych probek
    niczego nie wykrywa - wszystkie 30 przebiegow miesci sie wtedy w 6 mm.

    Cieciwy przez powietrze wystepuja w trzech miejscach (zmierzone):
      * dlon - biodro: SI przeskakuje o 120 mm z malego palca na biodro, obok
        ktorego zwisa reka (srodek 29 mm od skory), podobnie PC;
      * glowa: GB ma odcinki po 103-113 mm przecinajace krzywizne czaszki
        (srodek 13-14 mm od skory);
      * krocze: KI, LR, REN przeskakuja szczeline miedzy udami (10-14 mm).

    Sprawdzamy punkty 1/4, 1/2 i 3/4 kazdego odcinka. Zamiast cieciwy zostaje
    przerwa - widac wtedy, ze w tym miejscu przebieg jest wadliwy (zrodlo bledu:
    zbyt rzadkie probkowanie i dociaganie krzywych w model.py).
    """
    poly = np.asarray(poly, dtype=float)
    if len(poly) < 2:
        return []
    cuts = []
    for i in range(len(poly) - 1):
        a, b = poly[i], poly[i + 1]
        if any(distance_to_skin(a + (b - a) * t) > gap for t in (0.25, 0.5, 0.75)):
            cuts.append(i + 1)
    return [p for p in np.split(poly, cuts) if len(p) >= 2]
RING_OUTER = 0.009            # m - pierscien kursora
RING_INNER = 0.0068
POINT_RADIUS = 0.0032         # m - znacznik punktu do weryfikacji
NEAR_POINT = 0.015            # m - od tej odleglosci podajemy kod punktu pod kursorem
MIN_DISTANCE = 0.035          # m - blizej kamera nie podchodzi do powierzchni
ZOOM_STEP = 0.82              # ulamek odleglosci zostajacy po jednym kroku kolka


def _safe_print(text: str) -> None:
    enc = getattr(sys.stdout, "encoding", None) or "utf-8"
    sys.stdout.write(text.encode(enc, errors="replace").decode(enc, errors="replace") + "\n")
    sys.stdout.flush()


def _align_z_to(n: np.ndarray) -> np.ndarray:
    """Macierz obrotu przeprowadzajaca os Z na wektor n (wzor Rodriguesa)."""
    z = np.array([0.0, 0.0, 1.0])
    n = n / max(np.linalg.norm(n), 1e-12)
    v = np.cross(z, n)
    c = float(np.dot(z, n))
    if np.linalg.norm(v) < 1e-9:
        return np.eye(3) if c > 0 else np.diag([1.0, -1.0, -1.0])
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx * (1.0 / (1.0 + c))


class SurfaceViewer:
    def __init__(self, atlas: Atlas, off_screen: bool = False, window_size=(1500, 960)):
        self.atlas = atlas
        # Normalne liczone od nowa bez rozcinania na krawedziach. Siatka ma
        # 90 tys. trojkatow i przy domyslnym kacie ostrosci cieniowanie
        # lamalo sie na wypuklosciach (kolano, lopatki) w widoczne fasety.
        # Kopia, zeby nie zmieniac siatki, na ktorej atlas stawia punkty.
        # Normalne z kolejnosci wierzcholkow, BEZ auto_orient i bez obracania.
        # Rozstrzygniete pomiarem jasnosci renderu, nie geometria: ten wariant daje
        # cialo o jasnosci 176/255, a obrocony - 86/255 (ciemnoszara postac).
        # Test "normalna na tulowiu wskazuje od osi" wybieral wlasnie ten ciemny
        # wariant, wiec go usunieto. Przyczyny rozjazdu nie ustalono - kryterium
        # jest to, co widac na ekranie (tools: pomiar pikseli ciala).
        self.skin = atlas.layers["skin"].compute_normals(
            auto_orient_normals=False, consistent_normals=True,
            split_vertices=False, feature_angle=180.0, inplace=False)
        self.pl = pv.Plotter(off_screen=off_screen, window_size=window_size,
                             lighting="light kit",
                             title="Powierzchnia ciala - podglad")
        self.pl.set_background(THEME["bg_bottom"], top=THEME["bg_top"])
        self.show_labels = False
        self.channels_on: set = set()     # zaczynamy od czystej postaci
        self._mpos = np.array([m["pos"] for m in atlas.markers], dtype=float)
        self._mchan = np.array([m["channel"] for m in atlas.markers])
        self.channel_color = self._palette()
        self._last_move = 0.0
        self._last_code: str | None = None
        self._build()

    # ------------------------------------------------------------------ scena
    def _build(self) -> None:
        self.skin_actor = self.pl.add_mesh(
            self.skin, color=THEME["skin"], smooth_shading=True,
            ambient=0.22, diffuse=0.78, specular=0.18, specular_power=22,
            name="skin")

        # Kazdy kanal to osobne aktory: punkty, przebieg i kody. Zaznaczenie
        # w legendzie tylko przelacza widocznosc - nic sie nie przebudowuje.
        self.channel_actors: dict = {}
        sphere = pv.Sphere(radius=POINT_RADIUS, theta_resolution=12, phi_resolution=12)
        air_loc = vtk.vtkStaticCellLocator()
        air_loc.SetDataSet(self.skin)
        air_loc.BuildLocator()
        _c, _cid, _sub, _d2 = [0.0, 0.0, 0.0], vtk.reference(0), vtk.reference(0), vtk.reference(0.0)

        def to_skin(p) -> float:
            air_loc.FindClosestPoint([float(p[0]), float(p[1]), float(p[2])], _c, _cid, _sub, _d2)
            return float(_d2) ** 0.5

        self.cut_segments = 0
        for ch, info in self.atlas.channels.items():
            color = self.channel_color[ch]
            mask = self._mchan == ch
            actors = []
            if mask.any():
                pts = pv.PolyData(self._mpos[mask])
                actors.append(self.pl.add_mesh(
                    pts.glyph(geom=sphere, scale=False, orient=False),
                    color=color, smooth_shading=True, name=f"pts_{ch}"))
            tubes = []
            for side in ("R", "L"):
                for poly in self.atlas.meridian_lines(ch, side):
                    if len(poly) <= 3:
                        continue
                    pieces = split_through_air(poly, to_skin)
                    self.cut_segments += max(len(pieces) - 1, 0)
                    for piece in pieces:
                        tubes.append(pv.lines_from_points(piece).tube(
                            radius=LINE_RADIUS, n_sides=8))
                if info.get("clock", 0) >= 13:        # REN i DU sa pojedyncze
                    break
            if tubes:
                merged = tubes[0].merge(tubes[1:]) if len(tubes) > 1 else tubes[0]
                actors.append(self.pl.add_mesh(merged, color=color, smooth_shading=True,
                                               name=f"mer_{ch}"))
            labels = None
            if mask.any():
                labels = self.pl.add_point_labels(
                    self._mpos[mask], [m["key"] for m in self.atlas.markers
                                       if m["channel"] == ch],
                    font_size=10, text_color=THEME["text"], shape=None,
                    show_points=False, always_visible=False, name=f"lab_{ch}")
                labels.SetVisibility(False)
                labels.SetPickable(False)
            for act in actors:
                act.SetVisibility(False)
                act.SetPickable(False)
            self.channel_actors[ch] = {"body": actors, "labels": labels}
        self._legend()

        # kursor: pierscien w plaszczyznie XY, obracany macierza na normalna skory
        ring = pv.Disc(inner=RING_INNER, outer=RING_OUTER, normal=(0, 0, 1), c_res=48)
        self.ring_actor = self.pl.add_mesh(ring, color=THEME["ring"], lighting=False,
                                           name="ring")
        self.ring_actor.SetVisibility(False)
        self.ring_actor.SetPickable(False)

        # Trafianie tylko w skore, z lokatorem komorek. Lokator MUSI byc zbudowany
        # na siatce, ktora trzyma mapper aktora: add_mesh robi jej kopie, a picker
        # uzywa lokatora wylacznie dla tego samego obiektu. Zbudowany na naszej
        # siatce byl po cichu pomijany i kazde trafienie przegladalo 90 tys.
        # trojkatow - 11,5 ms zamiast 0,09 ms, a kursor nie nadazal za mysza.
        self.picker = vtk.vtkCellPicker()
        self.picker.SetTolerance(0.0005)
        self.picker.PickFromListOn()
        self.picker.AddPickList(self.skin_actor)
        self._locator = vtk.vtkCellLocator()
        self._locator.SetDataSet(self.skin_actor.GetMapper().GetInput())
        self._locator.BuildLocator()
        self.picker.AddLocator(self._locator)

        self.pl.add_text(
            "POWIERZCHNIA\n"
            "[lewy] obrot   [Shift+lewy] przesuw   [kolko] zblizenie do kursora\n"
            "[f] obracaj wokol miejsca pod kursorem   [r] widok calej postaci\n"
            "legenda po lewej: kanaly wg Wu Xing   [p] wszystkie / zadne   [l] kody\n"
            "[1] przod  [2] tyl  [3] bok prawy  [4] bok lewy  [5] gora",
            position="upper_right", font_size=8, color=THEME["dim"], name="help")
        self.info("")
        self._keys()
        self._view("xz", 180)

    def _keys(self) -> None:
        keys = {
            "p": self.toggle_points,
            "l": self.toggle_labels,
            "f": self.focus_under_cursor,
            "r": lambda: self._view("xz", 180),
            "1": lambda: self._view("xz", 180),
            "2": lambda: self._view("xz", 0),
            "3": lambda: self._view("yz", 0),
            "4": lambda: self._view("yz", 180),
            "5": lambda: self._view("xy", 0),
        }
        for k, fn in keys.items():
            self.pl.add_key_event(k, fn)

    def _interaction(self) -> None:
        """Wlasny styl kamery: kolko przybliza do punktu pod kursorem."""
        iren = self.pl.iren.interactor if hasattr(self.pl.iren, "interactor") else None
        if iren is None:
            return
        style = vtk.vtkInteractorStyleTrackballCamera()
        # obserwator na stylu ZASTEPUJE domyslna obsluge zdarzenia w VTK
        style.AddObserver("MouseWheelForwardEvent", lambda *_: self._wheel(ZOOM_STEP))
        style.AddObserver("MouseWheelBackwardEvent", lambda *_: self._wheel(1.0 / ZOOM_STEP))
        style.AddObserver("MouseMoveEvent", self._on_move)
        # Wbudowane skroty stylu VTK (f - lot kamery, r - reset, p - wybor,
        # w/s - siatka/powierzchnia, 3 - stereo, e/q - wyjscie) nakladalyby sie
        # na nasze klawisze. Pusty obserwator wylacza je; nasze skroty dzialaja
        # dalej, bo sa podpiete pod KeyPressEvent interaktora, nie pod styl.
        style.AddObserver("CharEvent", lambda *_: None)
        iren.SetInteractorStyle(style)
        self._style = style

    # --------------------------------------------------------------- pomocnicze
    def info(self, text: str) -> None:
        self.pl.add_text(text, position="lower_left", font_size=10,
                         color=THEME["text"], name="info")

    def _view(self, cam: str, azimuth: float) -> None:
        self.pl.camera_position = cam
        if azimuth:
            self.pl.camera.azimuth = azimuth
        self.pl.reset_camera()
        # 1.05, a nie 1.3: przy 1.3 widok z przodu obcinal glowe i stopy
        self.pl.camera.zoom(1.05)
        self.pl.render()

    def pick(self, x: int, y: int):
        """Punkt i normalna skory pod pikselem, albo (None, None)."""
        if not self.picker.Pick(x, y, 0, self.pl.renderer):
            return None, None
        p = np.array(self.picker.GetPickPosition())
        n = np.array(self.picker.GetPickNormal())
        # normalna ma wskazywac na zewnatrz, ku patrzacemu
        to_cam = np.array(self.pl.camera.position) - p
        if np.dot(n, to_cam) < 0:
            n = -n
        return p, n

    def _mouse(self):
        iren = self.pl.iren.interactor
        return iren.GetEventPosition()

    # ------------------------------------------------------------------ kursor
    def _on_move(self, obj, _event) -> None:
        obj.OnMouseMove()                    # obrot i przesuw dalej dzialaja
        now = time.perf_counter()
        if now - self._last_move < 0.012:    # ok. 80 razy na sekunde wystarczy
            return
        self._last_move = now
        x, y = self._mouse()
        p, n = self.pick(x, y)
        self.place_cursor(p, n)

    def place_cursor(self, p, n) -> None:
        """Kladzie pierscien na skorze. Tylko macierz aktora - zadnej geometrii."""
        if p is None:
            if self.ring_actor.GetVisibility():
                self.ring_actor.SetVisibility(False)
                self._set_code(None)
                self.pl.render()
            return
        R = _align_z_to(n)
        m = vtk.vtkMatrix4x4()
        pos = p + n * 0.0012                 # tuz nad skora, zeby nie migotal
        for i in range(3):
            for j in range(3):
                m.SetElement(i, j, float(R[i, j]))
            m.SetElement(i, 3, float(pos[i]))
        self.ring_actor.SetUserMatrix(m)
        self.ring_actor.SetVisibility(True)
        code = None
        if self.channels_on:
            # tylko punkty kanalow zaznaczonych w legendzie - ukryty punkt
            # nie moze podpowiadac kodu
            visible = np.isin(self._mchan, list(self.channels_on))
            d = np.linalg.norm(self._mpos - p, axis=1)
            d[~visible] = np.inf
            j = int(np.argmin(d))
            if d[j] <= NEAR_POINT:
                m = self.atlas.markers[j]
                el = self.atlas.channels[m["channel"]].get("element", "-")
                code = (f"{m['key']}  ({d[j] * 1000:.0f} mm od kursora)   "
                        f"{ELEMENTS.get(el, ELEMENTS['-'])[0].lower()}")
        self._set_code(code)
        self.pl.render()

    def _set_code(self, code) -> None:
        if code == self._last_code:
            return
        self._last_code = code
        self.info(code or "")

    # ------------------------------------------------------------------ kamera
    def _wheel(self, factor: float) -> None:
        x, y = self._mouse()
        p, _ = self.pick(x, y)
        self.zoom_toward(p, factor)

    def zoom_toward(self, target, factor: float) -> None:
        """Przybliza albo oddala, trzymajac miejsce pod kursorem w tym samym punkcie ekranu.

        Pozycja kamery i punkt ogniskowy przesuwaja sie o ten sam ulamek ku celowi,
        wiec kierunek patrzenia sie nie zmienia - obraz "rosnie" wokol kursora.
        """
        cam = self.pl.camera
        pos = np.array(cam.position)
        foc = np.array(cam.focal_point)
        if target is None:                   # kursor poza cialem: zwykle zblizenie
            cam.dolly(1.0 / factor)
        else:
            t = np.asarray(target, dtype=float)
            new_pos = t + (pos - t) * factor
            if np.linalg.norm(new_pos - t) < MIN_DISTANCE:
                return                       # nie wchodzimy pod skore
            cam.position = tuple(new_pos)
            cam.focal_point = tuple(t + (foc - t) * factor)
        self.pl.renderer.ResetCameraClippingRange()
        self.pl.render()

    def focus_under_cursor(self) -> None:
        """[f] - miejsce pod kursorem staje sie osia obrotu. Przejscie jest plynne."""
        x, y = self._mouse()
        p, _ = self.pick(x, y)
        if p is None:
            self.info("Kursor nie jest na ciele - najedz na miejsce, wokol ktorego obracac.")
            self.pl.render()
            return
        self.set_pivot(p, frames=10)

    def set_pivot(self, p, frames: int = 10) -> None:
        cam = self.pl.camera
        pos0 = np.array(cam.position)
        foc0 = np.array(cam.focal_point)
        shift = np.asarray(p, dtype=float) - foc0
        for k in range(1, frames + 1):
            s = k / frames
            s = s * s * (3 - 2 * s)          # lagodny start i koniec
            cam.focal_point = tuple(foc0 + shift * s)
            cam.position = tuple(pos0 + shift * s)
            self.pl.renderer.ResetCameraClippingRange()
            self.pl.render()
        self.info("Obrot wokol wskazanego miejsca.  [r] powrot do calej postaci")
        self.pl.render()

    # --------------------------------------------------------- kanaly i legenda
    def _palette(self) -> dict:
        """Kolor kazdego kanalu z jego przemiany (pole `element` w danych) i yin/yang."""
        out = {}
        for ch, info in self.atlas.channels.items():
            _, yin, yang = ELEMENTS.get(info.get("element", "-"), ELEMENTS["-"])
            out[ch] = yin if ch in YIN else yang
        return out

    def _legend(self) -> None:
        """Legenda-panel: przemiana zaznacza cala grupe, kanal - siebie."""
        self._channel_widgets: dict = {}
        self._element_widgets: dict = {}
        groups = {e: [] for e in ELEMENT_ORDER}
        for ch, info in self.atlas.channels.items():
            groups.get(info.get("element", "-"), groups["-"]).append(ch)
        x, y = 14, self.pl.window_size[1] - 44
        self.pl.add_text("KANALY  wg Wu Xing", position=(x, y + 8), font_size=9,
                         color=THEME["text"], name="leg_title")
        y -= 26
        for el in ELEMENT_ORDER:
            chans = groups[el]
            if not chans:
                continue
            name, yin, _ = ELEMENTS[el]
            self._element_widgets[el] = self.pl.add_checkbox_button_widget(
                self._element_cb(el), value=False, position=(x, y), size=18,
                border_size=1, color_on=yin, color_off=THEME["off"])
            self.pl.add_text(name, position=(x + 26, y + 1), font_size=9,
                             color=THEME["text"], name=f"leg_el_{el}")
            y -= 24
            for ch in chans:
                info = self.atlas.channels[ch]
                self._channel_widgets[ch] = self.pl.add_checkbox_button_widget(
                    self._channel_cb(ch), value=False, position=(x + 22, y), size=15,
                    border_size=1, color_on=self.channel_color[ch], color_off=THEME["off"])
                kind = "yin" if ch in YIN else "yang"
                self.pl.add_text(f"{ch}  {info['pl']}  ({kind})", position=(x + 44, y),
                                 font_size=8, color=THEME["dim"], name=f"leg_ch_{ch}")
                y -= 21
            y -= 6
        self._element_groups = groups

    def _channel_cb(self, ch: str):
        def cb(state: bool) -> None:
            self.set_channel(ch, bool(state))
        return cb

    def _element_cb(self, el: str):
        def cb(state: bool) -> None:
            for ch in self._element_groups[el]:
                self.set_channel(ch, bool(state), render=False)
            self.pl.render()
        return cb

    def set_channel(self, ch: str, on: bool, render: bool = True) -> None:
        (self.channels_on.add if on else self.channels_on.discard)(ch)
        acts = self.channel_actors.get(ch, {})
        for a in acts.get("body", []):
            a.SetVisibility(on)
        if acts.get("labels") is not None:
            acts["labels"].SetVisibility(on and self.show_labels)
        self._sync_legend()
        if on:
            info = self.atlas.channels[ch]
            name = ELEMENTS.get(info.get("element", "-"), ELEMENTS["-"])[0].lower()
            self.info(f"{ch} {info['pl']} - {name}, "
                      f"{'yin' if ch in YIN else 'yang'}")
        if render:
            self.pl.render()

    def _sync_legend(self) -> None:
        """Stan pol wyboru ma odpowiadac temu, co widac - takze po [p] i po grupie."""
        for ch, w in self._channel_widgets.items():
            w.GetRepresentation().SetState(ch in self.channels_on)
        for el, w in self._element_widgets.items():
            chans = self._element_groups[el]
            w.GetRepresentation().SetState(bool(chans) and all(c in self.channels_on
                                                               for c in chans))

    def toggle_points(self) -> None:
        """[p] - wszystkie kanaly naraz albo zadne."""
        on = len(self.channels_on) < len(self.channel_actors)
        for ch in self.channel_actors:
            self.set_channel(ch, on, render=False)
        self.info("Wszystkie kanaly" if on else "Czysta postac - kanaly ukryte")
        self.pl.render()

    def toggle_labels(self) -> None:
        self.show_labels = not self.show_labels
        for ch, acts in self.channel_actors.items():
            if acts.get("labels") is not None:
                acts["labels"].SetVisibility(self.show_labels and ch in self.channels_on)
        if self.show_labels and not self.channels_on:
            self.info("Kody pokaza sie po zaznaczeniu kanalu w legendzie.")
        self.pl.render()

    # ------------------------------------------------------------------- start
    def show(self) -> None:
        # WSLg potrafi zaparkowac nowe okno poza ekranem - jawna pozycja temu zapobiega
        try:
            rw = self.pl.ren_win
            rw.SetPosition(60, 60)
            rw.SetSize(*self.pl.window_size)
        except Exception:
            pass
        self.pl.show(auto_close=False, interactive_update=True)
        self._interaction()
        self.pl.iren.interactor.Start()
        self.pl.close()

    def screenshot(self, path: str) -> None:
        self.pl.screenshot(path)
