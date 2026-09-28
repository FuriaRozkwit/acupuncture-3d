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
import textwrap
import threading
import time
from pathlib import Path

import numpy as np
import pyvista as pv
import vtk

from .model import Atlas
from .translate_pl import (
    get_point_pl,
    polish_point_name,
    translate_action,
    translate_indications,
    translate_location,
    translate_needling,
    translate_palpation,
)
from .clinical_kb import (
    get_extraordinary_vessel_info,
    get_xi_cleft_info,
    get_yuan_luo_info,
    get_shu_mu_info,
    YUAN_LUO_PAIRS,
    SHU_MU_PAIRS,
    XI_CLEFT_POINTS,
)
from .extraordinary import (
    EXTRAORDINARY_VESSELS,
    EXTRAORDINARY_COLORS,
    get_vessel,
    get_point_extraordinary_memberships,
)
from .jing_jin import (
    JING_JIN_CHANNELS,
    MYERS_LINES,
    get_jing_jin,
    get_myers_line,
    analyze_stroke_fascia,
)
from .wire import Wire
from .klisza import Klisza, KliszaManager

THEME = {
    "skin": "#ece6dc",        # porcelanowa biel z cieplym odcieniem
    "bg_top": "#2e3440",      # elegancki kontrast studyjny
    "bg_bottom": "#121418",
    "cun_grid": "#be9a84",    # lagodny, dyskretny cieply nude/biszkopt (delikatna siatka cunowa)
    "ring": "#ffe36a",        # zloty akcent wyboru
    "text": "#e6e9ee",
    "dim": "#8e98a8",
    "off": "#282e38",
}

# Wu Xing: przemiana -> (nazwa, odcien yin, odcien yang).
# Wyrazista paleta o wysokiej chromie i nasyceniu, doskonale kontrastująca z porcelanową bielą.
ELEMENTS = {
    # Drewno: Żywy szmaragd (Yin - LR) i jasny promienny nefryt (Yang - GB)
    "Wood":  ("DREWNO", "#00a854", "#20bf6b"),
    # Ogień: Głęboki karmin (Yin - HT, PC) i płomienny koral (Yang - SI, SJ)
    "Fire":  ("OGIEN",  "#e02424", "#ff5252"),
    # Ziemia: Ciepły bursztyn (Yin - SP) i świetliste złoto (Yang - ST)
    "Earth": ("ZIEMIA", "#d97706", "#f59e0b"),
    # Metal: Elektryczny błękit kobaltowy (Yin - LU) i lodowy turkus/cyan (Yang - LI)
    "Metal": ("METAL",  "#0284c7", "#06b6d4"),
    # Woda: Głęboki szafirowy ultramaryna (Yin - KI) i żywy błękit lazurowy (Yang - BL)
    "Water": ("WODA",   "#1d4ed8", "#3b82f6"),
    # Poza Wu Xing: Królewski fiolet ametystowy (Yin - REN) i świetlista orchidea (Yang - DU)
    "-":     ("POZA WU XING", "#7c3aed", "#a855f7"),
}
ELEMENT_ORDER = ("Wood", "Fire", "Earth", "Metal", "Water")

# Kanoniczne rejony anatomiczne do selekcji i kadrowania punktów
REGIONS = {
    "ALL": {
        "label": "CALOSC",
        "full_name": "Całe ciało (Wszystkie punkty)",
        "focal": (0.0, 0.0, 0.90),
        "pos": (0.0, 2.45, 0.90),
    },
    "HEAD": {
        "label": "GLOWA",
        "full_name": "Głowa i szyja",
        "focal": (0.0, 0.0, 1.56),
        "pos": (0.0, 0.70, 1.56),
    },
    "CHEST": {
        "label": "KLATKA",
        "full_name": "Klatka piersiowa i brzuch",
        "focal": (0.0, 0.0, 1.12),
        "pos": (0.0, 0.90, 1.12),
    },
    "BACK": {
        "label": "PLECY",
        "full_name": "Plecy i lędźwie",
        "focal": (0.0, 0.0, 1.16),
        "pos": (0.0, -0.95, 1.16),
    },
    "ARMS": {
        "label": "RECE",
        "full_name": "Kończyny górne (Ręce)",
        "focal": (0.0, 0.0, 1.05),
        "pos": (0.0, 1.10, 1.05),
    },
    "LEGS": {
        "label": "NOGI",
        "full_name": "Kończyny dolne (Nogi)",
        "focal": (0.0, 0.0, 0.38),
        "pos": (0.0, 1.25, 0.38),
    },
}
REGION_ORDER = ("ALL", "HEAD", "CHEST", "BACK", "ARMS", "LEGS")
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
POINT_RADIUS = 0.0038         # m - wyrazniejszy znacznik punktu na posazku
NEAR_POINT = 0.022            # m - od tej odleglosci podajemy kod punktu pod kursorem
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
    def __init__(self, atlas: Atlas, off_screen: bool = False, window_size=(1500, 960),
                 channels: list[str] | set[str] | None = None):
        self.atlas = atlas
        self.wire = Wire()
        self.show_wire = False
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
        self.channels_on: set = set(channels) if channels is not None else set()
        self._mpos = np.array([m["pos"] for m in atlas.markers], dtype=float)
        self._mchan = np.array([m["channel"] for m in atlas.markers])
        self.channel_color = self._palette()
        self.active_region: str = "ALL"
        self.show_meridian_lines: bool = False
        self.show_points_mode: bool = False
        self.point_actors: dict = {}  # (ch, reg) -> actor

        # Szybka wektorowa klasyfikacja rejonów anatomicznych markerów
        from scipy.spatial import cKDTree
        _tree = cKDTree(self.wire.points[self.wire._valid])
        _, _n_valid = _tree.query(self._mpos)
        _wire_indices = self.wire._valid[_n_valid]

        self._mreg = []
        for i, m in enumerate(atlas.markers):
            w_idx = _wire_indices[i]
            seg = self.wire.segments[self.wire.seg_id[w_idx]]
            pos = m["pos"]
            if "head" in seg or "neck" in seg:
                reg = "HEAD"
            elif "arm" in seg or "hand" in seg or "forearm" in seg:
                reg = "ARMS"
            elif "thigh" in seg or "shank" in seg or "foot" in seg:
                reg = "LEGS"
            elif "trunk" in seg:
                reg = "CHEST" if pos[1] >= -0.01 else "BACK"
            else:
                reg = "HEAD"
            m["region"] = reg
            self._mreg.append(reg)
        self._mreg = np.array(self._mreg)

        self._last_move = 0.0
        self._last_pick_xy: tuple[int, int] | None = None
        self._last_code: str | None = None
        self.selected_point: str | None = None
        self._press_pos: tuple[int, int] | None = None
        self._press_time: float = 0.0

        # Lokator powierzchni skory do docinania rurek meridianow przez powietrze (lazy)
        self._air_loc = vtk.vtkStaticCellLocator()
        self._air_loc.SetDataSet(self.skin)
        self._air_loc.BuildLocator()
        self._c = [0.0, 0.0, 0.0]
        self._cid = vtk.reference(0)
        self._sub = vtk.reference(0)
        self._d2 = vtk.reference(0.0)
        # Narzedzia rysowania i oznaczania na glinianym posazku
        self.draw_mode: bool = False
        self.camera_frozen: bool = False
        self.draw_strokes: list = []
        self._current_stroke: list = []
        self.drawing_actors: list = []
        self._temp_stroke_actor = None

        # Baza i silnik analizy klisz diagnostycznych
        self.klisza_mgr = KliszaManager(self.atlas, wire=self.wire)
        self.active_klisza: Klisza | None = None
        self._browse_idx: int = -1
        self._browse_actors: list = []
        self._dock_buttons: list[dict] = []
        self._dock_actors: dict[str, Any] = {}
        self._dock_actor_map: dict[Any, dict] = {}
        self._dock_pressed_btn: dict | None = None
        self._dock_picker = vtk.vtkPropPicker()

        # 8 Naczyn Nadzwyczajnych (Qi Jing Ba Mai)
        self.active_vessel: str | None = None
        self.vessel_actors: dict[str, list] = {}
        # Jing Jin & Tasmy Myersa
        self.active_myers: str | None = None
        self.myers_actors: dict[str, list] = {}
        # Biao-Li (Yin-Yang) poziomy: 1 = Biao (Skora), 2 = Jing Jin (Powiez), 3 = Li (Wnetrze)
        self.biao_li_level: int = 1
        self.skeleton_actor = None
        self.muscles_actor = None
        self.console_listener_started: bool = False

        try:
            self.pl.ren_win.AddObserver("WindowResizeEvent", self._on_window_resize)
        except Exception:
            pass

        self._build()
        self._interaction()

    def _to_skin(self, p) -> float:
        self._air_loc.FindClosestPoint([float(p[0]), float(p[1]), float(p[2])], self._c, self._cid, self._sub, self._d2)
        return float(self._d2) ** 0.5

    # ------------------------------------------------------------------ scena
    def _build(self) -> None:
        # Skora jako porcelanowa biel o satynowym wykonczeniu
        self.skin_actor = self.pl.add_mesh(
            self.skin, color=THEME["skin"], smooth_shading=True,
            ambient=0.45, diffuse=0.62, specular=0.08,
            name="skin")

        # Wektorowa siatka cunowa - dyskretna, elegancka i czytelna na porcelanie
        self.wire_actor = self.pl.add_mesh(
            self.wire.mesh, color=THEME["cun_grid"], line_width=1.4, opacity=0.85,
            lighting=False,
            name="cun_grid")
        self.wire_actor.SetVisibility(self.show_wire)
        try:
            mapper = self.wire_actor.GetMapper()
            mapper.SetResolveCoincidentTopologyToPolygonOffset()
            mapper.SetRelativeCoincidentTopologyPolygonOffsetParameters(-1.0, -1.0)
        except Exception:
            pass

        # Kazdy kanal: znaczniki punktow z podzialem na rejony anatomiczne i wykonczeniem perlowym
        self.channel_actors: dict = {}
        self.point_actors: dict = {}
        sphere = pv.Sphere(radius=POINT_RADIUS, theta_resolution=12, phi_resolution=12)

        for ch in self.atlas.channels:
            color = self.channel_color[ch]
            mask_ch = (self._mchan == ch)
            is_on = ch in self.channels_on
            labels = None
            if mask_ch.any():
                for reg in ("HEAD", "CHEST", "BACK", "ARMS", "LEGS"):
                    mask_cr = mask_ch & (self._mreg == reg)
                    if mask_cr.any():
                        pts_cr = pv.PolyData(self._mpos[mask_cr])
                        act_cr = self.pl.add_mesh(
                            pts_cr.glyph(geom=sphere, scale=False, orient=False),
                            color=color, smooth_shading=True,
                            ambient=0.35, diffuse=0.75, specular=0.5, specular_power=25,
                            name=f"pts_{ch}_{reg}")
                        act_cr.SetVisibility(False)
                        act_cr.SetPickable(False)
                        self.point_actors[(ch, reg)] = act_cr

                labels = self.pl.add_point_labels(
                    self._mpos[mask_ch], [m["key"] for m in self.atlas.markers
                                           if m["channel"] == ch],
                    font_size=10, text_color=THEME["text"], shape=None,
                    show_points=False, always_visible=False, name=f"lab_{ch}")
                labels.SetVisibility(False)
                labels.SetPickable(False)

            self.channel_actors[ch] = {"meridian": None, "labels": labels}
            if is_on and self.show_meridian_lines:
                self._get_or_create_tubes(ch)

        self._legend()

        # Znacznik podswietlenia wybranego punktu (zloty, powiekszony)
        hl_sphere = pv.Sphere(radius=POINT_RADIUS * 1.6, theta_resolution=16, phi_resolution=16)
        self.highlight_actor = self.pl.add_mesh(
            hl_sphere, color="#ffe36a", emissive=True, lighting=False,
            name="selected_point_highlight")
        self.highlight_actor.SetVisibility(False)
        self.highlight_actor.SetPickable(False)

        # Znacznik punktu partnerskiego w parach TCM (turkusowy / cyan)
        hl_pair_sphere = pv.Sphere(radius=POINT_RADIUS * 1.5, theta_resolution=16, phi_resolution=16)
        self.highlight_actor_pair = self.pl.add_mesh(
            hl_pair_sphere, color="#38bdf8", emissive=True, lighting=False,
            name="selected_point_pair_highlight")
        self.highlight_actor_pair.SetVisibility(False)
        self.highlight_actor_pair.SetPickable(False)
        self.relation_tube_actor = None

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

        self._setup_right_dock()
        self._update_card(None)
        self.info("Gotowy. Wybierz punkt lub użyj panelu po prawej stronie.")
        self._keys()
        self._view("xz", 180)

    def _get_or_create_tubes(self, ch: str):
        """Tworzy rurki przebiegu kanalu na zadanie (lazy loading)."""
        entry = self.channel_actors.get(ch)
        if not entry:
            return None
        if entry.get("meridian") is not None:
            return entry["meridian"]
        info = self.atlas.channels.get(ch, {})
        color = self.channel_color[ch]
        tubes = []
        for side in ("R", "L"):
            for poly in self.atlas.meridian_lines(ch, side):
                if len(poly) <= 3:
                    continue
                pieces = split_through_air(poly, self._to_skin)
                for piece in pieces:
                    tubes.append(pv.lines_from_points(piece).tube(
                        radius=LINE_RADIUS, n_sides=8))
            if info.get("clock", 0) >= 13:        # REN i DU sa pojedyncze
                break
        if tubes:
            merged = tubes[0].merge(tubes[1:]) if len(tubes) > 1 else tubes[0]
            act = self.pl.add_mesh(merged, color=color, smooth_shading=True,
                                   name=f"mer_{ch}")
            act.SetPickable(False)
            entry["meridian"] = act
            return act
        return None

    def _keys(self) -> None:
        keys = {
            "p": self.toggle_points,
            "l": self.toggle_labels,
            "g": self.toggle_cun_grid,
            "m": self.toggle_all_main_meridians,
            "b": self.cycle_region,
            "e": self.cycle_extraordinary,
            "j": self.cycle_jing_jin,
            "y": self.cycle_biao_li,
            "u": self.cycle_yuan_luo,
            "o": self.cycle_shu_mu,
            "i": self.cycle_xi_cleft,
            "w": self.cycle_wu_xing,
            "d": self.toggle_draw_mode,
            "x": self.toggle_camera_freeze,
            "z": self.undo_draw_stroke,
            "c": self.clear_draw_strokes,
            "s": self.save_current_klisza,
            "k": self.browse_klisze,
            "a": self.analyze_current_strokes,
            "Escape": self.clear_selection,
            "f": self.toggle_camera_freeze,
            "r": self._view_full,
            "question": self.print_interactive_help,
            "slash": self.print_interactive_help,
            "1": self._view_front,
            "2": self._view_back,
            "3": self._view_right,
            "4": self._view_left,
            "5": lambda: self.set_region("ALL"),
            "6": lambda: self.set_region("HEAD"),
            "7": lambda: self.set_region("CHEST"),
            "8": lambda: self.set_region("BACK"),
            "9": lambda: self.set_region("ARMS"),
            "0": lambda: self.set_region("LEGS"),
        }
        for k, fn in keys.items():
            self.pl.add_key_event(k, fn)

    def _hex_to_rgb(self, h: str) -> tuple[float, float, float]:
        h = h.lstrip("#")
        return tuple(int(h[i:i+2], 16) / 255.0 for i in (0, 2, 4))

    def _is_over_dock(self, x: int, y: int) -> bool:
        """Sprawdza czy kursor znajduje sie w strefie prawego panelu przyciskow."""
        if hasattr(self, "_dock_picker"):
            hit = self._dock_picker.PickProp(x, y, self.pl.renderer)
            if hit and self._dock_picker.GetViewProp() in self._dock_actor_map:
                return True
        try:
            W, H = self.pl.ren_win.GetSize()
        except Exception:
            W, H = self.pl.window_size
        return (x >= W - 280 and y >= H - 720)

    def _on_window_resize(self, obj, _event) -> None:
        """Dynamicznie dopasowuje pozycje prawego panelu po zmianie rozmiaru okna."""
        try:
            W, H = obj.GetSize()
            self.pl.window_size = (W, H)
            bx = W - 262
            bw = 246
            bh = 26

            title = self.pl.actors.get("dock_title")
            if title:
                title.SetPosition(bx, H - 42)
            s1 = self.pl.actors.get("dock_sec_1")
            if s1:
                s1.SetPosition(bx, H - 68)
            s2 = self.pl.actors.get("dock_sec_2")
            if s2:
                s2.SetPosition(bx, H - 276)
            s3 = self.pl.actors.get("dock_sec_3")
            if s3:
                s3.SetPosition(bx, H - 458)
            s4 = self.pl.actors.get("dock_sec_4")
            if s4:
                s4.SetPosition(bx, H - 588)
            badge = self.pl.actors.get("draw_badge")
            if badge:
                badge.SetPosition(240, H - 34)
            fbadge = self.pl.actors.get("freeze_badge")
            if fbadge:
                fbadge.SetPosition(240, H - 34 if not self.draw_mode else H - 60)

            for btn in self._dock_buttons:
                y = H - btn.get("y_offset", 100)
                btn["y"] = y
                btn["rect"] = (bx, y - 4, bw, bh + 4)
                btn["actor"].SetPosition(bx, y)
        except Exception:
            pass

    def _view_front(self) -> None:
        _safe_print(">>> [KAMERA] Widok z przodu [1]")
        self._view("xz", 180)

    def _view_back(self) -> None:
        _safe_print(">>> [KAMERA] Widok z tylu [2]")
        self._view("xz", 0)

    def _view_right(self) -> None:
        _safe_print(">>> [KAMERA] Widok prawego boku [3]")
        self._view("yz", 0)

    def _view_left(self) -> None:
        _safe_print(">>> [KAMERA] Widok lewego boku [4]")
        self._view("yz", 180)

    def _view_full(self) -> None:
        _safe_print(">>> [KAMERA] Widok calej postaci [r]")
        self._view("xz", 180)

    def _setup_right_dock(self) -> None:
        """Buduje czytelny, zorganizowany panel przyciskow klinicznych po prawej stronie."""
        try:
            W, H = self.pl.ren_win.GetSize()
        except Exception:
            W, H = self.pl.window_size
        bx = W - 262
        bw = 246
        bh = 26

        self._dock_buttons.clear()
        self._dock_actors.clear()
        self._dock_actor_map.clear()

        # Naglowek panelu
        self.pl.add_text("PANEL KLINICZNY", position=(bx, H - 42), font_size=11,
                         color="#ffe36a", name="dock_title")

        def add_btn(key: str, label: str, sc: str, y_offset: int, action,
                    color="#1c2128", border="#485463", text_color="#f0f4f8"):
            y = H - y_offset
            gap = 26 - len(label) - len(sc)
            spaces = " " * max(gap, 2)
            txt = f" {label}{spaces}[{sc}] "
            act_name = f"dock_btn_{key}"
            actor = self.pl.add_text(txt, position=(bx, y), font_size=9,
                                    color=text_color, name=act_name)
            actor.SetPickable(True)
            tp = actor.prop
            r, g, b = self._hex_to_rgb(color)
            br, bg, bb = self._hex_to_rgb(border)
            tp.SetBackgroundColor(r, g, b)
            tp.SetBackgroundOpacity(0.92)
            tp.SetFrame(1)
            tp.SetFrameColor(br, bg, bb)
            tp.SetFrameWidth(1)

            btn_desc = {
                "key": key,
                "label": label,
                "shortcut": sc,
                "action": action,
                "y": y,
                "y_offset": y_offset,
                "rect": (bx, y - 4, bw, bh + 4),
                "actor": actor,
                "color": color,
                "border": border,
                "text_color": text_color,
                "active": False,
            }
            self._dock_buttons.append(btn_desc)
            self._dock_actors[key] = btn_desc
            self._dock_actor_map[actor] = btn_desc
            return btn_desc

        # ─── 1. SZKICOWNIK & KLISZE ───
        self.pl.add_text("── SZKICOWNIK & KLISZE ──", position=(bx, H - 68), font_size=8,
                         color="#8e98a8", name="dock_sec_1")
        add_btn("d", "Rysuj napiecia", "d", 94, self.toggle_draw_mode,
                color="#2c1717", border="#e63946", text_color="#ff7b72")
        add_btn("x", "Zamroz kamere",  "x", 120, self.toggle_camera_freeze,
                color="#16202c", border="#38bdf8", text_color="#7dd3fc")
        add_btn("z", "Cofnij linie",   "z", 146, self.undo_draw_stroke)
        add_btn("c", "Wyczysc szkic",  "c", 172, self.clear_draw_strokes)
        add_btn("s", "Zapisz klisze",  "s", 198, self.save_current_klisza,
                color="#15261e", border="#3fb950", text_color="#56d364")
        add_btn("k", "Przegladaj klisze", "k", 224, self.browse_klisze)
        add_btn("a", "Analizuj polaczenia", "a", 250, self.analyze_current_strokes,
                color="#2c2214", border="#d29922", text_color="#e3b341")

        # ─── 2. STRUKTURY MERIDIANOWE ───
        self.pl.add_text("── STRUKTURY MERIDIANOWE ──", position=(bx, H - 276), font_size=8,
                         color="#8e98a8", name="dock_sec_2")
        add_btn("g", "Siatka cun: OFF", "g", 302, self.toggle_cun_grid,
                color="#21262d", border="#485463", text_color="#8e98a8")
        add_btn("m", "12 Glownych: OFF", "m", 328, self.toggle_all_main_meridians,
                color="#16202c", border="#485463", text_color="#8e98a8")
        add_btn("p", "Punkty: OFF", "p", 354, self.toggle_points,
                color="#21262d", border="#485463", text_color="#8e98a8")
        add_btn("e", "8 Naczyn: OFF",  "e", 380, self.cycle_extraordinary,
                color="#261b2e", border="#485463", text_color="#8e98a8")
        add_btn("j", "Jing Jin: OFF",  "j", 406, self.cycle_jing_jin,
                color="#16222b", border="#485463", text_color="#8e98a8")
        add_btn("y", "Biao-Li: POWIERZCHNIA", "y", 432, self.cycle_biao_li,
                color="#232219", border="#eab308", text_color="#fef08a")

        # ─── 3. ZASADY & PARY TCM ───
        self.pl.add_text("── ZASADY & PARY TCM ──", position=(bx, H - 458), font_size=8,
                         color="#8e98a8", name="dock_sec_3")
        add_btn("u", "Yuan-Luo: CYKL", "u", 484, self.cycle_yuan_luo,
                color="#2a2415", border="#d29922", text_color="#fef08a")
        add_btn("o", "Shu-Mu: CYKL",   "o", 510, self.cycle_shu_mu,
                color="#15261e", border="#3fb950", text_color="#7ee787")
        add_btn("i", "Xi-Cleft: CYKL", "i", 536, self.cycle_xi_cleft,
                color="#2c1717", border="#f85149", text_color="#ff7b72")
        add_btn("w", "Wu Xing: CYKL",  "w", 562, self.cycle_wu_xing,
                color="#1c2538", border="#38bdf8", text_color="#7dd3fc")

        # ─── 4. RZUTY KAMERY & POMOC ───
        self.pl.add_text("── RZUTY KAMERY & POMOC ──", position=(bx, H - 588), font_size=8,
                         color="#8e98a8", name="dock_sec_4")
        add_btn("1", "Widok z przodu", "1", 614, self._view_front)
        add_btn("2", "Widok z tylu",   "2", 640, self._view_back)
        add_btn("r", "Cala postac",    "r", 666, self._view_full)
        add_btn("?", "Konsola & Lore TCM", "?", 692, self.print_interactive_help,
                color="#1a202c", border="#4a5568", text_color="#e2e8f0")

    def _check_dock_click(self, x: int, y: int) -> bool:
        """Sprawdza czy klikniecie trafilo w przycisk prawego panelu."""
        if hasattr(self, "_dock_picker"):
            hit = self._dock_picker.PickProp(x, y, self.pl.renderer)
            prop = self._dock_picker.GetViewProp()
            if hit and prop in self._dock_actor_map:
                btn = self._dock_actor_map[prop]
                _safe_print(f"\n>>> [KLIK PRZYCISKU] Wybrano: {btn['label']} [{btn['shortcut']}]")
                btn["action"]()
                self.pl.render()
                return True

        for btn in self._dock_buttons:
            bx, by, bw, bh = btn["rect"]
            if bx <= x <= bx + bw and by <= y <= by + bh:
                _safe_print(f"\n>>> [KLIK PRZYCISKU] Wybrano: {btn['label']} [{btn['shortcut']}]")
                btn["action"]()
                self.pl.render()
                return True
        return False

    def _update_dock_button(self, key: str, label: str | None = None,
                           color: str | None = None, border: str | None = None,
                           text_color: str | None = None) -> None:
        btn = self._dock_actors.get(key)
        if not btn:
            return
        if label is not None:
            btn["label"] = label
            sc = btn["shortcut"]
            gap = 26 - len(label) - len(sc)
            spaces = " " * max(gap, 2)
            txt = f" {label}{spaces}[{sc}] "
            btn["actor"].SetInput(txt)
        tp = btn["actor"].prop
        if color is not None:
            r, g, b = self._hex_to_rgb(color)
            tp.SetBackgroundColor(r, g, b)
        if border is not None:
            br, bg, bb = self._hex_to_rgb(border)
            tp.SetFrameColor(br, bg, bb)
        if text_color is not None:
            tr, tg, tb = self._hex_to_rgb(text_color)
            tp.SetColor(tr, tg, tb)

    def toggle_draw_mode(self) -> None:
        self.draw_mode = not self.draw_mode
        if self.draw_mode:
            self.camera_frozen = True
            try:
                self.pl.remove_actor("freeze_badge")
            except Exception:
                pass
            try:
                W, H = self.pl.ren_win.GetSize()
            except Exception:
                W, H = self.pl.window_size
            self.pl.add_text(
                "✏ TRYB RYSOWANIA [d] (KAMERA ZAMROŻONA): Przeciągaj myszą po ciele  |  [z] Cofnij  |  [s] Zapisz",
                position=(240, H - 34), font_size=9, color="#ff7b72", name="draw_badge"
            )
            self.info("TRYB RYSOWANIA (KAMERA ZAMROŻONA): Przeciągnij lewy przycisk po porcelanie. [d] zakończ, [s] zapisz.")
            self._update_dock_button("d", label="RYSOWANIE ON", color="#8b2626", border="#ff4d4d", text_color="#ffffff")
            self._update_dock_button("x", label="KAMERA: FREEZE", color="#0f2b38", border="#38bdf8", text_color="#7dd3fc")
            _safe_print("\n>>> [TRYB RYSOWANIA AKTYWNY - KAMERA ZAMROŻONA (FREEZE)]")
            _safe_print("    Przeciągaj lewym przyciskiem myszy po porcelanowym ciele, aby kreślić linie napięć.")
            _safe_print("    Kamera jest w 100% zablokowana - obraz NIE będzie się obracał podczas rysowania!")
            _safe_print("    [z] Cofnij linię  |  [c] Wyczyść szkic  |  [s] Zapisz kliszę  |  [x] Odmroź kamerę  |  [d] Zakończ\n")
        else:
            try:
                self.pl.remove_actor("draw_badge")
            except Exception:
                pass
            self.camera_frozen = False
            try:
                self.pl.remove_actor("freeze_badge")
            except Exception:
                pass
            self.info("Tryb eksploracji punktów (kamera odblokowana). [d] włącz rysowanie, [x] zamroź kamerę.")
            self._update_dock_button("d", label="Rysuj napiecia", color="#2c1717", border="#e63946", text_color="#ff7b72")
            self._update_dock_button("x", label="Zamroz kamere", color="#16202c", border="#38bdf8", text_color="#7dd3fc")
            _safe_print("\n>>> [TRYB EKSPLORACJI AKTYWNY] Kamera odblokowana. Obracaj model kamerą lub klikaj w punkty.")
            _safe_print("    Wciśnij [d] lub kliknij [Rysuj napięcia], aby wrócić do rysowania linii napięć.\n")
        self.pl.render()

    def toggle_camera_freeze(self, state: bool | None = None) -> None:
        """[x] / [f] Przełącza zamrożenie kamery (FREEZE) - uniemożliwia obrót podczas badania i szkicu."""
        if state is not None:
            self.camera_frozen = state
        else:
            self.camera_frozen = not self.camera_frozen

        try:
            W, H = self.pl.ren_win.GetSize()
        except Exception:
            W, H = self.pl.window_size

        if self.camera_frozen:
            if not self.draw_mode:
                self.pl.add_text(
                    "❄ KAMERA ZAMROŻONA (FREEZE) [x]: Obrót wyłączony  |  Możesz stabilnie rysować pasma [d]",
                    position=(240, H - 34),
                    font_size=9, color="#7dd3fc", name="freeze_badge"
                )
            self._update_dock_button("x", label="KAMERA: FREEZE",
                                    color="#0f2b38", border="#38bdf8", text_color="#7dd3fc")
            _safe_print(">>> [KAMERA] FREEZE: Kamera zamrożona [x]. Obrót i przesunięcie modelu zablokowane.")
            self.info("Kamera ZAMROŻONA (FREEZE) [x]. Obrót zablokowany. Możesz precyzyjnie rysować.")
        else:
            try:
                self.pl.remove_actor("freeze_badge")
            except Exception:
                pass
            self._update_dock_button("x", label="Zamroz kamere",
                                    color="#16202c", border="#38bdf8", text_color="#7dd3fc")
            _safe_print(">>> [KAMERA] FREEZE WYŁĄCZONE: Swobodny obrót i nawigacja 3D odblokowane [x].")
            self.info("Kamera ODBLOKOWANA [x]. Swobodny obrót 3D włączony.")

        self.pl.render()

    def save_current_klisza(self) -> None:
        """[s] Zapisuje biezace linie napiec jako klisze diagnostyczna do data/klisze."""
        if not self.draw_strokes:
            msg = "Brak narysowanych linii do zapisu kliszy. Włącz [d] i narysuj pasmo."
            self.info(msg)
            _safe_print(f">>> [ZAPISZ KLISZĘ] {msg}")
            self.pl.render()
            return

        k = self.klisza_mgr.create_klisza(self.draw_strokes)
        saved_path = self.klisza_mgr.save_klisza(k)
        self.active_klisza = k
        self.klisza_mgr.print_klisza_report(k)

        # Wizualna informacja o zapisie
        self._update_dock_button("s", label="ZAPISANO KLISZE!", color="#238636", border="#56d364", text_color="#ffffff")
        pts_count = len(k.all_matched_points)
        msg = f"Zapisano kliszę: {k.name} ({pts_count} pkt, {k.total_length_m*100:.1f} cm). Plik: data/klisze/{saved_path.name}"
        _safe_print(f">>> [ZAPISZ KLISZĘ SUKCES] {msg}")
        self.info(f"Zapisano klisze: {k.name} ({pts_count} pkt, {k.total_length_m*100:.1f} cm). Zapisano w data/klisze/{saved_path.name}")
        self.pl.render()

    def analyze_current_strokes(self) -> None:
        """[a] Analizuje biezace polaczenia napieciowe w modelu i drukuje raport."""
        if not self.draw_strokes:
            msg = "Brak linii napięć do analizy. Włącz [d] i narysuj linię na ciele."
            self.info(msg)
            _safe_print(f">>> [ANALIZA] {msg}")
            self.pl.render()
            return

        k = self.klisza_mgr.create_klisza(self.draw_strokes, name="Bieżąca analiza połączeń")
        self.klisza_mgr.print_klisza_report(k)
        dom = k.strokes[0].get("dominant_channel", "-") if k.strokes else "-"
        pts_count = len(k.all_matched_points)
        self.info(f"Analiza toru: {pts_count} punktow | {k.total_length_m*100:.1f} cm ({k.total_length_cun:.1f} cun) | Kanal: [{dom}] (szczegoly w konsoli)")
        self.pl.render()

    def browse_klisze(self) -> None:
        """[k] Przeglada i cyklicznie laduje zapisane klisze na model 3D."""
        all_meta = self.klisza_mgr.list_all_metadata()
        if not all_meta:
            msg = "Brak zapisanych klisz w bazie. Włącz [d], narysuj linię i naciśnij [s]."
            self.info(msg)
            _safe_print(f">>> [PRZEGLĄDAJ KLISZE] {msg}")
            self.pl.render()
            return

        self._browse_idx = (self._browse_idx + 1) % len(all_meta)
        target_meta = all_meta[self._browse_idx]
        k = self.klisza_mgr.load_klisza(target_meta["id"])
        if not k:
            self.info(f"Nie udało się wczytać kliszy: {target_meta['id']}")
            return

        # Usun poprzednie linie przegladania
        for act in self._browse_actors:
            try:
                self.pl.remove_actor(act)
            except Exception:
                pass
            self._browse_actors.clear()

        # Rysujemy pasma wczytanej kliszy w jasnym kolorze turkusowym (#38bdf8)
        for s_idx, s in enumerate(k.strokes):
            pts = np.array(s.get("points_3d", []))
            if len(pts) >= 2:
                lines = np.hstack([[len(pts)], np.arange(len(pts))])
                poly = pv.PolyData(pts, lines=lines)
                tube = poly.tube(radius=0.0032, n_sides=8)
                act = self.pl.add_mesh(tube, color="#38bdf8", smooth_shading=True,
                                       lighting=False, name=f"browse_stroke_{s_idx}")
                self._browse_actors.append(act)

        self.klisza_mgr.print_klisza_report(k)
        all_pts = k.all_matched_points
        if all_pts:
            self.select_point(all_pts[0]["key"])

        msg = f"Wczytano kliszę [{self._browse_idx+1}/{len(all_meta)}]: {k.name} ({len(all_pts)} pkt) [k: kolejna]"
        _safe_print(f">>> [PRZEGLĄDAJ KLISZE] {msg}")
        self.info(f"Wczytano klisze [{self._browse_idx+1}/{len(all_meta)}]: {k.name} ({len(all_pts)} pkt) [k: kolejna]")
        self.pl.render()

    def undo_draw_stroke(self) -> None:
        if self.drawing_actors:
            act = self.drawing_actors.pop()
            try:
                self.pl.remove_actor(act)
            except Exception:
                pass
            if self.draw_strokes:
                self.draw_strokes.pop()
            msg = f"Cofnięto linię szkicu. Pozostało pasm: {len(self.drawing_actors)}"
            _safe_print(f">>> [COFNIJ] {msg}")
            self.info(f"Cofnieto linie szkicu. Pozostalo: {len(self.drawing_actors)}")
            self.pl.render()
        else:
            msg = "Brak linii do cofnięcia w szkicowniku."
            _safe_print(f">>> [COFNIJ] {msg}")
            self.info("Brak linii do cofniecia.")
            self.pl.render()

    def clear_draw_strokes(self) -> None:
        for act in self.drawing_actors:
            try:
                self.pl.remove_actor(act)
            except Exception:
                pass
        if self._temp_stroke_actor is not None:
            try:
                self.pl.remove_actor(self._temp_stroke_actor)
            except Exception:
                pass
            self._temp_stroke_actor = None
        self.drawing_actors.clear()
        self.draw_strokes.clear()
        _safe_print(">>> [WYCZYŚĆ] Wyczyszczono wszystkie narysowane pasma napięć.")
        self.info("Wyczyszczono szkice napiec.")
        self.pl.render()

    def add_stroke_from_points(self, points: np.ndarray, color: str = "#e63946") -> None:
        """Dodaje linie napiecia bezposrednio ze wspolrzednych 3D (Nx3)."""
        pts = np.asarray(points, dtype=float)
        if len(pts) < 2:
            return
        self.draw_strokes.append(pts)
        lines = np.hstack([[len(pts)], np.arange(len(pts))])
        poly = pv.PolyData(pts, lines=lines)
        tube = poly.tube(radius=0.0028, n_sides=8)
        stroke_idx = len(self.draw_strokes)
        act = self.pl.add_mesh(
            tube, color=color, smooth_shading=True,
            lighting=False, name=f"draw_stroke_{stroke_idx}"
        )
        self.drawing_actors.append(act)
        self.pl.render()

    def _update_current_stroke_actor(self) -> None:
        if len(self._current_stroke) < 2:
            return
        pts = np.array(self._current_stroke)
        lines = np.hstack([[len(pts)], np.arange(len(pts))])
        poly = pv.PolyData(pts, lines=lines)
        tube = poly.tube(radius=0.0028, n_sides=8)
        if self._temp_stroke_actor is not None:
            try:
                self.pl.remove_actor(self._temp_stroke_actor)
            except Exception:
                pass
        self._temp_stroke_actor = self.pl.add_mesh(
            tube, color="#ff3344", smooth_shading=True,
            lighting=False, name="temp_stroke"
        )
        self.pl.render()

    def _interaction(self) -> None:
        """Wlasny styl kamery: blokada obrotu (FREEZE), rysowanie pasm, zoom pod kursorem."""
        if getattr(self, "_style", None) is not None:
            return
        iren = self.pl.iren.interactor if hasattr(self.pl.iren, "interactor") else None
        if iren is None:
            return
        style = vtk.vtkInteractorStyleTrackballCamera()
        style.AddObserver("MouseWheelForwardEvent", lambda *_: self._wheel(ZOOM_STEP))
        style.AddObserver("MouseWheelBackwardEvent", lambda *_: self._wheel(1.0 / ZOOM_STEP))
        style.AddObserver("MouseMoveEvent", self._on_move)
        style.AddObserver("LeftButtonPressEvent", self._on_left_down)
        style.AddObserver("LeftButtonReleaseEvent", self._on_left_up)
        style.AddObserver("RightButtonPressEvent", self._on_right_down)
        style.AddObserver("RightButtonReleaseEvent", self._on_right_up)
        style.AddObserver("MiddleButtonPressEvent", self._on_middle_down)
        style.AddObserver("MiddleButtonReleaseEvent", self._on_middle_up)
        # Wbudowane skroty stylu VTK wylaczamy, by nie kolidowaly z naszymi klawiszami
        style.AddObserver("CharEvent", lambda *_: None)

        self.pl.iren.style = style
        iren.SetInteractorStyle(style)
        self.pl.iren.update_style = lambda: iren.SetInteractorStyle(style)
        try:
            iren.Initialize()
        except Exception:
            pass
        self._style = style

    def _on_right_down(self, obj, _event) -> None:
        if not self.camera_frozen and not self.draw_mode:
            if hasattr(obj, "OnRightButtonDown"):
                obj.OnRightButtonDown()
            elif hasattr(self, "_style") and hasattr(self._style, "OnRightButtonDown"):
                self._style.OnRightButtonDown()

    def _on_right_up(self, obj, _event) -> None:
        if not self.camera_frozen and not self.draw_mode:
            if hasattr(obj, "OnRightButtonUp"):
                obj.OnRightButtonUp()
            elif hasattr(self, "_style") and hasattr(self._style, "OnRightButtonUp"):
                self._style.OnRightButtonUp()

    def _on_middle_down(self, obj, _event) -> None:
        if not self.camera_frozen and not self.draw_mode:
            if hasattr(obj, "OnMiddleButtonDown"):
                obj.OnMiddleButtonDown()
            elif hasattr(self, "_style") and hasattr(self._style, "OnMiddleButtonDown"):
                self._style.OnMiddleButtonDown()

    def _on_middle_up(self, obj, _event) -> None:
        if not self.camera_frozen and not self.draw_mode:
            if hasattr(obj, "OnMiddleButtonUp"):
                obj.OnMiddleButtonUp()
            elif hasattr(self, "_style") and hasattr(self._style, "OnMiddleButtonUp"):
                self._style.OnMiddleButtonUp()

    def _on_left_down(self, obj, _event) -> None:
        self._press_pos = self._mouse()
        self._press_time = time.perf_counter()
        x, y = self._press_pos
        self._dock_pressed_btn = None

        # 1. Sprawdzamy czy kliknięto w przycisk prawego docka przez vtkPropPicker
        if hasattr(self, "_dock_picker"):
            hit = self._dock_picker.PickProp(x, y, self.pl.renderer)
            prop = self._dock_picker.GetViewProp()
            if hit and prop in self._dock_actor_map:
                self._dock_pressed_btn = self._dock_actor_map[prop]
                try:
                    self._dock_pressed_btn["actor"].prop.SetFrameWidth(2)
                    self.pl.render()
                except Exception:
                    pass
                return

        # 2. Sprawdzamy prostokąty przycisków docka
        for btn in self._dock_buttons:
            bx, by, bw, bh = btn["rect"]
            if bx <= x <= bx + bw and by <= y <= by + bh:
                self._dock_pressed_btn = btn
                try:
                    btn["actor"].prop.SetFrameWidth(2)
                    self.pl.render()
                except Exception:
                    pass
                return

        # 3. Jeśli kliknięcie jest w strefie prawego panelu, nie rotujemy kamery
        if self._is_over_dock(x, y):
            return

        # 4. Tryb rysowania pasm napięć na ciele
        if self.draw_mode:
            self._current_stroke = []
            p, n = self.pick(x, y)
            if p is not None and n is not None:
                pt = p + n * 0.003
                self._current_stroke.append(pt)
            return

        # 5. Jeśli kamera jest zamrożona (FREEZE), nie wywołujemy rotacji VTK
        if self.camera_frozen:
            return

        # 6. Normalna rotacja kamery 3D
        if hasattr(obj, "OnLeftButtonDown"):
            obj.OnLeftButtonDown()
        elif hasattr(self, "_style") and hasattr(self._style, "OnLeftButtonDown"):
            self._style.OnLeftButtonDown()

    def _on_left_up(self, obj, _event) -> None:
        if not self.draw_mode and not self.camera_frozen and self._dock_pressed_btn is None:
            if hasattr(obj, "OnLeftButtonUp"):
                obj.OnLeftButtonUp()
            elif hasattr(self, "_style") and hasattr(self._style, "OnLeftButtonUp"):
                self._style.OnLeftButtonUp()
        if self._press_pos is None:
            return

        curr = self._mouse()
        dx = curr[0] - self._press_pos[0]
        dy = curr[1] - self._press_pos[1]
        dt = time.perf_counter() - self._press_time
        self._press_pos = None

        # 1. Obsługa przycisku docka (jeśli wciśnięto przy left_down)
        if self._dock_pressed_btn is not None:
            btn = self._dock_pressed_btn
            self._dock_pressed_btn = None
            try:
                btn["actor"].prop.SetFrameWidth(1)
            except Exception:
                pass
            _safe_print(f"\n>>> [KLIK PRZYCISKU] Wybrano: {btn['label']} [{btn['shortcut']}]")
            btn["action"]()
            self.pl.render()
            return

        # Fallback: hit w przycisk przy puszczeniu myszy
        if hasattr(self, "_dock_picker"):
            hit = self._dock_picker.PickProp(curr[0], curr[1], self.pl.renderer)
            prop = self._dock_picker.GetViewProp()
            if hit and prop in self._dock_actor_map:
                btn = self._dock_actor_map[prop]
                _safe_print(f"\n>>> [KLIK PRZYCISKU] Wybrano: {btn['label']} [{btn['shortcut']}]")
                btn["action"]()
                self.pl.render()
                return

        # Sprawdzenie prostokąta przycisków
        if self._check_dock_click(curr[0], curr[1]):
            return

        # Jeśli puszczenie myszy nastąpiło w strefie docka, nie wykonujemy akcji na ciele
        if self._is_over_dock(curr[0], curr[1]):
            return

        # 2. Tryb rysowania pasm napięć
        if self.draw_mode:
            if len(self._current_stroke) >= 2:
                pts = np.array(self._current_stroke)
                self.draw_strokes.append(pts)
                lines = np.hstack([[len(pts)], np.arange(len(pts))])
                poly = pv.PolyData(pts, lines=lines)
                tube = poly.tube(radius=0.0028, n_sides=8)
                if self._temp_stroke_actor is not None:
                    try:
                        self.pl.remove_actor(self._temp_stroke_actor)
                    except Exception:
                        pass
                    self._temp_stroke_actor = None
                stroke_idx = len(self.draw_strokes)
                act = self.pl.add_mesh(
                    tube, color="#e63946", smooth_shading=True,
                    lighting=False, name=f"draw_stroke_{stroke_idx}"
                )
                self.drawing_actors.append(act)
                msg = f"Narysowano pasmo napięcia #{stroke_idx} ({len(pts)} pkt). [z] Cofnij, [s] Zapisz kliszę."
                _safe_print(f"\n>>> [SZKICOWNIK] {msg}")
                self.info(msg)
            else:
                if self._temp_stroke_actor is not None:
                    try:
                        self.pl.remove_actor(self._temp_stroke_actor)
                    except Exception:
                        pass
                    self._temp_stroke_actor = None
            self._current_stroke = []
            self.pl.render()
            return

        # 3. Kliknięcie w model ciała (tolerancja przesunięcia do 25 px, czasu do 1.5 s)
        if dx * dx + dy * dy <= 625 and dt < 1.5:
            self.on_click(curr[0], curr[1])

    def _project_markers_to_screen(self):
        """Błyskawiczne rzutowanie wektorowe (NumPy) wszystkich punktów na współrzędne ekranu (0.04 ms)."""
        w, h = self.pl.window_size
        cam = self.pl.camera
        mat = cam.GetCompositeProjectionTransformMatrix(w / h, -1, 1)
        M = np.array([[mat.GetElement(r, c) for c in range(4)] for r in range(4)], dtype=np.float64)
        pts4 = np.hstack([self._mpos, np.ones((len(self._mpos), 1), dtype=np.float64)])
        clip = pts4 @ M.T
        w_clip = clip[:, 3]
        in_front = (w_clip > 0.01)
        ndc = clip[:, :3] / np.maximum(w_clip[:, None], 1e-9)
        disp_x = (ndc[:, 0] + 1.0) * 0.5 * w
        disp_y = (ndc[:, 1] + 1.0) * 0.5 * h
        depth = ndc[:, 2]
        return disp_x, disp_y, depth, in_front

    def on_click(self, x: int, y: int) -> None:
        """Obsługa kliknięcia: natychmiastowy wybór punktu i prezentacja karty Deadmana."""
        p, n = self.pick(x, y)
        cam_pos = np.array(self.pl.camera.position)

        disp_x, disp_y, depth, in_front = self._project_markers_to_screen()
        d_px = np.hypot(disp_x - x, disp_y - y)
        d_px[~in_front] = np.inf

        # Wykrywanie zasłonięcia (okluzja przez ciało z przodu)
        cam_to_m = np.linalg.norm(cam_pos - self._mpos, axis=1)
        occluded = np.zeros(len(self._mpos), dtype=bool)
        if p is not None:
            cam_to_p = np.linalg.norm(cam_pos - p)
            occluded = (cam_to_m > cam_to_p + 0.045)
            d_px[occluded] = np.inf

        # Sprawdzenie widoczności punktów w bieżącym trybie
        reg_mask = (self._mreg == self.active_region) if self.active_region != "ALL" else np.ones(len(self._mpos), dtype=bool)
        if self.channels_on:
            chan_mask = np.isin(self._mchan, list(self.channels_on))
        else:
            chan_mask = np.ones(len(self._mpos), dtype=bool) if self.show_points_mode else np.zeros(len(self._mpos), dtype=bool)
        visible = reg_mask & chan_mask

        picked_key = None

        # 1. Priorytet: bezpośrednie kliknięcie w WIDOCZNY punkt na ekranie (promień 36 px)
        vis_d_px = np.copy(d_px)
        vis_d_px[~visible] = np.inf
        if np.any(np.isfinite(vis_d_px)):
            best_vis = int(np.argmin(vis_d_px))
            if vis_d_px[best_vis] <= 36.0:
                picked_key = self.atlas.markers[best_vis]["key"]

        # 2. Jeśli brak widocznego w zasięgu, sprawdź dowolny punkt w pobliżu kursora (promień 26 px)
        if picked_key is None and np.any(np.isfinite(d_px)):
            best_any = int(np.argmin(d_px))
            if d_px[best_any] <= 26.0:
                picked_key = self.atlas.markers[best_any]["key"]

        # 3. Odległość 3D na powierzchni skóry (do 30 mm)
        if picked_key is None and p is not None:
            dists_3d = np.linalg.norm(self._mpos - p, axis=1)
            dists_3d[occluded] = np.inf
            closest_3d = int(np.argmin(dists_3d))
            if dists_3d[closest_3d] <= 0.030:
                picked_key = self.atlas.markers[closest_3d]["key"]


        if picked_key is not None:
            _safe_print(f"\n>>> [KLIK W PUNKT] Trafiono w punkt: {picked_key}")
            self.select_point(picked_key)
        else:
            if p is not None:
                # Klikniecie w puste miejsce na glinie: podaj wspolrzedne cunowe i palpacyjne
                try:
                    w_info = self.wire.describe(p)
                    seg = w_info["segment"]
                    cun = w_info["cun"]
                    flesh = w_info["flesh_mm"]
                    mus = w_info["muscle"]
                    mus_txt = f", brzusiec: {mus['label']}" if mus else ""
                    msg = f"Segment {seg} | Poziom {cun:.1f} cun | Tkanka {flesh:.1f} mm{mus_txt}"
                    _safe_print(f"\n>>> [KLIK NA SKÓRZE] {msg}")
                    self.info(f"Ciało: {msg}")
                except Exception:
                    pass
            else:
                _safe_print("\n>>> [KLIK POZA MODELEM] Odznaczenie punktu.")
                self.clear_selection()

    def select_point(self, key_or_code: str) -> bool:
        """Zaznacza punkt na ciele, podswietla go i wyswietla karte Deadmana."""
        normalized = key_or_code.replace(":", "·")
        # Obsluga powszechnych aliasow kanalow: TE <-> SJ, LIV <-> LR, CV <-> REN, GV <-> DU
        candidates = [normalized]
        if normalized.startswith("TE-"):
            candidates.append(normalized.replace("TE-", "SJ-", 1))
        elif normalized.startswith("SJ-"):
            candidates.append(normalized.replace("SJ-", "TE-", 1))
        elif normalized.startswith("LIV-"):
            candidates.append(normalized.replace("LIV-", "LR-", 1))
        elif normalized.startswith("LR-"):
            candidates.append(normalized.replace("LR-", "LIV-", 1))
        elif normalized.startswith("CV-"):
            candidates.append(normalized.replace("CV-", "REN-", 1))
        elif normalized.startswith("REN-"):
            candidates.append(normalized.replace("REN-", "CV-", 1))
        elif normalized.startswith("GV-"):
            candidates.append(normalized.replace("GV-", "DU-", 1))
        elif normalized.startswith("DU-"):
            candidates.append(normalized.replace("DU-", "GV-", 1))

        target = None
        for cand in candidates:
            if cand in self.atlas.by_key:
                target = self.atlas.by_key[cand]
                break
        if target is None:
            for cand in candidates:
                cand_clean = cand.upper().split("·")[0].split(":")[0].strip()
                for m in self.atlas.markers:
                    if m["code"].upper() == cand_clean or m["key"].upper() == cand.upper():
                        target = m
                        break
                if target is not None:
                    break
        if target is None:
            return False

        self.selected_point = target["key"]
        ch = target["channel"]
        if ch not in self.channels_on:
            self.set_channel(ch, True, render=False)

        if hasattr(self, "highlight_actor") and self.highlight_actor is not None:
            self.highlight_actor.SetPosition(target["pos"])
            self.highlight_actor.SetVisibility(True)

        d = self.atlas.describe(target["key"])
        self._update_card(d)
        self._print_deadman(d)
        code = d.get("code", target["key"])
        pl_name = polish_point_name(code, d.get("english", ""))
        pinyin = d.get("pinyin") or d.get("name_pinyin", "")
        if "(" in pl_name:
            self.info(f"Wybrano: {target['key']} — {pl_name} — [Esc] odznacz")
        else:
            disp_name = f"{pinyin} • {pl_name}" if pinyin else pl_name
            self.info(f"Wybrano: {target['key']} ({disp_name}) — [Esc] odznacz")
        self.pl.render()
        return True

    def clear_selection(self) -> None:
        """Czysci zaznaczenie punktu."""
        self.selected_point = None
        if hasattr(self, "highlight_actor") and self.highlight_actor is not None:
            self.highlight_actor.SetVisibility(False)
        if hasattr(self, "highlight_actor_pair") and self.highlight_actor_pair is not None:
            self.highlight_actor_pair.SetVisibility(False)
        if hasattr(self, "relation_tube_actor") and self.relation_tube_actor is not None:
            try:
                self.relation_tube_actor.SetVisibility(False)
            except Exception:
                pass
        self._update_card(None)
        self.info("Odznaczono punkt.")
        self.pl.render()

    def _get_point_pos(self, key_or_code: str) -> np.ndarray | None:
        normalized = key_or_code.replace(":", "·")
        candidates = [normalized]
        if normalized.startswith("TE-"):
            candidates.append(normalized.replace("TE-", "SJ-", 1))
        elif normalized.startswith("SJ-"):
            candidates.append(normalized.replace("SJ-", "TE-", 1))
        elif normalized.startswith("LIV-"):
            candidates.append(normalized.replace("LIV-", "LR-", 1))
        elif normalized.startswith("LR-"):
            candidates.append(normalized.replace("LR-", "LIV-", 1))
        elif normalized.startswith("CV-"):
            candidates.append(normalized.replace("CV-", "REN-", 1))
        elif normalized.startswith("REN-"):
            candidates.append(normalized.replace("REN-", "CV-", 1))
        elif normalized.startswith("GV-"):
            candidates.append(normalized.replace("GV-", "DU-", 1))
        elif normalized.startswith("DU-"):
            candidates.append(normalized.replace("DU-", "GV-", 1))

        for cand in candidates:
            if cand in self.atlas.by_key:
                return self.atlas.by_key[cand]["pos"]
            for side in ("R", "L"):
                with_side = f"{cand}·{side}"
                if with_side in self.atlas.by_key:
                    return self.atlas.by_key[with_side]["pos"]

        for cand in candidates:
            clean = cand.upper().split("·")[0].split(":")[0].strip()
            for m in self.atlas.markers:
                m_clean = m["code"].upper().split("·")[0].split(":")[0].strip()
                if m_clean == clean or m["key"].upper() == cand.upper():
                    return m["pos"]
        return None

    def _show_pair_relation(self, pos1: np.ndarray, pos2: np.ndarray, color: str = "#38bdf8") -> None:
        try:
            pts = np.array([pos1, pos2])
            lines = np.array([2, 0, 1])
            poly = pv.PolyData(pts, lines=lines)
            tube = poly.tube(radius=0.0022, n_sides=8)
            if hasattr(self, "relation_tube_actor") and self.relation_tube_actor is not None:
                try:
                    self.pl.remove_actor(self.relation_tube_actor)
                except Exception:
                    pass
            self.relation_tube_actor = self.pl.add_mesh(
                tube, color=color, smooth_shading=True, lighting=False,
                name="pair_relation_tube"
            )
            self.relation_tube_actor.SetVisibility(True)
            self.relation_tube_actor.SetPickable(False)
        except Exception:
            pass

    def _update_card(self, d: dict | None) -> None:
        """Prezentacja karty lokalizacji Deadmana w lewym dolnym rogu ekranu w jezyku polskim."""
        if d is None:
            txt = "[Kliknij punkt na ciele, aby wyswietlic opis Deadmana]"
            self.pl.add_text(txt, position=(14, 14), font_size=8,
                             color=THEME["dim"], name="deadman_card")
            return

        key = d.get("key", "")
        code = d.get("code", key)
        clean_code = code.replace("·L", "").replace("·R", "").replace(":L", "").replace(":R", "").strip()
        pl_rec = get_point_pl(clean_code)

        ch = d.get("channel", "")
        ch_info = self.atlas.channels.get(ch, {})
        ch_pl = (pl_rec.get("kanal_pl") if pl_rec else None) or ch_info.get("pl", ch)
        pinyin = (pl_rec.get("pinyin") if pl_rec else None) or d.get("pinyin") or d.get("name_pinyin", "")
        en = d.get("english") or d.get("name_en", "")
        meta = d.get("meta", {})

        side_str = " (P)" if "·R" in key or ":R" in key else (" (L)" if "·L" in key or ":L" in key else "")
        pl_name = (pl_rec.get("nazwa_pl") if pl_rec else None) or polish_point_name(code, en)
        disp_name = pl_name if "(" in pl_name else (f"{pinyin} ({pl_name})" if pinyin else pl_name)
        header = f"=== {code}{side_str}: {disp_name} - {ch_pl} [{ch}] ==="
        lines = [header]

        loc = (pl_rec.get("lokalizacja_pl") if pl_rec else None) or translate_location(clean_code, meta.get("location") or d.get("note_pl", ""), pinyin)
        if loc:
            lines.append("LOKALIZACJA ANATOMICZNA (Deadman):")
            for w in textwrap.wrap(loc, width=54):
                lines.append(f"  {w}")

        loc_note = (pl_rec.get("palpacja_pl") if pl_rec else None) or translate_palpation(clean_code, meta.get("locationNote") or "", pinyin)
        if loc_note:
            lines.append("WSKAZÓWKA PALPACYJNA:")
            for w in textwrap.wrap(loc_note, width=54):
                lines.append(f"  {w}")

        needling = (pl_rec.get("naklucie_pl") if pl_rec else None) or translate_needling(clean_code, meta.get("needling") or "")
        if needling:
            lines.append("NAKŁUCIE:")
            need_card = needling
            if len(need_card) > 170:
                cut = need_card[:170].rfind(" ")
                need_card = (need_card[:cut] if cut > 80 else need_card[:170]) + "..."
            for w in textwrap.wrap(need_card, width=54):
                lines.append(f"  {w}")

        dzialania = pl_rec.get("dzialania_pl") if pl_rec else None
        if dzialania:
            act_str = "; ".join(dzialania[:2])
            lines.append("DZIAŁANIE TCM:")
            for w in textwrap.wrap(act_str, width=54):
                lines.append(f"  {w}")
        else:
            actions = meta.get("actions")
            if actions:
                if isinstance(actions, list):
                    acts_list = [translate_action(a, clean_code) for a in actions if a and translate_action(a, clean_code)]
                    act_str = "; ".join(acts_list[:2])
                else:
                    act_str = translate_action(str(actions), clean_code)
                if act_str:
                    lines.append("DZIAŁANIE TCM:")
                    for w in textwrap.wrap(act_str, width=54):
                        lines.append(f"  {w}")

        # Relacje i zasady TCM dla punktu
        tcm_principles = []
        yl = get_yuan_luo_info(clean_code)
        if yl:
            p_part = yl[1].replace("Punkt Partnerski (Gość):", "Gość:").strip()
            tcm_principles.append(f"Yuan-Luo: {yl[0].split(':')[-1].strip()} <-> {p_part}")
        sm = get_shu_mu_info(clean_code)
        if sm:
            tcm_principles.append(f"Shu-Mu: {sm.get('organ', '')} (Mu: {sm.get('mu')}, Shu: {sm.get('shu')})")
        xi = get_xi_cleft_info(clean_code)
        if xi:
            tcm_principles.append(f"Szczelina Xi: {xi.split('(')[0].strip()}")
        bm = get_extraordinary_vessel_info(clean_code)
        if bm:
            tcm_principles.append(f"Ba Mai: {bm[0].split(':')[-1].strip()}")

        if tcm_principles:
            lines.append("ZASADY & POWIĄZANIA TCM:")
            for p in tcm_principles[:2]:
                for w in textwrap.wrap(p, width=54):
                    lines.append(f"  {w}")

        lines.append("[Esc] zamknij opis  --  pelne dane w konsoli")
        self.pl.add_text("\n".join(lines), position=(14, 14), font_size=8,
                         color="#f0f4f8", name="deadman_card")

    def _print_deadman(self, d: dict) -> None:
        """Pelny, nieobciety opis punktu wg Petera Deadmana w konsoli (PL)."""
        key = d.get("key", "")
        code = d.get("code", key)
        clean_code = code.replace("·L", "").replace("·R", "").replace(":L", "").replace(":R", "").strip()
        pl_rec = get_point_pl(clean_code)

        ch = d.get("channel", "")
        ch_info = self.atlas.channels.get(ch, {})
        ch_pl = (pl_rec.get("kanal_pl") if pl_rec else None) or ch_info.get("pl", ch)
        pinyin = (pl_rec.get("pinyin") if pl_rec else None) or d.get("pinyin") or d.get("name_pinyin", "")
        en = d.get("english") or d.get("name_en", "")
        meta = d.get("meta", {})

        side_str = " (Prawa strona)" if "·R" in key or ":R" in key else (" (Lewa strona)" if "·L" in key or ":L" in key else "")
        pl_name = (pl_rec.get("nazwa_pl") if pl_rec else None) or polish_point_name(code, en)
        disp_name = pl_name if "(" in pl_name else (f"{pinyin} ({pl_name})" if pinyin else pl_name)
        out = [
            "",
            "╔" + "═" * 76 + "╗",
            f"║  PUNKT: {key:<12} {disp_name}{side_str}".ljust(77) + "║",
            f"║  Kanał: {ch_pl} [{ch}]".ljust(77) + "║",
            "╠" + "═" * 76 + "╣",
        ]

        cats = (pl_rec.get("kategorie_pl") if pl_rec else None)
        if cats:
            cats_str = ", ".join(cats) if isinstance(cats, list) else str(cats)
            out.append(f"║  Kategorie: {cats_str}".ljust(77) + "║")
        elif meta.get("categories"):
            cats_raw = ", ".join(meta["categories"]) if isinstance(meta["categories"], list) else str(meta["categories"])
            out.append(f"║  Kategorie: {cats_raw}".ljust(77) + "║")

        loc = (pl_rec.get("lokalizacja_pl") if pl_rec else None) or translate_location(clean_code, meta.get("location") or d.get("note_pl", ""), pinyin)
        if loc:
            out.append("╟" + "─" * 76 + "╢")
            out.append("║  LOKALIZACJA ANATOMICZNA (Peter Deadman):".ljust(77) + "║")
            for w in textwrap.wrap(loc, width=72):
                out.append(f"║    {w}".ljust(77) + "║")

        loc_note = (pl_rec.get("palpacja_pl") if pl_rec else None) or translate_palpation(clean_code, meta.get("locationNote") or "", pinyin)
        if loc_note:
            out.append("╟" + "─" * 76 + "╢")
            out.append("║  WSKAZÓWKA PALPACYJNA I WERYFIKACJA:".ljust(77) + "║")
            for w in textwrap.wrap(loc_note, width=72):
                out.append(f"║    {w}".ljust(77) + "║")

        needling = (pl_rec.get("naklucie_pl") if pl_rec else None) or translate_needling(clean_code, meta.get("needling") or "")
        if needling:
            out.append("╟" + "─" * 76 + "╢")
            out.append("║  TECHNIKA NAKŁUCIA I BEZPIECZEŃSTWO:".ljust(77) + "║")
            for w in textwrap.wrap(needling, width=72):
                out.append(f"║    {w}".ljust(77) + "║")

        dzialania = (pl_rec.get("dzialania_pl") if pl_rec else None)
        if dzialania:
            out.append("╟" + "─" * 76 + "╢")
            out.append("║  DZIAŁANIE TRADYCYJNE (TCM):".ljust(77) + "║")
            for a in dzialania:
                out.append(f"║    * {a}".ljust(77) + "║")
        elif meta.get("actions"):
            out.append("╟" + "─" * 76 + "╢")
            out.append("║  DZIAŁANIE TRADYCYJNE (TCM):".ljust(77) + "║")
            acts = meta["actions"]
            if isinstance(acts, list):
                for a in acts:
                    if not a:
                        continue
                    act_text = translate_action(a, clean_code)
                    if act_text:
                        out.append(f"║    * {act_text}".ljust(77) + "║")
            else:
                act_text = translate_action(str(acts), clean_code)
                if act_text:
                    for w in textwrap.wrap(act_text, width=72):
                        out.append(f"║    {w}".ljust(77) + "║")

        ind = (pl_rec.get("wskazania_pl") if pl_rec else None) or translate_indications(clean_code, str(meta.get("indications") or ""))
        if ind:
            out.append("╟" + "─" * 76 + "╢")
            out.append("║  GŁÓWNE WSKAZANIA KLINICZNE:".ljust(77) + "║")
            for w in textwrap.wrap(ind, width=72):
                out.append(f"║    {w}".ljust(77) + "║")

        extraordinary_info = get_extraordinary_vessel_info(clean_code)
        if extraordinary_info:
            out.append("╟" + "─" * 76 + "╢")
            out.append("║  8 NACZYŃ NADZWYCZAJNYCH (Bā Mài Jiāo Huì Xué):".ljust(77) + "║")
            for line in extraordinary_info:
                for w in textwrap.wrap(line, width=72):
                    out.append(f"║    {w}".ljust(77) + "║")

        yuan_luo_info = get_yuan_luo_info(clean_code)
        if yuan_luo_info:
            out.append("╟" + "─" * 76 + "╢")
            out.append("║  POŁĄCZENIE MERIDIANOWE (Para Yuan-Luo):".ljust(77) + "║")
            for line in yuan_luo_info:
                for w in textwrap.wrap(line, width=72):
                    out.append(f"║    {w}".ljust(77) + "║")

        xi_cleft_info = get_xi_cleft_info(clean_code)
        if xi_cleft_info:
            out.append("╟" + "─" * 76 + "╢")
            out.append("║  PUNKT SZCZELINY XI (Ostry ból i stany nagłe):".ljust(77) + "║")
            for w in textwrap.wrap(xi_cleft_info, width=72):
                out.append(f"║    {w}".ljust(77) + "║")

        out.append("╚" + "═" * 76 + "╝")
        out.append("")
        _safe_print("\n".join(out))

    # --------------------------------------------------------------- pomocnicze
    def info(self, text: str) -> None:
        self.pl.add_text(text, position="lower_right", font_size=9,
                         color=THEME["text"], name="info")

    def toggle_cun_grid(self) -> None:
        """Przelacza widocznosc siatki cunowej [g]."""
        self.show_wire = not self.show_wire
        if hasattr(self, "wire_actor") and self.wire_actor is not None:
            self.wire_actor.SetVisibility(self.show_wire)
        self._sync_legend()
        self._update_dock_button("g", label="Siatka cun: ON" if self.show_wire else "Siatka cun: OFF",
                                border="#be9a84" if self.show_wire else "#485463",
                                text_color="#d9c2b2" if self.show_wire else "#8e98a8")
        self.info(f"Siatka cunowa: {'WLACZONA' if self.show_wire else 'WYLACZONA'} [g]")
        self.pl.render()

    def _wire_cb(self, state: bool) -> None:
        self.show_wire = bool(state)
        if hasattr(self, "wire_actor") and self.wire_actor is not None:
            self.wire_actor.SetVisibility(self.show_wire)
        self._update_dock_button("g", label="Siatka cun: ON" if self.show_wire else "Siatka cun: OFF",
                                border="#be9a84" if self.show_wire else "#485463",
                                text_color="#d9c2b2" if self.show_wire else "#8e98a8")
        self.info(f"Siatka cunowa: {'WLACZONA' if self.show_wire else 'WYLACZONA'} [g]")
        self.pl.render()

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
        if self.draw_mode and self._press_pos is not None:
            # Rysowanie na biezaco podczas przeciagania mysza po ciele
            x, y = self._mouse()
            p, n = self.pick(x, y)
            if p is not None and n is not None:
                pt = p + n * 0.003
                if not self._current_stroke:
                    self._current_stroke.append(pt)
                else:
                    last_pt = self._current_stroke[-1]
                    if np.linalg.norm(pt - last_pt) >= 0.0025:  # probkowanie co 2.5 mm
                        self._current_stroke.append(pt)
                        self._update_current_stroke_actor()
            return

        if not self.camera_frozen and not self.draw_mode:
            if hasattr(obj, "OnMouseMove"):
                obj.OnMouseMove()                    # obrot i przesuw dalej dzialaja gdy nie rysujemy i nie zamrozone
            elif hasattr(self, "_style") and hasattr(self._style, "OnMouseMove"):
                self._style.OnMouseMove()

        # Podczas przeciągania myszą (np. obracania kamery) NIE wykonujemy ciężkiego pickowania i podwójnego renderowania
        if self._press_pos is not None:
            return

        now = time.perf_counter()
        if now - self._last_move < 0.030:    # ok. 33 razy na sekunde wystarczy (plynne, 3x mniej obciazenia)
            return

        x, y = self._mouse()
        if hasattr(self, "_last_pick_xy") and self._last_pick_xy is not None:
            dx = x - self._last_pick_xy[0]
            dy = y - self._last_pick_xy[1]
            if dx * dx + dy * dy < 9:        # martwa strefa 3 pikseli
                return

        self._last_move = now
        self._last_pick_xy = (x, y)
        p, n = self.pick(x, y)
        self.place_cursor(p, n)

    def place_cursor(self, p, n) -> None:
        """Kladzie pierscien na skorze. Informuje o odleglosci do punktu i poziomie cun."""
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

        wire_txt = ""
        try:
            w_info = self.wire.describe(p)
            wire_txt = f"{w_info['segment']} • {w_info['cun']:.1f} cun"
        except Exception:
            pass

        code = None
        # Sprawdzamy zblizenie do widocznych punktow w aktywnym rejonie i kanale
        reg_mask = (self._mreg == self.active_region) if self.active_region != "ALL" else np.ones(len(self._mpos), dtype=bool)
        if self.channels_on:
            chan_mask = np.isin(self._mchan, list(self.channels_on))
        else:
            chan_mask = np.ones(len(self._mpos), dtype=bool) if self.show_points_mode else np.zeros(len(self._mpos), dtype=bool)
        visible = reg_mask & chan_mask
        d = np.linalg.norm(self._mpos - p, axis=1)
        d[~visible] = np.inf
        j = int(np.argmin(d))
        if d[j] <= NEAR_POINT:
            m_pt = self.atlas.markers[j]
            el = self.atlas.channels[m_pt["channel"]].get("element", "-")
            el_name = ELEMENTS.get(el, ELEMENTS['-'])[0].lower()
            pinyin = m_pt.get("pinyin") or m_pt.get("name_pinyin", "")
            pl_name = polish_point_name(m_pt["code"], m_pt.get("english", ""))
            name_str = f" ({pinyin} • {pl_name})" if pinyin else f" ({pl_name})"
            code = (f"{m_pt['key']}{name_str} - {d[j] * 1000:.0f} mm - "
                    f"{wire_txt} [{el_name}] (kliknij, aby wybrac)")
        if code is None and wire_txt:
            code = f"Siatka cunowa: {wire_txt}"
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
        """Legenda-panel: podział na Siatkę CUN, 12 Kanałów Głównych (Wu Xing) i 8 Naczyń Nadzwyczajnych."""
        self._channel_widgets: dict = {}
        self._element_widgets: dict = {}
        groups = {e: [] for e in ELEMENT_ORDER}
        for ch, info in self.atlas.channels.items():
            el = info.get("element", "-")
            if el in groups:
                groups[el].append(ch)
        self._element_groups = groups

        x, y = 14, self.pl.window_size[1] - 44

        # 1. Siatka cunowa
        self._wire_widget = self.pl.add_checkbox_button_widget(
            self._wire_cb, value=self.show_wire, position=(x, y), size=18,
            border_size=1, color_on=THEME["ring"], color_off=THEME["off"])
        self.pl.add_text("SIATKA CUNOWA [g]", position=(x + 26, y + 1), font_size=9,
                         color=THEME["ring"], name="leg_wire")
        y -= 26

        # 2. 12 Meridianów Głównych (Wu Xing)
        self.pl.add_text("12 MERIDIANOW GLOWNYCH [m]", position=(x, y + 6), font_size=9,
                         color=THEME["text"], name="leg_title")
        y -= 22
        for el in ELEMENT_ORDER:
            chans = groups[el]
            if not chans:
                continue
            name, yin, _ = ELEMENTS[el]
            el_on = bool(chans) and all(c in self.channels_on for c in chans)
            self._element_widgets[el] = self.pl.add_checkbox_button_widget(
                self._element_cb(el), value=el_on, position=(x, y), size=16,
                border_size=1, color_on=yin, color_off=THEME["off"])
            self.pl.add_text(name, position=(x + 24, y + 1), font_size=8,
                             color=THEME["text"], name=f"leg_el_{el}")
            y -= 20
            for ch in chans:
                info = self.atlas.channels[ch]
                ch_on = ch in self.channels_on
                self._channel_widgets[ch] = self.pl.add_checkbox_button_widget(
                    self._channel_cb(ch), value=ch_on, position=(x + 18, y), size=14,
                    border_size=1, color_on=self.channel_color[ch], color_off=THEME["off"])
                kind = "yin" if ch in YIN else "yang"
                self.pl.add_text(f"{ch}  {info['pl']}  ({kind})", position=(x + 38, y),
                                 font_size=8, color=THEME["dim"], name=f"leg_ch_{ch}")
                y -= 18
            y -= 4

        # 3. Naczynia Nadzwyczajne (Qi Jing Ba Mai)
        y -= 4
        self.pl.add_text("8 NACZYN NADZWYCZAJNYCH [e]", position=(x, y + 6), font_size=9,
                         color="#d8b4fe", name="leg_extra_title")
        y -= 22
        for ch in ("DU", "REN"):
            if ch in self.atlas.channels:
                info = self.atlas.channels[ch]
                ch_on = ch in self.channels_on
                self._channel_widgets[ch] = self.pl.add_checkbox_button_widget(
                    self._channel_cb(ch), value=ch_on, position=(x + 18, y), size=14,
                    border_size=1, color_on=self.channel_color[ch], color_off=THEME["off"])
                axis_txt = "os tylu / yang" if ch == "DU" else "os przodu / yin"
                code_disp = "RN" if ch == "REN" else ch
                self.pl.add_text(f"{code_disp}  {info['pl']}  ({axis_txt})", position=(x + 38, y),
                                 font_size=8, color="#e9d5ff", name=f"leg_ch_{ch}")
                y -= 19

    def _channel_cb(self, ch: str):
        def cb(state: bool) -> None:
            self.set_channel(ch, bool(state))
        return cb

    def _element_cb(self, el: str):
        def cb(state: bool) -> None:
            for ch in self._element_groups.get(el, []):
                self.set_channel(ch, bool(state), render=False)
            self.pl.render()
        return cb

    def _update_points_visibility(self) -> None:
        """Steruje widocznoscia punktow w zaleznosci od aktywnych kanalow i wybranego rejonu."""
        for (ch, reg), act in self.point_actors.items():
            if self.channels_on:
                ch_visible = (ch in self.channels_on)
            else:
                ch_visible = self.show_points_mode
            reg_visible = (self.active_region == "ALL" or reg == self.active_region)
            act.SetVisibility(ch_visible and reg_visible)

    def set_channel(self, ch: str, on: bool, render: bool = True) -> None:
        ch = ch.upper().strip()
        if ch == "RN":
            ch = "REN"
        elif ch == "TE":
            ch = "SJ"
        if ch not in self.atlas.channels:
            return
        (self.channels_on.add if on else self.channels_on.discard)(ch)
        self._update_points_visibility()
        if on:
            tubes = self._get_or_create_tubes(ch)
            if tubes is not None:
                tubes.SetVisibility(True)
            if ch in ("BL", "DU"):
                _safe_print(f">>> [MERIDIAN {ch}] Przebieg grzbietowy i nogi. Wcisnij [2], aby spojrzec na plecy.")
        else:
            entry = self.channel_actors.get(ch, {})
            if entry.get("meridian") is not None:
                entry["meridian"].SetVisibility(False)

        acts = self.channel_actors.get(ch, {})
        if acts.get("labels") is not None:
            acts["labels"].SetVisibility(on and self.show_labels)
        self._sync_legend()
        if on:
            info = self.atlas.channels.get(ch, {})
            el = info.get("element", "-")
            name = ELEMENTS.get(el, ELEMENTS["-"])[0].lower() if el != "-" else "naczynie nadzwyczajne"
            self.info(f"{ch} {info.get('pl', '')} - {name}, {'yin' if ch in YIN else 'yang'}")
        if render:
            self.pl.render()

    def toggle_all_main_meridians(self) -> None:
        """[m] - Przelacza 12 meridianow glownych (Zang-Fu) z punktami i liniami 3D."""
        main_12 = {"LU", "LI", "ST", "SP", "HT", "SI", "BL", "KI", "PC", "SJ", "GB", "LR"}
        all_on = main_12.issubset(self.channels_on)
        if all_on:
            for ch in main_12:
                self.set_channel(ch, False, render=False)
            self._update_dock_button("m", label="12 Glownych: OFF", border="#485463", text_color="#8e98a8")
            self.info("12 Meridianow Glownych: WYLACZONE [m]")
            _safe_print(">>> [12 MERIDIANOW GLOWNYCH] Wylaczono wszystkie kanaly Zang-Fu.")
        else:
            for ch in main_12:
                self.set_channel(ch, True, render=False)
            self._update_dock_button("m", label="12 Glownych: ON", border="#38bdf8", text_color="#7dd3fc")
            self.info("12 Meridianow Glownych: WLACZONE (punkty i linie 3D) [m]")
            _safe_print(">>> [12 MERIDIANOW GLOWNYCH] Wlaczono wszystkie 12 kanalow Zang-Fu (linie i punkty).")
        self.pl.render()

    def _sync_legend(self) -> None:
        """Stan pol wyboru ma odpowiadac temu, co widac - takze po [p], [g] i po grupie."""
        if hasattr(self, "_wire_widget") and self._wire_widget is not None:
            self._wire_widget.GetRepresentation().SetState(self.show_wire)
        for ch, w in self._channel_widgets.items():
            w.GetRepresentation().SetState(ch in self.channels_on)
        for el, w in self._element_widgets.items():
            chans = self._element_groups.get(el, [])
            w.GetRepresentation().SetState(bool(chans) and all(c in self.channels_on
                                                               for c in chans))

    def toggle_points(self) -> None:
        """[p] - wlacza / ukrywa punkty akupunkturowe."""
        self.show_points_mode = not self.show_points_mode
        self._update_points_visibility()
        self._update_dock_button("p", label="Punkty: WLACZONE" if self.show_points_mode else "Punkty: OFF",
                                border="#ffe36a" if self.show_points_mode else "#485463",
                                text_color="#ffe36a" if self.show_points_mode else "#8e98a8")
        self.info("Punkty akupunkturowe: WIDOCZNE [p]" if self.show_points_mode else "Punkty akupunkturowe: UKRYTE [p]")
        self.pl.render()

    def set_region(self, reg_id: str) -> None:
        """Wybiera aktywny rejon anatomiczny punktow i kadruje kamere."""
        if reg_id not in REGIONS:
            return
        self.active_region = reg_id
        self._update_points_visibility()

        r_info = REGIONS[reg_id]
        label = r_info["label"]
        self._update_dock_button(
            "b",
            label=f"Rejon: {label}",
            border="#ffe36a" if reg_id != "ALL" else "#485463",
            text_color="#ffe36a" if reg_id != "ALL" else "#f0f4f8"
        )

        # Kadrowanie kamery na wybrany rejon
        cam = self.pl.camera
        cam.focal_point = r_info["focal"]
        cam.position = r_info["pos"]
        cam.view_up = (0, 0, 1)
        self.pl.renderer.ResetCameraClippingRange()

        msg = f"Rejon punktow: {r_info['full_name']} [{reg_id}]"
        _safe_print(f"\n>>> [REJON ANATOMICZNY] {msg}")
        self.info(msg)
        self.pl.render()

    def cycle_region(self) -> None:
        """[b] Przelacza kolejno rejony anatomiczne."""
        cur_idx = REGION_ORDER.index(self.active_region) if self.active_region in REGION_ORDER else 0
        next_idx = (cur_idx + 1) % len(REGION_ORDER)
        self.set_region(REGION_ORDER[next_idx])

    def toggle_meridian_lines(self, state: bool | None = None) -> None:
        """[m] Przelacza ciagly przebieg linii/rurek meridianow w kolorach Wu Xing."""
        if state is not None:
            self.show_meridian_lines = state
        else:
            self.show_meridian_lines = not self.show_meridian_lines

        target_chans = self.channels_on if self.channels_on else set(self.atlas.channels.keys())
        for ch in self.atlas.channels:
            if self.show_meridian_lines and (ch in target_chans):
                tubes = self._get_or_create_tubes(ch)
                if tubes is not None:
                    tubes.SetVisibility(True)
            else:
                entry = self.channel_actors.get(ch, {})
                if entry.get("meridian") is not None:
                    entry["meridian"].SetVisibility(False)

        self._update_dock_button(
            "m",
            label="Meridiany: ON" if self.show_meridian_lines else "Linie meridianow",
            border="#38bdf8" if self.show_meridian_lines else "#485463",
            text_color="#7dd3fc" if self.show_meridian_lines else "#f0f4f8"
        )
        msg = f"Przebieg linii meridianow: {'WLACZONY' if self.show_meridian_lines else 'WYLACZONY'} [m]"
        _safe_print(f"\n>>> [LINIE MERIDIANOW] {msg}")
        self.info(msg)
        self.pl.render()


    def toggle_labels(self) -> None:
        self.show_labels = not self.show_labels
        for ch, acts in self.channel_actors.items():
            if acts.get("labels") is not None:
                acts["labels"].SetVisibility(self.show_labels and (ch in self.channels_on if self.channels_on else self.show_points_mode))
        self._update_dock_button("l", label="Etykiety: ON" if self.show_labels else "Etykiety kodow",
                                text_color="#ffe36a" if self.show_labels else "#f0f4f8")
        if self.show_labels and not self.channels_on and not self.show_points_mode:
            self.info("Kody punktow pokaza sie po wlaczeniu punktow lub kanalu.")
        self.pl.render()

    # ─── 8 NACZYŃ NADZWYCZAJNYCH (BA MAI) ───
    def _get_or_create_vessel_lines(self, vid: str) -> list:
        if vid in self.vessel_actors:
            return self.vessel_actors[vid]
        v = EXTRAORDINARY_VESSELS.get(vid)
        if not v:
            return []
        color = v["color"]
        seq = v["sequence"]
        actors = []
        sides = ("C",) if vid in ("REN_MAI", "DU_MAI") else ("L", "R")
        for side in sides:
            pts = []
            for code in seq:
                cand_keys = [f"{code}·{side}", f"{code}:{side}", code]
                m = None
                for k in cand_keys:
                    if k in self.atlas.by_key:
                        m = self.atlas.by_key[k]
                        break
                if m is None:
                    for mk in self.atlas.markers:
                        if mk["code"].upper() == code.upper() and (side == "C" or mk.get("side") == side or f"·{side}" in mk["key"]):
                            m = mk
                            break
                if m is not None:
                    pts.append(m["pos"])
            if len(pts) >= 2:
                pts_arr = np.array(pts, dtype=float)
                # Lekka wstęga/polilinia 3D (szybki rysunek - minimalny narzut geometryczny)
                lines = np.hstack([[len(pts_arr)], np.arange(len(pts_arr))])
                poly = pv.PolyData(pts_arr, lines=lines)
                tube = poly.tube(radius=0.0034, n_sides=6)
                act = self.pl.add_mesh(tube, color=color, smooth_shading=True, lighting=False, name=f"vessel_{vid}_{side}")
                act.SetPickable(False)
                actors.append(act)
        self.vessel_actors[vid] = actors
        return actors

    def cycle_extraordinary(self) -> None:
        order = [None, "DAI_MAI", "CHONG_MAI", "YIN_QIAO", "YANG_QIAO", "YIN_WEI", "YANG_WEI", "ALL"]
        curr_idx = -1
        if self.active_vessel is not None:
            try:
                curr_idx = order.index(self.active_vessel)
            except ValueError:
                curr_idx = -1
        nxt = order[(curr_idx + 1) % len(order)]
        self.set_extraordinary(nxt)

    def set_extraordinary(self, vid: str | None) -> None:
        for v_actors in self.vessel_actors.values():
            for act in v_actors:
                try:
                    act.SetVisibility(False)
                except Exception:
                    pass

        self.active_vessel = vid
        if vid is None:
            self._update_dock_button("e", label="8 Naczyn: OFF", color="#261b2e", border="#a855f7", text_color="#d8b4fe")
            self.info("8 Naczyn Nadzwyczajnych: WYLACZONE [e].")
            _safe_print(">>> [8 NACZYN] Wylaczono podglad naczyn nadzwyczajnych.")
        elif vid == "ALL":
            for v_key in ("DAI_MAI", "CHONG_MAI", "YIN_QIAO", "YANG_QIAO", "YIN_WEI", "YANG_WEI"):
                for act in self._get_or_create_vessel_lines(v_key):
                    act.SetVisibility(True)
            self._update_dock_button("e", label="8 Naczyn: WSZYSTKIE", color="#3b1d54", border="#c084fc", text_color="#f5d0fe")
            self.info("8 Naczyn Nadzwyczajnych: Wlaczone wszystkie naczynia (Ba Mai).")
            _safe_print(">>> [8 NACZYN] Wlaczono wszystkie Naczynia Nadzwyczajne (Ba Mai).")
        else:
            v = get_vessel(vid)
            if v:
                actual_id = v["id"]
                self.active_vessel = actual_id
                for act in self._get_or_create_vessel_lines(actual_id):
                    act.SetVisibility(True)
                lbl = v["name_pinyin"]
                self._update_dock_button("e", label=f"Ba Mai: {lbl}", color="#3b1d54", border="#a855f7", text_color="#f5d0fe")
                self.info(f"Aktywne Naczynie: {v['name_pl']} (punkt otwierajacy: {v['master_point']})")
                _safe_print(f"\n╔{'═'*72}╗")
                _safe_print(f"║  8 NACZYŃ NADZWYCZAJNYCH: {v['name_pl']}")
                _safe_print(f"║  Punkt Otwierający (Master): {v['master_point']}   |   Sprzężony: {v['coupled_point']}")
                _safe_print(f"║  Obszar: {v['area']}")
                _safe_print(f"║  Dynamika: {v['dynamics']}")
                _safe_print(f"║  Główne wskazania: {v['indications']}")
                _safe_print(f"╚{'═'*72}╝\n")
                if v.get("master_point"):
                    self.select_point(v["master_point"])
        self.pl.render()

    # ─── JING JIN & TAŚMY MYERSA ───
    def cycle_jing_jin(self) -> None:
        order = [None, "SBL", "SFL", "LL", "DFL", "DFAL", "SFAL", "DBAL", "SBAL"]
        curr_idx = -1
        if self.active_myers is not None:
            try:
                curr_idx = order.index(self.active_myers)
            except ValueError:
                curr_idx = -1
        nxt = order[(curr_idx + 1) % len(order)]
        self.set_jing_jin(nxt)

    def set_jing_jin(self, line_id: str | None) -> None:
        self.active_myers = line_id
        if line_id is None:
            self._update_dock_button("j", label="Jing Jin: OFF", color="#16222b", border="#38bdf8", text_color="#7dd3fc")
            self.info("Jing Jin & Myers: WYLACZONE [j].")
            _safe_print(">>> [JING JIN] Wylaczono podglad pasm powiezio-miesniowych Myersa.")
        else:
            m = get_myers_line(line_id)
            if m:
                jj_ch = m["jing_jin_channel"]
                self._update_dock_button("j", label=f"Jing Jin: {line_id}", color="#0f344d", border="#38bdf8", text_color="#bae6fd")
                self.info(f"Jing Jin / Myers: {m['name_pl']} [{jj_ch}]")
                _safe_print(f"\n╔{'═'*72}╗")
                _safe_print(f"║  TAŚMA MYERSA / JING JIN: {m['name_pl']}")
                _safe_print(f"║  Kanał Jing Jin: {m['jing_jin_name']} [{jj_ch}]")
                _safe_print(f"║  Biotensegracja i funkcja: {m['function']}")
                _safe_print(f"║  Wzorce dysfunkcji i objawy Lingshu: {m['dysfunction']}")
                _safe_print(f"║  Pasma składowe: {', '.join(m['tracks'][:3])}")
                _safe_print(f"╚{'═'*72}╝\n")
                self.set_channel(jj_ch, True, render=False)
                if line_id == "SBL":
                    self.set_region("BACK")
                elif line_id in ("SFL", "DFL"):
                    self.set_region("CHEST")
                elif line_id in ("DFAL", "SFAL", "DBAL", "SBAL"):
                    self.set_region("ARMS")
                elif line_id == "LL":
                    self._view_right()
        self.pl.render()

    # ─── BIAO-LI (YIN-YANG / GŁĘBOKOŚĆ) ───
    def cycle_biao_li(self) -> None:
        nxt = (self.biao_li_level % 3) + 1
        self.set_biao_li(nxt)

    def set_biao_li(self, level: int) -> None:
        self.biao_li_level = level
        READY_DIR = Path(__file__).resolve().parent.parent / "assets" / "ready"
        if level == 1:
            self.skin_actor.prop.SetOpacity(1.0)
            self.skin_actor.SetVisibility(True)
            if self.skeleton_actor:
                self.skeleton_actor.SetVisibility(False)
            self._update_dock_button("y", label="Biao-Li: POWIERZCHNIA", color="#232219", border="#eab308", text_color="#fef08a")
            self.info("Warstwa BIAO: Powierzchnia porcelanowa, skora, punkty i cun.")
            _safe_print(">>> [BIAO-LI] Poziom 1: BIAO (表 / Yang) - Powierzchnia ciala, skora i siatka cunowa.")
        elif level == 2:
            self.skin_actor.prop.SetOpacity(0.48)
            self.skin_actor.SetVisibility(True)
            if self.skeleton_actor:
                self.skeleton_actor.SetVisibility(False)
            self._update_dock_button("y", label="Biao-Li: POWIEZ (JING)", color="#362f13", border="#facc15", text_color="#fef08a")
            self.info("Warstwa JING JIN: Półprzezroczysta powłoka, pasma powięziowe.")
            _safe_print(">>> [BIAO-LI] Poziom 2: JING JIN (经筋) - Warstwa powiezio-miesniowa, biotensegracja.")
        elif level == 3:
            self.skin_actor.prop.SetOpacity(0.18)
            self.skin_actor.SetVisibility(True)
            if self.skeleton_actor is None:
                skel_path = READY_DIR / "skeleton.vtp"
                if skel_path.exists():
                    skel_mesh = pv.read(skel_path)
                    self.skeleton_actor = self.pl.add_mesh(
                        skel_mesh, color="#e2e8f0", smooth_shading=True,
                        ambient=0.4, diffuse=0.6, name="skeleton_biao_li"
                    )
            if self.skeleton_actor:
                self.skeleton_actor.SetVisibility(True)
            self._update_dock_button("y", label="Biao-Li: WNETRZE (LI)", color="#3f1e1e", border="#ef4444", text_color="#fca5a5")
            self.info("Warstwa LI: Wnętrze, głęboki szkielet Z-Anatomy i osie narządowe.")
            _safe_print(">>> [BIAO-LI] Poziom 3: LI (里 / Yin) - Wnetrze ciala, szkielet, osie Zang-Fu.")
        self.pl.render()

    # ─── INTERAKTYWNE SKOKI DO ZASAD TCM ───
    def cycle_yuan_luo(self) -> None:
        """Przelacza cyklicznie kolejne pary Yuan-Luo (Gospodarz-Gosc) [u]."""
        yuan_keys = list(YUAN_LUO_PAIRS.keys())
        if not hasattr(self, "_yuan_idx"):
            self._yuan_idx = -1
        self._yuan_idx = (self._yuan_idx + 1) % len(yuan_keys)
        y_code = yuan_keys[self._yuan_idx]
        self.jump_to_yuan_luo(y_code)

    def jump_to_yuan_luo(self, channel: str) -> None:
        ch = channel.upper().strip()
        target_pair = None
        for yuan_code, (role, partner) in YUAN_LUO_PAIRS.items():
            if yuan_code.startswith(ch + "-") or yuan_code == ch:
                target_pair = (yuan_code, partner.split()[0], role, partner)
                break
        if not target_pair:
            for yuan_code, (role, partner) in YUAN_LUO_PAIRS.items():
                if partner.startswith(ch + "-") or partner.split()[0] == ch:
                    target_pair = (yuan_code, partner.split()[0], role, partner)
                    break
        if target_pair:
            y_code, l_code, role, partner = target_pair
            y_ch = y_code.split("-")[0]
            l_ch = l_code.split("-")[0]
            self.set_channel(y_ch, True, render=False)
            self.set_channel(l_ch, True, render=False)
            self.select_point(y_code)

            pos_y = self._get_point_pos(y_code)
            pos_l = self._get_point_pos(l_code)
            if pos_l is not None and hasattr(self, "highlight_actor_pair") and self.highlight_actor_pair is not None:
                self.highlight_actor_pair.SetPosition(pos_l)
                self.highlight_actor_pair.SetVisibility(True)
            if pos_y is not None and pos_l is not None:
                self._show_pair_relation(pos_y, pos_l, color="#d29922")

            self._update_dock_button("u", label=f"Yuan-Luo: {y_code}+{l_code}",
                                    color="#2a2415", border="#d29922", text_color="#fef08a")
            _safe_print(f"\n╔{'═'*72}╗")
            _safe_print(f"║  ZASADA TCM: PARY GOSPODARZ - GOŚĆ (YUAN - LUO)")
            _safe_print(f"║  Główny punkt Yuan-Źródłowy (Gospodarz): {y_code} [{role}]")
            _safe_print(f"║  Sprzężony punkt Luo-Łączący (Gość): {partner}")
            _safe_print(f"║  Mechanizm: Przenosi nadmiar z kanału przeciążonego do pustego,")
            _safe_print(f"║  harmonizuje parę narządów Zang-Fu w relacji wnętrze-powierzchnia.")
            _safe_print(f"╚{'═'*72}╝\n")
            self.info(f"Yuan-Luo: {y_code} (Gospodarz) + {l_code} (Gość)")
            self.pl.render()
        else:
            _safe_print(f">>> Nie znaleziono pary Yuan-Luo dla kanału: {channel}. Dostępne: LU, LI, ST, SP, HT, SI, BL, KI, PC, TE, GB, LR")

    def cycle_shu_mu(self) -> None:
        """Przelacza cyklicznie kolejne pary Shu-Mu (Przod-Mu / Tyl-Shu) [o]."""
        channels = list(SHU_MU_PAIRS.keys())
        if not hasattr(self, "_shumu_idx"):
            self._shumu_idx = -1
        self._shumu_idx = (self._shumu_idx + 1) % len(channels)
        ch = channels[self._shumu_idx]
        self.jump_to_shu_mu(ch)

    def jump_to_shu_mu(self, channel: str) -> None:
        ch = channel.upper().strip()
        info = get_shu_mu_info(ch)
        if info:
            mu_pt = info["mu"]
            shu_pt = info["shu"]
            organ = info["organ"]
            desc = info["desc"]

            mu_ch = mu_pt.split("-")[0]
            shu_ch = shu_pt.split("-")[0]
            self.set_channel(mu_ch, True, render=False)
            self.set_channel(shu_ch, True, render=False)

            self.select_point(mu_pt)
            pos_mu = self._get_point_pos(mu_pt)
            pos_shu = self._get_point_pos(shu_pt)
            if pos_shu is not None and hasattr(self, "highlight_actor_pair") and self.highlight_actor_pair is not None:
                self.highlight_actor_pair.SetPosition(pos_shu)
                self.highlight_actor_pair.SetVisibility(True)
            if pos_mu is not None and pos_shu is not None:
                self._show_pair_relation(pos_mu, pos_shu, color="#3fb950")

            self._update_dock_button("o", label=f"Shu-Mu: {organ[:6]} ({mu_pt})",
                                    color="#15261e", border="#3fb950", text_color="#7ee787")
            _safe_print(f"\n╔{'═'*72}╗")
            _safe_print(f"║  ZASADA TCM: OŚ PRZÓD - TYŁ (PARY PRZEDNI-MU / TYLNY-SHU)")
            _safe_print(f"║  Narząd: {organ} [{ch}]")
            _safe_print(f"║  Punkt Przedni Mu (Alarmowy): {mu_pt} (brzuch/klatka - diagnoza ostrych stanów Yin)")
            _safe_print(f"║  Punkt Tylny Shu (Narządowy): {shu_pt} (plecy/Pęcherz - tonizacja narządu Yang)")
            _safe_print(f"║  Zastosowanie: {desc}")
            _safe_print(f"║  Kamera: Wciśnij [1] (przód - Mu) lub [2] (tył - Shu)")
            _safe_print(f"╚{'═'*72}╝\n")
            self.info(f"Shu-Mu: {organ} -> Mu: {mu_pt} (przód) | Shu: {shu_pt} (tył)")
            self.pl.render()
        else:
            _safe_print(f">>> Nie znaleziono pary Shu-Mu dla: {channel}. Dostępne: LU, LI, ST, SP, HT, SI, BL, KI, PC, TE, GB, LR")

    def cycle_xi_cleft(self) -> None:
        """Przelacza cyklicznie kolejne punkty Szczeliny Xi [i]."""
        xi_keys = list(XI_CLEFT_POINTS.keys())
        if not hasattr(self, "_xi_idx"):
            self._xi_idx = -1
        self._xi_idx = (self._xi_idx + 1) % len(xi_keys)
        x_code = xi_keys[self._xi_idx]
        self.jump_to_xi_cleft(x_code)

    def jump_to_xi_cleft(self, channel: str) -> None:
        ch = channel.upper().strip()
        target_code = None
        for code in XI_CLEFT_POINTS:
            if code.startswith(ch + "-") or code == ch:
                target_code = code
                break
        if target_code:
            desc = XI_CLEFT_POINTS[target_code]
            t_ch = target_code.split("-")[0]
            self.set_channel(t_ch, True, render=False)
            self.select_point(target_code)

            if hasattr(self, "highlight_actor_pair") and self.highlight_actor_pair is not None:
                self.highlight_actor_pair.SetVisibility(False)
            if hasattr(self, "relation_tube_actor") and self.relation_tube_actor is not None:
                try:
                    self.relation_tube_actor.SetVisibility(False)
                except Exception:
                    pass

            self._update_dock_button("i", label=f"Xi-Cleft: {target_code}",
                                    color="#2c1717", border="#f85149", text_color="#ff7b72")
            _safe_print(f"\n╔{'═'*72}╗")
            _safe_print(f"║  ZASADA TCM: PUNKT SZCZELINY XI (STANY NAGŁE I OSTRY BÓL)")
            _safe_print(f"║  Punkt Xi: {target_code}")
            _safe_print(f"║  Działanie kliniczne: {desc}")
            _safe_print(f"║  Mechanizm: Gwałtownie odblokowuje nagły zastój Qi i Krwi w kanale.")
            _safe_print(f"╚{'═'*72}╝\n")
            self.info(f"Punkt Xi: {target_code} ({desc})")
            self.pl.render()
        else:
            _safe_print(f">>> Nie znaleziono punktu Xi-Szczeliny dla: {channel}")

    def cycle_wu_xing(self) -> None:
        """Przelacza cyklicznie 5 Przemian (Drewno, Ogien, Ziemia, Metal, Woda) [w]."""
        order = list(ELEMENT_ORDER)
        if not hasattr(self, "_wuxing_idx"):
            self._wuxing_idx = -1
        self._wuxing_idx = (self._wuxing_idx + 1) % len(order)
        el = order[self._wuxing_idx]
        self.jump_to_wuxing(el)

    def jump_to_wuxing(self, element: str) -> None:
        el_map = {
            "DREWNO": "Wood", "DRZEWO": "Wood", "WOOD": "Wood",
            "OGIEN": "Fire", "OGIEŃ": "Fire", "FIRE": "Fire",
            "ZIEMIA": "Earth", "EARTH": "Earth",
            "METAL": "Metal",
            "WODA": "Water", "WATER": "Water"
        }
        raw = element.upper().strip()
        el_clean = el_map.get(raw, element.capitalize())
        if el_clean in ELEMENT_ORDER:
            chans = self._element_groups.get(el_clean, [])
            for ch in chans:
                self.set_channel(ch, True, render=False)
            name, yin, _ = ELEMENTS[el_clean]
            self._update_dock_button("w", label=f"Wu Xing: {name.capitalize()}",
                                    color="#1c2538", border="#38bdf8", text_color="#7dd3fc")
            _safe_print(f"\n╔{'═'*72}╗")
            _safe_print(f"║  ZASADA TCM: PIĘĆ PRZEMIAN (WU XING) - {name.upper()}")
            _safe_print(f"║  Powiązane kanały: {', '.join(chans)}")
            _safe_print(f"║  Włączono wizualizację meridianów i punktów tego żywiołu.")
            _safe_print(f"╚{'═'*72}╝\n")
            self.info(f"Wu Xing: Żywioł {name} ({', '.join(chans)})")
            self.pl.render()
        else:
            _safe_print(f">>> Nieznany żywioł Wu Xing: {element}. Dostępne: Drewno, Ogień, Ziemia, Metal, Woda.")

    def print_interactive_help(self) -> None:
        lines = [
            "",
            "╔══════════════════════════════════════════════════════════════════════════════╗",
            "║             KOMENDY TERMINALA I ZASADY MEDYCYNY CHIŃSKIEJ (TCM)              ║",
            "╠══════════════════════════════════════════════════════════════════════════════╣",
            "║  INTERAKTYWNE PRZYCISKI I SKRÓTY ZASAD TCM:                                  ║",
            "║    • [u] Yuan-Luo: CYKL        : przeskakuje przez 12 par Gospodarz-Gość     ║",
            "║    • [o] Shu-Mu: CYKL          : przeskakuje przez 12 osi Przód-Mu / Tył-Shu ║",
            "║    • [i] Xi-Cleft: CYKL        : przeskakuje przez 16 szczelin ostrego bólu  ║",
            "║    • [w] Wu Xing: CYKL         : cykl 5 Przemian (Drewno, Ogień, ...)        ║",
            "║    • [e] 8 Naczyń: CYKL        : 8 Naczyń Nadzwyczajnych (Du, Ren, Dai...)   ║",
            "║    • [j] Jing Jin: CYKL        : 12 taśm mięśniowo-powięziowych Mayersa      ║",
            "║    • [y] Biao-Li               : powłoki (Powierzchnia / Powięź / Wnętrze)   ║",
            "║    • [g] Siatka cun / [m] 12 Kanałów / [p] Punkty                            ║",
            "╟──────────────────────────────────────────────────────────────────────────────╢",
            "║  KOMENDY W KONSOLI (wpisz w oknie terminala i naciśnij Enter):               ║",
            "║    • yuan-luo [KANAL]          : np. 'yuan-luo LU' - para Gospodarz-Gość     ║",
            "║    • shumu [KANAL]             : np. 'shumu ST', 'shumu LU' - para Mu i Shu  ║",
            "║    • xicleft [KANAL]           : np. 'xicleft LU' - punkt ostrej blokady Xi  ║",
            "║    • wuxing [ŻYWIOŁ]           : np. 'wuxing drewno', 'wuxing woda'          ║",
            "║    • bamai [NAZWA]             : np. 'bamai dai', 'bamai chong', 'bamai qiao'║",
            "║    • punkt [KOD]               : np. 'punkt ST-36', 'punkt LU-7'             ║",
            "║    • rejon [REJON]             : np. 'rejon glowa', 'rejon plecy', 'klatka'  ║",
            "║    • biao-li [1|2|3]           : np. 'biao-li powiez', 'biao-li wnetrze'     ║",
            "║    • freeze / odmroz           : włącz / wyłącz blokadę obrotu do rysunku    ║",
            "║    • klisze                    : wyświetla listę zapisanych klisz            ║",
            "╚══════════════════════════════════════════════════════════════════════════════╝",
            "",
        ]
        _safe_print("\n".join(lines))

    def handle_console_command(self, cmd_line: str) -> None:
        cmd = cmd_line.strip()
        if not cmd:
            return
        parts = cmd.split()
        lead = parts[0].lower()

        if lead == "zasada" and len(parts) > 1:
            parts = parts[1:]
            lead = parts[0].lower()

        arg = parts[1] if len(parts) > 1 else ""

        if lead in ("pomoc", "help", "?"):
            self.print_interactive_help()
        elif lead in ("yuan-luo", "yuanluo", "yuan"):
            self.jump_to_yuan_luo(arg or "LU")
        elif lead in ("bamai", "ba-mai", "naczynie", "vessel"):
            self.set_extraordinary(arg or "DAI_MAI")
        elif lead in ("myers", "jingjin", "jing-jin"):
            self.set_jing_jin(arg or "SBL")
        elif lead in ("shumu", "shu-mu", "mu-shu"):
            self.jump_to_shu_mu(arg or "LU")
        elif lead in ("xicleft", "xi-cleft", "xi"):
            self.jump_to_xi_cleft(arg or "LU")
        elif lead in ("wuxing", "wu-xing", "element"):
            self.jump_to_wuxing(arg or "DREWNO")
        elif lead in ("punkt", "point", "p"):
            if arg:
                self.select_point(arg)
            else:
                _safe_print(">>> Podaj kod punktu, np. 'punkt ST-36'")
        elif lead in ("rejon", "region"):
            reg_map = {"glowa": "HEAD", "head": "HEAD", "klatka": "CHEST", "chest": "CHEST",
                       "plecy": "BACK", "back": "BACK", "rece": "ARMS", "arms": "ARMS",
                       "nogi": "LEGS", "legs": "LEGS", "calosc": "ALL", "all": "ALL"}
            reg_id = reg_map.get(arg.lower(), "ALL")
            self.set_region(reg_id)
        elif lead in ("biao-li", "biaoli", "warstwa"):
            lvl_map = {"skora": 1, "powierzchnia": 1, "biao": 1, "1": 1,
                       "powiez": 2, "jingjin": 2, "2": 2,
                       "wnetrze": 3, "szkielet": 3, "li": 3, "3": 3}
            lvl = lvl_map.get(arg.lower(), 1)
            self.set_biao_li(lvl)
        elif lead in ("kanal", "channel", "ch", "m", "meridian") or lead in (
            "bl", "pm", "ki", "lr", "gb", "ht", "si", "pc", "te", "sj", "sp", "st", "lu", "li", "du", "rn", "ren"
        ):
            target_raw = arg if lead in ("kanal", "channel", "ch", "m", "meridian") else lead
            if target_raw:
                if target_raw.lower() in ("all", "wszystkie", "glowne", "12"):
                    self.toggle_all_main_meridians()
                else:
                    ch_map = {
                        "PM": "BL", "PECHERZ": "BL", "BLADDER": "BL", "BL": "BL",
                        "N": "KI", "NERKI": "KI", "KIDNEY": "KI", "KI": "KI",
                        "W": "LR", "WATROBA": "LR", "WĄTROBA": "LR", "LIVER": "LR", "LR": "LR",
                        "PZ": "GB", "WZ": "GB", "GB": "GB", "GALLBLADDER": "GB",
                        "S": "HT", "SERCE": "HT", "HEART": "HT", "HT": "HT",
                        "JC": "SI", "JELITOCIENKIE": "SI", "SI": "SI",
                        "O": "PC", "OS": "PC", "OSIERDZIE": "PC", "PC": "PC",
                        "PO": "TE", "SJ": "TE", "TE": "TE", "SANJIAO": "TE",
                        "SL": "SP", "ŚL": "SP", "SLEDZIONA": "SP", "ŚLEDZIONA": "SP", "SP": "SP",
                        "Z": "ST", "Ż": "ST", "ZOLADEK": "ST", "ŻOŁĄDEK": "ST", "ST": "ST",
                        "P": "LU", "PL": "LU", "PLUCA": "LU", "PŁUCA": "LU", "LU": "LU",
                        "JG": "LI", "JELITOGRUBE": "LI", "LI": "LI",
                        "DU": "DU", "DUMAI": "DU", "GV": "DU",
                        "RN": "RN", "REN": "RN", "RENMAI": "RN", "CV": "RN",
                    }
                    target = ch_map.get(target_raw.upper(), target_raw.upper())
                    if target in self.atlas.channels:
                        self.set_channel(target, True)
                    else:
                        _safe_print(f">>> Nieznany kanal: '{target_raw}'. Dostepne: LU, LI, ST, SP, HT, SI, BL (PM), KI, PC, TE, GB, LR, DU, RN")
        elif lead in ("siatka", "cun", "grid"):
            self.toggle_cun_grid()
        elif lead in ("punkty", "points"):
            self.toggle_points()
        elif lead in ("glowne", "main"):
            self.toggle_all_main_meridians()
        elif lead in ("przod", "front"):
            self._view_front()
        elif lead in ("tyl", "tył"):
            self._view_back()
        elif lead in ("bok", "side"):
            self._view_right()
        elif lead in ("freeze", "zamroz"):
            self.toggle_camera_freeze(True)
        elif lead in ("odmroz", "unfreeze"):
            self.toggle_camera_freeze(False)
        elif lead in ("klisze", "klisza"):
            self.browse_klisze()
        elif lead in ("analizuj", "analiza"):
            self.analyze_current_strokes()
        elif lead in ("czysc", "clear"):
            self.clear_draw_strokes()
        else:
            if not self.select_point(cmd):
                _safe_print(f">>> Nieznane polecenie: '{cmd}'. Wpisz 'pomoc', aby zobaczyć dostępne zasady i komendy.")

    def _start_console_listener(self) -> None:
        if self.console_listener_started or self.pl.off_screen:
            return
        self.console_listener_started = True

        def _worker():
            time.sleep(0.6)
            _safe_print("\n" + "=" * 76)
            _safe_print("  INTERAKTYWNY TERMINAL KLINICZNY I ZASAD TCM")
            _safe_print("  Wpisz komendę bezpośrednio tutaj i wciśnij Enter:")
            _safe_print("    • 'pomoc'                 -> pełna lista komend i zasad TCM")
            _safe_print("    • 'zasada yuan-luo LU'    -> para Gospodarz-Gość kanału Płuc")
            _safe_print("    • 'zasada bamai dai'      -> aktywacja naczynia opasującego Dai Mai")
            _safe_print("    • 'zasada myers sbl'      -> Taśma Powierzchowna Tylna (Zu Taiyang)")
            _safe_print("    • 'zasada shumu ST'       -> para Przedni Mu / Tylny Shu Żołądka")
            _safe_print("    • 'punkt ST-36'           -> natychmiastowe zbliżenie do Zusanli")
            _safe_print("=" * 76 + "\n")
            while True:
                try:
                    line = sys.stdin.readline()
                    if not line:
                        break
                    cmd = line.strip()
                    if cmd:
                        self.handle_console_command(cmd)
                except Exception:
                    break

        t = threading.Thread(target=_worker, daemon=True, name="acu3d_console_listener")
        t.start()

    def _print_welcome_guide(self) -> None:
        """Drukuje czytelny przewodnik startowy w oknie konsoli / terminala."""
        lines = [
            "",
            "╔══════════════════════════════════════════════════════════════════════════════╗",
            "║             ATLAS AKUPUNKTURY 3D & SZKICOWNIK KLINICZNY (ATLAS-3D)           ║",
            "╠══════════════════════════════════════════════════════════════════════════════╣",
            "║  Porcelanowy posążek ciała z siatką cunową i bazą wiedzy Petera Deadmana     ║",
            "╟──────────────────────────────────────────────────────────────────────────────╢",
            "║  INTERAKCJA I STEROWANIE MYSZĄ:                                              ║",
            "║    • Lewy klik w punkt na ciele  -> Wybór punktu i opis Deadmana w konsoli  ║",
            "║    • Lewy klik w skórę (pomiędzy)-> Pomiar cunowy, segment i tkanka podskórna║",
            "║    • Przeciąganie myszą          -> Płynny obrót 3D modelu wokół osi         ║",
            "║    • Prawy panel (przyciski)     -> Kliknij myszą dowolny przycisk           ║",
            "║    • Rolka myszy                 -> Przybliżanie / oddalanie pod kursor      ║",
            "╟──────────────────────────────────────────────────────────────────────────────╢",
            "║  PRZYCISKI I SKRÓTY KLAWISZOWE:                                              ║",
            "║    [b] Rejon punktów    : Przełącz rejon (Głowa, Klatka, Plecy, Ręce, Nogi)  ║",
            "║    [5]..[0] Rejony      : [5] Całość, [6] Głowa, [7] Klatka, [8] Plecy,     ║",
            "║                           [9] Ręce, [0] Nogi                                 ║",
            "║    [m] Linie meridianów : Przełącz ciągły przebieg 3D meridianów             ║",
            "║    [e] 8 Naczyń (Ba Mai): Przełącz naczynia nadzwyczajne (Dai, Chong, Qiao..)║",
            "║    [j] Jing Jin (Myers) : Przełącz taśmy powięziowo-mięśniowe Myersa (SBL...)║",
            "║    [y] Biao-Li (Yin-Yang: Przełącz głębokość (Powierzchnia -> Powięź -> Wnętrze)║",
            "║    [d] Rysuj napięcia   : Włącz/wyłącz tryb rysowania pasm na porcelanie     ║",
            "║    [x] Zamroź kamerę    : FREEZE - blokada obrotu kamery do precyzji rysunku ║",
            "║    [z] Cofnij linię     : Cofa ostatnio narysowaną linię napięcia            ║",
            "║    [c] Wyczyść szkic    : Usuwa wszystkie narysowane linie                   ║",
            "║    [s] Zapisz kliszę    : Zapisuje narysowany wzorzec do pliku kliszy (JSON) ║",
            "║    [k] Przeglądaj klisze: Wczytuje i wyświetla zapisane klisze na modelu     ║",
            "║    [a] Analizuj         : Raport korelacji pasm z meridianami w konsoli      ║",
            "║    [g] Siatka cunowa    : Przełącz widoczność wektorowej siatki cunowej      ║",
            "║    [p] Wszystkie punkty : Włącz / ukryj wszystkie punkty akupunkturowe       ║",
            "║    [l] Etykiety kodów   : Włącz / wyłącz kody punktów                        ║",
            "║    [?] / [/] Pomoc      : Pełna lista reguł TCM i komend konsoli             ║",
            "║    [1]/[2]/[3]/[4]/[r]  : Widok Przód / Tył / Bok Prawy / Lewy / Całość      ║",
            "║    [Esc]                : Zamknij opis / odznacz wybrany punkt               ║",
            "╟──────────────────────────────────────────────────────────────────────────────╢",
            "║  KOMENDY INTERAKTYWNE TERMINALA:                                             ║",
            "║    Wpisz bezpośrednio w tym oknie terminala i naciśnij Enter:                ║",
            "║    • 'zasada yuan-luo LU'    -> skok do pary Gospodarz-Gość kanału Płuc      ║",
            "║    • 'zasada bamai dai'      -> aktywacja naczynia opasującego Dai Mai       ║",
            "║    • 'zasada myers sbl'      -> Taśma Powierzchowna Tylna (Zu Taiyang)       ║",
            "║    • 'zasada shumu ST'       -> para Przedni Mu / Tylny Shu Żołądka          ║",
            "║    • 'zasada wuxing drewno'  -> żywioł Drewna (Wątroba i Pęcherzyk)          ║",
            "║    • 'punkt ST-36'           -> natychmiastowe zbliżenie do Zusanli          ║",
            "║    • 'pomoc'                 -> pełna lista komend i zasad klinicznych       ║",
            "╚══════════════════════════════════════════════════════════════════════════════╝",
            "",
        ]
        _safe_print("\n".join(lines))

    # ------------------------------------------------------------------- start
    def show(self) -> None:
        # WSLg potrafi zaparkowac nowe okno poza ekranem - jawna pozycja temu zapobiega
        try:
            rw = self.pl.ren_win
            rw.SetPosition(60, 60)
            rw.SetSize(*self.pl.window_size)
            rw.SetDesiredUpdateRate(60.0)
        except Exception:
            pass
        self._interaction()
        self._print_welcome_guide()
        self._start_console_listener()
        self.pl.show()

    def screenshot(self, path: str) -> None:
        self.pl.screenshot(path)

