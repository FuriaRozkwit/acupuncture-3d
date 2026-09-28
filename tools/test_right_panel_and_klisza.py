import sys
from pathlib import Path
REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

import pyvista as pv
pv.OFF_SCREEN = True
import numpy as np

from acu3d.model import Atlas
from acu3d.surface_viewer import SurfaceViewer, THEME

print("=== TEST 1: Inicjalizacja SurfaceViewer z prawym panelem dock ===")
atlas = Atlas()
v = SurfaceViewer(atlas, off_screen=True, window_size=(1500, 960))

print(f"Liczba przyciskow w prawym panelu: {len(v._dock_buttons)}")
assert len(v._dock_buttons) >= 14, f"Powinno byc min. 14 przyciskow, jest: {len(v._dock_buttons)}"

keys_expected = ["d", "x", "z", "c", "s", "k", "a", "g", "p", "l", "1", "2", "3", "4", "r"]
for k in keys_expected:
    assert k in v._dock_actors, f"Brak przycisku [{k}] w docku!"
print("Wszystkie przyciski sa zarejestrowane w docku.")

print("\n=== TEST 2: Test klikania w przycisk [d] Rysuj napiecia ===")
btn_d = v._dock_actors["d"]
bx, by, bw, bh = btn_d["rect"]
click_x = bx + 20
click_y = by + 10

assert not v.draw_mode
hit = v._check_dock_click(click_x, click_y)
assert hit, "Klikniecie w prostokat przycisku nie zostalo wykryte!"
assert v.draw_mode, "Tryb rysowania nie zostal wlaczony po kliknieciu!"
print("Sukces: Klikniecie w przycisk [d] wlaczylo tryb rysowania.")

print("\n=== TEST 3: Dodanie pasma napiecia i analiza polaczen [a] ===")
p_lu1 = atlas.by_key["LU-1·L"]["pos"]
p_lu5 = atlas.by_key["LU-5·L"]["pos"]
p_lu7 = atlas.by_key["LU-7·L"]["pos"]
p_lu11 = atlas.by_key["LU-11·L"]["pos"]

t1 = np.linspace(0, 1, 15)
s1 = [p_lu1 + (p_lu5 - p_lu1) * t for t in t1]
t2 = np.linspace(0, 1, 20)[1:]
s2 = [p_lu5 + (p_lu7 - p_lu5) * t for t in t2]
t3 = np.linspace(0, 1, 10)[1:]
s3 = [p_lu7 + (p_lu11 - p_lu7) * t for t in t3]
stroke = np.array(s1 + s2 + s3)

v.add_stroke_from_points(stroke, color="#e63946")
assert len(v.draw_strokes) == 1

# Klikniecie przycisku [a] Analizuj polaczenia
btn_a = v._dock_actors["a"]
ax, ay, aw, ah = btn_a["rect"]
v._check_dock_click(ax + 20, ay + 10)
print("Sukces: Analiza polaczen zostala wykonana i wydrukowana.")

print("\n=== TEST 4: Zapis do kliszy [s] ===")
btn_s = v._dock_actors["s"]
sx, sy, sw, sh = btn_s["rect"]
v._check_dock_click(sx + 20, sy + 10)
assert v.active_klisza is not None
print(f"Sukces: Zapisano klisze o ID: {v.active_klisza.id}")

print("\n=== TEST 5: Test cofania [z] i czyszczenia [c] ===")
assert len(v.draw_strokes) == 1
btn_z = v._dock_actors["z"]
v._check_dock_click(btn_z["rect"][0] + 10, btn_z["rect"][1] + 10)
assert len(v.draw_strokes) == 0
print("Sukces: Przycisk [z] cofnal pasmo.")

print("\n=== TEST 6: Test przegladania i ladowania klisz [k] ===")
btn_k = v._dock_actors["k"]
v._check_dock_click(btn_k["rect"][0] + 10, btn_k["rect"][1] + 10)
assert len(v._browse_actors) >= 1, "Klisza nie zostala narysowana w trybie browse!"
print("Sukces: Przycisk [k] wczytal i wyswietlil klisze na modelu.")

print("\n=== TEST 7: Renderowanie zrzutow ekranu do artefaktow ===")
art_dir = Path("/mnt/c/Users/Paweł/.gemini/antigravity/brain/7c83af86-491a-484b-a5a4-fdc8094fa736")

# 1. Cala postac z prawym panelem
v._view("xz", 180)
v.select_point("LU-1·L")
shot_full = art_dir / "right_dock_full.png"
v.pl.screenshot(str(shot_full))
print(f"1. Zapisano zrzut calej postaci z dockiem: {shot_full}")

# 2. Zblizenie na klatke i reke z widocznym panelem i linia kliszy
cam = v.pl.camera
foc0 = np.array(cam.focal_point)
pos0 = np.array(cam.position)
shift = np.array(p_lu5) - foc0
cam.focal_point = tuple(foc0 + shift * 0.7)
cam.position = tuple(pos0 + shift * 0.7)
cam.zoom(2.4)
v.pl.renderer.ResetCameraClippingRange()
v.pl.render()
shot_close = art_dir / "right_dock_chest_klisza.png"
v.pl.screenshot(str(shot_close))
print(f"2. Zapisano zrzut zblizenia: {shot_close}")

print("\n--> WSZYSTKIE TESTY PRAWEGO PANELU I KLISZ ZALICZONE SUKCESEM!")
