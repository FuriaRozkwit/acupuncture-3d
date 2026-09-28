import sys
from pathlib import Path
REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

import pyvista as pv
pv.OFF_SCREEN = True
import numpy as np

from acu3d.model import Atlas
from acu3d.surface_viewer import SurfaceViewer, THEME

print("=== TEST 1: Inicjalizacja SurfaceViewer (Porcelana + Cielista siatka) ===")
atlas = Atlas()
v = SurfaceViewer(atlas, off_screen=True, window_size=(1600, 1100))

print("THEME skin:", THEME["skin"])
print("THEME cun_grid:", THEME["cun_grid"])
assert THEME["skin"] == "#ece6dc", "Kolor skory powinien byc porcelanowy (#ece6dc)!"
assert THEME["cun_grid"] == "#b37356", "Siatka powinna byc cielista (#b37356)!"

print("=== TEST 2: Wlaczenie trybu rysowania [d] ===")
assert not v.draw_mode
v.toggle_draw_mode()
assert v.draw_mode, "Tryb rysowania nie zostal wlaczony!"
print("Tryb rysowania aktywny:", v.draw_mode)

print("=== TEST 3: Rysowanie anatomicznych linii napiecia powiezii ===")
lu1 = atlas.by_key["LU-1·L"]["pos"]
lu5 = atlas.by_key["LU-5·L"]["pos"]
lu7 = atlas.by_key["LU-7·L"]["pos"]
lu11 = atlas.by_key["LU-11·L"]["pos"]
ren17 = atlas.by_key["REN-17"]["pos"]

# 1. Linia ramienna (Arm Line - LU): od zgiecia lokciowego do nadgarstka i kciuka
t1 = np.linspace(0, 1, 15)
stroke_forearm = [lu5 + (lu7 - lu5) * t for t in t1]
t2 = np.linspace(0, 1, 15)[1:]
stroke_hand = [lu7 + (lu11 - lu7) * t for t in t2]
arm_line = np.array(stroke_forearm + stroke_hand)
v.add_stroke_from_points(arm_line, color="#e63946")

# 2. Linia piersiowa powiezii: od mostka (REN-17) do kata pachowego i dolu podobojczykowego (LU-1)
chest_line = np.array([ren17 + (lu1 - ren17) * t for t in np.linspace(0, 1, 20)])
v.add_stroke_from_points(chest_line, color="#e63946")

art_dir = Path("/mnt/c/Users/Paweł/.gemini/antigravity/brain/7c83af86-491a-484b-a5a4-fdc8094fa736")

# 1. Zrzut calej postaci z przodu
v.select_point("ST-36·L")
v._view("xz", 180)
v.pl.renderer.ResetCameraClippingRange()
shot_full = art_dir / "porcelain_cielisty_full.png"
v.pl.screenshot(str(shot_full))
print(f"1. Zapisano zrzut calej postaci: {shot_full}")

# 2. Zblizenie na klatke piersiowa i LU-1
v.select_point("LU-1·L")
v._view("xz", 180)
cam = v.pl.camera
foc0 = np.array(cam.focal_point)
pos0 = np.array(cam.position)
shift = np.array(lu1) - foc0
cam.focal_point = tuple(foc0 + shift)
cam.position = tuple(pos0 + shift)
cam.zoom(3.5)
v.pl.renderer.ResetCameraClippingRange()
v.pl.render()
shot_chest = art_dir / "porcelain_cielisty_chest_draw.png"
v.pl.screenshot(str(shot_chest))
print(f"2. Zapisano zrzut klatki piersiowej z linia napiecia: {shot_chest}")

# 3. Zblizenie na przedramie i LU-7
v.select_point("LU-7·L")
v._view("xz", 180)
cam = v.pl.camera
foc0 = np.array(cam.focal_point)
pos0 = np.array(cam.position)
shift = np.array(lu7) - foc0
cam.focal_point = tuple(foc0 + shift)
cam.position = tuple(pos0 + shift)
cam.zoom(4.2)
v.pl.renderer.ResetCameraClippingRange()
v.pl.render()
shot_arm = art_dir / "porcelain_cielisty_arm_draw.png"
v.pl.screenshot(str(shot_arm))
print(f"3. Zapisano zrzut przedramienia z linia napiecia: {shot_arm}")

print("\n--> WSZYSTKIE TESTY ZALICZONE SUKCESEM!")
