import sys
from pathlib import Path
REPO = Path("/mnt/c/Users/Paweł/repos/acupuncture-3d")
sys.path.insert(0, str(REPO))

import pyvista as pv
pv.OFF_SCREEN = True
import numpy as np

from acu3d.model import Atlas
from acu3d.surface_viewer import SurfaceViewer, ELEMENTS, REGIONS

print(">>> [TEST] Ładowanie atlasu i inicjalizacja SurfaceViewer...")
atlas = Atlas()
viewer = SurfaceViewer(atlas, off_screen=True, window_size=(1500, 960))

# 1. Sprawdzenie palety Wu Xing
print("\n1. Sprawdzanie palety Wu Xing...")
assert ELEMENTS["Wood"][1] == "#00a854", f"Błędny odcień drewna yin: {ELEMENTS['Wood']}"
assert ELEMENTS["Fire"][1] == "#e02424", f"Błędny odcień ognia yin: {ELEMENTS['Fire']}"
assert ELEMENTS["Earth"][1] == "#d97706", f"Błędny odcień ziemi yin: {ELEMENTS['Earth']}"
assert ELEMENTS["Metal"][1] == "#0284c7", f"Błędny odcień metalu yin: {ELEMENTS['Metal']}"
assert ELEMENTS["Water"][1] == "#1d4ed8", f"Błędny odcień wody yin: {ELEMENTS['Water']}"
assert ELEMENTS["-"][1] == "#7c3aed", f"Błędny odcień extra yin: {ELEMENTS['-']}"
print("   Paleta Wu Xing zweryfikowana pomyślnie!")

# 2. Sprawdzenie podziału aktorów na (ch, reg)
print("\n2. Sprawdzanie aktorów punktów (ch, reg)...")
assert len(viewer.point_actors) > 25, f"Zbyt mało aktorów punktów: {len(viewer.point_actors)}"
total_pts = sum(m["region"] in REGIONS for m in viewer.atlas.markers)
assert total_pts == 675, f"Oczekiwano 675 punktów z przypisanym rejonem, otrzymano: {total_pts}"
print(f"   Znaleziono {len(viewer.point_actors)} grup aktorów dla 675 punktów.")

ART_DIR = Path("/mnt/c/Users/Paweł/.gemini/antigravity/brain/7c83af86-491a-484b-a5a4-fdc8094fa736")

# 3. Test widoku calej postaci z wyrazistymi punktami Wu Xing
viewer.set_region("ALL")
viewer.pl.render()
viewer.pl.screenshot(str(ART_DIR / "wuxing_points_all_body.png"))
print("   Zapisano: wuxing_points_all_body.png")

# 4. Testy selekcji rejonami
for reg in ("HEAD", "CHEST", "BACK", "ARMS", "LEGS"):
    viewer.set_region(reg)
    viewer.pl.render()
    
    # Weryfikacja widoczności aktorów
    for (ch, r), act in viewer.point_actors.items():
        expected = (r == reg)
        actual = bool(act.GetVisibility())
        assert actual == expected, f"Błąd widoczności dla ({ch}, {r}) w rejonie {reg}: oczekiwano {expected}, jest {actual}"
    
    shot_name = f"region_{reg.lower()}_points.png"
    viewer.pl.screenshot(str(ART_DIR / shot_name))
    print(f"   Rejon {reg}: widoczność zweryfikowana. Zapisano: {shot_name}")

# 5. Test włącznika ciągłych linii meridianów [m]
print("\n5. Sprawdzanie włącznika linii meridianów [m]...")
viewer.set_region("ALL")
viewer.toggle_meridian_lines(True)
assert viewer.show_meridian_lines is True
visible_tubes = 0
for ch in atlas.channels:
    tube = viewer._get_or_create_tubes(ch)
    if tube and tube.GetVisibility():
        visible_tubes += 1
assert visible_tubes == 14, f"Oczekiwano 14 widocznych rurek meridianów, jest: {visible_tubes}"
viewer.pl.render()
viewer.pl.screenshot(str(ART_DIR / "meridian_lines_wuxing_active.png"))
print("   Włącznik meridianów [m] ON: 14 meridianów aktywnych. Zapisano: meridian_lines_wuxing_active.png")

# Zbliżenie na klatkę i twarz z meridianami i punktami perłowymi
viewer.pl.camera.focal_point = (0.0, 0.0, 1.25)
viewer.pl.camera.position = (0.0, -0.75, 1.25)
viewer.pl.renderer.ResetCameraClippingRange()
viewer.pl.render()
viewer.pl.screenshot(str(ART_DIR / "wuxing_meridians_and_pearls_closeup.png"))
print("   Zapisano zbliżenie: wuxing_meridians_and_pearls_closeup.png")

# Wyłączenie linii meridianów
viewer.toggle_meridian_lines(False)
assert viewer.show_meridian_lines is False
for ch in atlas.channels:
    tube = viewer._get_or_create_tubes(ch)
    if tube:
        assert not tube.GetVisibility(), f"Rurka meridianu {ch} powinna być ukryta!"
print("   Włącznik meridianów [m] OFF: rurki ukryte poprawnie.")

# 6. Test klikania przycisków w docku
print("\n6. Sprawdzanie przycisków w docku...")
m_btn = viewer._dock_actors.get("m")
assert m_btn is not None, "Brak przycisku 'm' w docku!"
m_btn["action"]() # toggle on
assert viewer.show_meridian_lines is True
m_btn["action"]() # toggle off
assert viewer.show_meridian_lines is False

b_btn = viewer._dock_actors.get("b")
assert b_btn is not None, "Brak przycisku 'b' w docku!"
init_reg = viewer.active_region
b_btn["action"]() # cycle
assert viewer.active_region != init_reg, "Przycisk 'b' nie zmienił rejonu!"
print(f"   Przycisk 'b' zmienił rejon z {init_reg} na {viewer.active_region}")

print("\n>>> WSZYSTKIE TESTY ZAKOŃCZONE PEŁNYM SUKCESEM! <<<")
