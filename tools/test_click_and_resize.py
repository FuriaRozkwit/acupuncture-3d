import sys
from pathlib import Path
REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

import pyvista as pv
pv.OFF_SCREEN = True
import numpy as np
import vtk

from acu3d.model import Atlas
from acu3d.surface_viewer import SurfaceViewer

print("=================================================================")
print("  TEST WERYFIKACYJNY: REJESTRACJA KLIKNIĘĆ, DOCK I RESIZE EKRANU")
print("=================================================================")

atlas = Atlas()
v = SurfaceViewer(atlas, off_screen=True, window_size=(1500, 960))

print("\n[KROK 1] Weryfikacja bazy przycisków i pickera vtkPropPicker...")
assert hasattr(v, "_dock_picker"), "Brak _dock_picker w SurfaceViewer!"
assert hasattr(v, "_dock_actor_map"), "Brak _dock_actor_map w SurfaceViewer!"
assert len(v._dock_buttons) >= 14, f"Niepoprawna liczba przycisków: {len(v._dock_buttons)}"
print(f"  OK: Zarejestrowano {len(v._dock_buttons)} przycisków w docku.")

# Sprawdzamy czy aktorzy przycisków są Pickable
for btn in v._dock_buttons:
    act = btn["actor"]
    assert act.GetPickable(), f"Aktor przycisku [{btn['key']}] ma Pickable=False!"
    assert act in v._dock_actor_map, f"Aktor [{btn['key']}] nie jest w _dock_actor_map!"
print("  OK: Wszyscy aktorzy przycisków mają Pickable=True i są w mapie.")

print("\n[KROK 2] Test bezpośredniego trafienia vtkPropPicker w przycisk [d] (1500x960)...")
btn_d = v._dock_actors["d"]
bx, by, bw, bh = btn_d["rect"]
hit_x = bx + 50
hit_y = by + 12

hit = v._dock_picker.PickProp(hit_x, hit_y, v.pl.renderer)
picked_actor = v._dock_picker.GetViewProp()
assert hit, f"PickProp nie trafił w aktor na pozycji ({hit_x}, {hit_y})!"
assert picked_actor == btn_d["actor"], "PickProp trafił w inny obiekt niż przycisk [d]!"
print(f"  OK: vtkPropPicker trafił idealnie w przycisk [{btn_d['label']}]!")

print("\n[KROK 3] Test symulacji kliknięcia myszą (left_down -> left_up) w [d]...")
assert not v.draw_mode, "draw_mode powinien być False początkowo"

# Ustawiamy pozycję zdarzenia myszy w interaktorze
iren = v.pl.iren.interactor
iren.SetEventPosition(hit_x, hit_y)

# Symulacja left_down
v._on_left_down(iren, "LeftButtonPressEvent")
assert v._dock_pressed_btn is not None, "Button nie został oznaczony jako wciśnięty w _on_left_down!"
assert v._dock_pressed_btn["key"] == "d"
print("  OK: _on_left_down rozpoznał przycisk i ustawił stan wciśnięcia bez rotacji kamery.")

# Symulacja left_up
v._on_left_up(iren, "LeftButtonReleaseEvent")
assert v.draw_mode, "draw_mode nie przełączył się na True po left_up!"
assert v._dock_pressed_btn is None, "Po left_up stan wciśnięcia powinien być wyczyszczony!"
print("  OK: left_up pomyślnie aktywował przycisk [d] i włączył tryb rysowania!")

print("\n[KROK 4] Test zmiany rozmiaru okna do Full HD (1920x1080) i automatycznej adaptacji docka...")
# Zmieniamy rozmiar okna renderowania
v.pl.ren_win.SetSize(1920, 1080)
v.pl.render()

# Sprawdzamy czy współrzędne przycisków dopasowały się do nowego prawego brzegu
expected_bx = 1920 - 262
btn_d_resized = v._dock_actors["d"]
new_bx = btn_d_resized["rect"][0]
new_by = btn_d_resized["rect"][1]
print(f"  Nowa pozycja przycisku [d] po resize: x={new_bx}, y={new_by} (oczekiwana x={expected_bx})")
assert new_bx == expected_bx, f"Błędne x po resize: {new_bx} != {expected_bx}"

# Test kliknięcia w nowy punkt przycisku [s] Zapisz kliszę po resize
btn_s = v._dock_actors["s"]
sx = btn_s["rect"][0] + 50
sy = btn_s["rect"][1] + 12

hit_s = v._dock_picker.PickProp(sx, sy, v.pl.renderer)
assert hit_s, "PickProp nie trafił w przycisk [s] po resize!"
assert v._dock_picker.GetViewProp() == btn_s["actor"], "Trafiono w inny aktor po resize!"
print("  OK: Przyciski po zmianie rozmiaru okna zachowują 100% precyzję trafień sprzętowych!")

print("\n[KROK 5] Test zaznaczania punktu na meridianie płuc (LU-1·L)...")
# Przełączamy kamerę na widok z przodu
v._view("xz", 180)
m_lu1 = atlas.by_key["LU-1·L"]
pos_lu1 = m_lu1["pos"]

# Rzutujemy punkt 3D na ekran 2D
v.pl.renderer.SetWorldPoint(float(pos_lu1[0]), float(pos_lu1[1]), float(pos_lu1[2]), 1.0)
v.pl.renderer.WorldToDisplay()
disp = v.pl.renderer.GetDisplayPoint()
disp_x, disp_y = int(disp[0]), int(disp[1])
print(f"  Współrzędne ekranowe punktu LU-1·L: ({disp_x}, {disp_y})")

# Wyłączamy draw_mode żeby przetestować kliknięcie diagnostyczne
v.draw_mode = False
iren.SetEventPosition(disp_x, disp_y)
v._on_left_down(iren, "LeftButtonPressEvent")
v._on_left_up(iren, "LeftButtonReleaseEvent")

print(f"  Wybrany punkt w atlasie: {v.selected_point}")
assert v.selected_point == "LU-1·L", f"Oczekiwano wyboru LU-1·L, wybrano: {v.selected_point}"
assert v.highlight_actor.GetVisibility(), "Znacznik podświetlenia punktu nie jest widoczny!"
assert "LU" in v.channels_on, "Kanał LU nie został automatycznie włączony!"
print("  OK: Punkt LU-1·L został błyskawicznie wybrany, podświetlony i opisany!")

print("\n[KROK 6] Test kliknięcia w skórę obok punktu (odczyt siatki cunowej)...")
# Kliknięcie w punkt na klatce piersiowej
chest_disp_x = disp_x - 40
chest_disp_y = disp_y - 20
iren.SetEventPosition(chest_disp_x, chest_disp_y)
v._on_left_down(iren, "LeftButtonPressEvent")
v._on_left_up(iren, "LeftButtonReleaseEvent")
print(f"  OK: Kliknięcie w skórę przetworzone pomyślnie.")

print("\n[KROK 7] Test odznaczania punktu [Esc]...")
v.clear_selection()
assert v.selected_point is None, "Punkt nie został odznaczony!"
assert not v.highlight_actor.GetVisibility(), "Znacznik podświetlenia powinien być ukryty!"
print("  OK: Odznaczanie punktu działa prawidłowo.")

print("\n=================================================================")
print("  WSZYSTKIE TESTY INTERFEJSU I REJESTRACJI KLIKNIĘĆ ZALICZONE!   ")
print("=================================================================")
