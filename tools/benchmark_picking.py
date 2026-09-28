import time
import sys
from pathlib import Path
import numpy as np
import vtk

HERE = Path(__file__).resolve().parent.parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from acu3d.model import Atlas
from acu3d.surface_viewer import SurfaceViewer

atlas = Atlas()
sv = SurfaceViewer(atlas, off_screen=True)

# Test 1: czas pick(x, y)
t0 = time.perf_counter()
for _ in range(50):
    sv.pick(400, 300)
t_pick = (time.perf_counter() - t0) / 50 * 1000
print(f"1. Czas jednego sv.pick(x, y): {t_pick:.2f} ms")

# Test 2: czas place_cursor
p, n = sv.pick(400, 300)
t0 = time.perf_counter()
for _ in range(50):
    sv.place_cursor(p, n)
t_cursor = (time.perf_counter() - t0) / 50 * 1000
print(f"2. Czas jednego sv.place_cursor: {t_cursor:.2f} ms")

# Test 3: czas starej petli WorldToDisplay (675 punktow)
t0 = time.perf_counter()
for _ in range(10):
    for m in sv.atlas.markers:
        sv.pl.renderer.SetWorldPoint(float(m["pos"][0]), float(m["pos"][1]), float(m["pos"][2]), 1.0)
        sv.pl.renderer.WorldToDisplay()
        d = sv.pl.renderer.GetDisplayPoint()
t_old_loop = (time.perf_counter() - t0) / 10 * 1000
print(f"3. Czas starej petli WorldToDisplay: {t_old_loop:.2f} ms")

# Test 4: czas wektorowego rzutowania macierza kamery (675 punktow w numpy)
cam = sv.pl.camera
w, h = sv.pl.window_size
aspect = w / h
mat = cam.GetCompositeProjectionTransformMatrix(aspect, -1, 1)
# pobierz 16 elementow macierzy 4x4
M = np.zeros((4, 4), dtype=np.float64)
for r in range(4):
    for c in range(4):
        M[r, c] = mat.GetElement(r, c)

t0 = time.perf_counter()
for _ in range(50):
    pts4 = np.hstack([sv._mpos, np.ones((len(sv._mpos), 1), dtype=np.float64)])
    clip = pts4 @ M.T
    w_clip = clip[:, 3:4]
    w_clip[w_clip == 0] = 1e-9
    ndc = clip[:, :3] / w_clip
    disp_x = (ndc[:, 0] + 1.0) * 0.5 * w
    disp_y = (ndc[:, 1] + 1.0) * 0.5 * h
t_vector = (time.perf_counter() - t0) / 50 * 1000
print(f"4. Czas wektorowego rzutowania (NumPy): {t_vector:.3f} ms")
print(f"   -> Przyspieszenie: {t_old_loop / t_vector:.1f}x!")
