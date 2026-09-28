import sys
import os
sys.path.insert(0, "/mnt/c/Users/Paweł/repos/acupuncture-3d")

import numpy as np
import pyvista as pv
from acu3d.model import Atlas
from acu3d.surface_viewer import SurfaceViewer

print(">>> Initializing Atlas and SurfaceViewer off_screen...")
atlas = Atlas()
viewer = SurfaceViewer(atlas, off_screen=True)

# 1. Test initial camera position
cam_pos_0 = np.array(viewer.pl.camera.position)
print("Initial cam pos:", cam_pos_0)
print("Initial camera_frozen:", getattr(viewer, "camera_frozen", None))

# 2. Test toggle_camera_freeze
viewer.toggle_camera_freeze()
print("After toggle, camera_frozen:", viewer.camera_frozen)
assert viewer.camera_frozen is True, "Camera should be frozen!"

# 3. Test simulating mouse drag while frozen
iren = viewer.pl.iren.interactor
iren.SetEventPosition(500, 500)
iren.LeftButtonPressEvent()
iren.SetEventPosition(600, 600)
iren.MouseMoveEvent()
iren.LeftButtonReleaseEvent()

cam_pos_1 = np.array(viewer.pl.camera.position)
print("Cam pos after drag while frozen:", cam_pos_1)
assert np.allclose(cam_pos_0, cam_pos_1), "Camera MUST NOT move while frozen!"
print("PASSED: Camera did not rotate while frozen!")

# 4. Test draw_mode activation and automatic freeze
viewer.toggle_camera_freeze(False)
assert viewer.camera_frozen is False
viewer.toggle_draw_mode()
assert viewer.draw_mode is True, "Draw mode should be active!"
assert viewer.camera_frozen is True, "Entering draw mode MUST freeze camera!"
print("PASSED: Entering draw mode automatically freezes camera!")

# 5. Simulate drawing on the body:
# Find a 2D screen coordinate that hits the skin
W, H = viewer.pl.window_size
hit_coords = []
for test_y in range(300, 700, 20):
    p, n = viewer.pick(W // 2, test_y)
    if p is not None:
        hit_coords.append((W // 2, test_y))

print(f"Found {len(hit_coords)} valid points on torso/skin for drawing.")
assert len(hit_coords) >= 5, "Should find at least 5 points on skin"

# Start drag
start_x, start_y = hit_coords[0]
iren.SetEventPosition(start_x, start_y)
iren.LeftButtonPressEvent()

for pt_x, pt_y in hit_coords[1:6]:
    iren.SetEventPosition(pt_x, pt_y)
    iren.MouseMoveEvent()

end_x, end_y = hit_coords[5]
iren.SetEventPosition(end_x, end_y)
iren.LeftButtonReleaseEvent()

print(f"Draw strokes count: {len(viewer.draw_strokes)}")
assert len(viewer.draw_strokes) == 1, f"Expected 1 stroke drawn, got {len(viewer.draw_strokes)}"
print(f"Stroke #1 points: {len(viewer.draw_strokes[0])}")
assert len(viewer.draw_strokes[0]) >= 2, "Stroke should have >= 2 points"

# Camera still must NOT have moved during drawing
cam_pos_2 = np.array(viewer.pl.camera.position)
assert np.allclose(cam_pos_0, cam_pos_2), "Camera MUST NOT move during drawing!"
print("PASSED: Drawing created stroke without moving camera!")

print("\n>>> ALL TESTS PASSED SUCCESSFULLY! <<<")
