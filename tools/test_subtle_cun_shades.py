import sys
from pathlib import Path
REPO = Path("/mnt/c/Users/Paweł/repos/acupuncture-3d")
sys.path.insert(0, str(REPO))

import pyvista as pv
pv.OFF_SCREEN = True
import numpy as np
from acu3d.model import Atlas
from acu3d.surface_viewer import SurfaceViewer, THEME

atlas = Atlas()

candidates = [
    ("baseline", "#b37356", 2.2, 0.95, "Obecny: Cielista sienna (2.2px, 0.95)"),
    ("cand_A_cashmere", "#c29b88", 1.4, 0.85, "Kandydat A: Jasny kaszmirowy cielisty (#c29b88, 1.4px, 0.85)"),
    ("cand_B_biscuit", "#b8977e", 1.4, 0.85, "Kandydat B: Cieply piaskowy biszkopt (#b8977e, 1.4px, 0.85)"),
    ("cand_C_sand_stone", "#baa38b", 1.3, 0.82, "Kandydat C: Miekki piaskowy taupe (#baa38b, 1.3px, 0.82)"),
    ("cand_D_warm_nude", "#c49a82", 1.4, 0.85, "Kandydat D: Lagodny cieply nude (#c49a82, 1.4px, 0.85)"),
    ("cand_E_soft_clay", "#ad8c75", 1.4, 0.85, "Kandydat E: Stonowana glinka (#ad8c75, 1.4px, 0.85)"),
]

art_dir = Path("/mnt/c/Users/Paweł/.gemini/antigravity/brain/7c83af86-491a-484b-a5a4-fdc8094fa736")

for name, color, width, opac, desc in candidates:
    THEME["cun_grid"] = color
    v = SurfaceViewer(atlas, off_screen=True, window_size=(1400, 960))
    # Adjust wire actor properties directly
    v.wire_actor.prop.SetColor(*v._hex_to_rgb(color))
    v.wire_actor.prop.SetLineWidth(width)
    v.wire_actor.prop.SetOpacity(opac)
    v.pl.render()
    
    # Save full postac
    out_full = art_dir / f"cun_shade_{name}_full.png"
    v.screenshot(str(out_full))
    
    # Zoom on chest & arm
    cam = v.pl.camera
    foc = np.array([0.0, 0.1, 0.85])
    cam.focal_point = tuple(foc)
    cam.position = tuple(foc + np.array([0.0, 1.4, 0.1]))
    v.pl.renderer.ResetCameraClippingRange()
    v.pl.render()
    out_zoom = art_dir / f"cun_shade_{name}_zoom.png"
    v.screenshot(str(out_zoom))
    print(f"Rendered {name}: {desc}")

print("All candidates rendered!")
