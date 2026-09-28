import sys
from acu3d.model import Atlas
from acu3d.surface_viewer import SurfaceViewer
import numpy as np

def main():
    atlas = Atlas()
    sv = SurfaceViewer(atlas, off_screen=True)
    
    print("=== TESTING RIGHT DOCK BUTTON CLICKS ===")
    print(f"Total buttons in _dock_buttons: {len(sv._dock_buttons)}")
    
    for btn in sv._dock_buttons:
        bx, by, bw, bh = btn["rect"]
        center_x = bx + bw // 2
        center_y = by + bh // 2
        
        # Test direct _check_dock_click at button center
        hit = sv._check_dock_click(center_x, center_y)
        print(f"Button [{btn['shortcut']}] '{btn['label']}': rect=({bx}, {by}, {bw}, {bh}), click at ({center_x}, {center_y}) -> hit={hit}")
        if not hit:
            print(f"  FAILED to hit button {btn['label']}!")

if __name__ == "__main__":
    main()
