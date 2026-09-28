import os
import sys
from acu3d.model import Atlas
from acu3d.wire import Wire
import numpy as np

def main():
    atlas = Atlas()
    print("=== TEST ATLAS MERIDIAN LINES FOR BL, DU, RN ===")
    for ch in ("BL", "DU", "RN"):
        for side in ("R", "L"):
            lines = atlas.meridian_lines(ch, side)
            print(f"Channel {ch} side {side}: {len(lines)} paths, lens: {[len(p) for p in lines]}")
            for i, p in enumerate(lines):
                if len(p) > 0:
                    P = np.array(p)
                    # Y < -0.01 is back
                    back_pts = np.sum(P[:, 1] < -0.01)
                    print(f"  Path {i}: {len(p)} points, Y range: [{P[:,1].min():.3f}, {P[:,1].max():.3f}], Back points: {back_pts}/{len(p)}")
            if ch in ("DU", "RN"):
                break

if __name__ == "__main__":
    main()
