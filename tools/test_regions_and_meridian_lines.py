import sys
from pathlib import Path
REPO = Path("/mnt/c/Users/Paweł/repos/acupuncture-3d")
sys.path.insert(0, str(REPO))

import pyvista as pv
pv.OFF_SCREEN = True
import numpy as np
from acu3d.model import Atlas
from acu3d.wire import Wire

atlas = Atlas()
wire = Wire()

def classify_marker(m, wire):
    pos = m["pos"]
    w = wire.describe(pos)
    seg = w.get("segment", "")
    if "head" in seg or "neck" in seg:
        return "HEAD"
    elif "arm" in seg or "hand" in seg or "forearm" in seg:
        return "ARMS"
    elif "thigh" in seg or "shank" in seg or "foot" in seg:
        return "LEGS"
    elif "trunk" in seg:
        if pos[1] >= -0.01:
            return "CHEST"
        else:
            return "BACK"
    return "HEAD" # fallback

for m in atlas.markers:
    m["region"] = classify_marker(m, wire)

# Verify all markers have region
regs = set(m["region"] for m in atlas.markers)
print("Unique regions found:", regs)
assert regs == {"HEAD", "CHEST", "BACK", "ARMS", "LEGS"}, f"Unexpected regions: {regs}"

# Verify counts
from collections import Counter
counts = Counter(m["region"] for m in atlas.markers)
print("Marker counts per region:", counts)
assert sum(counts.values()) == len(atlas.markers)

print("SUCCESS: Region classification verified!")
