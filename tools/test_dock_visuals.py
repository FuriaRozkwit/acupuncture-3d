import sys
from pathlib import Path
REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

import pyvista as pv
pv.OFF_SCREEN = True
import numpy as np

p = pv.Plotter(off_screen=True, window_size=(1500, 960))
p.set_background("#121418", top="#2e3440")

# Draw a sphere as mockup of center body
sphere = pv.Sphere(radius=0.3)
p.add_mesh(sphere, color="#ece6dc")

W, H = 1500, 960
bx = W - 255
bw = 240
bh = 26

def make_btn(label, y, color="#21262d", border="#ffe36a", text_color="#f0f4f8"):
    act = p.add_text(f" {label} ", position=(bx, y), font_size=9, color=text_color)
    tp = act.prop
    # parse hex to rgb 0-1
    def hex_to_rgb(h):
        h = h.lstrip("#")
        return tuple(int(h[i:i+2], 16) / 255.0 for i in (0, 2, 4))
    r, g, b = hex_to_rgb(color)
    br, bg, bb = hex_to_rgb(border)
    tp.SetBackgroundColor(r, g, b)
    tp.SetBackgroundOpacity(0.85)
    tp.SetFrame(1)
    tp.SetFrameColor(br, bg, bb)
    tp.SetFrameWidth(1)
    return act

# Header
p.add_text("PANEL KLINICZNY", position=(bx, H - 42), font_size=11, color="#ffe36a")

# Section 1
p.add_text("── SZKICOWNIK & KLISZE ──", position=(bx, H - 74), font_size=8, color="#8e98a8")
make_btn("✏  Rysuj napięcia   [d]", H - 108, color="#2a1b1b", border="#e63946", text_color="#ff7b72")
make_btn("↺  Cofnij linię    [z]", H - 142)
make_btn("🗑  Wyczyść szkic   [c]", H - 176)
make_btn("💾  Zapisz kliszę   [s]", H - 210, color="#1c2d27", border="#4cb574", text_color="#7ee787")
make_btn("🎞  Klisze / Historia [k]", H - 244)
make_btn("⚡  Analizuj tor    [a]", H - 278, color="#2d261c", border="#dcab3f", text_color="#f2cc60")

# Section 2
p.add_text("── WIDOK & ANATOMIA ──", position=(bx, H - 316), font_size=8, color="#8e98a8")
make_btn("🌐  Siatka cunowa   [g]", H - 350)
make_btn("●   Kanały i punkty [p]", H - 384)
make_btn("🏷  Etykiety kodów  [l]", H - 418)

# Section 3
p.add_text("── RZUTY KAMERY ──", position=(bx, H - 456), font_size=8, color="#8e98a8")
make_btn("👁  Przód          [1]", H - 490)
make_btn("👁  Tył            [2]", H - 524)
make_btn("👁  Bok prawy / lewy [3/4]", H - 558)
make_btn("🔄  Cała postać    [r]", H - 592)

p.screenshot("/tmp/test_dock_visuals.png")
print("Screenshot saved to /tmp/test_dock_visuals.png successfully!")
