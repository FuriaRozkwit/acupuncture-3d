import sys
from pathlib import Path
import pyvista as pv
pv.OFF_SCREEN = True

p = pv.Plotter(off_screen=True, window_size=(1500, 960))
p.set_background("#121418", top="#2e3440")

W, H = 1500, 960
bx = W - 260
bw = 244

def hex_to_rgb(h):
    h = h.lstrip("#")
    return tuple(int(h[i:i+2], 16) / 255.0 for i in (0, 2, 4))

def make_btn(icon, name, sc, y, color="#1c2128", border="#485463", text_color="#f0f4f8"):
    # Consistent 24 character button line
    label = f" {icon}  {name:<14} [{sc:>1}] "
    act = p.add_text(label, position=(bx, y), font_size=9, color=text_color, font="courier")
    tp = act.prop
    r, g, b = hex_to_rgb(color)
    br, bg, bb = hex_to_rgb(border)
    tp.SetBackgroundColor(r, g, b)
    tp.SetBackgroundOpacity(0.92)
    tp.SetFrame(1)
    tp.SetFrameColor(br, bg, bb)
    tp.SetFrameWidth(1)
    return act

p.add_text("PANEL KLINICZNY", position=(bx, H - 42), font_size=11, color="#ffe36a")
p.add_text("── SZKICOWNIK & KLISZE ──", position=(bx, H - 72), font_size=8, color="#8e98a8")

make_btn("✏", "Rysuj napięcia", "d", H - 104, color="#2c1717", border="#e63946", text_color="#ff7b72")
make_btn("↺", "Cofnij linię",   "z", H - 136)
make_btn("🗑", "Wyczyść szkic",  "c", H - 168)
make_btn("💾", "Zapisz kliszę",  "s", H - 200, color="#15261e", border="#3fb950", text_color="#56d364")
make_btn("🎞", "Przeglądaj",     "k", H - 232)
make_btn("⚡", "Analizuj tor",   "a", H - 264, color="#2c2214", border="#d29922", text_color="#e3b341")

p.add_text("── WIDOK & ANATOMIA ──", position=(bx, H - 302), font_size=8, color="#8e98a8")
make_btn("🌐", "Siatka cun",     "g", H - 334, color="#21262d", border="#b37356", text_color="#d2a89b")
make_btn("●", "Kanały i punkty","p", H - 366)
make_btn("🏷", "Etykiety kodów", "l", H - 398)

p.add_text("── RZUTY KAMERY ──", position=(bx, H - 436), font_size=8, color="#8e98a8")
make_btn("👁", "Widok z przodu", "1", H - 468)
make_btn("👁", "Widok z tyłu",   "2", H - 500)
make_btn("👁", "Bok prawy",      "3", H - 532)
make_btn("👁", "Bok lewy",       "4", H - 564)
make_btn("🔄", "Reset kamery",   "r", H - 596)

p.screenshot("/tmp/test_dock_align.png")
print("Saved align screenshot")
