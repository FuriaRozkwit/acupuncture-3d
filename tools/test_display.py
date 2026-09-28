import os
import sys
from pathlib import Path

# Fix WSLg sockets
if not os.path.exists("/run/user/0"):
    os.makedirs("/run/user/0", mode=0o700, exist_ok=True)
if not os.path.exists("/run/user/0/wayland-0") and os.path.exists("/mnt/wslg/runtime-dir/wayland-0"):
    try:
        os.symlink("/mnt/wslg/runtime-dir/wayland-0", "/run/user/0/wayland-0")
    except Exception as e:
        print("wayland symlink:", e)

if not os.path.exists("/tmp/.X11-unix") and os.path.exists("/mnt/wslg/.X11-unix"):
    try:
        os.symlink("/mnt/wslg/.X11-unix", "/tmp/.X11-unix")
    except Exception as e:
        print("x11 symlink:", e)

os.environ["DISPLAY"] = os.environ.get("DISPLAY", ":0")
os.environ["WAYLAND_DISPLAY"] = os.environ.get("WAYLAND_DISPLAY", "wayland-0")
os.environ["XDG_RUNTIME_DIR"] = os.environ.get("XDG_RUNTIME_DIR", "/run/user/0")

print(f"DISPLAY={os.environ['DISPLAY']}")
print(f"WAYLAND_DISPLAY={os.environ['WAYLAND_DISPLAY']}")
print(f"XDG_RUNTIME_DIR={os.environ['XDG_RUNTIME_DIR']}")

import pyvista as pv
print("PyVista version:", pv.__version__)

try:
    p = pv.Plotter(off_screen=False)
    p.add_text("Test Display", position='upper_left')
    print("Supports OpenGL:", p.ren_win.SupportsOpenGL())
    # show briefly
    p.show(auto_close=True, interactive=False)
    print("Window opened and closed successfully!")
except Exception as exc:
    print("DISPLAY TEST ERROR:", exc)
    sys.exit(1)
