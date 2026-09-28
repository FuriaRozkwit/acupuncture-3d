import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from acu3d.model import Atlas
from acu3d.surface_viewer import SurfaceViewer

print("1. Inicjalizacja Atlas...")
atlas = Atlas()

print("2. Inicjalizacja SurfaceViewer(off_screen=False)...")
try:
    sv = SurfaceViewer(atlas, off_screen=False)
    print("   SurfaceViewer zainicjalizowany pomyslnie.")
    
    print("3. Ustawienie kamery i render...")
    sv.pl.render()
    
    print("4. Zrzut ekranu...")
    sv.pl.screenshot("test_onscreen_shot.png")
    print("   Zrzut test_onscreen_shot.png zapisany.")
    
    print("5. Test przelaczenia kanalu BL...")
    sv.set_channel("BL", True)
    sv._view_back()
    sv.pl.render()
    sv.pl.screenshot("test_onscreen_bl_back.png")
    print("   Zrzut test_onscreen_bl_back.png zapisany.")
    
    print("6. Test przelaczenia Yuan-Luo...")
    sv.cycle_yuan_luo()
    sv.pl.render()
    sv.pl.screenshot("test_onscreen_yuan_luo.png")
    print("   Zrzut test_onscreen_yuan_luo.png zapisany.")

    print("7. Zamykanie...")
    sv.pl.close()
    print("=== ONSCREEN TEST PASSED! ===")
except Exception as e:
    import traceback
    print("BŁĄD W ONSCREEN:")
    traceback.print_exc()
    sys.exit(1)
