import sys
from acu3d.model import Atlas
from acu3d.surface_viewer import SurfaceViewer

def main():
    atlas = Atlas()
    sv = SurfaceViewer(atlas, off_screen=True)
    
    # 1. Czysty porcelanowy ludzik (stan startowy)
    sv.screenshot("test_clean_porcelain.png")
    print("1. Saved test_clean_porcelain.png (pure porcelain, no points, no wire)")
    
    # 2. Włączenie PM (BL) - automatycznie włącza punkty i linie 3D!
    sv.set_channel("BL", True)
    sv._view_back()
    sv.screenshot("test_pm_selected_back.png")
    print("2. Saved test_pm_selected_back.png (BL lines + points on back)")
    
    # 3. Włączenie Du Mai (DU) obok BL - oba widoczne na plecach!
    sv.set_channel("DU", True)
    sv.screenshot("test_pm_and_du_back.png")
    print("3. Saved test_pm_and_du_back.png (both BL and DU active together on spine)")
    
    # 4. Włączenie siatki CUN
    sv.toggle_cun_grid()
    sv.screenshot("test_pm_du_and_cun.png")
    print("4. Saved test_pm_du_and_cun.png (BL + DU + CUN grid)")

if __name__ == "__main__":
    main()
