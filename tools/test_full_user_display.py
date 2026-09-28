"""Kompleksowy test uzytkownika weryfikujacy wyswietlanie wszystkich struktur i zasad TCM:
1. Czysty porcelanowy posazek (przod i tyl)
2. Siatka cunowa proporcjonalna Fēn-Cùn (przod i tyl)
3. Punkty akupunkturowe (widocznosc perel)
4. Przebieg grzbietowy: kanal Pecherza (BL) i Du Mai (DU) na plecach (brak znikania)
5. 12 glownych meridianow Zang-Fu
6. Zasada Yuan-Luo (Gospodarz-Gosc z przestrzenna tuba relacji 3D)
7. Zasada Shu-Mu (Przedni Mu i Tylny Shu)
8. Punkty szczelin Xi-Cleft
9. Piec Przemian Wu Xing
10. Szkicownik powiezio-napieciowy i analiza polaczen
"""
import sys
import os
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from acu3d.model import Atlas
from acu3d.surface_viewer import SurfaceViewer

ARTIFACT_DIR = Path("/mnt/c/Users/Paweł/.gemini/antigravity/brain/7c83af86-491a-484b-a5a4-fdc8094fa736")
if not ARTIFACT_DIR.exists():
    ARTIFACT_DIR = HERE

def run_display_tests():
    print("======================================================================")
    print("   ROZPOCZECIE KOMPLEKSOWYCH TESTOW WYSWIETLANIA DLA UZYTKOWNIKA     ")
    print("======================================================================")

    atlas = Atlas()
    sv = SurfaceViewer(atlas, off_screen=True)
    assert sv is not None, "Inicjalizacja SurfaceViewer nie powiodla sie"

    # TEST 1: Czysty porcelanowy posazek
    print("\n[TEST 1] Czysty porcelanowy posazek (Start bez szumu informacyjnego)...")
    assert sv.show_wire is False
    assert sv.show_points_mode is False
    assert len(sv.channels_on) == 0
    sv._view_front()
    p1_front = str(ARTIFACT_DIR / "user_test_1_porcelain_front.png")
    sv.pl.screenshot(p1_front)
    sv._view_back()
    p1_back = str(ARTIFACT_DIR / "user_test_1_porcelain_back.png")
    sv.pl.screenshot(p1_back)
    print(f"  -> Zapisano zrzuty: {p1_front}, {p1_back}")

    # TEST 2: Siatka CUN
    print("\n[TEST 2] Siatka proporcjonalna Cun [g]...")
    sv.toggle_cun_grid()
    assert sv.show_wire is True, "Siatka cun powinna byc wlaczona!"
    sv._view_front()
    p2_front = str(ARTIFACT_DIR / "user_test_2_cun_grid_front.png")
    sv.pl.screenshot(p2_front)
    sv._view_back()
    p2_back = str(ARTIFACT_DIR / "user_test_2_cun_grid_back.png")
    sv.pl.screenshot(p2_back)
    print(f"  -> Zapisano zrzuty: {p2_front}, {p2_back}")

    # TEST 3: Wszystkie punkty akupunkturowe
    print("\n[TEST 3] Wszystkie punkty akupunkturowe [p]...")
    sv.toggle_points()
    assert sv.show_points_mode is True, "Punkty powinny byc wlaczone!"
    p3 = str(ARTIFACT_DIR / "user_test_3_points_all.png")
    sv.pl.screenshot(p3)
    print(f"  -> Zapisano zrzut: {p3}")

    # Wylaczamy punkty i siatke dla czytelnosci kolejnych testow
    sv.toggle_points()
    sv.toggle_cun_grid()

    # TEST 4: Trakt grzbietowy Tai Yang (BL) + Du Mai (DU)
    print("\n[TEST 4] Trakt grzbietowy: Pecherz (BL) i Du Mai (DU)...")
    sv.set_channel("BL", True)
    sv.set_channel("DU", True)
    assert "BL" in sv.channels_on and "DU" in sv.channels_on
    sv._view_back()
    p4 = str(ARTIFACT_DIR / "user_test_4_meridians_bl_and_du_back.png")
    sv.pl.screenshot(p4)
    print(f"  -> Zapisano zrzut: {p4}")

    # TEST 5: 12 Glownych Meridianow Zang-Fu
    print("\n[TEST 5] 12 Glownych Meridianow Zang-Fu [m]...")
    sv.toggle_all_main_meridians()
    sv._view_front()
    p5 = str(ARTIFACT_DIR / "user_test_5_all_12_main_meridians.png")
    sv.pl.screenshot(p5)
    print(f"  -> Zapisano zrzut: {p5}")

    # Wylaczamy 12 glownych dla czytelnosci zasad TCM
    sv.toggle_all_main_meridians()
    sv.set_channel("DU", False)

    # TEST 6: Zasada Yuan-Luo (Gospodarz - Gosc)
    print("\n[TEST 6] Zasada TCM: Yuan-Luo (LU-9 + LI-6 oraz BL-64 + KI-4)...")
    sv.jump_to_yuan_luo("LU")
    assert sv.selected_point is not None
    assert bool(sv.highlight_actor_pair.GetVisibility()) is True
    assert sv.relation_tube_actor is not None and bool(sv.relation_tube_actor.GetVisibility()) is True
    p6_lu = str(ARTIFACT_DIR / "user_test_6_yuan_luo_lu_li.png")
    sv.pl.screenshot(p6_lu)

    sv.jump_to_yuan_luo("BL")
    sv._view_back()
    p6_bl = str(ARTIFACT_DIR / "user_test_6_yuan_luo_bl_ki.png")
    sv.pl.screenshot(p6_bl)
    print(f"  -> Zapisano zrzuty: {p6_lu}, {p6_bl}")

    # TEST 7: Zasada Shu-Mu (Przedni Mu i Tylny Shu)
    print("\n[TEST 7] Zasada TCM: Shu-Mu (GB-24 + BL-19)...")
    sv.jump_to_shu_mu("GB")
    p7 = str(ARTIFACT_DIR / "user_test_7_shumu_gb.png")
    sv.pl.screenshot(p7)
    print(f"  -> Zapisano zrzut: {p7}")

    # TEST 8: Punkty Szczelin Xi-Cleft
    print("\n[TEST 8] Zasada TCM: Xi-Cleft (LU-6)...")
    sv.jump_to_xi_cleft("LU")
    p8 = str(ARTIFACT_DIR / "user_test_8_xicleft_lu6.png")
    sv.pl.screenshot(p8)
    print(f"  -> Zapisano zrzut: {p8}")

    # TEST 9: Piec Przemian Wu Xing
    print("\n[TEST 9] Zasada TCM: Wu Xing (Drewno - GB i LR)...")
    sv.jump_to_wuxing("drewno")
    p9 = str(ARTIFACT_DIR / "user_test_9_wuxing_wood.png")
    sv.pl.screenshot(p9)
    print(f"  -> Zapisano zrzut: {p9}")

    # TEST 10: Rysownik i Analiza Polaczen
    print("\n[TEST 10] Szkicownik powiezio-napieciowy i analiza korelacji...")
    sv.toggle_camera_freeze(True)
    assert sv.camera_frozen is True
    # Symulacja narysowania pasma napiecia w strefie klatki/rak
    import numpy as np
    dummy_stroke = np.array([
        [0.05, 0.10, 0.35],
        [0.08, 0.12, 0.30],
        [0.12, 0.14, 0.25],
        [0.16, 0.15, 0.20],
    ])
    sv.draw_strokes.append(dummy_stroke)
    sv.analyze_current_strokes()
    p10 = str(ARTIFACT_DIR / "user_test_10_sketch_and_analysis.png")
    sv.pl.screenshot(p10)
    print(f"  -> Zapisano zrzut: {p10}")

    print("\n======================================================================")
    print("   WSZYSTKIE TESTY WYSWIETLANIA DLA UZYTKOWNIKA ZAKONCZONE SUKCESEM!  ")
    print("======================================================================")

if __name__ == "__main__":
    run_display_tests()
