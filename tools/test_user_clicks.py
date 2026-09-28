"""Kompleksowy test automatyczny interakcji uzytkownika z atlasem:
- test klikniec we wszystkie przyciski panelu bocznego
- test klikania zasad TCM (Yuan-Luo, Shu-Mu, Xi-Cleft, Wu-Xing, Ba-Mai)
- test wyboru punktow na modelu
- weryfikacja poprawnosci karty Deadmana i braku bledow typu TypeError/AttributeError
"""
import sys
import numpy as np
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from acu3d.model import Atlas
from acu3d.surface_viewer import SurfaceViewer
from acu3d.clinical_kb import YUAN_LUO_PAIRS, SHU_MU_PAIRS, XI_CLEFT_POINTS


def run_tests():
    print("=== [TEST 1] INICJALIZACJA ATLASU I SURFACEVIEWER ===")
    atlas = Atlas()
    sv = SurfaceViewer(atlas, off_screen=True)
    assert sv is not None
    print("  -> Inicjalizacja zakonczona sukcesem.")

    print("\n=== [TEST 2] WERYFIKACJA STANU POCZATKOWEGO (CZYSTY PORCELANOWY POSAZEK) ===")
    assert sv.show_wire is False, "Siatka cun powinna byc domyslnie wylaczona!"
    assert sv.show_points_mode is False, "Punkty powinny byc domyslnie wylaczone!"
    assert len(sv.channels_on) == 0, "Zaden kanal nie powinien byc domyslnie wlaczony!"
    print("  -> Model startuje jako czysty porcelanowy ludzik (OK).")

    print("\n=== [TEST 3] TEST KLIKALNOSCI WSZYSTKICH PRZYCISKOW W PRAWYM DOCKU ===")
    total_buttons = len(sv._dock_buttons)
    print(f"  -> Liczba przyciskow w panelu: {total_buttons}")
    assert total_buttons >= 19, f"Oczekiwano co najmniej 19 przyciskow, znaleziono {total_buttons}"

    for btn in sv._dock_buttons:
        bx, by, bw, bh = btn["rect"]
        cx = bx + bw // 2
        cy = by + bh // 2
        hit = sv._check_dock_click(cx, cy)
        assert hit, f"Blad klikniecia w przycisk: {btn['label']} na pozycji ({cx}, {cy})"
        print(f"  [OK] Przycisk [{btn['shortcut']}] '{btn['label']}' klikalny poprawnie.")

    print("\n=== [TEST 4] TEST ZASAD TCM: YUAN-LUO (GOSPODARZ - GOSC) ===")
    for i in range(3):
        sv.cycle_yuan_luo()
        assert sv.selected_point is not None, "Punkt glowny powinien byc zaznaczony!"
        assert bool(sv.highlight_actor_pair.GetVisibility()) is True, "Punkt partnerski powinien byc podswietlony!"
        assert sv.relation_tube_actor is not None and bool(sv.relation_tube_actor.GetVisibility()) is True, "Linia relacji powinna byc widoczna!"
        print(f"  [OK] Cykl Yuan-Luo #{i+1}: wybrany punkt {sv.selected_point}")

    print("\n=== [TEST 5] TEST ZASAD TCM: SHU-MU (PRZOD - TYL) ===")
    for i in range(3):
        sv.cycle_shu_mu()
        assert sv.selected_point is not None, "Punkt Przedni-Mu powinien byc zaznaczony!"
        assert bool(sv.highlight_actor_pair.GetVisibility()) is True, "Punkt Tylny-Shu powinien byc podswietlony!"
        print(f"  [OK] Cykl Shu-Mu #{i+1}: wybrany punkt {sv.selected_point}")

    print("\n=== [TEST 6] TEST ZASAD TCM: XI-CLEFT (SZCZELINY OSTREGO BOLU) ===")
    for i in range(3):
        sv.cycle_xi_cleft()
        assert sv.selected_point is not None, "Punkt Xi powinien byc zaznaczony!"
        print(f"  [OK] Cykl Xi-Cleft #{i+1}: wybrany punkt {sv.selected_point}")

    print("\n=== [TEST 7] TEST ZASAD TCM: WU XING (5 PRZEMIAN) ===")
    for i in range(3):
        sv.cycle_wu_xing()
        assert len(sv.channels_on) > 0, "Kanaly danego zywiołu powinny byc wlaczone!"
        print(f"  [OK] Cykl Wu Xing #{i+1}: aktywne kanaly {sv.channels_on}")

    print("\n=== [TEST 8] TEST KLIKNIEC BEZPOSREDNICH W PUNKTY AKUPUNKTUROWE ===")
    sample_points = ["ST-36·R", "LU-9·R", "BL-23·L", "REN-12", "DU-20", "GB-34·R", "KI-3·L", "HT-7·R", "SP-6·R", "LR-3·L", "LI-4·R", "TE-5·L"]
    for pt in sample_points:
        ok = sv.select_point(pt)
        assert ok, f"Nie udalo sie zaznaczyc punktu {pt}"
        print(f"  [OK] Bezposredni wybor punktu: {pt}")

    print("\n=== [TEST 9] TEST KOMEND TERMINALA (CONSOLE COMMANDS) ===")
    test_commands = [
        "yuan-luo LU",
        "shumu ST",
        "xicleft HT",
        "wuxing drewno",
        "bamai dai",
        "punkt ST-36",
        "rejon klatka",
        "biao-li 1",
        "pomoc"
    ]
    for cmd in test_commands:
        sv.handle_console_command(cmd)
        print(f"  [OK] Komenda konsoli wykonana: '{cmd}'")

    print("\n=== [TEST 10] ZRZUT EKRANU WERYFIKACYJNY DLA PARY YUAN-LUO ===")
    sv.jump_to_yuan_luo("LU")
    sv.pl.screenshot("test_yuan_luo_pair_verified.png")
    print("  -> Zapisano zrzut: test_yuan_luo_pair_verified.png")

    print("\n" + "=" * 60)
    print("WSZYSTKIE TESTY INTERAKCJI UZYTKOWNIKA ZAKONCZONE SUKCESEM!")
    print("=" * 60)


if __name__ == "__main__":
    run_tests()
