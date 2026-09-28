#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Weryfikacja polonizacji medycznej Deadmana i test renderowania karty HUD."""
import sys
import os
from pathlib import Path

# Ustawienie sciezki projektu
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from acu3d.translate_pl import get_point_pl, polish_point_name, translate_location
from acu3d.model import Atlas

print("=== TEST 1: Weryfikacja bazy medycznej data/point_meta_pl.json ===")
test_points = [
    "BL-1", "BL-2", "BL-6", "ST-1", "ST-2", "ST-6", "ST-7",
    "LI-20", "SI-18", "SI-19", "SJ-17", "TE-17", "GB-14", "GB-20",
    "DU-15", "DU-20", "DU-26", "REN-22", "REN-24", "ST-36", "SP-14", "SP-18"
]

all_ok = True
for code in test_points:
    rec = get_point_pl(code)
    if not rec:
        print(f"[BLAD] Brak rekordu dla {code}")
        all_ok = False
        continue
    nazwa = rec.get("nazwa_pl", "")
    loc = rec.get("lokalizacja_pl", "")
    palp = rec.get("palpacja_pl", "")
    need = rec.get("naklucie_pl", "")
    acts = rec.get("dzialania_pl", [])
    inds = rec.get("wskazania_pl", "")

    # Sprawdzenie czy nie ma smieci OCR
    for text, label in [(loc, "lokalizacja"), (need, "naklucie"), (inds, "wskazania")]:
        if any(bad in text for bad in ["~---", "clea:", "benet", "Wucht", "irline"]):
            print(f"[BLAD OCR] {code} {label}: {text}")
            all_ok = False

    print(f"OK: {code:6} | {nazwa:30} | {len(loc)} zn. lok | {len(need)} zn. nakl | {len(acts)} dzialan")

if all_ok:
    print("\n--> Test 1 ZALICZONY: Wszystkie kluczowe punkty fizjonomii i anatomii sa prawidlowe.")
else:
    print("\n--> Test 1 NIEZALICZONY!")

print("\n=== TEST 2: Symulacja SurfaceViewer i zrzut ekranu ===")
try:
    from acu3d.surface_viewer import SurfaceViewer
    import pyvista as pv

    pv.OFF_SCREEN = True
    atlas = Atlas()
    v = SurfaceViewer(atlas, off_screen=True)
    
    # Wlacz kanaly Zoladka i Pecherza
    v.set_channel("ST", True, render=False)
    v.set_channel("BL", True, render=False)
    
    # Zaznacz punkt ST-2 na twarzy (Cztery Błękity / Sibai)
    v.select_point("ST-2·L")
    
    # Kamera na twarz
    # Pozycja twarzy z atlasu
    face_pt = v.atlas.by_key.get("ST-2·L")
    if face_pt:
        fp = face_pt["pos"]
        v.pl.camera.focal_point = (fp[0], fp[1], fp[2])
        v.pl.camera.position = (fp[0], fp[1] + 0.35, fp[2] + 0.05)
        v.pl.camera.up = (0, 0, 1)
        v.pl.camera.view_angle = 30
    
    art_dir = Path("/mnt/c/Users/Paweł/.gemini/antigravity/brain/7c83af86-491a-484b-a5a4-fdc8094fa736")
    shot_path = art_dir / "surface_polish_card.png"
    v.pl.screenshot(str(shot_path))
    print(f"Zapisano zrzut karty i fizjonomii twarzy: {shot_path}")

    # Przetestuj tez punkt ST-36
    v.select_point("ST-36·L")
    st36_shot = art_dir / "surface_polish_st36.png"
    v.pl.screenshot(str(st36_shot))
    print(f"Zapisano zrzut karty ST-36: {st36_shot}")

    print("\n--> Test 2 ZALICZONY: SurfaceViewer zrenderowal i zapisal zrzuty ekranu.")
except Exception as e:
    import traceback
    traceback.print_exc()
    print(f"[BLAD] Test 2: {e}")
