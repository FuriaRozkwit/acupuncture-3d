"""Pobranie i przygotowanie zasobow anatomicznych (setki MB, jednorazowo).

Siatek nie trzymamy w repozytorium - ten skrypt odtwarza katalog `assets/` od zera:

    1. Z-Anatomy.zip          (~87 MB)  -> Startup.blend z cala anatomia
    2. eksport z Blendera                -> PLY: 684 miesnie, szkielet, kosci wzorcowe
    3. BodyParts3D isa OBJ    (~140 MB) -> skora (FMA7163 = element FJ2810)
    4. obrot, decymacja, dopasowanie     -> assets/ready/*.vtp + data/*.json

Wymaga zainstalowanego Blendera (>= 4.0). Sciezke mozna podac przez --blender.

    python tools/fetch_assets.py [--blender "C:/Program Files/Blender Foundation/Blender 5.2/blender.exe"]

Licencje pobieranych danych: patrz ATTRIBUTION.md (CC BY-SA).
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import urllib.request
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DL = ROOT / "assets" / "download"
BLEND_DIR = ROOT / "assets" / "blend"

Z_ANATOMY = "https://github.com/Z-Anatomy/Models-of-human-anatomy/raw/master/Z-Anatomy.zip"
TA2 = "https://raw.githubusercontent.com/Z-Anatomy/Models-of-human-anatomy/master/TA2.csv"
BP3D = "https://dbarchive.biosciencedbc.jp/data/bodyparts3d/LATEST/isa_BP3D_4.0_obj_99.zip"
SKIN_ELEMENT = "FJ2810.obj"          # FMA7163 "skin" wg isa_element_parts.txt

BLENDER_CANDIDATES = [
    r"C:\Program Files\Blender Foundation\Blender 5.2\blender.exe",
    r"C:\Program Files\Blender Foundation\Blender 4.2\blender.exe",
    "/usr/bin/blender", "/Applications/Blender.app/Contents/MacOS/Blender",
]

REF_BONES = ",".join(
    ["Atlas (C1)", "Axis (C2)"] + [f"Vertebra C{i}" for i in range(3, 8)]
    + [f"Vertebra T{i}" for i in range(1, 13)] + [f"Vertebra L{i}" for i in range(1, 6)]
    + ["Sacrum", "Coccyx", "Manubrium of sternum", "Body of sternum", "Xiphoid process",
       "Hip bone.r", "Femur.r", "Patella.r", "Tibia.r", "Fibula.r", "Calcaneus.r", "Talus.r",
       "Humerus.r", "Radius.r", "Ulna.r", "Scapula.r", "Clavicle.r", "Mandible",
       "Third metacarpal bone.r", "First metacarpal bone.r"])


def download(url: str, dst: Path) -> Path:
    if dst.exists() and dst.stat().st_size > 1000:
        print(f"  jest juz: {dst.name} ({dst.stat().st_size / 1e6:.0f} MB)")
        return dst
    dst.parent.mkdir(parents=True, exist_ok=True)
    print(f"  pobieram {url}")
    with urllib.request.urlopen(url) as r, open(dst, "wb") as f:
        shutil.copyfileobj(r, f)
    print(f"  zapisano {dst.name} ({dst.stat().st_size / 1e6:.0f} MB)")
    return dst


def find_blender(explicit: str | None) -> str:
    for c in ([explicit] if explicit else []) + BLENDER_CANDIDATES:
        if c and Path(c).exists():
            return c
    found = shutil.which("blender")
    if found:
        return found
    raise SystemExit("Nie znalazlem Blendera. Podaj sciezke: --blender <plik>")


def run(cmd: list) -> None:
    print("  $ " + " ".join(str(c) for c in cmd[:3]) + " ...")
    p = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8", errors="replace")
    for line in (p.stdout or "").splitlines():
        if line.startswith("[") or "blad" in line.lower():
            print("    " + line)
    if p.returncode != 0:
        print(p.stderr[-2000:], file=sys.stderr)
        raise SystemExit(f"krok zakonczony bledem ({p.returncode})")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--blender", default=None)
    ap.add_argument("--skip-download", action="store_true")
    args = ap.parse_args()
    blender = find_blender(args.blender)
    py = sys.executable

    print("1/5  pobieranie Z-Anatomy")
    zanat = DL / "Z-Anatomy.zip"
    if not args.skip_download:
        download(Z_ANATOMY, zanat)
        download(TA2, DL / "TA2.csv")
    blend = BLEND_DIR / "Z-Anatomy" / "Startup.blend"
    if not blend.exists():
        print("  rozpakowuje Startup.blend (~300 MB)")
        BLEND_DIR.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(zanat) as z:
            z.extract("Z-Anatomy/Startup.blend", BLEND_DIR)

    print("2/5  eksport siatek z Blendera")
    run([blender, "--background", "--factory-startup", "--python",
         str(ROOT / "tools" / "blender_export.py"), "--", str(blend),
         str(ROOT / "assets" / "mesh"), "muscles,skeleton"])
    run([blender, "--background", "--factory-startup", "--python",
         str(ROOT / "tools" / "blender_export_named.py"), "--", str(blend),
         str(ROOT / "assets" / "mesh" / "ref"), REF_BONES])

    print("3/5  pobieranie skory z BodyParts3D")
    bp = DL / "bp3d_isa.zip"
    if not args.skip_download:
        download(BP3D, bp)
    skin_obj = DL / "isa_BP3D_4.0_obj_99" / SKIN_ELEMENT
    if not skin_obj.exists():
        with zipfile.ZipFile(bp) as z:
            name = [n for n in z.namelist() if n.endswith(SKIN_ELEMENT)][0]
            z.extract(name, DL)
        print(f"  wypakowano {SKIN_ELEMENT}")

    print("4/5  obrot i decymacja siatek")
    run([py, str(ROOT / "tools" / "prepare_meshes.py")])
    run([py, str(ROOT / "tools" / "prepare_skin.py")])

    print("5/5  kalibracja rigu z kosci")
    run([py, str(ROOT / "tools" / "calibrate_rig.py")])

    size = sum(f.stat().st_size for f in (ROOT / "assets" / "ready").rglob("*")) / 1e6
    print(f"\ngotowe. assets/ready: {size:.0f} MB")
    print("uruchom:  python -m acu3d")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
