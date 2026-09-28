"""Przegladarka powierzchni z trybem edycji punktow (gizmo + panel wlasciwosci).

Osobny punkt wejscia, zeby tryb edycji nie wymagal zmian w __main__.py i
surface_viewer.py - te pliki maja teraz niezacommitowana prace z innej galezi.

    python -m acu3d.edit                    wszystkie kanaly
    python -m acu3d.edit -c LU -c LI        wybrane kanaly
    python -m acu3d.edit --point LU-7       od razu edycja tego punktu
    python -m acu3d.edit --point LU-7 --screenshot out.png
    python -m acu3d.edit --fabryczna        konfiguracja fabryczna zamiast domyslnej
    python -m acu3d.edit --profil NAZWA     konkretny profil

W oknie: [t] wlacza tryb edycji, klik w punkt pokazuje gizmo i panel,
[v] pokazuje konfiguracje (fabryczna i profile) i pozwala wybrac domyslna.
"""
from __future__ import annotations

import argparse
import sys

from .model import Atlas


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m acu3d.edit")
    ap.add_argument("-c", "--channel", action="append", help="kanal do pokazania (mozna kilka)")
    ap.add_argument("--point", help="punkt do edycji od razu, np. LU-7 albo LU-7·L")
    ap.add_argument("--screenshot", help="zapisz obraz zamiast otwierac okno")
    ap.add_argument("--profiles", action="store_true", help="otwórz panel wyboru konfiguracji [v]")
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--profil", help="wczytaj ten profil zamiast domyslnego")
    g.add_argument("--fabryczna", action="store_true", help="wczytaj konfiguracje fabryczna")
    args = ap.parse_args(argv)

    from .surface_viewer import SurfaceViewer
    from .edit_mode import attach

    atlas = Atlas(profile="fabryczna" if args.fabryczna else (args.profil or "auto"))
    channels = [c.upper() for c in args.channel] if args.channel else None
    sv = SurfaceViewer(atlas, off_screen=bool(args.screenshot), channels=channels)
    em = attach(sv)
    em.toggle()
    if args.profiles:
        em.toggle_profiles()
    if args.point:
        sv.select_point(args.point)
    if args.screenshot:
        sv.screenshot(args.screenshot)
        print(f"zapisano {args.screenshot}")
        return 0
    sv.show()
    return 0


if __name__ == "__main__":
    sys.exit(main())
