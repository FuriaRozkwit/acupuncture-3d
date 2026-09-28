"""Zarzadzanie konfiguracjami punktow z wiersza polecen.

    python -m acu3d.profil lista                 profile, domyslny oznaczony *
    python -m acu3d.profil pokaz [NAZWA]         co profil zmienia wzgledem fabrycznej
    python -m acu3d.profil uzyj NAZWA            ustaw domyslny (laduje sie przy starcie)
    python -m acu3d.profil fabryczna             domyslna = sama fabryka
    python -m acu3d.profil nowy NAZWA [--z INNY] nowy profil, pusty albo kopia
    python -m acu3d.profil przywroc NAZWA KOD    jeden punkt z powrotem na fabryczny
    python -m acu3d.profil usun NAZWA

Fabryczny data/points.json nie jest tu nigdy zapisywany.
"""
from __future__ import annotations

import argparse
import sys

from . import profiles as P


def _fmt(a: dict) -> str:
    return ", ".join(f"{k} {v}" for k, v in a.items() if k != "region")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m acu3d.profil")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("lista")
    s = sub.add_parser("pokaz"); s.add_argument("nazwa", nargs="?")
    s = sub.add_parser("uzyj"); s.add_argument("nazwa")
    sub.add_parser("fabryczna")
    s = sub.add_parser("nowy"); s.add_argument("nazwa"); s.add_argument("--z", dest="z")
    s = sub.add_parser("przywroc"); s.add_argument("nazwa"); s.add_argument("kod")
    s = sub.add_parser("usun"); s.add_argument("nazwa")
    a = ap.parse_args(argv)

    try:
        if a.cmd == "lista":
            default = P.active_name()
            for name in [P.FACTORY] + P.names():
                n = "" if name == P.FACTORY else f"{len(P.load(name)['points']):3d} pkt"
                print(f" {'*' if name == default else ' '} {name:<24s}{n}")
            print("\n * = laduje sie domyslnie")
        elif a.cmd == "pokaz":
            name = a.nazwa or P.active_name()
            if name == P.FACTORY:
                print("fabryczna: brak zmian (to jest punkt odniesienia)")
                return 0
            fac = {p["code"]: p["anchor"] for p in P.factory_points()}
            doc = P.load(name)
            print(f"profil '{name}': {len(doc['points'])} punktow zmienionych")
            for code, rec in doc["points"].items():
                stale = rec.get("factory_anchor") not in (None, fac.get(code))
                print(f"\n {code}{'   ! fabryczna zmienila sie od edycji' if stale else ''}")
                print(f"   fabryczna: {_fmt(fac.get(code, {}))}")
                print(f"   profil:    {_fmt(rec['anchor'])}")
                if rec.get("source"):
                    print(f"   zrodlo:    {rec['source'][:110]}")
        elif a.cmd == "uzyj":
            P.set_active(a.nazwa)
            print(f"domyslny: {P.active_name()}")
        elif a.cmd == "fabryczna":
            P.set_active(P.FACTORY)
            print("domyslna: fabryczna")
        elif a.cmd == "nowy":
            print(f"utworzono '{P.create(a.nazwa, a.z)}'" + (f" (kopia '{a.z}')" if a.z else ""))
        elif a.cmd == "przywroc":
            ok = P.remove_point(a.nazwa, a.kod.upper())
            print(f"{a.kod.upper()}: {'fabryczny' if ok else 'nie byl zmieniony w tym profilu'}")
        elif a.cmd == "usun":
            P.delete(a.nazwa)
            print(f"usunieto '{a.nazwa}'; domyslny: {P.active_name()}")
    except (KeyError, ValueError, FileExistsError, FileNotFoundError) as e:
        print(f"blad: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
