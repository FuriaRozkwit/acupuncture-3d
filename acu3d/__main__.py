"""Punkt wejscia: python -m acu3d [opcje]"""
from __future__ import annotations

import argparse
import sys

from .model import Atlas
from .viewer import Viewer


def main(argv=None) -> int:
    for stream in (sys.stdout, sys.stderr):          # opisy sa w UTF-8, konsola bywa cp1250
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass

    ap = argparse.ArgumentParser(
        prog="acu3d",
        description="Przegladarka 3D punktow akupunkturowych, meridianow i miesni")
    ap.add_argument("--surface", action="store_true",
                    help="sama powloka ciala do ogladania: bez warstw i bez oznaczen")
    ap.add_argument("--wire", action="store_true",
                    help="silnik napiec: sama siatka cunowa, punkty i linie sil")
    ap.add_argument("--bare", action="store_true",
                    help="sam silnik: ruch, kanaly i zaznaczanie punktow, bez miesni")
    ap.add_argument("--procedural", action="store_true",
                    help="model proceduralny zamiast siatek anatomicznych")
    ap.add_argument("--height", type=float, default=None,
                    help="wzrost modelu proceduralnego w metrach")
    ap.add_argument("--channel", "-c", action="append", default=None,
                    help="pokaz tylko wskazany kanal, np. -c ST (mozna wielokrotnie)")
    ap.add_argument("--point", "-p", default=None,
                    help="podswietl punkt na starcie, np. ST-36 albo ST-36:L")
    ap.add_argument("--muscle", "-M", default=None,
                    help="wybierz miesien po fragmencie nazwy, np. gastrocnemius")
    ap.add_argument("--activation", type=float, default=None,
                    help="napiecie wybranego miesnia 0..1")
    ap.add_argument("--hypertrophy", type=float, default=None,
                    help="przerost wybranego miesnia, 1.0 = bez zmian")
    ap.add_argument("--list-muscles", action="store_true",
                    help="wypisz katalog miesni i zakoncz")
    ap.add_argument("--layer", "-l", action="append", default=None,
                    choices=["skin", "muscle", "skeleton"], help="warstwy widoczne na starcie")
    ap.add_argument("--skin-opacity", type=float, default=None, help="krycie skory 0..1")
    ap.add_argument("--no-meridians", action="store_true", help="ukryj linie kanalow")
    ap.add_argument("--labels", action="store_true", help="pokaz kody punktow")
    ap.add_argument("--tension", default=None,
                    help="wczytaj zapisana mape napiec (data/tension/*.json)")
    ap.add_argument("--tension-report", action="store_true",
                    help="wypisz zestawienie wczytanej mapy i zakoncz")
    ap.add_argument("--view", default=None, choices=["front", "back", "side", "top"])
    ap.add_argument("--screenshot", default=None, help="zapisz zrzut i zakoncz")
    args = ap.parse_args(argv)

    if args.list_muscles:
        from .muscles import MuscleSet
        for m in MuscleSet().muscles:
            print(f"{m.id:46s} {m.label}")
        return 0

    atlas = Atlas(args.height, procedural=args.procedural)
    if args.surface:
        from .surface_viewer import SurfaceViewer
        channels = [c.upper() for c in args.channel] if args.channel else None
        sv = SurfaceViewer(atlas, off_screen=bool(args.screenshot), channels=channels)
        if args.labels:
            sv.toggle_labels()
        if args.point:
            sv.select_point(args.point)
        sv._sync_legend()
        if args.screenshot:
            sv.screenshot(args.screenshot)
            print(f"zapisano {args.screenshot}")
            return 0
        sv.show()
        return 0
    if args.wire:
        from .wire_viewer import WireViewer
        w = WireViewer(atlas, off_screen=bool(args.screenshot))
        if args.tension:
            from pathlib import Path as _P
            print(f"wczytano {w.tension.load(_P(args.tension))} oznaczen")
            w._refresh()
        if args.tension_report:
            print(w.tension.line_report())
            return 0
        if args.view:
            cam, az = {"front": ("xz", 180), "back": ("xz", 0),
                       "side": ("yz", 0), "top": ("xy", 0)}[args.view]
            w._view(cam, az)
        if args.labels:
            w.toggle_labels()
        if args.screenshot:
            w.screenshot(args.screenshot)
            print(f"zapisano {args.screenshot}")
            return 0
        w.show()
        return 0

    v = Viewer(atlas, off_screen=bool(args.screenshot), bare=args.bare)

    if args.channel:
        chans = {c.upper() for c in args.channel}
        unknown = chans - set(atlas.channels)
        if unknown:
            print(f"nieznane kanaly: {', '.join(sorted(unknown))}", file=sys.stderr)
            print(f"dostepne: {', '.join(atlas.channels)}", file=sys.stderr)
            return 2
        v.set_channels(chans)
    if args.layer:
        for name, actor in v.layer_actors.items():
            actor.SetVisibility(name in args.layer)
    if args.bare and (args.muscle or args.tension or args.tension_report):
        print("--bare wyklucza miesnie i mape napiec", file=sys.stderr)
        return 2
    if args.skin_opacity is not None:
        v._set_skin_opacity(args.skin_opacity)
    if args.no_meridians:
        v.toggle_meridians()
    if args.labels:
        v.toggle_labels()
    if args.view:
        cam, az = {"front": ("xz", 180), "back": ("xz", 0),
                   "side": ("yz", 0), "top": ("xy", 0)}[args.view]
        v._view(cam, az)
    if args.muscle:
        if v.muscles is None:
            print("miesnie sa dostepne tylko w trybie anatomicznym", file=sys.stderr)
            return 2
        hits = v.muscles.find(args.muscle)
        if not hits:
            print(f"nie znaleziono miesnia: {args.muscle}", file=sys.stderr)
            return 2
        m = hits[0]
        if args.activation is not None:
            m.activation = args.activation
        if args.hypertrophy is not None:
            m.hypertrophy = args.hypertrophy
        v.muscles.apply(m)
        v.muscles.refresh()
        v.select_muscle(m)
        v.layer_actors["muscle"].SetVisibility(True)
    if args.tension:
        from pathlib import Path as _P
        v._ensure_tension()
        n = v.tension.load(_P(args.tension))
        v.tension_mode = True
        v._refresh_tension()
        print(f"wczytano {n} oznaczen z {args.tension}")
    if args.tension_report:
        if v.tension is None or not v.tension.marks:
            print("brak oznaczen - podaj --tension <plik>", file=sys.stderr)
            return 2
        print(v.tension.report())
        return 0
    if args.point:
        key = args.point.upper().replace(":", "·")
        if key not in atlas.by_key:
            key = key + "·R"
        if key not in atlas.by_key:
            print(f"nieznany punkt: {args.point}", file=sys.stderr)
            return 2
        v.focus(key)

    if args.screenshot:
        v.screenshot(args.screenshot)
        print(f"zapisano {args.screenshot}")
        return 0
    v.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
