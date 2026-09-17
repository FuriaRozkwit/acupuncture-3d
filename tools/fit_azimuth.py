"""Dopasowanie azymutow punktow konczynowych do anatomii -> data/azimuth_fit.json.

Azymut (obrot punktu wokol konczyny) byl najslabszym miejscem danych: dlugosc
wzdluz osi wynika z kosci, ale kat wpisywalem na oko. Tymczasem opisy zrodlowe
wymieniaja struktury: "miedzy sciegnami", "na miesniu piszczelowym przednim",
"promieniowo od sciegna dwuglowego". Skoro mamy te miesnie jako osobne siatki,
kat da sie ZMIERZYC: skanujemy pelny obrot i wybieramy ten, przy ktorym punkt
lezy najblizej struktury nazwanej w podreczniku.

Nazwa struktury pochodzi z ksiazki, geometria z Z-Anatomy, a wynik nie zalezy od
mojej pierwotnej wartosci - poza tym, ze odrzucamy skoki wieksze niz `--max-shift`,
zeby literowka w OCR nie przerzucila punktu na druga strone konczyny.

    python tools/fit_azimuth.py [--max-shift 50] [--apply]
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
warnings.simplefilter("ignore")

from acu3d.model import Atlas  # noqa: E402
from acu3d.muscles import MuscleSet  # noqa: E402

LIMB_SEGS = {"upper_arm", "forearm", "thigh", "shank", "neck"}
STOP = {"muscle", "tendon", "tendons", "border", "belly", "head", "process", "bone"}


def aliases(name: str) -> set:
    """Warianty nazwy, ktorymi podrecznik moze sie poslugiwac.

    Katalog ma "Brachioradialis muscle" i "Long head of biceps brachii",
    a ksiazka pisze "brachioradialis" i "biceps brachii" - bez tego dopasowania
    przepadala wiekszosc opisow.
    """
    n = name.lower().strip()
    out = {n}
    for suffix in (" muscle", " muscle of hand", " muscle of foot"):
        if n.endswith(suffix):
            out.add(n[: -len(suffix)])
    if " of " in n:                       # "long head of biceps brachii"
        tail = n.split(" of ", 1)[1]
        out.add(tail)
        for suffix in (" muscle",):
            if tail.endswith(suffix):
                out.add(tail[: -len(suffix)])
    return {a for a in out if len(a) >= 8}


def named_structures(text: str, catalog: dict) -> list:
    """Struktury wymienione w opisie punktu (miesnie i kosci)."""
    low = " " + re.sub(r"[^a-z ]", " ", text.lower()) + " "
    hits = []
    for name, pts in catalog.items():
        if any(f" {a} " in low for a in aliases(name)):
            hits.append((name, pts))
    return hits


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-shift", type=float, default=50.0, help="maksymalna zmiana kata")
    ap.add_argument("--step", type=float, default=4.0)
    ap.add_argument("--max-dist", type=float, default=12.0,
                    help="ile mm na strukture wolno zostac, zeby uznac dopasowanie")
    args = ap.parse_args()

    atlas = Atlas()
    if not atlas.meta:
        print("brak data/point_meta.json - nie ma opisow do dopasowania", file=sys.stderr)
        return 2
    ms = MuscleSet()
    body = atlas.anatomy

    catalog: dict = {}
    for m in ms.muscles:                       # miesnie: nazwa -> chmura wierzcholkow
        if m.side in ("R", "C"):
            catalog.setdefault(m.name.lower(), []).append(ms.rest[m.slice])
    for bone in ("Tibia", "Fibula", "Radius", "Ulna", "Humerus", "Femur",
                 "Scaphoid bone", "Navicular bone"):
        try:
            catalog.setdefault(bone.lower(), []).append(
                np.asarray(body.struct(bone).mesh.points))
        except KeyError:
            pass
    catalog = {k: np.vstack(v) for k, v in catalog.items()}

    fits, skipped, rejected = {}, [], []
    for p in atlas.points:
        a = p["anchor"]
        if a.get("seg") not in LIMB_SEGS or "az" not in a:
            continue
        # "po stronie promieniowej sciegna X" - minimum odleglosci kladzie punkt NA X,
        # wiec takie azymuty ustala sie pomiarem brzegu struktury i tu sie ich nie rusza
        if p.get("az_source", "").startswith("pomiar:"):
            skipped.append(p["code"])
            continue
        text = atlas.meta.get(p["code"], {}).get("location", "")
        named = named_structures(text, catalog)
        if not named:
            skipped.append(p["code"])
            continue
        clouds = [pts for _, pts in named]
        seg = body.seg(a["seg"])
        t = __import__("acu3d.anchors", fromlist=["_axial_t"])._axial_t(body, a)
        az0 = float(a["az"])

        def score(pos):
            """Jedna struktura: punkt ma na niej lezec. Kilka: ma lezec MIEDZY nimi,
            wiec sumujemy odleglosci do kazdej z osobna - minimum sumy wypada
            pomiedzy, a nie na najblizszej."""
            return float(sum(np.min(np.linalg.norm(c - pos, axis=1)) for c in clouds))

        best, best_d = az0, None
        for az in np.arange(-180.0, 180.0, args.step):
            if abs(((az - az0 + 180) % 360) - 180) > args.max_shift:
                continue
            d = score(seg.surface_by_azimuth(t, float(az)))
            if best_d is None or d < best_d:
                best, best_d = float(az), d
        d0 = score(seg.surface_by_azimuth(t, az0))
        if best_d is None:
            continue
        # Bramka jakosci: jesli punkt i tak nie dosiega nazwanej struktury, to znaczy,
        # ze nazwa trafila do opisu mimochodem (podpisy rycin w OCR) - odrzucamy.
        limit = args.max_dist * len(clouds) / 1000.0
        if best_d > limit:
            rejected.append((p["code"], round(best_d * 1000, 1), [n for n, _ in named][:1]))
            continue
        fits[p["code"]] = {"az_old": az0, "az_new": round(best, 1),
                           "dist_old_mm": round(d0 * 1000, 1),
                           "dist_new_mm": round(best_d * 1000, 1),
                           "structures": [n for n, _ in named], "n_struct": len(named)}

    out = {"note": "azymuty dopasowane do struktur wymienionych w opisach zrodlowych",
           "tool": "tools/fit_azimuth.py", "fits": fits}
    (Path(__file__).resolve().parent.parent / "data" / "azimuth_fit.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")

    moved = {k: v for k, v in fits.items() if abs(v["az_new"] - v["az_old"]) >= 5}
    gain = [v["dist_old_mm"] - v["dist_new_mm"] for v in fits.values()]
    print(f"przyjete dopasowania: {len(fits)}   odrzucone (nazwa mimochodem): {len(rejected)}"
          f"   bez nazwanej struktury: {len(skipped)}")
    if rejected:
        print("  odrzucone: " + ", ".join(f"{c} ({d:.0f} mm)" for c, d, _ in rejected[:10]))
    print(f"kat zmieniony o >=5 stopni: {len(moved)}")
    if gain:
        print(f"poprawa odleglosci od struktury: mediana {np.median(gain):.1f} mm, "
              f"srednia {np.mean(gain):.1f} mm")
    print("\nnajwieksze korekty:")
    for code, v in sorted(moved.items(), key=lambda kv: -abs(kv[1]["az_new"] - kv[1]["az_old"]))[:14]:
        print(f"  {code:8s} {v['az_old']:+7.1f} -> {v['az_new']:+7.1f}   "
              f"odleglosc {v['dist_old_mm']:5.1f} -> {v['dist_new_mm']:5.1f} mm   "
              f"({v['structures'][0][:30]})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
