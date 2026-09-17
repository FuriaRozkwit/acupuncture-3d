"""Weryfikacja polozenia punktow wzgledem opisow zrodlowych.

Kotwice w data/points.json sa autorskie - wpisalem je recznie na podstawie opisow
z podrecznika. Ten skrypt sprawdza je NIEZALEZNIE: wyciaga z tekstu zrodlowego
liczby, ktore da sie zmierzyc ("2 cun bocznie od pepka", "1,5 cun od wyrostka L2",
"3 cun powyzej kostki przysrodkowej") i porownuje je z tym, co faktycznie wyszlo
w modelu 3D.

To nie jest test tautologiczny: tekst pochodzi z ksiazki (data/point_meta.json),
a pomiar - z geometrii anatomicznej i linijki cun wyliczonej z kosci.

    python tools/verify_points.py [--tol 0.5] [--all]

Wymaga lokalnego data/point_meta.json (patrz ATTRIBUTION.md).
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

ALIAS = {"L.I.": "LI", "L.L": "LI", "LL": "LI", "KID": "KI", "HE": "HT", "HT": "HT",
         "P": "PC", "PC": "PC", "LIV": "LR", "LR": "LR", "TB": "SJ", "SJ": "SJ"}
CODE_RE = r"(L\.I\.|L\.L|LL|KID|LIV|REN|DU|LU|ST|SP|SI|BL|GB|SJ|TB|HE|PC|P)[-‑\s.]?(\d{1,2})"
NUM = r"(\d+(?:\.\d+)?)"

# ile cun ma 1 metr w danym rejonie - do przeliczania zmierzonych odleglosci
TRUNK_SEGS = {"trunk", "neck"}


def norm_code(ch: str, num: str) -> str:
    return f"{ALIAS.get(ch.upper(), ch.upper())}-{int(num)}"


class Verifier:
    def __init__(self, atlas: Atlas):
        self.a = atlas
        self.an = atlas.anatomy
        self.anchors = {p["code"]: p["anchor"] for p in atlas.points}
        self.claims = []
        try:
            self.ms = MuscleSet()
            self.muscle_names = {m.name.lower() for m in self.ms.muscles
                                 if len(m.name) > 7}
        except Exception:
            self.ms, self.muscle_names = None, set()

    # --- pomiary w modelu ------------------------------------------------
    def pos(self, code: str, side: str = "R"):
        m = self.a.by_key.get(code) or self.a.by_key.get(f"{code}·{side}")
        return None if m is None else np.asarray(m["pos"])

    def cun_scale(self, code: str, axis: str) -> float:
        """Metry na 1 cun w rejonie danego punktu."""
        seg = self.anchors.get(code, {}).get("seg", "trunk")
        if seg in TRUNK_SEGS or seg == "head":
            return self.an.cun_trunk_lat if axis == "lat" else self.an.cun_trunk_vert
        return self.an.seg(seg).cun_m

    def landmark(self, name: str):
        n = name.lower()
        if "umbilic" in n:
            return self.pos("REN-8")
        if "medial malleolus" in n:
            return np.array(self.an.landmarks["medial_malleolus"])
        if "lateral malleolus" in n:
            return np.array(self.an.landmarks["lateral_malleolus"])
        if "wrist" in n:
            return np.array(self.an.landmarks["wrist"])
        if "anterior hairline" in n:
            z = self.an.levels["hairline_ant"]
            return np.array([0.0, 0.0, z])
        if "pubic symphysis" in n:
            return np.array(self.an.landmarks["pubis"])
        if "sternocostal angle" in n or "xiphoid" in n:
            return np.array(self.an.landmarks["xiphoid"])
        return None

    def add(self, code, kind, claimed, measured, detail):
        if measured is None:
            return
        self.claims.append({"code": code, "kind": kind, "claimed": claimed,
                            "measured": measured, "err": measured - claimed,
                            "detail": detail})

    # --- parsery opisow --------------------------------------------------
    def check_lateral(self, code: str, text: str) -> None:
        m = re.search(NUM + r"\s*cun lateral to the (anterior |posterior )?(midline|"
                      r"umbilicus|" + CODE_RE + ")", text, re.I)
        if not m:
            return
        claimed = float(m.group(1))
        p = self.pos(code)
        if p is None:
            return
        ref_txt = m.group(3) or ""
        if re.match(CODE_RE, ref_txt, re.I):
            g = re.match(CODE_RE, ref_txt, re.I)
            ref = self.pos(norm_code(g.group(1), g.group(2)))
            if ref is None:
                return
            dist = abs(abs(p[0]) - abs(ref[0]))
        else:
            dist = abs(p[0])                       # od linii posrodkowej
        self.add(code, "bocznie", claimed, dist / self.cun_scale(code, "lat"),
                 f"{claimed} cun bocznie")

    def check_axial(self, code: str, text: str) -> None:
        m = re.search(NUM + r"\s*cun (superior|inferior|proximal|distal|above|below)"
                      r"(?: to| from)? (?:the )?([A-Za-z .‑-]{2,28})", text, re.I)
        if not m:
            return
        claimed, ref_txt = float(m.group(1)), m.group(3)
        p = self.pos(code)
        if p is None:
            return
        g = re.search(CODE_RE, ref_txt, re.I)
        ref = self.pos(norm_code(g.group(1), g.group(2))) if g else self.landmark(ref_txt)
        if ref is None:
            return
        seg = self.anchors.get(code, {}).get("seg", "trunk")
        if seg in TRUNK_SEGS:
            # linijka tulowia jest ODCINKOWA (klatka, nadbrzusze, podbrzusze maja
            # wlasne podzialy), wiec mierzymy roznice wspolrzednych w cun, a nie
            # metry przeliczone jedna skala
            measured = abs(self.an.cun_of_z(float(p[2])) - self.an.cun_of_z(float(ref[2])))
        else:
            measured = float(np.linalg.norm(p - ref)) / self.cun_scale(code, "ax")
        self.add(code, "wzdluz", claimed, measured,
                 f"{claimed} cun od {ref_txt.strip()[:20]}")

    def check_vertebra(self, code: str, text: str) -> None:
        m = re.search(r"spinous process of the .*?\((T|L|C)(\d{1,2})\)", text, re.I)
        if not m:
            m = re.search(r"(?:^|\s)(T|L|C)(\d{1,2})(?:\s|\)|\.)", text)
        if not m:
            return
        key = f"{m.group(1).upper()}{int(m.group(2))}"
        z_ref = self.an.vertebra_z.get(key)
        p = self.pos(code)
        if z_ref is None or p is None:
            return
        # punkt ma lezec przy DOLNYM brzegu wyrostka, czyli ciut nizej
        err_cm = (p[2] - z_ref) * 100
        self.claims.append({"code": code, "kind": "krag", "claimed": 0.0,
                            "measured": err_cm, "err": err_cm,
                            "detail": f"poziom {key}"})

    def check_muscle(self, code: str, text: str) -> None:
        """Test NIEZALEZNY: czy miesien wymieniony w opisie jest tym, ktory model
        ma pod punktem. Nazwa pochodzi z ksiazki, geometria z Z-Anatomy - dwa
        rozne zrodla, wiec zgodnosc nie moze wyjsc z zalozenia."""
        if self.ms is None:
            return
        low = " " + re.sub(r"[^a-z ]", " ", text.lower()) + " "
        named = {n for n in self.muscle_names if f" {n} " in low}
        if not named:
            return
        p = self.pos(code)
        if p is None:
            return
        near = self.nearest_muscles(self.in_muscle_frame(p, code), k=3)
        hit = any(any(n in m.name.lower() for n in named) for m, _ in near)
        self.claims.append({
            "code": code, "kind": "miesien", "claimed": 1.0,
            "measured": 1.0 if hit else 0.0, "err": 0.0 if hit else 1.0,
            "detail": f"opis: {sorted(named)[0][:22]} | model: {near[0][0].name[:26]}"
            f" ({near[0][1] * 1000:.0f} mm)"})

    def in_muscle_frame(self, p: np.ndarray, code: str) -> np.ndarray:
        """Punkt przeniesiony do ukladu siatek Z-Anatomy.

        Punkty na ramieniu i przedramieniu liczymy od srodka przekroju SKORY
        (BodyParts3D), a miesnie sa z Z-Anatomy - oba zbiory dzieli na konczynach
        12-17 mm. Bez tej poprawki test "czy pod punktem jest miesien z opisu"
        mierzylby blad dopasowania siatek, a nie polozenie punktu.
        """
        seg_name = self.anchors.get(code, {}).get("seg")
        shift = getattr(self.an, "section_shift", None)
        if shift is None or seg_name not in getattr(self.an, "LIMBS", ()):
            return p
        seg = self.an.seg(seg_name)
        t = float(np.dot(p - seg.p0, seg.e_ax)) / seg.length
        return p - shift(seg_name, float(np.clip(t, 0.0, 1.0)))

    def nearest_muscles(self, p: np.ndarray, k: int = 3) -> list:
        d = np.linalg.norm(self.ms.rest - p, axis=1)
        out, taken = [], set()
        for j in np.argsort(d)[:4000]:
            mi = int(self.ms._mid_points[j])
            if mi in taken:
                continue
            taken.add(mi)
            out.append((self.ms.muscles[mi], float(d[j])))
            if len(out) >= k:
                break
        return out

    def run(self) -> None:
        for code, meta in self.a.meta.items():
            if code not in self.anchors:
                continue
            text = meta.get("location", "")
            if not text:
                continue
            self.check_lateral(code, text)
            self.check_axial(code, text)
            self.check_vertebra(code, text)
            self.check_muscle(code, text)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tol", type=float, default=0.5, help="tolerancja w cun")
    ap.add_argument("--all", action="store_true", help="wypisz wszystkie porownania")
    args = ap.parse_args()

    atlas = Atlas()
    if not atlas.meta:
        print("brak data/point_meta.json - nie ma z czym porownywac", file=sys.stderr)
        return 2
    v = Verifier(atlas)
    v.run()

    by_kind = {}
    for c in v.claims:
        by_kind.setdefault(c["kind"], []).append(c)

    print(f"model: {atlas.mode}, 1 cun tulowia: pion {atlas.anatomy.cun_trunk_vert * 100:.2f} cm, "
          f"poziom {atlas.anatomy.cun_trunk_lat * 100:.2f} cm")
    print(f"sprawdzonych twierdzen z tekstu zrodlowego: {len(v.claims)}\n")

    for kind, rows in sorted(by_kind.items()):
        errs = np.array([r["err"] for r in rows])
        unit = "cm" if kind == "krag" else "cun"
        ok = np.abs(errs) <= (1.5 if kind == "krag" else args.tol)
        print(f"--- {kind} ({len(rows)} punktow, jednostka: {unit}) ---")
        print(f"    zgodnych w tolerancji: {ok.sum()}/{len(rows)} ({100 * ok.mean():.0f}%)")
        print(f"    blad: mediana {np.median(np.abs(errs)):.2f}, "
              f"srednia {np.abs(errs).mean():.2f}, max {np.abs(errs).max():.2f} {unit}")
        worst = sorted(rows, key=lambda r: -abs(r["err"]))[:6]
        for w in worst:
            print(f"      {w['code']:8s} deklarowane {w['claimed']:5.1f} -> "
                  f"zmierzone {w['measured']:6.2f}  ({w['err']:+.2f} {unit})  {w['detail']}")
        print()

    if args.all:
        for c in sorted(v.claims, key=lambda r: r["code"]):
            print(f"{c['code']:8s} {c['kind']:8s} {c['claimed']:5.1f} -> {c['measured']:6.2f} "
                  f"({c['err']:+.2f})  {c['detail']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
