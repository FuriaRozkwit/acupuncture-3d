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


# "na poziomie wyrostka kolczystego X" znaczy przy jego DOLNYM brzegu - o tyle
# cun nizej niz punkt odniesienia kregu w `rig_fitted`. To konwencja atlasu,
# spisana raz tutaj, zeby test mial z czym porownywac.
VERT_DZ_CUN = -0.4
REF = Path(__file__).resolve().parent.parent / "assets" / "mesh" / "ref"


def in_gloss(text: str, at: int) -> bool:
    """Czy trafienie siedzi w nawiasie objasniajacym INNY punkt.

    "mniej wiecej na poziomie REN-4 (3 cun ponizej pepka)" nie mowi, gdzie lezy
    opisywany punkt - mowi, gdzie lezy REN-4. Bez tego GB-27 wychodzil o 0.8 cun
    obok, choc siedzi 2 mm od kolca biodrowego, czyli tam, gdzie ma siedziec.
    """
    o = text.rfind("(", 0, at)
    if o < 0 or 0 <= text.rfind(")", 0, at) > o:
        return False
    if text.find(")", at) < 0:
        return False
    return re.search(CODE_RE + r"\s*$", text[:o].strip(), re.I) is not None


class Verifier:
    def __init__(self, atlas: Atlas):
        self.a = atlas
        self.an = atlas.anatomy
        self.anchors = {p["code"]: p["anchor"] for p in atlas.points}
        self.claims = []
        try:
            self.ms = MuscleSet()
            # Ksiazka pisze "brachioradialis", katalog "Brachioradialis muscle" -
            # bez odciecia koncowki test gubil polowe wzmianek (LU-7 widzial tylko
            # jeden z dwoch sciegien z opisu i nie mogl sprawdzic "w szczelinie miedzy").
            self.muscle_names = set()
            for m in self.ms.muscles:
                n = m.name.lower()
                for form in (n, re.sub(r"\s+muscle$", "", n)):
                    if len(form) > 7:
                        self.muscle_names.add(form)
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
        if in_gloss(text, m.start()):
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
        # "w polowie miedzy SI-10 a wyrostkiem T2" nie jest twierdzeniem o poziomie
        # tego wyrostka - punkt ma lezec WYZEJ. Bez tego SI-13 wychodzil o 1.9 cm
        # za nisko, choc opis mowi co innego.
        if re.search(r"(midway|halfway)\s+between", text, re.I):
            return
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
        # "na poziomie wyrostka kolczystego T5" znaczy przy jego DOLNYM brzegu,
        # a nie na jego srodku - kotwice maja na to `dz_cun` (domyslnie -0.4 cun).
        # Porownywanie z samym `z_ref` mierzylo wiec wlasna konwencje atlasu:
        # 31 punktow na 31 wychodzilo "za nisko" o -0.90 cm przy odchyleniu 0.37,
        # czyli test nie mogl wykryc zadnej pomylki.
        # UWAGA: oczekiwanie bierzemy ze stalej konwencji, a NIE z kotwicy punktu -
        # inaczej test porownywalby kotwice sama ze soba i zawsze przechodzil.
        dz_cun = self.spinous_low(key)
        want_cm = dz_cun * self.cun_scale(code, "vert") * 100
        got_cm = (p[2] - z_ref) * 100
        self.claims.append({"code": code, "kind": "krag", "claimed": round(want_cm, 2),
                            "measured": round(got_cm, 2), "err": got_cm - want_cm,
                            "detail": f"poziom {key} ({dz_cun:+.1f} cun od wyrostka)"})

    def spinous_low(self, key: str) -> float:
        """Zmierzony dolny brzeg wyrostka kolczystego, w cun pod jego czubkiem.

        `rig_fitted` trzyma czubek (sprawdzone: 0.00 cun od niego). Opis mowi
        "na poziomie DOLNEGO brzegu", a to nie jest jedna liczba dla calego
        kregoslupa: mediana -0.30 cun, ale T1 -0.08 a L1 -0.56 (odchylenie 0.14).
        Stala VERT_DZ_CUN zostaje tylko jako awaryjna, gdy brak siatki kregu.
        """
        cache = getattr(self, "_low", None)
        if cache is None:
            cache = self._low = {}
        if key in cache:
            return cache[key]
        f = REF / f"vertebra-{key.lower()}.ply"
        v = self.an.fit["vertebra"].get(key)
        out = VERT_DZ_CUN
        if f.exists() and v is not None:
            import pyvista as pv
            Q = np.asarray(pv.read(f).points).copy()
            Q[:, 0] *= -1.0
            Q[:, 1] *= -1.0
            sp = np.asarray(v["spinous"])
            tail = Q[(Q[:, 1] < np.percentile(Q[:, 1], 15))
                     & (np.abs(Q[:, 0] - sp[0]) < 0.010)]
            if len(tail) >= 20:
                out = float((tail[:, 2].min() - sp[2]) / self.an.cun_trunk_vert)
        cache[key] = out
        return out

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
        q = self.in_muscle_frame(p, code)

        # "w szczelinie MIEDZY sciegnami A i B" to inne twierdzenie niz "nad A":
        # punkt ma lezec POSRODKU, a najblizszy moze byc trzeci miesien, ktory
        # opis pomija. Tak jest z LU-9: 22.8 mm do APL i 22.7 do FCR, czyli
        # dokladnie tam, gdzie ma byc - a regula "nazwany ma byc najblizszy"
        # zglaszala blad, bo najblizej jest ramienno-promieniowy.
        pair = self.between_pair(text, named)
        if pair:
            (na, da), (nb, db) = [(n, self.muscle_gap(q, n)) for n in pair]
            if da is not None and db is not None:
                # Mierzymy odleglosc od DALSZEJ z dwoch struktur. Sama symetria
                # nie wystarcza: punkt po drugiej stronie przedramienia lezy
                # 48 mm od obu sciegien, czyli "idealnie posrodku" i zupelnie
                # nie tam, gdzie trzeba.
                u = self.cun_scale(code, "ax")
                got = max(da, db) / u
                self.claims.append({
                    "code": code, "kind": "miesien", "claimed": 0.0,
                    "measured": round(got, 2), "err": got,
                    "detail": f"miedzy {na[:12]} ({da*1000:.0f} mm) a {nb[:12]}"
                              f" ({db*1000:.0f} mm)"})
                return

        near = self.nearest_muscles(q, k=3)
        hit = any(any(n in m.name.lower() for n in named) for m, _ in near)
        self.claims.append({
            "code": code, "kind": "miesien", "claimed": 1.0,
            "measured": 1.0 if hit else 0.0, "err": 0.0 if hit else 1.0,
            "detail": f"opis: {sorted(named)[0][:22]} | model: {near[0][0].name[:26]}"
            f" ({near[0][1] * 1000:.0f} mm)"})

    def between_pair(self, text: str, named: set):
        """Dwa miesnie z opisu 'between A and B', jesli oba znamy z siatki."""
        low = " " + re.sub(r"[^a-z ]", " ", text.lower()) + " "
        if " between " not in low:
            return None
        tail = low.split(" between ", 1)[1]
        hits = [n for n in named if f" {n} " in f" {tail} "]
        # "miedzy glowa mostkowa a obojczykowa MOS" to jeden miesien w dwoch
        # slowach - obie nazwy trafiaja w te sama siatke i pomiar nic nie znaczy.
        seen, out = [], []
        for n in sorted(hits):
            idx = frozenset(self.muscle_idx(n))
            if idx and not any(idx & s for s in seen):
                seen.append(idx)
                out.append(n)
        return tuple(out[:2]) if len(out) >= 2 else None

    def muscle_idx(self, name: str) -> list:
        idx = [m.index for m in self.ms.muscles
               if name in m.name.lower() and m.side in ("R", "r")]
        return idx or [m.index for m in self.ms.muscles if name in m.name.lower()]

    def muscle_gap(self, q: np.ndarray, name: str):
        """Najmniejsza odleglosc punktu od siatki miesnia o tej nazwie."""
        idx = self.muscle_idx(name)
        if not idx:
            return None
        sel = self.ms.rest[np.isin(self.ms._mid_points, idx)]
        return None if not len(sel) else float(np.linalg.norm(sel - q, axis=1).min())

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
