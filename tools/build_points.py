"""Zrodlo prawdy o polozeniu punktow -> data/points.json.

Kazdy punkt opisany jest kotwica proporcjonalna (patrz acu3d/anchors.py), a nie
sztywnymi wspolrzednymi. Tabele nizej sa recznie zakotwiczone wg opisow z
"A Manual of Acupuncture" (Deadman) - jednostka jest cun, punktami odniesienia
faldy, wyrostki kolczyste i linie wlosow.

Znacznik `accuracy`:
    anchored  - kotwica odwzorowuje opis podrecznikowy wprost (tulow, konczyny)
    approx    - model nie ma tej struktury (palce, malzowina, oczodol) -> polozenie
                pogladowe, do weryfikacji wzrokowej z opisem w panelu
Uruchomienie:  python tools/build_points.py
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
APPROX_SEGS = {"hand", "foot", "head"}

POINTS: list = []


# --- helpery kotwic ---------------------------------------------------------
def _pos(cun=None, t=None, frm="p0", z_frac=None, z_from=None, z_cun=0.0,
         vertebra=None, dz_cun=-0.4) -> dict:
    if cun is not None:
        return {"cun": cun, "from": frm}
    if t is not None:
        return {"t": t}
    if z_frac is not None:
        return {"z_frac": z_frac}
    if z_from is not None:
        return {"z_from": z_from, "z_cun": z_cun}
    if vertebra is not None:
        return {"vertebra": vertebra, "dz_cun": dz_cun}
    raise ValueError("brak polozenia osiowego")


def A(code: str, seg: str, az: float = None, lat: float = None, face: str = "anterior",
      out: float = 0.0, side: str = "bilateral", accuracy: str = None, note: str = "",
      **pos) -> None:
    """Dodaje punkt: segment + polozenie osiowe + azymut albo odleglosc boczna w cun."""
    anchor = {"seg": seg}
    anchor.update(_pos(**pos))
    if az is not None:
        anchor["az"] = az
    else:
        anchor["lat"] = lat
        anchor["face"] = face
    if out:
        anchor["out"] = out
    acc = accuracy or ("approx" if seg in APPROX_SEGS else "anchored")
    POINTS.append({"code": code, "anchor": anchor, "side": side,
                   "accuracy": acc, "note": note})


def SC(code: str, arc: float, lat: float = 0.0, side: str = "bilateral",
       note: str = "") -> None:
    """Punkt na sklepieniu czaszki: `arc` cun od przedniej linii wlosow po luku."""
    anchor = {"seg": "head", "arc_cun": arc}
    if lat:
        anchor["lat"] = lat
    POINTS.append({"code": code, "anchor": anchor,
                   "side": "midline" if lat == 0 and side == "midline" else side,
                   "accuracy": "anchored", "note": note})


# ============================================================================
# PLUCA (LU) - taiyin reki
# ============================================================================
A("LU-1", "trunk", lat=6, z_frac=0.772, note="1. przestrzen miedzyzebrowa, 6 cun od linii posrodkowej")
A("LU-2", "trunk", lat=6, z_frac=0.795, note="pod obojczykiem, w zaglebieniu podobojczykowym")
A("LU-3", "upper_arm", az=35, cun=3, note="3 cun ponizej przedniego faldu pachowego")
A("LU-4", "upper_arm", az=35, cun=4)
A("LU-5", "forearm", az=45, cun=0, note="fald lokciowy, promieniowo od sciegna dwuglowego")
A("LU-6", "forearm", az=55, cun=5, note="7 cun powyzej faldu nadgarstka")
A("LU-7", "forearm", az=70, cun=1.5, frm="p1")
A("LU-8", "forearm", az=72, cun=1, frm="p1")
A("LU-9", "forearm", az=68, cun=0.3, frm="p1", note="fald nadgarstka, przy tetnicy promieniowej")
A("LU-10", "hand", az=78, t=0.30, note="srodek I kosci srodrecza, strona dloniowa")
A("LU-11", "hand", az=95, t=0.52, note="kacik paznokcia kciuka od strony promieniowej")

# ============================================================================
# JELITO GRUBE (LI) - yangming reki
# ============================================================================
A("LI-1", "hand", az=100, t=0.62, note="kacik paznokcia palca wskazujacego")
A("LI-2", "hand", az=105, t=0.55)
A("LI-3", "hand", az=108, t=0.45)
A("LI-4", "hand", az=120, t=0.30, note="grzbiet dloni, miedzy I a II koscia srodrecza")
A("LI-5", "hand", az=125, t=0.06, note="tabakierka anatomiczna")
A("LI-6", "forearm", az=115, cun=3, frm="p1")
A("LI-7", "forearm", az=115, cun=5, frm="p1")
A("LI-8", "forearm", az=110, cun=4)
A("LI-9", "forearm", az=110, cun=3)
A("LI-10", "forearm", az=110, cun=2, note="2 cun ponizej LI-11")
A("LI-11", "forearm", az=118, cun=0, note="boczny koniec faldu lokciowego przy zgietym lokciu")
A("LI-12", "upper_arm", az=120, cun=1, frm="p1")
A("LI-13", "upper_arm", az=115, cun=3, frm="p1")
A("LI-14", "upper_arm", az=110, cun=6, frm="p1", note="przyczep miesnia naramiennego")
A("LI-15", "upper_arm", az=105, t=0.02, note="zaglebienie przednio-dolne wyrostka barkowego")
A("LI-16", "trunk", lat=4.2, z_frac=0.828, out=0.2, note="miedzy wyrostkiem barkowym a obojczykiem")
A("LI-17", "neck", az=62, z_frac=0.840)
A("LI-18", "neck", az=62, z_frac=0.845, note="na wysokosci chrzastki tarczowatej, miedzy glowami SCM")
A("LI-19", "head", az=14, z_frac=0.898)
A("LI-20", "head", az=18, z_frac=0.893, note="w bruzdzie nosowo-wargowej, przy skrzydelku nosa")

# ============================================================================
# ZOLADEK (ST) - yangming nogi
# ============================================================================
A("ST-1", "head", az=15, z_frac=0.918, note="miedzy galka oczna a dolnym brzegiem oczodolu")
A("ST-2", "head", az=15, z_frac=0.910, note="otwor podoczodolowy")
A("ST-3", "head", az=17, z_frac=0.899)
A("ST-4", "head", az=13, z_frac=0.889, note="kacik ust")
A("ST-5", "head", az=45, z_frac=0.878, note="przed miesniem zwaczem, na tetnicy twarzowej")
A("ST-6", "head", az=60, z_frac=0.881, note="wypuklosc zwacza przy zacisnietych zebach")
A("ST-7", "head", az=68, z_frac=0.900, note="pod lukiem jarzmowym")
A("ST-8", "head", az=40, z_frac=0.957, note="kacik czolowy linii wlosow")
A("ST-9", "neck", lat=1.5, z_frac=0.845, note="na wysokosci chrzastki tarczowatej, przy tetnicy szyjnej")
A("ST-10", "neck", lat=1.8, z_frac=0.828)
A("ST-11", "trunk", lat=2, z_frac=0.806)
A("ST-12", "trunk", lat=4, z_frac=0.812, note="srodek dolu nadobojczykowego")
A("ST-13", "trunk", lat=4, z_frac=0.800, note="pod obojczykiem")
A("ST-14", "trunk", lat=4, z_frac=0.775)
A("ST-15", "trunk", lat=4, z_frac=0.755)
A("ST-16", "trunk", lat=4, z_frac=0.737)
A("ST-17", "trunk", lat=4, z_frac=0.720, note="srodek brodawki sutkowej")
A("ST-18", "trunk", lat=4, z_frac=0.702)
for i, cun in zip(range(19, 31), [6, 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5]):
    A(f"ST-{i}", "trunk", lat=2, z_from="umbilicus", z_cun=cun,
      note=f"2 cun bocznie od linii posrodkowej, {abs(cun)} cun {'powyzej' if cun > 0 else 'ponizej'} pepka"
      if cun else "2 cun bocznie od pepka")
A("ST-31", "thigh", az=30, cun=0, note="na wysokosci dolnego brzegu spojenia lonowego")
A("ST-32", "thigh", az=20, cun=6, frm="p1", note="6 cun powyzej gornego brzegu rzepki")
A("ST-33", "thigh", az=25, cun=3, frm="p1")
A("ST-34", "thigh", az=22, cun=2, frm="p1")
A("ST-35", "shank", az=25, cun=0, note="zaglebienie bocznie od wiezadla rzepki")
A("ST-36", "shank", az=22, cun=3, note="3 cun ponizej ST-35, szerokosc palca bocznie od grzebienia piszczeli")
A("ST-37", "shank", az=22, cun=6)
A("ST-38", "shank", az=22, cun=8)
A("ST-39", "shank", az=22, cun=9)
A("ST-40", "shank", az=38, cun=8, note="2 szerokosci palca bocznie od ST-38")
A("ST-41", "foot", az=0, t=0.20, accuracy="approx",
  note="srodek faldu stawu skokowego; w modelu przesuniete nieco na grzbiet stopy, bo kostka i stopa zlewaja sie w jedna brylę")
A("ST-42", "foot", az=5, t=0.35)
A("ST-43", "foot", az=10, t=0.52)
A("ST-44", "foot", az=12, t=0.72)
A("ST-45", "foot", az=15, t=0.90, note="kacik paznokcia II palca stopy")

# ============================================================================
# SLEDZIONA (SP) - taiyin nogi
# ============================================================================
A("SP-1", "foot", az=-70, t=0.90, note="kacik paznokcia palucha od strony przysrodkowej")
A("SP-2", "foot", az=-75, t=0.80)
A("SP-3", "foot", az=-80, t=0.62)
A("SP-4", "foot", az=-80, t=0.45)
A("SP-5", "foot", az=-70, t=0.12, note="przednio-dolnie od kostki przysrodkowej")
A("SP-6", "shank", az=-80, cun=3, frm="p1", note="3 cun powyzej kostki przysrodkowej, przy brzegu piszczeli")
A("SP-7", "shank", az=-80, cun=6, frm="p1")
A("SP-8", "shank", az=-75, cun=3, note="3 cun ponizej SP-9")
A("SP-9", "shank", az=-80, cun=0.5, note="pod klykciem przysrodkowym piszczeli")
A("SP-10", "thigh", az=-45, cun=2, frm="p1", note="2 cun powyzej rzepki, wypuklosc miesnia obszernego przysrodkowego")
A("SP-11", "thigh", az=-55, cun=6, frm="p1")
A("SP-12", "trunk", lat=3.5, z_frac=0.492, note="fald pachwinowy, bocznie od tetnicy udowej")
A("SP-13", "trunk", lat=4, z_frac=0.512)
A("SP-14", "trunk", lat=4, z_from="umbilicus", z_cun=-1.3)
A("SP-15", "trunk", lat=4, z_from="umbilicus", z_cun=0, note="4 cun bocznie od pepka")
A("SP-16", "trunk", lat=4, z_from="umbilicus", z_cun=3)
A("SP-17", "trunk", lat=6, z_frac=0.702)
A("SP-18", "trunk", lat=6, z_frac=0.720)
A("SP-19", "trunk", lat=6, z_frac=0.737)
A("SP-20", "trunk", lat=6, z_frac=0.755)
A("SP-21", "trunk", lat=6, face="lateral", z_frac=0.735, note="linia pachowa srodkowa, 6. przestrzen miedzyzebrowa")

# ============================================================================
# SERCE (HT) - shaoyin reki
# ============================================================================
A("HT-1", "trunk", lat=6, face="lateral", z_frac=0.790, accuracy="approx",
  note="srodek dolu pachowego - model nie ma zachylka pachowego, punkt pokazany na bocznej scianie klatki")
A("HT-2", "upper_arm", az=-35, cun=3, frm="p1")
A("HT-3", "forearm", az=-55, cun=0, note="przysrodkowy koniec faldu lokciowego")
A("HT-4", "forearm", az=-60, cun=1.5, frm="p1")
A("HT-5", "forearm", az=-62, cun=1, frm="p1")
A("HT-6", "forearm", az=-62, cun=0.5, frm="p1")
A("HT-7", "forearm", az=-60, cun=0.3, frm="p1", note="fald nadgarstka, przy kosci grochowatej")
A("HT-8", "hand", az=-35, t=0.35, note="dlon, miedzy IV a V koscia srodrecza")
A("HT-9", "hand", az=-60, t=0.60, note="kacik paznokcia V palca")

# ============================================================================
# JELITO CIENKIE (SI) - taiyang reki
# ============================================================================
A("SI-1", "hand", az=-110, t=0.62, note="kacik paznokcia V palca od strony lokciowej")
A("SI-2", "hand", az=-115, t=0.50)
A("SI-3", "hand", az=-120, t=0.35, note="koniec faldu dloni przy zacisnietej piesci")
A("SI-4", "hand", az=-125, t=0.15)
A("SI-5", "forearm", az=-125, cun=0, frm="p1")
A("SI-6", "forearm", az=-130, cun=1, frm="p1")
A("SI-7", "forearm", az=-140, cun=5, frm="p1")
A("SI-8", "forearm", az=-150, cun=0, note="miedzy wyrostkiem lokciowym a nadklykciem przysrodkowym")
A("SI-9", "trunk", lat=5.5, face="posterior", z_frac=0.772, note="1 cun powyzej tylnego faldu pachowego")
A("SI-10", "trunk", lat=5, face="posterior", z_frac=0.800)
A("SI-11", "trunk", lat=4, face="posterior", z_frac=0.760, note="srodek lopatki")
A("SI-12", "trunk", lat=4, face="posterior", z_frac=0.795, note="srodek dolu nadgrzebieniowego")
A("SI-13", "trunk", lat=3, face="posterior", z_frac=0.800)
A("SI-14", "trunk", lat=3, face="posterior", vertebra="T1")
A("SI-15", "trunk", lat=2, face="posterior", vertebra="C7")
A("SI-16", "neck", az=115, z_frac=0.845)
A("SI-17", "neck", az=135, z_frac=0.862, note="za katem zuchwy, przy przednim brzegu SCM")
A("SI-18", "head", az=55, z_frac=0.912, note="pod kascia jarzmowa, pionowo pod bocznym katem oka")
A("SI-19", "head", az=78, z_frac=0.918, note="przed skrawkiem malzowiny, przy otwartych ustach")

# ============================================================================
# PECHERZ MOCZOWY (BL) - taiyang nogi
# ============================================================================
A("BL-1", "head", az=12, z_frac=0.926, note="przy przysrodkowym kaciku oka")
A("BL-2", "head", az=18, z_frac=0.933, note="przysrodkowy koniec brwi")
SC("BL-3", 0.5, lat=0.75, note="0,5 cun w obrebie linii wlosow, ok. 0,6 cun bocznie")
SC("BL-4", 0.5, lat=1.5, note="0,5 cun powyzej linii wlosow, 1,5 cun bocznie")
SC("BL-5", 1.0, lat=1.5)
SC("BL-6", 2.5, lat=1.5)
SC("BL-7", 4.0, lat=1.5)
SC("BL-8", 5.5, lat=1.5)
SC("BL-9", 9.5, lat=1.3, note="2,5 cun powyzej tylnej linii wlosow")
A("BL-10", "neck", az=165, z_frac=0.862, note="1,3 cun bocznie od DU-15, na miesniu czworobocznym")
# linia przykregoslupowa 1,5 cun
for code, vert in [("BL-11", "T1"), ("BL-12", "T2"), ("BL-13", "T3"), ("BL-14", "T4"),
                   ("BL-15", "T5"), ("BL-16", "T6"), ("BL-17", "T7"), ("BL-18", "T9"),
                   ("BL-19", "T10"), ("BL-20", "T11"), ("BL-21", "T12"), ("BL-22", "L1"),
                   ("BL-23", "L2"), ("BL-24", "L3"), ("BL-25", "L4"), ("BL-26", "L5"),
                   ("BL-27", "S1"), ("BL-28", "S2"), ("BL-29", "S3"), ("BL-30", "S4")]:
    A(code, "trunk", lat=1.5, face="posterior", vertebra=vert,
      note=f"1,5 cun bocznie od dolnego brzegu wyrostka kolczystego {vert}")
for code, vert in [("BL-31", "S1"), ("BL-32", "S2"), ("BL-33", "S3"), ("BL-34", "S4")]:
    A(code, "trunk", lat=0.8, face="posterior", vertebra=vert,
      note=f"otwor krzyzowy {vert[1]}")
A("BL-35", "trunk", lat=0.5, face="posterior", vertebra="coccyx", dz_cun=0.0,
  note="0,5 cun bocznie od szczytu kosci guzicznej")
A("BL-36", "thigh", az=180, cun=0, note="srodek faldu posladkowego")
A("BL-37", "thigh", az=180, cun=6)
A("BL-38", "thigh", az=175, cun=1, frm="p1", note="1 cun powyzej faldu podkolanowego")
A("BL-39", "shank", az=160, cun=0, note="boczny koniec faldu podkolanowego")
A("BL-40", "shank", az=180, cun=0, note="srodek faldu podkolanowego")
# linia przykregoslupowa 3 cun
for code, vert in [("BL-41", "T2"), ("BL-42", "T3"), ("BL-43", "T4"), ("BL-44", "T5"),
                   ("BL-45", "T6"), ("BL-46", "T7"), ("BL-47", "T9"), ("BL-48", "T10"),
                   ("BL-49", "T11"), ("BL-50", "T12"), ("BL-51", "L1"), ("BL-52", "L2"),
                   ("BL-53", "S2"), ("BL-54", "S4")]:
    A(code, "trunk", lat=3, face="posterior", vertebra=vert,
      note=f"3 cun bocznie od wyrostka kolczystego {vert}")
A("BL-55", "shank", az=180, cun=2)
A("BL-56", "shank", az=180, cun=5)
A("BL-57", "shank", az=180, cun=8, note="szczyt brzuscow miesnia brzuchatego lydki")
A("BL-58", "shank", az=168, cun=7, frm="p1")
A("BL-59", "shank", az=170, cun=3, frm="p1")
A("BL-60", "shank", az=150, cun=0, frm="p1", note="miedzy kostka boczna a sciegnem Achillesa")
A("BL-61", "foot", az=150, t=0.08)
A("BL-62", "foot", az=110, t=0.10, note="tuz ponizej kostki bocznej")
A("BL-63", "foot", az=115, t=0.20)
A("BL-64", "foot", az=115, t=0.35)
A("BL-65", "foot", az=115, t=0.55)
A("BL-66", "foot", az=118, t=0.75)
A("BL-67", "foot", az=120, t=0.92, note="kacik paznokcia V palca stopy")

# ============================================================================
# NERKI (KI) - shaoyin nogi
# ============================================================================
A("KI-1", "foot", az=180, t=0.42, note="podeszwa, na granicy 1/3 przedniej i 2/3 tylnych")
A("KI-2", "foot", az=-140, t=0.25)
A("KI-3", "shank", az=-140, cun=0, frm="p1", note="miedzy kostka przysrodkowa a sciegnem Achillesa")
A("KI-4", "foot", az=-160, t=0.05)
A("KI-5", "foot", az=-145, t=0.08)
A("KI-6", "foot", az=-110, t=0.03, note="tuz ponizej kostki przysrodkowej")
A("KI-7", "shank", az=-160, cun=2, frm="p1")
A("KI-8", "shank", az=-130, cun=2, frm="p1")
A("KI-9", "shank", az=-150, cun=5, frm="p1")
A("KI-10", "shank", az=-170, cun=0, note="przysrodkowy koniec faldu podkolanowego")
for i, cun in zip(range(11, 22), [-5, -4, -3, -2, -1, 0, 2, 3, 4, 5, 6]):
    A(f"KI-{i}", "trunk", lat=0.5, z_from="umbilicus", z_cun=cun,
      note="0,5 cun bocznie od linii posrodkowej brzucha")
A("KI-22", "trunk", lat=2, z_frac=0.702)
A("KI-23", "trunk", lat=2, z_frac=0.720)
A("KI-24", "trunk", lat=2, z_frac=0.737)
A("KI-25", "trunk", lat=2, z_frac=0.755)
A("KI-26", "trunk", lat=2, z_frac=0.775)
A("KI-27", "trunk", lat=2, z_frac=0.795, note="pod obojczykiem, 2 cun bocznie")

# ============================================================================
# OSIERDZIE (PC) - jueyin reki
# ============================================================================
A("PC-1", "trunk", lat=5, face="lateral", z_frac=0.745, note="1 cun bocznie od brodawki sutkowej")
A("PC-2", "upper_arm", az=-10, cun=2)
A("PC-3", "forearm", az=-8, cun=0, note="fald lokciowy, lokciowo od sciegna dwuglowego")
A("PC-4", "forearm", az=0, cun=5, frm="p1")
A("PC-5", "forearm", az=0, cun=3, frm="p1")
A("PC-6", "forearm", az=0, cun=2, frm="p1", note="2 cun powyzej faldu nadgarstka, miedzy sciegnami")
A("PC-7", "forearm", az=0, cun=0.3, frm="p1", note="srodek faldu nadgarstka")
A("PC-8", "hand", az=0, t=0.33, note="srodek dloni, miedzy II a III koscia srodrecza")
A("PC-9", "hand", az=0, t=0.95, note="szczyt III palca")

# ============================================================================
# POTROJNY OGRZEWACZ (SJ) - shaoyang reki
# ============================================================================
A("SJ-1", "hand", az=-160, t=0.68, note="kacik paznokcia IV palca od strony lokciowej")
A("SJ-2", "hand", az=-170, t=0.58)
A("SJ-3", "hand", az=185, t=0.42)
A("SJ-4", "forearm", az=180, cun=0.3, frm="p1", note="srodek faldu nadgarstka od strony grzbietowej")
A("SJ-5", "forearm", az=180, cun=2, frm="p1", note="2 cun powyzej SJ-4, miedzy koscia promieniowa a lokciowa")
A("SJ-6", "forearm", az=180, cun=3, frm="p1")
A("SJ-7", "forearm", az=195, cun=3, frm="p1")
A("SJ-8", "forearm", az=180, cun=4, frm="p1")
A("SJ-9", "forearm", az=180, cun=5)
A("SJ-10", "upper_arm", az=180, cun=1, frm="p1", note="1 cun powyzej wyrostka lokciowego")
A("SJ-11", "upper_arm", az=180, cun=2, frm="p1")
A("SJ-12", "upper_arm", az=180, cun=4, frm="p1")
A("SJ-13", "upper_arm", az=170, cun=5, frm="p1")
A("SJ-14", "upper_arm", az=150, t=0.02, note="tylno-dolnie od wyrostka barkowego")
A("SJ-15", "trunk", lat=3.2, face="posterior", z_frac=0.812)
A("SJ-16", "neck", az=150, z_frac=0.855)
A("SJ-17", "head", az=95, z_frac=0.895, note="za platkiem ucha")
A("SJ-18", "head", az=105, z_frac=0.905)
A("SJ-19", "head", az=110, z_frac=0.918)
A("SJ-20", "head", az=100, z_frac=0.938, note="nad szczytem malzowiny")
A("SJ-21", "head", az=82, z_frac=0.923, note="powyzej SI-19, przed wcięciem malzowiny")
A("SJ-22", "head", az=72, z_frac=0.928)
A("SJ-23", "head", az=38, z_frac=0.930, note="boczny koniec brwi")

# ============================================================================
# PECHERZYK ZOLCIOWY (GB) - shaoyang nogi
# ============================================================================
A("GB-1", "head", az=45, z_frac=0.925, note="0,5 cun bocznie od bocznego kacika oka")
A("GB-2", "head", az=80, z_frac=0.908)
A("GB-3", "head", az=75, z_frac=0.918, note="nad lukiem jarzmowym")
A("GB-4", "head", az=62, z_frac=0.942)
A("GB-5", "head", az=70, z_frac=0.945)
A("GB-6", "head", az=78, z_frac=0.942)
A("GB-7", "head", az=88, z_frac=0.938)
A("GB-8", "head", az=95, z_frac=0.950, note="1,5 cun powyzej szczytu malzowiny")
A("GB-9", "head", az=108, z_frac=0.950)
A("GB-10", "head", az=118, z_frac=0.942)
A("GB-11", "head", az=128, z_frac=0.932)
A("GB-12", "head", az=130, z_frac=0.902, note="w zaglebieniu za wyrostkiem sutkowatym")
SC("GB-13", 0.5, lat=3.0)
A("GB-14", "head", az=25, z_frac=0.947, note="1 cun powyzej srodka brwi")
SC("GB-15", 0.5, lat=2.25)
SC("GB-16", 1.5, lat=2.25)
SC("GB-17", 2.5, lat=2.25)
SC("GB-18", 4.0, lat=2.25)
SC("GB-19", 8.5, lat=2.25)
A("GB-20", "neck", az=145, z_frac=0.868, note="pod potylica, miedzy SCM a czworobocznym")
A("GB-21", "trunk", lat=5.5, face="posterior", z_frac=0.827, out=0.3,
  note="szczyt barku, w polowie miedzy DU-14 a wyrostkiem barkowym")
A("GB-22", "trunk", lat=6.5, face="lateral", z_frac=0.752, note="3 cun ponizej dolu pachowego")
A("GB-23", "trunk", lat=6.0, face="lateral", z_frac=0.745)
A("GB-24", "trunk", lat=4, z_frac=0.700, note="7. przestrzen miedzyzebrowa")
A("GB-25", "trunk", lat=5.5, face="lateral", z_frac=0.660, note="wolny koniec XII zebra")
A("GB-26", "trunk", lat=6, face="lateral", z_from="umbilicus", z_cun=0)
A("GB-27", "trunk", lat=5, z_frac=0.560, note="przy kolcu biodrowym przednim gornym")
A("GB-28", "trunk", lat=5.5, z_frac=0.548)
A("GB-29", "trunk", lat=6.5, face="lateral", z_frac=0.530,
  note="miedzy kolcem biodrowym przednim gornym a krętarzem wiekszym")
A("GB-30", "thigh", az=160, t=0.02, note="posladek, 1/3 odleglosci od krętarza do rozworu krzyzowego")
A("GB-31", "thigh", az=90, cun=7, frm="p1", note="tam, gdzie siega czubek srodkowego palca zwisajacej reki")
A("GB-32", "thigh", az=90, cun=5, frm="p1")
A("GB-33", "thigh", az=95, cun=0.5, frm="p1")
A("GB-34", "shank", az=80, cun=1, note="w zaglebieniu przednio-dolnie od glowy strzalki")
A("GB-35", "shank", az=100, cun=7, frm="p1")
A("GB-36", "shank", az=85, cun=7, frm="p1")
A("GB-37", "shank", az=88, cun=5, frm="p1")
A("GB-38", "shank", az=88, cun=4, frm="p1")
A("GB-39", "shank", az=90, cun=3, frm="p1", note="3 cun powyzej kostki bocznej")
A("GB-40", "foot", az=95, t=0.12, note="przednio-dolnie od kostki bocznej")
A("GB-41", "foot", az=60, t=0.50)
A("GB-42", "foot", az=60, t=0.62)
A("GB-43", "foot", az=65, t=0.80)
A("GB-44", "foot", az=80, t=0.92, note="kacik paznokcia IV palca stopy")

# ============================================================================
# WATROBA (LR) - jueyin nogi
# ============================================================================
A("LR-1", "foot", az=-20, t=0.90, note="kacik paznokcia palucha od strony bocznej")
A("LR-2", "foot", az=-40, t=0.78)
A("LR-3", "foot", az=-30, t=0.60, note="zaglebienie miedzy I a II koscia srodstopia")
A("LR-4", "foot", az=-55, t=0.10, note="przed kostka przysrodkowa")
A("LR-5", "shank", az=-95, cun=5, frm="p1")
A("LR-6", "shank", az=-95, cun=7, frm="p1")
A("LR-7", "shank", az=-85, cun=5)
A("LR-8", "shank", az=-140, cun=0, note="przysrodkowy koniec faldu podkolanowego przy zgietym kolanie")
A("LR-9", "thigh", az=-110, cun=4, frm="p1")
A("LR-10", "thigh", az=-95, cun=3)
A("LR-11", "thigh", az=-100, cun=3.5)
A("LR-12", "trunk", lat=2.5, z_frac=0.495, note="fald pachwinowy")
A("LR-13", "trunk", lat=5.5, face="lateral", z_frac=0.640, note="wolny koniec XI zebra")
A("LR-14", "trunk", lat=4, z_frac=0.700, note="6. przestrzen miedzyzebrowa, pionowo pod brodawka")

# ============================================================================
# NACZYNIE POCZECIA (REN) - linia posrodkowa przednia
# ============================================================================
A("REN-1", "trunk", az=0, z_frac=0.452, side="midline", note="krocze")
for i, (frm, cun) in enumerate([("pubis", 0), ("pubis", 1), ("pubis", 2), ("pubis", 3)], start=2):
    A(f"REN-{i}", "trunk", lat=0, z_from=frm, z_cun=cun, side="midline",
      note=f"{cun} cun powyzej spojenia lonowego" if cun else "gorny brzeg spojenia lonowego")
for i, cun in zip(range(6, 16), [-1.5, -1, 0, 1, 2, 3, 4, 5, 6, 7]):
    A(f"REN-{i}", "trunk", lat=0, z_from="umbilicus", z_cun=cun, side="midline",
      note=("srodek pepka" if cun == 0 else
            f"{abs(cun)} cun {'powyzej' if cun > 0 else 'ponizej'} pepka"))
A("REN-16", "trunk", lat=0, z_frac=0.700, side="midline")
A("REN-17", "trunk", lat=0, z_frac=0.720, side="midline", note="na mostku, na wysokosci brodawek")
A("REN-18", "trunk", lat=0, z_frac=0.737, side="midline")
A("REN-19", "trunk", lat=0, z_frac=0.755, side="midline")
A("REN-20", "trunk", lat=0, z_frac=0.772, side="midline")
A("REN-21", "trunk", lat=0, z_frac=0.788, side="midline")
A("REN-22", "trunk", lat=0, z_frac=0.802, side="midline", note="wciecie szyjne mostka")
A("REN-23", "neck", az=0, z_frac=0.850, side="midline", note="powyzej kosci gnykowej")
A("REN-24", "head", az=0, z_frac=0.879, side="midline", note="bruzda bródkowo-wargowa")

# ============================================================================
# NACZYNIE GUBERNATORA (DU) - linia posrodkowa tylna
# ============================================================================
A("DU-1", "trunk", lat=0, face="posterior", z_frac=0.472, side="midline",
  note="miedzy szczytem kosci guzicznej a odbytem")
A("DU-2", "trunk", lat=0, face="posterior", z_frac=0.512, side="midline", note="rozwor krzyzowy")
for code, vert in [("DU-3", "L4"), ("DU-4", "L2"), ("DU-5", "L1"), ("DU-6", "T11"),
                   ("DU-7", "T10"), ("DU-8", "T9"), ("DU-9", "T7"), ("DU-10", "T6"),
                   ("DU-11", "T5"), ("DU-12", "T3"), ("DU-13", "T1"), ("DU-14", "C7")]:
    A(code, "trunk", lat=0, face="posterior", vertebra=vert, side="midline",
      note=f"ponizej wyrostka kolczystego {vert}")
A("DU-15", "neck", az=180, z_frac=0.858, side="midline", note="0,5 cun powyzej tylnej linii wlosow")
A("DU-16", "neck", az=180, z_frac=0.870, side="midline", note="pod guzowatoscia potyliczna zewnetrzna")
SC("DU-17", 9.5, side="midline")
SC("DU-18", 8.0, side="midline")
SC("DU-19", 6.5, side="midline")
SC("DU-20", 5.0, side="midline", note="5 cun za przednia linia wlosow, na szczycie glowy")
SC("DU-21", 3.5, side="midline")
SC("DU-22", 2.0, side="midline")
SC("DU-23", 1.0, side="midline")
SC("DU-24", 0.5, side="midline")
A("DU-25", "head", az=0, z_frac=0.906, side="midline", note="czubek nosa")
A("DU-26", "head", az=0, z_frac=0.897, side="midline", note="gorna 1/3 rowka nosowo-wargowego")
A("DU-27", "head", az=0, z_frac=0.891, side="midline")
A("DU-28", "head", az=0, z_frac=0.887, side="midline", note="wedzidelko wargi gornej")


# ============================================================================
# Kanaly: nazwy, zywiol, kolor, kolejnosc w obiegu dobowym
# ============================================================================
CHANNELS = {
    "LU":  {"pl": "Pluca",              "en": "Lung",             "zh": "手太陰肺经", "element": "Metal", "color": "#e3e8f0", "clock": 1},
    "LI":  {"pl": "Jelito grube",       "en": "Large Intestine",  "zh": "手阳明大肠经", "element": "Metal", "color": "#b9c4d4", "clock": 2},
    "ST":  {"pl": "Zoladek",            "en": "Stomach",          "zh": "足阳明胃经", "element": "Earth", "color": "#e0a83c", "clock": 3},
    "SP":  {"pl": "Sledziona",          "en": "Spleen",           "zh": "足太陰脾经", "element": "Earth", "color": "#c98b2e", "clock": 4},
    "HT":  {"pl": "Serce",              "en": "Heart",            "zh": "手少陰心经", "element": "Fire",  "color": "#e04a44", "clock": 5},
    "SI":  {"pl": "Jelito cienkie",     "en": "Small Intestine",  "zh": "手太阳小肠经", "element": "Fire",  "color": "#b83b36", "clock": 6},
    "BL":  {"pl": "Pecherz moczowy",    "en": "Bladder",          "zh": "足太阳膡胱经", "element": "Water", "color": "#4d8fd6", "clock": 7},
    "KI":  {"pl": "Nerki",              "en": "Kidney",           "zh": "足少陰肾经", "element": "Water", "color": "#2f5f9e", "clock": 8},
    "PC":  {"pl": "Osierdzie",          "en": "Pericardium",      "zh": "手厥陰心包经", "element": "Fire",  "color": "#d977b8", "clock": 9},
    "SJ":  {"pl": "Potrojny ogrzewacz", "en": "San Jiao",         "zh": "手少阳三焦经", "element": "Fire",  "color": "#a95a9c", "clock": 10},
    "GB":  {"pl": "Pecherzyk zolciowy", "en": "Gall Bladder",     "zh": "足少阳胆经", "element": "Wood",  "color": "#57c08a", "clock": 11},
    "LR":  {"pl": "Watroba",            "en": "Liver",            "zh": "足厥陰肝经", "element": "Wood",  "color": "#2f9463", "clock": 12},
    "REN": {"pl": "Naczynie poczecia",  "en": "Conception vessel", "zh": "任脉", "element": "-", "color": "#f0c674", "clock": 13},
    "DU":  {"pl": "Naczynie gubernatora", "en": "Governing vessel", "zh": "督脉", "element": "-", "color": "#c6d3e4", "clock": 14},
}

# Przebieg kanalu do rysowania - domyslnie kolejnosc numeryczna, tu wyjatki
PATH_OVERRIDES = {
    "BL": [
        [f"BL-{i}" for i in range(1, 36)],
        [f"BL-{i}" for i in range(41, 55)],
        ["BL-36", "BL-37", "BL-38", "BL-40"] + [f"BL-{i}" for i in range(55, 68)],
    ],
}


def build() -> None:
    meta = json.loads((ROOT / "data" / "point_meta.json").read_text(encoding="utf-8"))
    by_code = {}
    out = []
    for p in POINTS:
        code = p["code"]
        if code in by_code:
            raise SystemExit(f"duplikat kotwicy: {code}")
        ch, idx = code.split("-")
        m = meta.get(code, {})
        rec = {
            "code": code,
            "channel": ch,
            "index": int(idx),
            "pinyin": m.get("pinyin", ""),
            "english": m.get("english", ""),
            "categories": m.get("categories", []),
            "side": p["side"],
            "accuracy": p["accuracy"],
            "anchor": p["anchor"],
            "note_pl": p["note"],
        }
        by_code[code] = rec
        out.append(rec)

    missing = sorted(set(meta) - set(by_code), key=lambda c: (c.split("-")[0], int(c.split("-")[1])))
    extra = sorted(set(by_code) - set(meta))
    doc = {
        "version": 1,
        "coordinate_system": "+X = prawa strona, +Y = przod, +Z = gora; kotwice w cun",
        "channels": CHANNELS,
        "points": sorted(out, key=lambda r: (CHANNELS[r["channel"]]["clock"], r["index"])),
    }
    (ROOT / "data" / "points.json").write_text(
        json.dumps(doc, ensure_ascii=False, indent=1), encoding="utf-8")

    paths = {}
    for ch in CHANNELS:
        if ch in PATH_OVERRIDES:
            paths[ch] = PATH_OVERRIDES[ch]
        else:
            codes = [r["code"] for r in sorted(
                (x for x in out if x["channel"] == ch), key=lambda r: r["index"])]
            paths[ch] = [codes]
    (ROOT / "data" / "meridians.json").write_text(
        json.dumps({"paths": paths}, ensure_ascii=False, indent=1), encoding="utf-8")

    n_appr = sum(1 for r in out if r["accuracy"] == "approx")
    print(f"punktow: {len(out)}  (anchored: {len(out) - n_appr}, approx: {n_appr})")
    if missing:
        print(f"BRAK kotwic dla {len(missing)}: {', '.join(missing)}")
    if extra:
        print(f"kody spoza zrodla: {', '.join(extra)}")


if __name__ == "__main__":
    build()
