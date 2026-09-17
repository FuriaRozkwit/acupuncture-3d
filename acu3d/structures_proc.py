"""Proceduralne odpowiedniki drobnych struktur anatomicznych.

Tryb `anatomy` bierze paliczki, kosci srodrecza czy kosc jarzmowa wprost z siatek
Z-Anatomy (`Anatomy.struct`). Tryb proceduralny takich siatek nie ma, a bez nich
89 z 361 punktow nie da sie rozwiazac - to punkty opisane nie w cun, tylko przez
anatomie: kacik paznokcia, przestrzen miedzy I a II koscia srodrecza, okolica
przed skrawkiem ucha.

Ten modul odtwarza te struktury z samego rigu. Kazda opisana jest w lokalnej ramce
swojego segmentu, w wielkosciach wzglednych, wiec skaluje sie razem z wzrostem:

    t       zakres wzdluz osi segmentu (0 = koniec blizszy, 1 = dalszy);
            wartosci spoza [0, 1] sa dozwolone - pieta lezy za poczatkiem stopy
    lat     przesuniecie w bok, w ulamkach polosi rx przekroju (+ = na zewnatrz)
    ant     przesuniecie do przodu, w ulamkach polosi ry (+ = do przodu)
    axis    zwrot osi struktury; domyslnie os segmentu

Zwrot osi ma znaczenie tylko tam, gdzie kotwica uzywa `along` innego niz 0.5, a nie
podaje `region` (dla dloni i stopy zwrot koryguje `frame()["long"]`). Dotyczy to
czterech struktur i przy kazdej stoi komentarz, z czego wynika kierunek.

Model jest stylizowany: proceduralna dlon jest jedna bryla, a glowa nie ma malzowiny.
Struktury palcow i ucha sa wiec umownymi miejscami na powierzchni, nie ksztaltami
kosci - maja dawac poprawna topografie punktu, nie wierna anatomie.
"""
from __future__ import annotations

import numpy as np

SPEC: dict[str, dict] = {}


def _hand_digits() -> None:
    """Kosci srodrecza i paliczki dloni. Kciuk jest krotszy i lezy bardziej dloniowo."""
    lat = {1: 0.80, 2: 0.42, 3: 0.02, 4: -0.36, 5: -0.72}
    mc = {1: (0.05, 0.36), 2: (0.10, 0.50), 3: (0.10, 0.51), 4: (0.10, 0.49), 5: (0.10, 0.46)}
    pp = {1: (0.38, 0.60), 2: (0.52, 0.72), 3: (0.53, 0.75), 4: (0.52, 0.72), 5: (0.48, 0.66)}
    dp = {1: (0.62, 0.76), 2: (0.86, 0.98), 3: (0.88, 1.00), 4: (0.85, 0.96), 5: (0.80, 0.90)}
    names = {1: "first", 2: "second", 3: "third", 4: "fourth", 5: "fifth"}
    ordn = {1: "First", 2: "Second", 3: "Third", 4: "Fourth", 5: "Fifth"}
    for i in range(1, 6):
        SPEC[ordn[i] + " metacarpal bone"] = dict(
            seg="hand", t=mc[i], lat=lat[i], ant=(0.30 if i == 1 else 0.05))
        SPEC["Proximal phalanx of " + names[i] + " finger of hand"] = dict(
            seg="hand", t=pp[i], lat=lat[i] * 1.05, ant=(0.25 if i == 1 else 0.0))
        SPEC["Distal phalanx of " + names[i] + " finger of hand"] = dict(
            seg="hand", t=dp[i], lat=lat[i] * 1.10, ant=(0.20 if i == 1 else 0.0))


def _foot_digits() -> None:
    """Kosci srodstopia i paliczki stopy. Paluch (I) jest przysrodkowy i dluzszy."""
    lat = {1: -0.72, 2: -0.34, 3: 0.02, 4: 0.36, 5: 0.70}
    mt = {1: (0.30, 0.66), 2: (0.32, 0.70), 3: (0.32, 0.70), 4: (0.32, 0.69), 5: (0.30, 0.66)}
    pp = {1: (0.68, 0.88), 2: (0.72, 0.86), 3: (0.72, 0.86), 4: (0.71, 0.84), 5: (0.70, 0.82)}
    dp = {1: (0.88, 1.00), 2: (0.88, 0.98), 3: (0.87, 0.97), 4: (0.86, 0.95), 5: (0.84, 0.93)}
    names = {1: "first", 2: "second", 3: "third", 4: "fourth", 5: "fifth"}
    ordn = {1: "First", 2: "Second", 3: "Third", 4: "Fourth", 5: "Fifth"}
    for i in range(1, 6):
        SPEC[ordn[i] + " metatarsal bone"] = dict(
            seg="foot", t=mt[i], lat=lat[i], ant=-0.10)
        SPEC["Proximal phalanx of " + names[i] + " finger of foot"] = dict(
            seg="foot", t=pp[i], lat=lat[i] * 1.05, ant=-0.15)
        SPEC["Distal phalanx of " + names[i] + " finger of foot"] = dict(
            seg="foot", t=dp[i], lat=lat[i] * 1.10, ant=-0.15)


_hand_digits()
_foot_digits()

SPEC.update({
    # --- nadgarstek i step
    "Scaphoid bone":  dict(seg="hand", t=(0.00, 0.08), lat=0.50, ant=-0.15),
    "Talus":          dict(seg="foot", t=(0.00, 0.12), lat=0.00, ant=0.10),
    "Calcaneus":      dict(seg="foot", t=(-0.34, -0.02), lat=0.00, ant=-0.10),
    "Navicular bone": dict(seg="foot", t=(0.16, 0.28), lat=-0.42, ant=0.00),
    "Cuboid bone":    dict(seg="foot", t=(0.20, 0.34), lat=0.46, ant=-0.05),

    # --- podudzie; obie kosci koncza sie kostka, wiec zakres siega konca segmentu
    "Tibia":  dict(seg="shank", t=(0.02, 1.00), lat=-0.34, ant=0.34),
    "Fibula": dict(seg="shank", t=(0.08, 1.00), lat=0.56, ant=-0.05),

    # --- czaszka
    # Os w dol: kotwice traktuja along=1.0 jako dolny koniec kosci nosowej,
    # bo stamtad schodza w dol do czubka nosa (DU-25) i rynienki (DU-26..28).
    "Nasal bone": dict(seg="head", t=(0.30, 0.52), lat=0.00, ant=0.62,
                       axis=(0.0, 0.0, -1.0)),
    "Anterior segment of eyeball": dict(seg="head", t=(0.41, 0.45), lat=0.42, ant=0.60,
                                        axis=(0.0, 1.0, 0.0)),
    "Frontal bone": dict(seg="head", t=(0.50, 0.62), lat=0.00, ant=0.55),
    # Os w lewo: along < 0.5 ma wypadac po prawej stronie (ST-5, ST-6),
    # a 0.5 w linii posrodkowej (REN-24).
    "Mandible": dict(seg="head", t=(0.06, 0.14), lat=0.00, ant=0.50,
                     axis=(-1.0, 0.0, 0.0), length=0.075),
    # Os do przodu: along=0.25 ma byc tylna czescia luku, przy uchu (ST-7).
    "Zygomatic bone": dict(seg="head", t=(0.34, 0.40), lat=0.72, ant=0.30,
                           axis=(0.0, 1.0, 0.0), length=0.045),
    "Temporal bone": dict(seg="head", t=(0.35, 0.41), lat=0.90, ant=-0.05),

    # --- malzowina; model proceduralny jej nie ma, wiec sa to miejsca na bocznej
    #     powierzchni glowy na wysokosci ucha
    # Os w gore: SJ-18 (along 0.25) lezy nizej niz SJ-19 (along 0.55).
    "Helix": dict(seg="head", t=(0.33, 0.50), lat=1.00, ant=-0.15, axis=(0.0, 0.0, 1.0)),
    "Apex of auricle":           dict(seg="head", t=(0.46, 0.50), lat=1.00, ant=-0.10),
    "Tragus":                    dict(seg="head", t=(0.36, 0.39), lat=0.95, ant=0.05),
    "Antitragus":                dict(seg="head", t=(0.34, 0.37), lat=0.95, ant=-0.03),
    "Lobule of auricle":         dict(seg="head", t=(0.28, 0.33), lat=0.90, ant=-0.05),
    "Anterior notch of auricle": dict(seg="head", t=(0.37, 0.40), lat=0.95, ant=0.08),
})


def _unit(v: np.ndarray) -> np.ndarray:
    n = float(np.linalg.norm(v))
    return v / n if n > 1e-12 else v


class ProcStruct:
    """Struktura wyliczona z rigu. Interfejs zgodny z `anatomy.StructRef`."""

    def __init__(self, name: str, centroid: np.ndarray, axis: np.ndarray, extent: tuple):
        self.name = name
        self.centroid = np.asarray(centroid, dtype=float)
        self.axis = _unit(np.asarray(axis, dtype=float))
        self.extent = (float(extent[0]), float(extent[1]))

    def at(self, along: float, ref_dir: np.ndarray = None) -> np.ndarray:
        axis, (lo, hi) = self.axis, self.extent
        if ref_dir is not None and float(np.dot(axis, ref_dir)) < 0:
            axis, lo, hi = -axis, -hi, -lo
        return self.centroid + axis * (lo + (hi - lo) * float(np.clip(along, 0.0, 1.0)))

    def extreme(self, direction: np.ndarray) -> np.ndarray:
        d = _unit(np.asarray(direction, dtype=float))
        ends = [self.at(0.0), self.at(1.0)]
        return ends[int(np.argmax([float(e @ d) for e in ends]))]


def build(rig, name: str) -> ProcStruct:
    """Buduje strukture z tabeli SPEC w ukladzie swiata (prawa strona ciala)."""
    spec = SPEC.get(name)
    if spec is None:
        raise KeyError("nie znam struktury proceduralnej: " + name)
    seg = rig.seg(spec["seg"])
    t0, t1 = spec["t"]
    rx, ry, dy = seg.radii(0.5 * (t0 + t1))

    def place(t: float) -> np.ndarray:
        base = seg.p0 + (seg.p1 - seg.p0) * t          # bez obcinania do [0, 1]
        return (base + seg.e_lat * (spec["lat"] * rx)
                + seg.e_ant * (spec["ant"] * ry + dy))

    a, b = place(t0), place(t1)
    centroid = 0.5 * (a + b)
    axis = np.array(spec["axis"], dtype=float) if "axis" in spec else seg.e_ax
    if "length" in spec:
        length = float(spec["length"]) * rig.H / 1.72   # tabela w metrach dla H = 1.72
    else:
        length = float(np.linalg.norm(b - a))
    return ProcStruct(name, centroid, axis, (-0.5 * length, 0.5 * length))
