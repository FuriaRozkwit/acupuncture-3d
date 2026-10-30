"""Czy `english` i `pinyin` w data/points.json zgadzaja sie z ksiazka (SPRAWA-007).

Wzor poprawnosci lezy poza repozytorium, w `ksiazka/` projektu fiszki: naglowki
rozdzialow 05-18 oraz indeksy 25, 26 i 29. Sciezke zmienia ACU3D_KSIAZKA. Brak
tych plikow to blad, nie zaliczenie: arbiter bez wzoru niczego nie dowodzi.

Nazwa jest potwierdzona, gdy:
  1. zaczyna sie od niej naglowek punktu w rozdziale (to jedyne zrodlo, ktore
     odroznia nazwe glowna od alternatywnej: indeksy podaja obie, a kursywa
     zniknela przy OCR), albo
  2. naglowek jest nieczytelny, a oba indeksy podaja pod tym kodem te jedna nazwe
     (gdy nazw jest wiecej, indeksy nie wskazuja glownej i musi rozstrzygnac naglowek), albo
  3. jest w WYJATKI z podanym powodem.

    python tools/check_names.py
    python tools/check_names.py --mutacje    # samokontrola: psuje punkty, liczy, co przeszlo
"""
from __future__ import annotations

import collections
import difflib
import json
import os
import random
import re
import sys
import unicodedata
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
KSIAZKA = Path(os.environ.get("ACU3D_KSIAZKA", ROOT.parent / "fiszki" / "ksiazka"))
BUNDLES = ("meridian-bundle.json", "meridian-bundle-anatomy.json")
ROZDZIALY = {"05": "LU", "06": "LI", "07": "ST", "08": "SP", "09": "HT", "10": "SI", "11": "BL",
             "12": "KI", "13": "PC", "14": "SJ", "15": "GB", "16": "LR", "17": "REN", "18": "DU"}
ALIAS = {"KID": "KI", "HE": "HT", "P": "PC", "LIV": "LR", "LL": "LI"}
KANAL = r"(L\.?[IL1]\.?|[A-Z]{1,3})"
NAGLOWEK = re.compile(rf"^([A-Z][A-Z' ]{{2,}})\s+{KANAL}\s*-\s*(\d+)\b")
INDEKS_25 = re.compile(rf"^(\S[^\n(]*?)\s+{KANAL}\s*-\s*(\d+)\s*\(([^)\n]+)\)", re.M)
INDEKS_26 = re.compile(rf"^(.+?)\s*\(([^()\n]*?)\s+{KANAL}\s*-\s*(\d+)\)", re.M)
INDEKS_29 = re.compile(rf"^{KANAL}\s*-\s*(\d+)\s+(\S+)", re.M)
ZASIEG = 5   # ile niepustych linii pod naglowkiem moze nalezec do nazwy
PODOBIENSTWO = 0.85   # pinyin w naglowku zepsuty przez OCR (ZHONGJf, BAIHUi)

# Ksiazka nie rozstrzyga tych pol sama z powodu OCR; kazdy wpis mowi dlaczego.
WYJATKI_ENGLISH: dict = {}
# Tu ksiazka przeczy sama sobie, wiec zgoda indeksow niczego nie dowodzi.
TYLKO_NAGLOWEK = {"BL-39": "naglowek: Outside of the Crook, oba indeksy: Lateral to the Crook"}
WYJATKI_PINYIN = {
    "GB-15": ("Toulinqi", "OCR: indeksy maja Toulingi i Tulfnqi"),
    "GB-41": ("Zulinqi", "OCR: indeksy maja Zulingi i Ztltnqi"),
}


def kod(kanal: str, numer: str) -> str:
    k = kanal.replace(".", "").replace("1", "I").strip()
    return f"{ALIAS.get(k, k)}-{int(numer)}"


def klucz(text: str) -> str:
    text = unicodedata.normalize("NFD", text or "")
    return re.sub(r"[^a-z]", "", "".join(c for c in text if not unicodedata.combining(c)).lower())


def nazwa_z_linii(linia: str, nazwa: str) -> bool:
    """Linia zaczyna sie od calej nazwy, a to, co dalej, jest szumem OCR."""
    start = klucz(nazwa)
    return bool(start) and klucz(" ".join(linia.split()[:len(nazwa.split())])) == start


def linie_naglowka(linie: list, pinyin: str) -> list:
    """Numery linii, od ktorych zaczyna sie punkt. OCR psuje kod kanalu (`SjJ-2`, `GB`
    urwane na koncu wiersza), wiec punkt szukany jest po pinyinie z indeksu 29: pierwsze
    slowo wiersza pisane wielkimi literami. Dopasowanie scisle ma pierwszenstwo."""
    wzor = klucz(pinyin)
    scisle, bliskie = [], []
    for i, linia in enumerate(linie):
        slowa = linia.split()
        if not slowa or sum(ch.isupper() for ch in slowa[0]) < 0.7 * len(slowa[0]):
            continue
        k = klucz(slowa[0])
        if k == wzor:
            scisle.append(i)
        elif difflib.SequenceMatcher(None, k, wzor).ratio() >= PODOBIENSTWO:
            bliskie.append(i)
    return scisle or bliskie


def wczytaj_ksiazke() -> dict:
    naglowki = collections.defaultdict(list)
    pinyin = collections.defaultdict(set)
    tekst = (KSIAZKA / "29-index-point-numbers.md").read_text(encoding="utf-8-sig")
    numery = {kod(m.group(1), m.group(2)): m.group(3) for m in INDEKS_29.finditer(tekst)}
    for numer, kanal in ROZDZIALY.items():
        plik = next(KSIAZKA.glob(f"{numer}-*.md"))
        linie = [x.strip() for x in plik.read_text(encoding="utf-8-sig").splitlines()]
        starty = collections.defaultdict(set)
        for i, linia in enumerate(linie):
            m = NAGLOWEK.match(linia)
            if m and kod(m.group(2), m.group(3)).startswith(kanal + "-"):
                starty[kod(m.group(2), m.group(3))].add(i)
                pinyin[kod(m.group(2), m.group(3))].add(klucz(m.group(1)))
        for c, py in numery.items():
            if c.startswith(kanal + "-"):
                starty[c].update(linie_naglowka(linie, py))
        for c, wiersze in starty.items():
            for i in sorted(wiersze):
                pod = [x for x in linie[i + 1:i + 12] if x and not x.startswith("#")]
                naglowki[c].append(pod[:ZASIEG])

    i25 = collections.defaultdict(set)
    tekst = (KSIAZKA / "25-index-chinese-point-names.md").read_text(encoding="utf-8-sig")
    for m in INDEKS_25.finditer(tekst):
        c = kod(m.group(2), m.group(3))
        i25[c].add(klucz(m.group(4)))
        pinyin[c].add(klucz(m.group(1)))
    i26 = collections.defaultdict(set)
    tekst = (KSIAZKA / "26-index-english-point-names.md").read_text(encoding="utf-8-sig")
    for m in INDEKS_26.finditer(tekst):
        i26[kod(m.group(3), m.group(4))].add(klucz(m.group(1)))
    tekst = (KSIAZKA / "29-index-point-numbers.md").read_text(encoding="utf-8-sig")
    for m in INDEKS_29.finditer(tekst):
        pinyin[kod(m.group(1), m.group(2))].add(klucz(m.group(3)))
    return {"naglowki": naglowki, "i25": i25, "i26": i26, "pinyin": pinyin}


def ocen_nazwe(point: dict, wzor: dict) -> str:
    c, nazwa = point["code"], point.get("english") or ""
    if c in WYJATKI_ENGLISH and WYJATKI_ENGLISH[c][0] == nazwa:
        return "wyjatek"
    if any(nazwa_z_linii(x, nazwa) for pod in wzor["naglowki"].get(c, []) for x in pod):
        return "naglowek"
    if c not in TYLKO_NAGLOWEK and {klucz(nazwa)} == wzor["i25"].get(c) == wzor["i26"].get(c):
        return "oba_indeksy"
    return "BLAD"


def sprawdz(points: list, wzor: dict) -> tuple:
    bledy, poziomy = [], collections.Counter()
    for p in points:
        poziom = ocen_nazwe(p, wzor)
        poziomy[poziom] += 1
        if poziom == "BLAD":
            bledy.append(f"data/points.json {p['code']}: english = {p.get('english')!r} bez potwierdzenia")
        c, py = p["code"], p.get("pinyin") or ""
        if c in WYJATKI_PINYIN and WYJATKI_PINYIN[c][0] == py:
            continue
        if klucz(py) not in wzor["pinyin"].get(c, ()):
            bledy.append(f"data/points.json {c}: pinyin = {py!r} bez potwierdzenia")
    return bledy, poziomy


def sprawdz_paczki(points: list) -> list:
    """Paczka dla Terrapki niesie te same nazwy co data/points.json, dla obu stron ciala."""
    wzor = {p["code"]: p.get("english") for p in points}
    bledy = []
    for name in BUNDLES:
        path = ROOT / "dist" / name
        if not path.exists():
            continue
        for p in json.loads(path.read_text(encoding="utf-8"))["points"]:
            if p["code"] not in wzor:
                bledy.append(f"dist/{name} {p['code']}: kodu nie ma w data/points.json")
            elif p.get("english") != wzor[p["code"]]:
                bledy.append(f"dist/{name} {p['code']}: english = {p.get('english')!r}")
    return bledy


def mutacje(points: list, wzor: dict) -> int:
    """Czy arbiter widzi blad: psuje kazdy punkt i liczy, co przeszlo. Dane zostaja nietkniete."""
    def wykryte(point: dict, pole: str, wartosc: str) -> bool:
        return bool(sprawdz([{**point, pole: wartosc}], wzor)[0])

    los = random.Random(30)
    tekst = (KSIAZKA / "26-index-english-point-names.md").read_text(encoding="utf-8-sig")
    inne = collections.defaultdict(set)
    for m in INDEKS_26.finditer(tekst):
        inne[kod(m.group(3), m.group(4))].add(m.group(1).strip())
    proby = collections.defaultdict(list)
    for p in points:
        e = p["english"]
        proby["zamiana nazwy"].append((p, "english", los.choice([x for x in points if x["english"] != e])["english"]))
        proby["zamiana pinyinu"].append((p, "pinyin", los.choice([x for x in points if x["pinyin"] != p["pinyin"]])["pinyin"]))
        for v in (e[:-1], e + " M", e[:2] + "q" + e[3:], e + "s"):
            proby["literowka"].append((p, "english", v))
        proby["nazwa alternatywna"] += [(p, "english", x) for x in inne[p["code"]] if klucz(x) != klucz(e)]
    przeszlo = 0
    for nazwa, lista in proby.items():
        zle = [(p["code"], v) for p, pole, v in lista if not wykryte(p, pole, v)]
        przeszlo += len(zle)
        print(f"mutacje - {nazwa}: {len(lista)}, przeszlo {len(zle)} {zle[:5] if zle else ''}")
    return przeszlo


def main() -> int:
    if not (KSIAZKA / "26-index-english-point-names.md").exists():
        print(f"brak wzoru w {KSIAZKA} - nie ma z czym porownywac", file=sys.stderr)
        return 2
    points = json.loads((ROOT / "data" / "points.json").read_text(encoding="utf-8"))["points"]
    wzor = wczytaj_ksiazke()
    if "--mutacje" in sys.argv:
        return 1 if mutacje(points, wzor) else 0
    bledy, poziomy = sprawdz(points, wzor)
    bledy += sprawdz_paczki(points)
    for line in bledy:
        print(f"NAZWA: {line}")
    print(f"potwierdzenia nazw: {dict(poziomy)}; bledow: {len(bledy)}")
    return 1 if bledy else 0


if __name__ == "__main__":
    raise SystemExit(main())
