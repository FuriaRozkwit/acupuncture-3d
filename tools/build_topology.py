"""Struktura kanalow: role punktow i miejsca, w ktorych kanaly sie spotykaja.

Po co
-----
Zeby polaczenie dwoch punktow bieglo tak, jak biegnie w ciele, a nie na skrot,
trzeba wiedziec nie tylko GDZIE punkty leza, ale JAK kanaly sie ze soba lacza.
Ten skrypt wyciaga te wiedze z opisow kategorii w `data/points.json` i zapisuje
ja jako czysta strukture: kody, role, listy kanalow. Zadnej prozy.

To ma dwa skutki. Model zyskuje topologie kanalow - wie, ze LU-7 jest punktem
luo Pluc, a GB-20 miejscem spotkania kilku kanalow. A dane traca smieci OCR,
ktore siedza w surowych opisach ("Lung Is", "Kz", 'Eb.'), bo z kazdego zdania
zostaje tylko wydobyty fakt.

    python tools/build_topology.py
"""
from __future__ import annotations

import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"

# nazwa kanalu w opisie -> kod. Kolejnosc ma znaczenie: "Large Intestine" musi
# byc sprawdzone przed "Intestine", a "Small Intestine" przed "Large".
CHANNELS = (
    ("large intestine", "LI"), ("small intestine", "SI"),
    ("gall bladder", "GB"), ("urinary bladder", "BL"), ("bladder", "BL"),
    ("pericardium", "PC"), ("sanjiao", "SJ"), ("triple burner", "SJ"),
    ("stomach", "ST"), ("spleen", "SP"), ("heart", "HT"), ("kidney", "KI"),
    ("liver", "LR"), ("lung", "LU"),
    ("conception vessel", "REN"), ("directing vessel", "REN"),
    ("governing vessel", "DU"),
)
# naczynia nadzwyczajne - lacza sie z kanalami glownymi przez punkty otwierajace
EXTRA = (
    ("penetrating vessel", "CHONG"), ("girdle vessel", "DAI"), ("girdling vessel", "DAI"),
    ("yin linking vessel", "YINWEI"), ("yang linking vessel", "YANGWEI"),
    ("yin motility vessel", "YINQIAO"), ("yang motility vessel", "YANGQIAO"),
    ("yin heel vessel", "YINQIAO"), ("yang heel vessel", "YANGQIAO"),
)
# role, ktore chcemy wydobyc: wzorzec w opisie -> nazwa roli
ROLES = (
    ("luo-connecting", "luo"), ("yuan-source", "yuan"), ("xi-cleft", "xi"),
    ("front-mu", "front_mu"), ("back-shu", "back_shu"),
    ("hui-meeting", "hui"), ("confluent", "confluent"),
    ("jing-well", "jing_well"), ("ying-spring", "ying_spring"),
    ("shu-stream", "shu_stream"), ("jing-river", "jing_river"), ("he-sea", "he_sea"),
    ("window of heaven", "window"), ("command point", "command"),
    ("sea of", "sea"), ("lower he-sea", "lower_he_sea"),
)

# Pary kanalow. To sa fakty tradycji, nie cytat: kazdy podrecznik podaje je
# identycznie, bo wynikaja z ukladu piecu przemian i szesciu podzialow.
INTERIOR_EXTERIOR = (("LU", "LI"), ("ST", "SP"), ("HT", "SI"),
                     ("BL", "KI"), ("PC", "SJ"), ("GB", "LR"))
SAME_NAME = {"taiyin": ("LU", "SP"), "yangming": ("LI", "ST"),
             "shaoyin": ("HT", "KI"), "taiyang": ("SI", "BL"),
             "jueyin": ("PC", "LR"), "shaoyang": ("SJ", "GB")}


def channels_in(text: str) -> list:
    """Kody kanalow wymienione w opisie, w kolejnosci wystapienia."""
    low = text.lower()
    found = []
    for name, code in CHANNELS:
        i = low.find(name)
        if i >= 0 and code not in [c for _, c in found]:
            found.append((i, code))
    return [c for _, c in sorted(found)]


def extras_in(text: str) -> list:
    low = text.lower()
    return sorted({code for name, code in EXTRA if name in low})


def main() -> int:
    doc = json.loads((DATA / "points.json").read_text(encoding="utf-8"))
    roles: dict = {}
    crossings: dict = {}
    opens: dict = {}
    for p in doc["points"]:
        code, own = p["code"], p["channel"]
        for cat in (p.get("categories") or []):
            low = cat.lower()
            ch = channels_in(cat)
            ex = extras_in(cat)
            hit = False
            if "great luo" in low:          # wielki luo sledziony - osobna rola
                roles.setdefault(code, []).append(
                    {"role": "great_luo", "of": (ch[0] if ch else own)})
                continue
            # Jeden opis potrafi nadac kilka rol naraz: LU-9 to "Shu-Stream,
            # Yuan-Source and Earth point". Zatrzymanie sie na pierwszym trafieniu
            # gubilo punkty shu kanalow yin i psulo caly obraz piatki transportu.
            # Dolny he-sea Jelita Grubego to ST-37, a NIE punkt he-sea tego kanalu
            # (tym jest LI-11). Wzorzec "he-sea" zawiera sie w "lower he-sea",
            # wiec bez tego wyjatku ST-37 wchodzil do piatki shu
            # transportu obcego kanalu i psul caly pomiar samopodobienstwa.
            lower = "lower he-sea" in low
            for pat, role in ROLES:
                if pat not in low:
                    continue
                if lower and role == "he_sea":
                    continue
                target = (ex[0] if ex else (ch[0] if ch else own))
                if not any(r["role"] == role and r["of"] == target
                           for r in roles.get(code, [])):
                    roles.setdefault(code, []).append({"role": role, "of": target})
                if role == "confluent":
                    opens[code] = sorted(set(opens.get(code, []) + [target]))
                hit = True
            # "Meeting point of the Bladder, Gall Bladder and Sanjiao channels"
            if "meeting point of" in low and "hui-meeting" not in low:
                meet = sorted(set(ch + ex) - {own})
                if meet:
                    crossings[code] = sorted(set(crossings.get(code, []) + meet))
                hit = True
            if ex and ("confluent" in low or "opening" in low):
                opens[code] = sorted(set(opens.get(code, []) + ex))
                hit = True
            if not hit and ch and ch[0] != own:
                crossings[code] = sorted(set(crossings.get(code, []) + [ch[0]]))

    # skrot: dla kazdego kanalu jego punkt luo i yuan
    luo = {}
    for c, rs in sorted(roles.items()):
        for r in rs:
            if r["role"] == "luo" and r["of"] not in luo:
                luo[r["of"]] = c
    yuan = {r["of"]: c for c, rs in roles.items() for r in rs if r["role"] == "yuan"}

    out = {
        "source": "wydobyte z pol categories w data/points.json (tools/build_topology.py); "
                  "same fakty strukturalne - kody, role, listy kanalow",
        "interior_exterior": [list(x) for x in INTERIOR_EXTERIOR],
        "same_name": {k: list(v) for k, v in SAME_NAME.items()},
        "roles": {c: rs for c, rs in sorted(roles.items())},
        "crossings": dict(sorted(crossings.items())),
        "opening_points": dict(sorted(opens.items())),
        "luo_of": dict(sorted(luo.items())),
        "yuan_of": dict(sorted(yuan.items())),
    }
    (DATA / "channel_topology.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")

    n_roles = sum(len(v) for v in roles.values())
    print(f"punktow z rola: {len(roles)}  (rol lacznie {n_roles})")
    print(f"punktow skrzyzowan: {len(crossings)}")
    print(f"punktow otwierajacych naczynia nadzwyczajne: {len(opens)}")
    print(f"kanalow z punktem luo: {len(luo)} / 14,  z yuan: {len(yuan)} / 14")
    miss = [c for c in doc["channels"] if c not in luo and c not in ("REN", "DU")]
    if miss:
        print(f"bez punktu luo w danych: {', '.join(miss)}")
    print("zapisano data/channel_topology.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
