"""Dzienniczek kombinacji punktow: 2-3 punkty wybrane razem, z tym co miedzy nimi.

Kombinacja to nie to samo co linia napiecia. Linia napiecia bierze sie z ciala -
stawiasz ja tam, gdzie palec cos zastal. Kombinacja bierze sie z wiedzy: wybierasz
punkty, o ktorych wiadomo, ze dzialaja razem, i patrzysz, co ich polaczenie
przechodzi w tym konkretnym ciele.

Zapis jest lokalny (`data/combos.json`, poza repozytorium). Powod jest ten sam,
dla ktorego poza repozytorium jest `point_meta.json`: opisy dzialania kombinacji
w podrecznikach sa materialem chronionym. Kody punktow i tradycyjne nazwy to
fakty i mozna je zapisywac swobodnie - przepisane akapity juz nie.
"""
from __future__ import annotations

import json

import numpy as np
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BOOK = ROOT / "data" / "combos.json"


class Combination:
    """Jedna kombinacja: punkty w kolejnosci wyboru plus wlasna notatka."""

    def __init__(self, keys: list, name: str = "", note: str = "", ts: str = None):
        self.keys = list(keys)                 # np. ["LI-4·R", "LR-3·R"]
        self.name = name or " + ".join(k.split("·")[0] for k in self.keys)
        self.note = note
        self.ts = ts or datetime.now(timezone.utc).isoformat(timespec="seconds")

    @property
    def codes(self) -> list:
        return [k.split("·")[0] for k in self.keys]

    def to_dict(self) -> dict:
        return {"keys": self.keys, "name": self.name, "note": self.note, "ts": self.ts}

    @classmethod
    def from_dict(cls, d: dict) -> "Combination":
        return cls(d["keys"], d.get("name", ""), d.get("note", ""), d.get("ts"))


class CombinationBook:
    """Dzienniczek na dysku. Startuje pusty - wypelnia go uzytkownik."""

    def __init__(self, path: Path | None = None):
        self.path = path or BOOK
        self.items: list = []
        self.load()

    def load(self) -> int:
        if not self.path.exists():
            return 0
        doc = json.loads(self.path.read_text(encoding="utf-8"))
        self.items = [Combination.from_dict(d) for d in doc.get("combinations", [])]
        return len(self.items)

    def save(self) -> Path:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        doc = {"note": "wlasny dzienniczek kombinacji; opisow z podrecznikow tu nie wpisywac",
               "saved": datetime.now(timezone.utc).isoformat(timespec="seconds"),
               "combinations": [c.to_dict() for c in self.items]}
        self.path.write_text(json.dumps(doc, ensure_ascii=False, indent=1),
                             encoding="utf-8")
        return self.path

    def add(self, keys: list, name: str = "", note: str = "") -> Combination:
        c = Combination(keys, name, note)
        self.items.append(c)
        return c

    def find(self, keys: list):
        want = sorted(keys)
        for c in self.items:
            if sorted(c.keys) == want:
                return c
        return None

    def __len__(self) -> int:
        return len(self.items)


def channel_route(tension, atlas, chan, key_a: str, key_b: str):
    """Trasa wzdluz kanalow miedzy dwoma punktami - albo None, gdy jej nie ma.

    Plan podaje `channels.ChannelMap`: kody, przez ktore trzeba przejsc. Miedzy
    kolejnymi kodami idziemy po powierzchni ciala, bo odcinki sa krotkie i wtedy
    trasa trzyma sie skory zamiast ciac przez wnetrze.
    """
    ca, cb = key_a.split("·")[0], key_b.split("·")[0]
    side = key_a.split("·")[1] if "·" in key_a else "R"
    plan = chan.connect(ca, cb)
    codes = plan.get("codes") or []
    if len(codes) < 2 or tension.wire is None:
        return None, plan
    poly = []
    prev = None
    for c in codes:
        m = (atlas.by_key.get(f"{c}·{side}") or atlas.by_key.get(c)
             or atlas.by_key.get(f"{c}·R"))
        if m is None:
            continue
        if prev is not None:
            seg, _ = tension.wire.route(prev, m["pos"])
            poly.extend(seg.tolist() if not poly else seg.tolist()[1:])
        prev = m["pos"]
    return (np.array(poly) if len(poly) >= 2 else None), plan


def describe(tension, markers_by_key: dict, keys: list, corridor: float = 0.03,
             atlas=None, chan=None) -> dict:
    """Co lezy miedzy wybranymi punktami - dla kazdej pary osobno.

    Przy dwoch punktach jest jedna para, przy trzech - trzy, wiec od razu widac,
    czy trojka tworzy lancuch (A-B, B-C krotkie, A-C dlugie) czy trojkat.
    """
    pts = [markers_by_key[k] for k in keys if k in markers_by_key]
    out = {"keys": keys, "pairs": [], "channels": {}, "segments": []}
    for m in pts:
        out["channels"][m["channel"]] = out["channels"].get(m["channel"], 0) + 1
    for i in range(len(pts)):
        for j in range(i + 1, len(pts)):
            a, b = pts[i], pts[j]
            poly, plan = (None, None)
            if chan is not None and atlas is not None:
                poly, plan = channel_route(tension, atlas, chan, keys[i], keys[j])
            r = tension.analyse_span(a["pos"], b["pos"], corridor, path=poly)
            r["from"], r["to"], r["plan"] = a, b, plan
            # punkty koncowe nie sa "po drodze" - zostawiamy to, co miedzy nimi
            ends = {a["code"], b["code"]}
            r["between"] = [q for q in r["points"] if q["code"] not in ends]
            out["pairs"].append(r)
    return out


def report(tension, markers_by_key: dict, keys: list, corridor: float = 0.03,
           atlas=None, chan=None) -> str:
    d = describe(tension, markers_by_key, keys, corridor, atlas, chan)
    names = " + ".join(k.split("·")[0] for k in keys)
    lines = [f"KOMBINACJA: {names}", "=" * 78]
    ch = ", ".join(f"{c} ({n})" for c, n in sorted(d["channels"].items(),
                                                   key=lambda kv: -kv[1]))
    lines.append(f"kanaly: {ch}")
    lines.append("")
    for r in d["pairs"]:
        a, b = r["from"], r["to"]
        lines.append(f"{a['code']} -> {b['code']}   po ciele {r['length_mm']:.0f} mm"
                     f"   (w linii prostej {r['chord_mm']:.0f} mm)")
        plan = r.get("plan")
        if plan:
            lines.append(f"   polaczenie: {plan['why']}")
        if r["between"]:
            lines.append("   po drodze: " + ", ".join(
                f"{q['code']} ({q['off_mm']:.0f} mm od osi)" for q in r["between"][:10]))
        else:
            lines.append("   po drodze: nic w korytarzu")
        for m in r["muscles"][:3]:
            side = {"L": "lewy", "R": "prawy", "C": ""}.get(m["side"], "")
            lines.append(f"   wzdluz: {m['name']}{(' (' + side + ')') if side else ''}"
                         f"  os pod {m['angle_deg']} st., najplycej {m['depth_mm']:.0f} mm")
        lines.append("")
    lines.append("Trasa biegnie po powierzchni ciala (graf siatki cunowej), a nie")
    lines.append("na skrot przez przestrzen. To jest geometria tego modelu,")
    lines.append("nie uzasadnienie kombinacji.")
    lines.append("Mowi, ktore punkty i miesnie lezy na odcinku miedzy wybranymi -")
    lines.append("czy to cokolwiek znaczy klinicznie, rozstrzyga wiedza, nie odleglosc.")
    return "\n".join(lines)
