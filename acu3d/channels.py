"""Meridianologia: jak dwa punkty lacza sie ze soba przez kanaly.

Prosta miedzy punktami nie jest polaczeniem - to tylko cieciwa. W tradycji
punkty lacza sie WZDLUZ kanalu, a kanaly ze soba w okreslonych miejscach:

  * punkt luo laczy kanal z jego para wewnetrzno-zewnetrzna (LU-7 z Jelitem
    Grubym, ST-40 ze Sledziona);
  * punkt skrzyzowania (jiaohui) to miejsce, gdzie kilka kanalow faktycznie
    przechodzi przez to samo miejsce w ciele - GB-20, DU-14 i tak dalej;
  * osiem punktow otwierajacych laczy kanal glowny z naczyniem nadzwyczajnym.

Ten modul szuka takiego polaczenia i mowi, KTOREDY ono idzie. Struktura kanalow
pochodzi z `data/channel_topology.json` (`tools/build_topology.py`).

Czego tu nie ma: przebiegow wewnetrznych, kanalow rozchodzacych sie (jing bie),
sciegnisto-miesniowych (jing jin) i podluznych luo. Model zna tylko przebieg
zewnetrzny - ten, ktory da sie pokazac na skorze.
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"

PL_ROLE = {"luo": "punkt luo", "yuan": "punkt yuan", "xi": "punkt xi",
           "front_mu": "punkt mu", "back_shu": "punkt shu grzbietu",
           "confluent": "punkt otwierajacy", "hui": "punkt hui",
           "great_luo": "wielki luo"}
PL_VESSEL = {"CHONG": "Naczynie Przenikajace", "DAI": "Naczynie Opasujace",
             "YINWEI": "Naczynie Yin Wei", "YANGWEI": "Naczynie Yang Wei",
             "YINQIAO": "Naczynie Yin Qiao", "YANGQIAO": "Naczynie Yang Qiao",
             "REN": "Naczynie Poczecia", "DU": "Naczynie Gubernatora"}


class ChannelMap:
    def __init__(self, atlas):
        self.atlas = atlas
        t = json.loads((DATA / "channel_topology.json").read_text(encoding="utf-8"))
        self.luo_of = t["luo_of"]
        self.yuan_of = t["yuan_of"]
        self.crossings = t["crossings"]
        self.opening = t["opening_points"]
        self.roles = t["roles"]
        self.pairs = {}
        for a, b in t["interior_exterior"]:
            self.pairs[a] = b
            self.pairs[b] = a
        self.same_name = {}
        for div, (a, b) in t["same_name"].items():
            self.same_name[a] = (b, div)
            self.same_name[b] = (a, div)
        paths = json.loads((DATA / "meridians.json").read_text(encoding="utf-8"))["paths"]
        self.order = {ch: [c for branch in br for c in branch] for ch, br in paths.items()}
        # kanaly przechodzace przez dany punkt skrzyzowania
        self.meets: dict = {}
        for code, chans in self.crossings.items():
            own = code.split("-")[0]
            self.meets[code] = sorted(set(chans) | {own})

    # --- odcinek kanalu ---------------------------------------------------
    def segment(self, channel: str, a: str, b: str) -> list:
        """Kody punktow od `a` do `b` wzdluz kanalu (pusta lista, gdy sie nie da)."""
        seq = self.order.get(channel, [])
        if a not in seq or b not in seq:
            return []
        i, j = seq.index(a), seq.index(b)
        return seq[i:j + 1] if i <= j else seq[j:i + 1][::-1]

    def role_of(self, code: str) -> str:
        rs = self.roles.get(code, [])
        names = [PL_ROLE.get(r["role"], r["role"]) for r in rs if r["role"] in PL_ROLE]
        return ", ".join(dict.fromkeys(names))

    # --- polaczenie dwoch punktow ----------------------------------------
    def connect(self, code_a: str, code_b: str) -> dict:
        """Jak dojsc z jednego punktu do drugiego trzymajac sie kanalow.

        Kolejnosc szukania nie jest dowolna - idzie od polaczenia najsilniejszego
        do najslabszego:

        1. ten sam kanal - nie ma czego szukac, idziemy po nim;
        2. oba punkty otwieraja naczynia nadzwyczajne - to sa klasyczne pary
           (SI-3 z BL-62, LU-7 z KI-6), a nie przypadkowa zbieznosc;
        3. para wewnetrzno-zewnetrzna - przez punkt luo, bo po to on jest;
        4. punkt skrzyzowania - kanaly faktycznie przechodza przez to samo miejsce;
        5. ten sam podzial szesciu kanalow - pokrewienstwo bez styku.

        Bez tej kolejnosci skrzyzowanie wygrywalo z luo i GB-34 laczylo sie
        z LR-3 objazdem przez podzebrze zamiast przez GB-37.
        """
        ca, cb = code_a.split("-")[0], code_b.split("-")[0]

        if ca == cb:
            seq = self.segment(ca, code_a, code_b)
            if seq:
                return {"codes": seq, "via": "channel",
                        "why": f"ten sam kanal {ca} - odcinek {len(seq)} punktow"}

        va = set(self.opening.get(code_a, []))
        vb = set(self.opening.get(code_b, []))
        if va and vb:
            both = va & vb
            if both:
                v = sorted(both)[0]
                return {"codes": [code_a, code_b], "via": "extraordinary",
                        "vessel": v,
                        "why": f"oba punkty otwieraja {PL_VESSEL.get(v, v)}"}
            return {"codes": [code_a, code_b], "via": "extraordinary_pair",
                    "why": f"{code_a} otwiera {PL_VESSEL.get(sorted(va)[0], '?')}, "
                           f"{code_b} otwiera {PL_VESSEL.get(sorted(vb)[0], '?')} - "
                           f"klasyczna para punktow otwierajacych"}

        if self.pairs.get(ca) == cb:
            lo = self.luo_of.get(ca)
            lo_b = self.luo_of.get(cb)
            first = self.segment(ca, code_a, lo) if lo else []
            second = self.segment(cb, lo_b, code_b) if lo_b else []
            if first and second:
                # z punktu luo kanalu A wchodzimy w kanal B przy jego wlasnym luo:
                # to jest ta para, ktora oba punkty luo lacza
                return {"codes": first + second, "via": "luo", "meeting": lo,
                        "why": f"{ca} i {cb} to para wewnetrzno-zewnetrzna; "
                               f"laczy je punkt luo {lo} (po stronie {cb}: {lo_b})"}
            if first:
                return {"codes": first, "via": "luo", "meeting": lo,
                        "why": f"{ca} i {cb} to para wewnetrzno-zewnetrzna; "
                               f"laczy je punkt luo {lo}"}

        best = None
        for code, chans in self.meets.items():
            if ca in chans and cb in chans:
                s1 = self._to_crossing(ca, code_a, code)
                s2 = self._to_crossing(cb, code_b, code)
                if s1 and s2:
                    n = len(s1) + len(s2)
                    if best is None or n < best[0]:
                        best = (n, code, s1, s2)
        if best is not None:
            _, code, s1, s2 = best
            return {"codes": s1 + s2[::-1], "via": "crossing", "meeting": code,
                    "why": f"kanaly {ca} i {cb} spotykaja sie w {code}"}

        nb = self.same_name.get(ca)
        if nb and nb[0] == cb:
            return {"codes": [], "via": "same_name", "division": nb[1],
                    "why": f"{ca} i {cb} to ten sam podzial ({nb[1]}), "
                           f"ale zadnym punktem sie w tym modelu nie stykaja"}

        return {"codes": [], "via": None,
                "why": f"{ca} i {cb} nie lacza sie w strukturze, ktora model zna - "
                       f"trasa pojdzie po powierzchni ciala"}

    def _to_crossing(self, channel: str, code: str, crossing: str) -> list:
        """Odcinek kanalu od punktu do miejsca skrzyzowania.

        Skrzyzowanie ma wlasny kod (np. GB-20), ale lezy takze na innych kanalach.
        Dla kanalu, ktory tylko przez nie przechodzi, bierzemy najblizszy jego
        wlasny punkt - dalej trasa idzie juz po powierzchni ciala.
        """
        seq = self.order.get(channel, [])
        if crossing in seq:
            return self.segment(channel, code, crossing)
        if not seq or code not in seq:
            return []
        by_key = self.atlas.by_key
        tgt = by_key.get(crossing) or by_key.get(crossing + "·R")
        if tgt is None:
            return []
        import numpy as np
        best, bd = None, 1e9
        for c in seq:
            m = by_key.get(c) or by_key.get(c + "·R")
            if m is None:
                continue
            d = float(np.linalg.norm(m["pos"] - tgt["pos"]))
            if d < bd:
                best, bd = c, d
        return self.segment(channel, code, best) if best else []
