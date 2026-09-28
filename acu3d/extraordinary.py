"""Topologiczne mapowanie 8 Naczyn Nadzwyczajnych (Qi Jing Ba Mai - 奇經八脈).

Naczynia nadzwyczajne (z wyjatkiem Ren Mai i Du Mai) nie posiadaja wlasnych punktow,
lecz spajaja punkty 12 meridianow glownych w nadrzedne magistrale topologiczne
(Punkty Skrzyzowan / Jiao Hui Xue - 交會穴).

Modul definiuje sekwencje wezlow, punkty otwierajace, punkty Xi, dynamike
kliniczna oraz kolorystyke rysunkowa dla kazdego z 8 naczyn.
"""
from __future__ import annotations

from typing import Any

# Kolorystyka rysunkowa 8 Naczyn Nadzwyczajnych (odcinajaca sie od Wu Xing)
EXTRAORDINARY_COLORS: dict[str, str] = {
    "DAI_MAI": "#eab308",      # Złoty bursztyn (Pas / Jedyny meridian poziomy)
    "CHONG_MAI": "#dc2626",    # Głęboki cynober / szkarłat (Morze Krwi)
    "YIN_QIAO": "#06b6d4",     # Jasny turkus (Wznoszenie Yin od pięty do oka)
    "YANG_QIAO": "#f97316",    # Żywy oranż (Ruch Yang od stopy przez bark do potylicy)
    "YIN_WEI": "#8b5cf6",      # Fiolet orchidei (Wiązanie kanałów Yin tułowia)
    "YANG_WEI": "#ec4899",     # Magenta / rubin (Wiązanie kanałów Yang głowy i boków)
    "REN_MAI": "#a855f7",      # Królewski ametyst (Oś przednia Yin)
    "DU_MAI": "#6366f1",       # Indygo (Oś tylna Yang)
}

EXTRAORDINARY_VESSELS: dict[str, dict[str, Any]] = {
    "DAI_MAI": {
        "id": "DAI_MAI",
        "name_pl": "Dai Mai (Naczynie Opasujące / Pas)",
        "name_zh": "帶脈",
        "name_pinyin": "Dài Mài",
        "master_point": "GB-41",
        "coupled_point": "TE-5",
        "xi_point": None,
        "is_horizontal": True,
        "color": EXTRAORDINARY_COLORS["DAI_MAI"],
        # Sekwencja wezlow tworzaca pas wokol talii (od wolnego zebra do podbrzusza)
        "sequence": ["LR-13", "GB-26", "GB-27", "GB-28"],
        "area": "Talia, pas, podbrzusze, miednica, biodra, narządy rodne",
        "dynamics": "Jedyny meridian o przebiegu poprzecznym. Działa jak pas ściągający wszystkie pionowe meridiany, stabilizuje miednicę, kontroluje upławy i tonus powłok brzusznych.",
        "indications": "Bóle lędźwi jak przy siedzeniu w wodzie, osłabienie nóg, wypadanie narządów, zaburzenia ginekologiczne, przepukliny, blokady pasa biodrowego.",
    },
    "CHONG_MAI": {
        "id": "CHONG_MAI",
        "name_pl": "Chong Mai (Naczynie Przenikające / Morze Krwi)",
        "name_zh": "衝脈",
        "name_pinyin": "Chōng Mài",
        "master_point": "SP-4",
        "coupled_point": "PC-6",
        "xi_point": None,
        "is_horizontal": False,
        "color": EXTRAORDINARY_COLORS["CHONG_MAI"],
        # Sekwencja brzuszno-piersiowa wzdluz linii Nerek (0.5c od linii posrodkowej)
        "sequence": [
            "REN-1", "ST-30",
            "KI-11", "KI-12", "KI-13", "KI-14", "KI-15",
            "KI-16", "KI-17", "KI-18", "KI-19", "KI-20", "KI-21"
        ],
        "area": "Brzuch, klatka piersiowa, serce, żołądek, macica, krew, dół biodrowy",
        "dynamics": "Morze Krwi (Xue Hai) i 12 Meridianów. Reguluje cykle krwi, odżywia esencję Jing, harmonizuje Qi żołądka skierowane podprąd (nudności, ucisk w klatce).",
        "indications": "Nieregularne miesiączki, bolesne skurcze podbrzusza, dławienie w klatce piersiowej, kołatanie serca, refluks i wymioty, uczucie uderzenia gorąca w górę (Ni Qi).",
    },
    "YIN_QIAO": {
        "id": "YIN_QIAO",
        "name_pl": "Yin Qiao Mai (Naczynie Piętowe Yin)",
        "name_zh": "陰蹻脈",
        "name_pinyin": "Yīn Qiāo Mài",
        "master_point": "KI-6",
        "coupled_point": "LU-7",
        "xi_point": "KI-8",
        "is_horizontal": False,
        "color": EXTRAORDINARY_COLORS["YIN_QIAO"],
        # Od piety wewnetrznej wzdluz przysrodkowego uda, gardla do oka
        "sequence": ["KI-6", "KI-8", "ST-9", "BL-1"],
        "area": "Wewnętrzna strona nóg, krocze, podbrzusze, gardło, oczy (kącik wewnętrzny)",
        "dynamics": "Prowadzi energię Yin w górę od pięty do oczu. Równoważy napięcie zginaczy i przywodzicieli. Zapewnia zamykanie oczu i sen w nocy.",
        "indications": "Nadmierna senność w dzień lub bezsenność nocna, suchość gardła, kurcze przywodzicieli nóg (noga w rotacji wewnętrznej), bolesne miesiączkowanie, trudności w oddawaniu moczu.",
    },
    "YANG_QIAO": {
        "id": "YANG_QIAO",
        "name_pl": "Yang Qiao Mai (Naczynie Piętowe Yang)",
        "name_zh": "陽蹻脈",
        "name_pinyin": "Yáng Qiāo Mài",
        "master_point": "BL-62",
        "coupled_point": "SI-3",
        "xi_point": "BL-59",
        "is_horizontal": False,
        "color": EXTRAORDINARY_COLORS["YANG_QIAO"],
        # Od piety zewnetrznej, przez biodro, bark, twarz do oka i potylicy
        "sequence": [
            "BL-62", "BL-61", "BL-59", "GB-29", "SI-10",
            "LI-15", "LI-16", "ST-4", "ST-3", "ST-1", "BL-1", "GB-20"
        ],
        "area": "Boczna strona nóg, biodro, tylny bark, policzek, kącik oka, potylica",
        "dynamics": "Wznosi energię Yang od stóp ku głowie. Odpowiada za czuwanie, otwarcie oczu w dzień, napięcie mięśni prostowników i taśmy tylno-bocznej.",
        "indications": "Bezsenność, oczy otwierające się mimowolnie, sztywność i kurcze mięśni bocznej strony nogi (noga w rotacji zewnętrznej), bóle barku i karku, porażenie nerwu twarzowego.",
    },
    "YIN_WEI": {
        "id": "YIN_WEI",
        "name_pl": "Yin Wei Mai (Naczynie Łączące Yin)",
        "name_zh": "陰維脈",
        "name_pinyin": "Yīn Wéi Mài",
        "master_point": "PC-6",
        "coupled_point": "SP-4",
        "xi_point": "KI-9",
        "is_horizontal": False,
        "color": EXTRAORDINARY_COLORS["YIN_WEI"],
        # Od podudzia wewnetrznego przez brzuch i zebra do gardla
        "sequence": [
            "KI-9", "SP-12", "SP-13", "SP-15", "SP-16",
            "LR-14", "REN-22", "REN-23"
        ],
        "area": "Wewnętrzne podudzie, podbrzusze, boki brzucha, żebra, klatka piersiowa, gardło, serce",
        "dynamics": "Wiąże i harmonizuje wszystkie meridiany Yin wnętrza ciała. Odpowiada za dobrostan emocjonalny, łagodzi lęki, ucisk w klatce piersiowej i ból serca.",
        "indications": "Bóle serca, ucisk za mostkiem, duszność, melancholia, lęki i koszmary senne, ból w podżebrzu, trudności w połykaniu.",
    },
    "YANG_WEI": {
        "id": "YANG_WEI",
        "name_pl": "Yang Wei Mai (Naczynie Łączące Yang)",
        "name_zh": "陽維脈",
        "name_pinyin": "Yáng Wéi Mài",
        "master_point": "TE-5",
        "coupled_point": "GB-41",
        "xi_point": "GB-35",
        "is_horizontal": False,
        "color": EXTRAORDINARY_COLORS["YANG_WEI"],
        # Od stopy zewnetrznej, przez biodro, bark, skron i czolo do karku
        "sequence": [
            "BL-63", "GB-35", "GB-29", "SI-10", "SJ-15", "GB-21",
            "ST-8", "GB-13", "GB-14", "GB-15", "GB-20", "DU-16", "DU-15"
        ],
        "area": "Zewnętrzna stopa i łydka, biodro, bark, skroń, czoło, kark, strefa Shaoyang",
        "dynamics": "Wiąże i koordynuje wszystkie kanały Yang powierzchni ciała. Reguluje Wei Qi (odporność obronną), chroni przed czynnikami patogennymi zewnętrznymi (Wiatr-Zimno).",
        "indications": "Dreszcze i gorączka naprzemienna, bóle głowy skroniowe i potyliczne, sztywność karku, bóle w strefie Shaoyang, osłabienie odporności.",
    },
    "REN_MAI": {
        "id": "REN_MAI",
        "name_pl": "Ren Mai (Naczynie Poczęcia / Morze Yin)",
        "name_zh": "任脈",
        "name_pinyin": "Rèn Mài",
        "master_point": "LU-7",
        "coupled_point": "KI-6",
        "xi_point": None,
        "is_horizontal": False,
        "color": EXTRAORDINARY_COLORS["REN_MAI"],
        "sequence": [f"REN-{i}" for i in range(1, 25)],
        "area": "Brzuszna linia pośrodkowa od krocza do wargi dolnej",
        "dynamics": "Morze wszystkich kanałów Yin (Yin Mai Zhi Hai). Kieruje układem rozrodczym, płodnością, ciążą i odżywianiem narządów wewnętrznych.",
        "indications": "Zaburzenia płodności, upławy, bóle brzucha, duszność, obrzęki, zaburzenia oddawania moczu.",
    },
    "DU_MAI": {
        "id": "DU_MAI",
        "name_pl": "Du Mai (Naczynie Zarządzające / Morze Yang)",
        "name_zh": "督脈",
        "name_pinyin": "Dū Mài",
        "master_point": "SI-3",
        "coupled_point": "BL-62",
        "xi_point": None,
        "is_horizontal": False,
        "color": EXTRAORDINARY_COLORS["DU_MAI"],
        "sequence": [f"DU-{i}" for i in range(1, 29)],
        "area": "Grzbietowa linia pośrodkowa wzdłuż kręgosłupa, głowy do wargi górnej",
        "dynamics": "Morze wszystkich kanałów Yang (Yang Mai Zhi Hai). Wznosi Yang Qi, odżywia mózg i szpik, stabilizuje oś kręgosłupa i postawę pionową.",
        "indications": "Sztywność kręgosłupa, ból pleców, zawroty głowy, opadanie narządów, osłabienie Yang, udary, drgawki.",
    },
}


def get_vessel(key: str) -> dict[str, Any] | None:
    """Zwraca dane naczynia po kluczu lub nazwie (np. 'DAI_MAI', 'dai', 'chong', 'qiao')."""
    k = key.upper().replace(" ", "_").replace("-", "_")
    if k in EXTRAORDINARY_VESSELS:
        return EXTRAORDINARY_VESSELS[k]
    # Przeszukiwanie czesciowe
    for vk, v in EXTRAORDINARY_VESSELS.items():
        if k in vk or k in v["name_pinyin"].upper() or k in v["name_pl"].upper():
            return v
    return None


def get_point_extraordinary_memberships(point_code: str) -> list[dict[str, str]]:
    """Zwraca liste naczyn nadzwyczajnych, w ktorych dany punkt bierze udzial."""
    code = point_code.upper().strip()
    clean = code.replace("·L", "").replace("·R", "").replace(":L", "").replace(":R", "")
    memberships = []
    for vid, v in EXTRAORDINARY_VESSELS.items():
        seq = [p.upper() for p in v["sequence"]]
        role = []
        if clean in seq:
            idx = seq.index(clean) + 1
            role.append(f"węzeł #{idx}/{len(seq)} na trasie")
        if v.get("master_point") == clean:
            role.append("GŁÓWNY PUNKT OTWIERAJĄCY (Bā Mài)")
        if v.get("coupled_point") == clean:
            role.append("PUNKT SPRZĘŻONY (Para Bā Mài)")
        if v.get("xi_point") == clean:
            role.append("PUNKT SZCZELINY XI (Ostry ból)")
        if role:
            memberships.append({
                "vessel_id": vid,
                "vessel_name": v["name_pl"],
                "role": "; ".join(role),
                "color": v["color"],
            })
    return memberships
