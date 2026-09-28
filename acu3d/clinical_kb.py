"""Baza wiedzy klinicznej: 8 Naczyn Nadzwyczajnych, pary Yuan-Luo, lancuchy napiec i synergie."""
from __future__ import annotations

# 8 Naczyn Nadzwyczajnych (Ba Mai Jiao Hui Xue - 八脈交會穴)
EXTRAORDINARY_PAIRS: dict[str, dict[str, str]] = {
    "LU-7": {
        "vessel": "Ren Mai (Naczynie Poczęcia / Ren Channel)",
        "pair": "KI-6 (Zhaohai)",
        "pair_vessel": "Yin Qiao Mai (Naczynie Piętowe Yin)",
        "area": "Brzuszna linia pośrodkowa, klatka piersiowa, gardło, płuca, układ rozrodczy",
        "dynamics": "Otwiera obieg przedni ciała, odżywia Yin, chłodzi gardło, reguluje opadanie Qi Płuc i mocz.",
    },
    "KI-6": {
        "vessel": "Yin Qiao Mai (Naczynie Piętowe Yin)",
        "pair": "LU-7 (Lieque)",
        "pair_vessel": "Ren Mai (Naczynie Poczęcia)",
        "area": "Przyśrodkowa strona kończyn dolnych, gardło, oczy, pęcherz, narządy miednicy mniejszej",
        "dynamics": "Prowadzi Yin w górę do oczu i krtani, znosi bezsenność nocną, suchość i kurcze przywodzicieli.",
    },
    "SI-3": {
        "vessel": "Du Mai (Naczynie Zarządzające / Du Channel)",
        "pair": "BL-62 (Shenmai)",
        "pair_vessel": "Yang Qiao Mai (Naczynie Piętowe Yang)",
        "area": "Kręgosłup, kark, potylica, mózg, oś grzbietowa, oczy",
        "dynamics": "Wznosi Yang, prostuje kręgosłup, gasi Wiatr wewnętrzny/zewnętrzny, rozluźnia taśmę tylną.",
    },
    "BL-62": {
        "vessel": "Yang Qiao Mai (Naczynie Piętowe Yang)",
        "pair": "SI-3 (Houxi)",
        "pair_vessel": "Du Mai (Naczynie Zarządzające)",
        "area": "Boczna strona kończyn dolnych, plecy, kark, głowa, oczy",
        "dynamics": "Reguluje tonus taśmy boczno-tylnej, znosi nadmierną senność dzienną i asymetrię postawy.",
    },
    "SP-4": {
        "vessel": "Chong Mai (Naczynie Przenikające / Morze Krwi)",
        "pair": "PC-6 (Neiguan)",
        "pair_vessel": "Yin Wei Mai (Naczynie Łączące Yin)",
        "area": "Środkowy i dolny ogrzewacz, żołądek, jelita, serce, macica, krew",
        "dynamics": "Harmonizuje trawienie, porusza zastój Krwi w miednicy, łagodzi ból brzucha i nudności.",
    },
    "PC-6": {
        "vessel": "Yin Wei Mai (Naczynie Łączące Yin)",
        "pair": "SP-4 (Gongsun)",
        "pair_vessel": "Chong Mai (Naczynie Przenikające)",
        "area": "Klatka piersiowa, serce, nadbrzusze, gardło, Shen (Duch)",
        "dynamics": "Udrażnia klatkę piersiową, uspokaja Shen, znosi mdłości, kołatanie serca i ucisk podmostkowy.",
    },
    "GB-41": {
        "vessel": "Dai Mai (Naczynie Opasujące / Pas)",
        "pair": "TE-5 (Waiguan) / SJ-5",
        "pair_vessel": "Yang Wei Mai (Naczynie Łączące Yang)",
        "area": "Obręcz miednicy, pas, biodra, boki tułowia, skronie, oczy",
        "dynamics": "Wiąże wszystkie meridiany pionowe, stabilizuje miednicę, łagodzi migreny skroniowe i ból bioder.",
    },
    "TE-5": {
        "vessel": "Yang Wei Mai (Naczynie Łączące Yang)",
        "pair": "GB-41 (Zulinqi)",
        "pair_vessel": "Dai Mai (Naczynie Opasujące)",
        "area": "Boki ciała, barki, szyja, uszy, skronie, strefa Shaoyang",
        "dynamics": "Łączy wszystkie meridiany Yang, rozprasza Wiatr-Gorąco, wzmacnia Wei Qi (odporność zewnętrzną).",
    },
    "SJ-5": {
        "vessel": "Yang Wei Mai (Naczynie Łączące Yang)",
        "pair": "GB-41 (Zulinqi)",
        "pair_vessel": "Dai Mai (Naczynie Opasujące)",
        "area": "Boki ciała, barki, szyja, uszy, skronie, strefa Shaoyang",
        "dynamics": "Łączy wszystkie meridiany Yang, rozprasza Wiatr-Gorąco, wzmacnia Wei Qi (odporność zewnętrzną).",
    },
}

# Pary Gospodarz - Gość (Yuan-Luo / Host-Guest combinations)
YUAN_LUO_PAIRS: dict[str, tuple[str, str]] = {
    "LU-9": ("Punkt Yuan-Źródłowy Kanału Płuc", "LI-6 (Pianli) — Punkt Luo Jelita Grubego"),
    "LI-4": ("Punkt Yuan-Źródłowy Kanału Jelita Grubego", "LU-7 (Lieque) — Punkt Luo Płuc"),
    "ST-42": ("Punkt Yuan-Źródłowy Kanału Żołądka", "SP-4 (Gongsun) — Punkt Luo Śledziony"),
    "SP-3": ("Punkt Yuan-Źródłowy Kanału Śledziony", "ST-40 (Fenglong) — Punkt Luo Żołądka"),
    "HT-7": ("Punkt Yuan-Źródłowy Kanału Serca", "SI-7 (Zhizheng) — Punkt Luo Jelita Cienkiego"),
    "SI-4": ("Punkt Yuan-Źródłowy Kanału Jelita Cienkiego", "HT-5 (Tongli) — Punkt Luo Serca"),
    "BL-64": ("Punkt Yuan-Źródłowy Kanału Pęcherza Moczowego", "KI-4 (Dazhong) — Punkt Luo Nerek"),
    "KI-3": ("Punkt Yuan-Źródłowy Kanału Nerek", "BL-58 (Feiyang) — Punkt Luo Pęcherza Moczowego"),
    "PC-7": ("Punkt Yuan-Źródłowy Kanału Osierdzia", "TE-5 / SJ-5 (Waiguan) — Punkt Luo Potrójnego Ogrzewacza"),
    "TE-4": ("Punkt Yuan-Źródłowy Kanału Potrójnego Ogrzewacza", "PC-6 (Neiguan) — Punkt Luo Osierdzia"),
    "SJ-4": ("Punkt Yuan-Źródłowy Kanału Potrójnego Ogrzewacza", "PC-6 (Neiguan) — Punkt Luo Osierdzia"),
    "GB-40": ("Punkt Yuan-Źródłowy Kanału Pęcherzyka Żółciowego", "LR-5 (Ligou) — Punkt Luo Wątroby"),
    "LR-3": ("Punkt Yuan-Źródłowy Kanału Wątroby", "GB-37 (Guangming) — Punkt Luo Pęcherzyka Żółciowego"),
}

# Punkty Xi-Szczeliny (ostre zastoje, nagły ból)
XI_CLEFT_POINTS: dict[str, str] = {
    "LU-6": "Płuca (ostry kaszel, krwioplucie, ból klatki piersiowej)",
    "LI-7": "Jelito Grube (ostry ból gardła, obrzęk twarzy, ból zębów)",
    "ST-34": "Żołądek (ostry ból żołądka, zapalenie żołądka, ostry ból kolana)",
    "SP-8": "Śledziona (ostre bolesne miesiączkowanie, ból brzucha, zastoje krwi)",
    "HT-6": "Serce (ostry ból serca, dławica piersiowa, krwawienie z nosa)",
    "SI-6": "Jelito Cienkie (ostry ból barku, karku i łopatki)",
    "BL-63": "Pęcherz Moczowy (ostre zapalenie pęcherza, nagły lumbago)",
    "KI-5": "Nerki (ostre bolesne miesiączkowanie, trudności w oddawaniu moczu)",
    "PC-4": "Osierdzie (ostry ból klatki piersiowej, kołatanie serca, wymioty krwią)",
    "TE-7": "Potrójny Ogrzewacz (ostry ból ucha, ból głowy)",
    "SJ-7": "Potrójny Ogrzewacz (ostry ból ucha, ból głowy)",
    "GB-36": "Pęcherzyk Żółciowy (ostry ból w klatce i podżebrzu, wściekłość)",
    "LR-6": "Wątroba (ostry ból podbrzusza, krwotoki maciczne, zastoje)",
}


def get_extraordinary_vessel_info(code: str) -> list[str] | None:
    """Zwraca relacje z 8 Naczyniami Nadzwyczajnymi jesli punkt jest punktem zbiegu."""
    clean = code.split("·")[0].split(":")[0].strip().upper()
    if clean in EXTRAORDINARY_PAIRS:
        data = EXTRAORDINARY_PAIRS[clean]
        return [
            f"Punkt Nadrzędny (Master Point): {data['vessel']}",
            f"Punkt Sprzężony (Couple Point): {data['pair']} -> {data['pair_vessel']}",
            f"Strefa anatomiczna: {data['area']}",
            f"Dynamika terapeutyczna: {data['dynamics']}",
        ]
    return None


def get_yuan_luo_info(code: str) -> list[str] | None:
    """Zwraca relacje pary Gospodarz-Gosc (Yuan-Luo) dla danego punktu."""
    clean = code.split("·")[0].split(":")[0].strip().upper()
    if clean in YUAN_LUO_PAIRS:
        role, partner = YUAN_LUO_PAIRS[clean]
        return [
            f"Rola w parze: {role}",
            f"Punkt Partnerski (Gość): {partner}",
            "Zastosowanie: Jednoczesne wzmocnienie narządu macierzystego i udrożnienie meridianu sprzężonego.",
        ]
    # Sprawdz czy to punkt Luo w parze
    for y_code, (role, partner) in YUAN_LUO_PAIRS.items():
        if partner.startswith(clean):
            return [
                f"Rola w parze: Punkt Luo-Łączący sprzężony z punktem Yuan {y_code}",
                "Zastosowanie: Otwiera przepływ do kanału Yin/Yang o przeciwnej biegunowości.",
            ]
    return None


def get_xi_cleft_info(code: str) -> str | None:
    """Zwraca informacje o dzialaniu punktu Xi-Szczeliny."""
    clean = code.split("·")[0].split(":")[0].strip().upper()
    return XI_CLEFT_POINTS.get(clean)


# Pary Przednie Mu (Alarmowe) - Tylne Shu (Wzmacniajace narzadowo)
SHU_MU_PAIRS: dict[str, dict[str, str]] = {
    "LU": {"organ": "Płuca", "mu": "LU-1", "shu": "BL-13", "desc": "Przedni Mu (Zhongfu) diagnozuje ostry stan Płuc; Tylny Shu (Feishu) odżywia Qi Płuc i leczy przewlekły kaszel."},
    "LI": {"organ": "Jelito Grube", "mu": "ST-25", "shu": "BL-25", "desc": "Przedni Mu (Tianshu) reguluje perystaltykę i biegunkę; Tylny Shu (Dachangshu) tonizuje jelita i leczy zaparcia."},
    "ST": {"organ": "Żołądek", "mu": "REN-12", "shu": "BL-21", "desc": "Przedni Mu (Zhongwan) harmonizuje trawienie i żołądek; Tylny Shu (Weishu) wzmacnia Yang żołądka."},
    "SP": {"organ": "Śledziona", "mu": "LR-13", "shu": "BL-20", "desc": "Przedni Mu (Zhangmen) usuwa zastoje w brzuchu; Tylny Shu (Pishu) wytwarza Qi i Krew ze Śledziony."},
    "HT": {"organ": "Serce", "mu": "REN-14", "shu": "BL-15", "desc": "Przedni Mu (Juque) uspokaja Ducha Shen; Tylny Shu (Xinshu) odżywia Krew i Yang Serca."},
    "SI": {"organ": "Jelito Cienkie", "mu": "REN-4", "shu": "BL-27", "desc": "Przedni Mu (Guanyuan) wzmacnia Yuan Qi; Tylny Shu (Xiaochangshu) oddziela czyste od mętnego."},
    "BL": {"organ": "Pęcherz Moczowy", "mu": "REN-3", "shu": "BL-28", "desc": "Przedni Mu (Zhongji) usuwa Wilgoć-Gorąco z pęcherza; Tylny Shu (Pangguangshu) otwiera drogi moczowe."},
    "KI": {"organ": "Nerki", "mu": "GB-25", "shu": "BL-23", "desc": "Przedni Mu (Jingmen) wzmacnia esencję Nerek; Tylny Shu (Shenshu) tonizuje Yin i Yang Nerek."},
    "PC": {"organ": "Osierdzie", "mu": "REN-17", "shu": "BL-14", "desc": "Przedni Mu (Danzhong) otwiera klatkę piersiową i porusza Qi; Tylny Shu (Jueyinshu) reguluje rytm serca."},
    "TE": {"organ": "Potrójny Ogrzewacz", "mu": "REN-5", "shu": "BL-22", "desc": "Przedni Mu (Shimen) reguluje dolne drogi wodne; Tylny Shu (Sanjiaoshu) udrażnia transformację płynów."},
    "SJ": {"organ": "Potrójny Ogrzewacz", "mu": "REN-5", "shu": "BL-22", "desc": "Przedni Mu (Shimen) reguluje dolne drogi wodne; Tylny Shu (Sanjiaoshu) udrażnia transformację płynów."},
    "GB": {"organ": "Pęcherzyk Żółciowy", "mu": "GB-24", "shu": "BL-19", "desc": "Przedni Mu (Riyue) oczyszcza Gorąco z pęcherzyka; Tylny Shu (Danshu) wzmacnia odwagę i łagodzi żółtaczkę."},
    "LR": {"organ": "Wątroba", "mu": "LR-14", "shu": "BL-18", "desc": "Przedni Mu (Qimen) rozprasza zastoje Qi Wątroby; Tylny Shu (Ganshu) odżywia Krew Wątroby i wzmacnia oczy."},
}


def get_shu_mu_info(channel_or_point: str) -> dict[str, Any] | None:
    """Zwraca pare Shu-Mu dla danego kanalu lub sprawdza czy punkt nalezy do pary."""
    key = channel_or_point.split("·")[0].split(":")[0].strip().upper()
    if key in SHU_MU_PAIRS:
        return SHU_MU_PAIRS[key]
    for ch, d in SHU_MU_PAIRS.items():
        if d["mu"] == key:
            return {"channel": ch, "role": "Punkt Przedni-Mu (Alarmowy)", **d}
        if d["shu"] == key:
            return {"channel": ch, "role": "Punkt Tylny-Shu (Narządowy)", **d}
    return None

