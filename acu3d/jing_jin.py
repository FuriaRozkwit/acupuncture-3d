"""Mapowanie 12 Meridianow Sciegnisto-Miesniowych (Jing Jin - 经筋) oraz Tasm Anatomicznych Myersa (Anatomy Trains).

W medycynie chinskiej Jing Jin nie wnikaja do narzadow Zang-Fu, lecz tworza
szerokie pasma powiezio-miesniowe odpowiedzialne za motoryke, postawe pionowa
i przenoszenie naprezen (biotensegracja).

Koresponduja bezposrednio z tasmami powieziowymi Thomasa Myersa:
- Zu Taiyang (BL)  <-> Superficial Back Line (SBL - Tasma Powierzchowna Tylna)
- Zu Yangming (ST) <-> Superficial Front Line (SFL - Tasma Powierzchowna Przednia)
- Zu Shaoyang (GB) <-> Lateral Line (LL - Tasma Boczna)
- Zu Taiyin (SP) / Zu Shaoyin (KI) / Zu Jueyin (LR) <-> Deep Front Line (DFL - Tasma Gleboka Przednia)
- Kanaly reki (LU, LI, HT, SI, PC, SJ) <-> Cztery Tasmy Ramienne (Arm Lines)
"""
from __future__ import annotations

from typing import Any

# Tasmy Myersa z pelnym profilem biomechanicznym
MYERS_LINES: dict[str, dict[str, Any]] = {
    "SBL": {
        "id": "SBL",
        "name_pl": "Taśma Powierzchowna Tylna (Superficial Back Line)",
        "name_en": "Superficial Back Line",
        "jing_jin_channel": "BL",
        "jing_jin_name": "Zu Taiyang Jing Jin (Pęcherza Moczowego)",
        "color": "#1d4ed8",
        "tracks": [
            "Rozcięgno podeszwowe i zginacze krótkie palców",
            "Ścięgno Achillesa i m. brzuchaty łydki / płaszczkowaty",
            "Grupa kulszowo-goleniowa (biceps femoris, semitendinosus, semimembranosus)",
            "Więzadło krzyżowo-guzowe i powięź piersiowo-lędźwiowa",
            "Mięsień prostownik grzbietu (erector spinae)",
            "Czepiec ścięgnisty (galea aponeurotica) i łuk brwiowy",
        ],
        "muscle_patterns": [
            "gastrocnemius", "soleus", "biceps-femoris", "semitendinosus",
            "semimembranosus", "longissimus", "iliocostalis", "spinalis"
        ],
        "fascia_patterns": [
            "plantar-aponeurosis", "epicranial-aponeurosis"
        ],
        "function": "Utrzymanie wyprostu postawy pionowej, ochrona całej osi grzbietowej, amortyzacja zgięcia ku przodowi.",
        "dysfunction": "Ograniczenie skłonu w przód, przykurcz łydki i ścięgna Achillesa, ból lędźwi, bóle głowy od potylicy ku oczodołom.",
    },
    "SFL": {
        "id": "SFL",
        "name_pl": "Taśma Powierzchowna Przednia (Superficial Front Line)",
        "name_en": "Superficial Front Line",
        "jing_jin_channel": "ST",
        "jing_jin_name": "Zu Yangming Jing Jin (Żołądka)",
        "color": "#d97706",
        "tracks": [
            "Grzbietowe ścięgna prostowników palców i m. piszczelowy przedni",
            "Więzadło rzepki i guzowatość piszczeli",
            "Mięsień prosty uda (rectus femoris)",
            "Mięsień prosty brzucha (rectus abdominis) i powięź mostkowa",
            "Mięsień mostkowo-obojczykowo-sutkowy (MOS) i powięź czaszki",
        ],
        "muscle_patterns": [
            "tibialis-anterior", "rectus-femoris", "rectus-abdominis",
            "sternocleidomastoid"
        ],
        "fascia_patterns": [
            "sheath"
        ],
        "function": "Równoważenie taśmy tylnej, ochrona trzewi brzusznych, wyprost kolana i zgięcie tułowia.",
        "dysfunction": "Zapadnięcie klatki piersiowej, protrakcja głowy, przykurcz zginaczy biodra, przeprost w kolanach, refluks.",
    },
    "LL": {
        "id": "LL",
        "name_pl": "Taśma Boczna (Lateral Line)",
        "name_en": "Lateral Line",
        "jing_jin_channel": "GB",
        "jing_jin_name": "Zu Shaoyang Jing Jin (Pęcherzyka Żółciowego)",
        "color": "#16a34a",
        "tracks": [
            "Mięśnie strzałkowe (fibularis longus / brevis)",
            "Więzadło poboczne strzałkowe i pasmo biodrowo-piszczelowe (ITB)",
            "Mięsień napinacz powięzi szerokiej (TFL) i m. pośladkowy wielki / średni",
            "Mięśnie skośne brzucha (obliquus externus / internus)",
            "Mięśnie międzyżebrowe zewnętrzne i wewnętrzne",
            "Mięsień płatowaty głowy i MOS",
        ],
        "muscle_patterns": [
            "fibularis", "gluteus-medius", "gluteus-maximus", "tensor-fasciae-latae",
            "obliquus-externus", "obliquus-internus", "splenius"
        ],
        "fascia_patterns": [
            "fascia"
        ],
        "function": "Stabilizacja boczna miednicy i tułowia, kontrola rotacji osiowej, zapobieganie koślawieniu kolana i pronacji stopy.",
        "dysfunction": "Asymetria miednicy (skrót pozorny nogi), ból pasma ITB, ból krętarza, migreny skroniowe, skolioza.",
    },
    "DFL": {
        "id": "DFL",
        "name_pl": "Taśma Głęboka Przednia (Deep Front Line) – Rdzeń Biomechaniczny",
        "name_en": "Deep Front Line",
        "jing_jin_channel": "KI",
        "jing_jin_name": "Zu Shaoyin / Taiyin / Jueyin (Nerki / Śledziona / Wątroba)",
        "color": "#7c3aed",
        "tracks": [
            "Podeszwowa strona stopy i m. piszczelowy tylny / zginacz długi palców",
            "Przegroda międzymięśniowa przyśrodkowa i mięśnie przywodziciele uda",
            "Dno miednicy, powięź zasłonowa i kość guziczna",
            "Powięź przedkręgowa, mięsień lędźwiowy większy (psoas) i m. biodrowy",
            "Przepona oddechowa (crura) i osierdzie / powięź wewnątrzpiersiowa",
            "Mięśnie podgnykowe, nagłośnia i mięśnie dna jamy ustnej",
        ],
        "muscle_patterns": [
            "tibialis-posterior", "adductor-magnus", "adductor-longus",
            "psoas-major", "iliacus", "gracilis"
        ],
        "fascia_patterns": [
            "tendinous-arch-of-levator-ani"
        ],
        "function": "Osiowy rdzeń tensegracyjny ciała (core), stabilizacja oddechu i trzewi, podparcie łuku podłużnego stopy.",
        "dysfunction": "Zaburzenia oddechu, ucisk w klatce, płaskostopie poprzeczne, opadanie narządów miednicy, lęki, ból pachwiny.",
    },
    "DFAL": {
        "id": "DFAL",
        "name_pl": "Taśma Przednia Głęboka Ramienia (Deep Front Arm Line)",
        "name_en": "Deep Front Arm Line",
        "jing_jin_channel": "LU",
        "jing_jin_name": "Shou Taiyin Jing Jin (Płuc)",
        "color": "#0284c7",
        "tracks": [
            "Mięsień kłębu kciuka i powięź kłębu",
            "Powięź przednia przedramienia i m. ramienno-promieniowy",
            "Mięsień dwugłowy ramienia (biceps brachii)",
            "Mięsień piersiowy mniejszy (pectoralis minor) i wyrostek kruczy",
        ],
        "muscle_patterns": [
            "biceps-brachii", "brachioradialis", "pectoralis-minor", "abductor-pollicis"
        ],
        "fascia_patterns": ["palmar-aponeurosis"],
        "function": "Kontrola kciuka, chwyt precyzyjny, supinacja przedramienia.",
        "dysfunction": "Przykurcz klatki, bark protrakcyjny, zespół cieśni nadgarstka, ból kłębu.",
    },
    "SFAL": {
        "id": "SFAL",
        "name_pl": "Taśma Przednia Powierzchowna Ramienia (Superficial Front Arm Line)",
        "name_en": "Superficial Front Arm Line",
        "jing_jin_channel": "LI",
        "jing_jin_name": "Shou Yangming Jing Jin (Jelita Grubego)",
        "color": "#06b6d4",
        "tracks": [
            "Troczek zginaczy i mięśnie zginacze nadgarstka / palców",
            "Nadkłykieć przyśrodkowy kości ramiennej",
            "Mięsień piersiowy większy (pectoralis major) i m. najszerszy grzbietu",
        ],
        "muscle_patterns": [
            "pectoralis-major", "pronator-teres", "flexor-carpi"
        ],
        "fascia_patterns": ["palmar-aponeurosis"],
        "function": "Mocny chwyt dłoni, przywiedzenie i rotacja wewnętrzna ramienia.",
        "dysfunction": "Ból przedniego brzegu barku, łokieć golfisty, ograniczenie odwiedzenia ramienia.",
    },
    "DBAL": {
        "id": "DBAL",
        "name_pl": "Taśma Tylna Głęboka Ramienia (Deep Back Arm Line)",
        "name_en": "Deep Back Arm Line",
        "jing_jin_channel": "SJ",
        "jing_jin_name": "Shou Shaoyang Jing Jin (Potrójnego Ogrzewacza)",
        "color": "#f43f5e",
        "tracks": [
            "Kłębik małego palca i prostowniki palców",
            "Mięsień trójgłowy ramienia (triceps brachii)",
            "Mięsień podgrzebieniowy, obły mniejszy i m. równoległoboczne",
        ],
        "muscle_patterns": [
            "triceps-brachii", "infraspinatus", "rhomboideus", "extensor-digitorum"
        ],
        "fascia_patterns": [],
        "function": "Stabilizacja łopatki na klatce piersiowej, wyprost łokcia.",
        "dysfunction": "Sztywność międzyłopatkowa, ból tylnej torebki barku, łokieć tenisisty.",
    },
    "SBAL": {
        "id": "SBAL",
        "name_pl": "Taśma Tylna Powierzchowna Ramienia (Superficial Back Arm Line)",
        "name_en": "Superficial Back Arm Line",
        "jing_jin_channel": "SI",
        "jing_jin_name": "Shou Taiyang Jing Jin (Jelita Cienkiego)",
        "color": "#ef4444",
        "tracks": [
            "Brzeg łokciowy dłoni i prostownik łokciowy nadgarstka",
            "Powięź łopatkowa i tylna część m. naramiennego",
            "Mięsień czworoboczny grzbietu (trapezius)",
        ],
        "muscle_patterns": [
            "trapezius", "deltoid", "extensor-carpi-ulnaris"
        ],
        "fascia_patterns": [],
        "function": "Odwiedzenie i uniesienie obręczy barkowej, wyprost karku.",
        "dysfunction": "Sztywność karku, ból kąta łopatki, ból promieniujący do potylicy.",
    },
}

# 12 Meridianow Sciegnisto-Miesniowych (Jing Jin)
JING_JIN_CHANNELS: dict[str, dict[str, Any]] = {
    "BL": {
        "code": "BL",
        "name_pl": "Zu Taiyang Jing Jin (Ścięgnisto-Mięśniowy Pęcherza Moczowego)",
        "name_zh": "足太陽經筋",
        "myers_line": "SBL",
        "jing_well_point": "BL-67",
        "binding_points": ["BL-40 (kolano)", "BL-10 (potylica)", "BL-1 (kącik oka)"],
        "clinical_sign": "Skurcz ścięgna Achillesa, sztywność prostowników kręgosłupa, niemożność zgięcia ku przodowi, odgięcie głowy w tył (opis z Lingshu rozdz. 13).",
    },
    "GB": {
        "code": "GB",
        "name_pl": "Zu Shaoyang Jing Jin (Ścięgnisto-Mięśniowy Pęcherzyka Żółciowego)",
        "name_zh": "足少陽經筋",
        "myers_line": "LL",
        "jing_well_point": "GB-44",
        "binding_points": ["GB-34 (główka strzałki)", "GB-29 (biodro)", "GB-20 (kąt żuchwy/potylica)"],
        "clinical_sign": "Kurcze mięśni strzałkowych, ból zewnętrznej strony kolana, sztywność bocznej strony tułowia, niemożność otwarcia oka lub skrętu szyi.",
    },
    "ST": {
        "code": "ST",
        "name_pl": "Zu Yangming Jing Jin (Ścięgnisto-Mięśniowy Żołądka)",
        "name_zh": "足陽明經筋",
        "myers_line": "SFL",
        "jing_well_point": "ST-45",
        "binding_points": ["ST-41 (staw skokowy)", "ST-35 / ST-36 (rzepka)", "ST-12 (nadobojczyk)", "ST-5 / ST-6 (żuchwa)"],
        "clinical_sign": "Przykurcz mięśni goleni, przeprost lub ból z przodu kolana, skurcz mięśni uda, skurcz powłok brzusznych, opadanie kącika ust (paraliż twarzy).",
    },
    "SP": {
        "code": "SP",
        "name_pl": "Zu Taiyin Jing Jin (Ścięgnisto-Mięśniowy Śledziony)",
        "name_zh": "足太陰經筋",
        "myers_line": "DFL",
        "jing_well_point": "SP-1",
        "binding_points": ["SP-5 (kostka przyśrodkowa)", "SP-9 (kłykieć piszczeli)", "SP-12 (pachwina)", "REN-3 (spojenie łonowe)"],
        "clinical_sign": "Skurcz palucha, ból wewnętrznej strony kostki i kolana, ból przywodzicieli, skurcz podbrzusza promieniujący do kręgosłupa.",
    },
    "KI": {
        "code": "KI",
        "name_pl": "Zu Shaoyin Jing Jin (Ścięgnisto-Mięśniowy Nerek)",
        "name_zh": "足少陰經筋",
        "myers_line": "DFL",
        "jing_well_point": "KI-1",
        "binding_points": ["KI-3 / KI-6 (pod kostką)", "KI-10 (dół podkolanowy)", "REN-1 (krocze)", "kręgosłup lędźwiowy"],
        "clinical_sign": "Kurcze podeszwy stopy, sztywność i ból wewnętrznej strony podudzia i uda, skurcz mięśni lędźwiowych uniemożliwiający wyprost.",
    },
    "LR": {
        "code": "LR",
        "name_pl": "Zu Jueyin Jing Jin (Ścięgnisto-Mięśniowy Wątroby)",
        "name_zh": "足厥陰經筋",
        "myers_line": "DFL",
        "jing_well_point": "LR-1",
        "binding_points": ["LR-4 (przed kostką)", "LR-8 (przyśrodkowe kolano)", "REN-2 (narządy rodne)"],
        "clinical_sign": "Ból palucha, ból przywodzicieli pachwiny, skurcz mięśni dna miednicy, zaburzenia funkcji narządów rodnych.",
    },
    "LU": {
        "code": "LU",
        "name_pl": "Shou Taiyin Jing Jin (Ścięgnisto-Mięśniowy Płuc)",
        "name_zh": "手太陰經筋",
        "myers_line": "DFAL",
        "jing_well_point": "LU-11",
        "binding_points": ["LU-9 (nadgarstek)", "LU-5 (zgięcie łokcia)", "LU-1 (dół podobojczykowy)"],
        "clinical_sign": "Kurcz kłębu kciuka, ból ścięgien przedramienia, ból z przodu barku i pod obojczykiem, ból klatki piersiowej przy głębokim wdechu.",
    },
    "LI": {
        "code": "LI",
        "name_pl": "Shou Yangming Jing Jin (Ścięgnisto-Mięśniowy Jelita Grubego)",
        "name_zh": "手陽明經筋",
        "myers_line": "SFAL",
        "jing_well_point": "LI-1",
        "binding_points": ["LI-5 (tabakierka anatomiczna)", "LI-11 (zgięcie łokcia)", "LI-15 (wyrostek barkowy)", "LI-20 (skrzydełko nosa)"],
        "clinical_sign": "Niemożność uniesienia ramienia, ból wyrostka barkowego, sztywność karku, ból kącika ust i zębów.",
    },
    "HT": {
        "code": "HT",
        "name_pl": "Shou Shaoyin Jing Jin (Ścięgnisto-Mięśniowy Serca)",
        "name_zh": "手少陰經筋",
        "myers_line": "DBAL",
        "jing_well_point": "HT-9",
        "binding_points": ["HT-7 (kość grochowata)", "HT-3 (przyśrodkowy łokieć)", "pacha", "mostek"],
        "clinical_sign": "Kurcze mięśni dłoni po stronie małego palca, ból pachy, uczucie ucisku i sztywności w klatce piersiowej.",
    },
    "SI": {
        "code": "SI",
        "name_pl": "Shou Taiyang Jing Jin (Ścięgnisto-Mięśniowy Jelita Cienkiego)",
        "name_zh": "手太陽經筋",
        "myers_line": "SBAL",
        "jing_well_point": "SI-1",
        "binding_points": ["SI-5 (nadgarstek)", "SI-8 (wyrostek łokciowy)", "SI-9..SI-11 (łopatka)", "ząb żuchwy i ucho"],
        "clinical_sign": "Ból małego palca, ból tylnej torebki stawu ramiennego, sztywność łopatki, ból karku z promieniowaniem do ucha.",
    },
    "PC": {
        "code": "PC",
        "name_pl": "Shou Jueyin Jing Jin (Ścięgnisto-Mięśniowy Osierdzia)",
        "name_zh": "手厥陰經筋",
        "myers_line": "DFAL",
        "jing_well_point": "PC-9",
        "binding_points": ["PC-7 (środek nadgarstka)", "PC-3 (dół łokciowy)", "pacha", "żebra"],
        "clinical_sign": "Przykurcz palca środkowego, ból ścięgien zginaczy przedramienia, ból pachy i ucisk podżebrza.",
    },
    "SJ": {
        "code": "SJ",
        "name_pl": "Shou Shaoyang Jing Jin (Ścięgnisto-Mięśniowy Potrójnego Ogrzewacza)",
        "name_zh": "手少陽經筋",
        "myers_line": "DBAL",
        "jing_well_point": "SJ-1",
        "binding_points": ["SJ-4 (grzbiet nadgarstka)", "SJ-10 (łokieć)", "SJ-14 (bark)", "skroń i kąt żuchwy"],
        "clinical_sign": "Ból palca serdecznego, ból bocznego łokcia, skurcz mięśni szyi i boku głowy, trudności w połykaniu.",
    },
}


def get_jing_jin(channel_code: str) -> dict[str, Any] | None:
    """Zwraca dane meridianu Jing Jin dla danego kanalu (np. 'BL', 'GB', 'ST')."""
    return JING_JIN_CHANNELS.get(channel_code.upper().strip())


def get_myers_line(line_code: str) -> dict[str, Any] | None:
    """Zwraca dane tasmy Myersa (np. 'SBL', 'SFL', 'LL', 'DFL')."""
    return MYERS_LINES.get(line_code.upper().strip())


def analyze_stroke_fascia(matched_points: list[dict[str, Any]], dominant_channel: str) -> dict[str, Any]:
    """Analizuje narysowane pasmo pod katem naciagniecia tasmy Myersa i Jing Jin."""
    jj = get_jing_jin(dominant_channel)
    myers_id = jj["myers_line"] if jj else None
    myers = get_myers_line(myers_id) if myers_id else None

    return {
        "channel": dominant_channel,
        "jing_jin_name": jj["name_pl"] if jj else f"Pasmo {dominant_channel}",
        "myers_id": myers_id,
        "myers_name": myers["name_pl"] if myers else "Brak bezpośredniego odpowiednika",
        "myers_color": myers["color"] if myers else "#94a3b8",
        "tracks": myers["tracks"] if myers else [],
        "function": myers["function"] if myers else "",
        "dysfunction": myers["dysfunction"] if myers else (jj["clinical_sign"] if jj else ""),
        "suggested_muscles": myers["muscle_patterns"] if myers else [],
    }
