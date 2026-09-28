"""Tłumaczenie i polonizacja opisów anatomicznych oraz klinicznych punktów (Deadman/TCM).

Zapewnia spójne, profesjonalne polskie słownictwo anatomiczne, terminologię
fizjonomii twarzy i czaszki oraz bezpieczne wskazania kliniczne.
Wykorzystuje prekompilowaną bazę data/point_meta_pl.json z natychmiastowym dostępem O(1).
"""
from __future__ import annotations

import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PL_META_FILE = ROOT / "data" / "point_meta_pl.json"

PL_META: dict[str, dict] = {}
if PL_META_FILE.exists():
    try:
        with open(PL_META_FILE, "r", encoding="utf-8") as f:
            PL_META = json.load(f)
    except Exception:
        PL_META = {}


def normalize_code(code: str) -> str:
    """Oczyszcza kod z oznaczeń stron i formatowania (np. 'SP-14·L' -> 'SP-14')."""
    if not code:
        return ""
    c = code.replace("·L", "").replace("·R", "").replace(":L", "").replace(":R", "").strip()
    return c


def get_point_pl(code: str) -> dict | None:
    """Pobiera kompletny polski rekord punktu z bazy data/point_meta_pl.json."""
    c = normalize_code(code)
    if c in PL_META:
        return PL_META[c]
    # Sprawdz aliasy TE <-> SJ oraz LR <-> LIV
    alt = c.replace("TE-", "SJ-") if c.startswith("TE-") else (c.replace("SJ-", "TE-") if c.startswith("SJ-") else c)
    if alt in PL_META:
        return PL_META[alt]
    alt2 = c.replace("LIV-", "LR-") if c.startswith("LIV-") else (c.replace("LR-", "LIV-") if c.startswith("LR-") else c)
    if alt2 in PL_META:
        return PL_META[alt2]
    return None


def polish_point_name(code: str, default_en: str = "") -> str:
    """Zwraca kanoniczną polską nazwę punktu z transkrypcją pinyin."""
    rec = get_point_pl(code)
    if rec and rec.get("nazwa_pl"):
        return rec["nazwa_pl"]
    return default_en or code


def _find_code_and_fallback(arg1: str, arg2: str = "") -> tuple[str, str]:
    """Rozpoznaje, który argument jest kodem punktu (np. 'LU-1', 'ST-36'), a który tekstem opisu."""
    s1 = str(arg1 or "").strip()
    s2 = str(arg2 or "").strip()
    is_code1 = bool(re.match(r"^[A-Za-z]{1,4}[-·:]\d+", s1))
    is_code2 = bool(re.match(r"^[A-Za-z]{1,4}[-·:]\d+", s2))
    if is_code1:
        return normalize_code(s1), s2
    elif is_code2:
        return normalize_code(s2), s1
    return normalize_code(s1), s2


def translate_location(code_or_text: str, fallback_text: str = "", *args, **kwargs) -> str:
    """Zwraca polski opis lokalizacji anatomicznej punktu."""
    c, raw = _find_code_and_fallback(code_or_text, fallback_text)
    if c:
        rec = get_point_pl(c)
        if rec and rec.get("lokalizacja_pl"):
            return rec["lokalizacja_pl"]

    out = raw or code_or_text or ""
    # Oczyszczenie smieci OCR
    out = re.sub(r"[~一—]{2,}", " ", out)
    out = re.sub(r"\b[A-Z]{1,3}-\d+\s+[A-Z]{1,3}-\d+\b", "", out)
    out = re.sub(r"\s{2,}", " ", out)
    return out.strip()


def translate_palpation(code_or_text: str, fallback_text: str = "", *args, **kwargs) -> str:
    """Zwraca polską wskazówkę palpacyjną i weryfikację położenia."""
    c, raw = _find_code_and_fallback(code_or_text, fallback_text)
    if c:
        rec = get_point_pl(c)
        if rec and rec.get("palpacja_pl"):
            return rec["palpacja_pl"]
    return raw.strip() if raw else ""


def translate_needling(code_or_text: str, fallback_text: str = "", *args, **kwargs) -> str:
    """Zwraca polski opis techniki nakłucia, głębokości i ostrzeżenia kliniczne."""
    c, raw = _find_code_and_fallback(code_or_text, fallback_text)
    if c:
        rec = get_point_pl(c)
        if rec and rec.get("naklucie_pl"):
            return rec["naklucie_pl"]
    return raw.strip() if raw else ""


def translate_action(action: str | None, code: str = "", *args, **kwargs) -> str:
    """Zwraca polskie działanie TCM punktu."""
    if not action:
        return ""
    action_str = str(action).strip()
    c, _ = _find_code_and_fallback(code, action_str)
    if c:
        rec = get_point_pl(c)
        if rec and rec.get("dzialania_pl"):
            for act in rec["dzialania_pl"]:
                if action_str.lower() in act.lower() or act.lower() in action_str.lower():
                    return act
    return action_str


def translate_indications(code_or_text: str, fallback_text: str = "", *args, **kwargs) -> str:
    """Zwraca polskie wskazania kliniczne i patologie leczone danym punktem."""
    c, raw = _find_code_and_fallback(code_or_text, fallback_text)
    if c:
        rec = get_point_pl(c)
        if rec and rec.get("wskazania_pl"):
            return rec["wskazania_pl"]
    return raw.strip() if raw else ""


def translate_channel(channel_eng: str) -> str:
    """Zwraca polską nazwę kanału (meridianu)."""
    mapping = {
        "Lung": "Kanał Płuc (Shou Tai Yin)",
        "Large Intestine": "Kanał Jelita Grubego (Shou Yang Ming)",
        "Stomach": "Kanał Żołądka (Zu Yang Ming)",
        "Spleen": "Kanał Śledziony (Zu Tai Yin)",
        "Heart": "Kanał Serca (Shou Shao Yin)",
        "Small Intestine": "Kanał Jelita Cienkiego (Shou Tai Yang)",
        "Bladder": "Kanał Pęcherza Moczowego (Zu Tai Yang)",
        "Kidney": "Kanał Nerek (Zu Shao Yin)",
        "Pericardium": "Kanał Osierdzia (Shou Jue Yin)",
        "Sanjiao": "Kanał Potrójnego Ogrzewacza (Shou Shao Yang)",
        "Gall Bladder": "Kanał Pęcherzyka Żółciowego (Zu Shao Yang)",
        "Liver": "Kanał Wątroby (Zu Jue Yin)",
        "Conception vessel": "Naczynie Poczęcia (Ren Mai)",
        "Governing vessel": "Naczynie Rządzące (Du Mai)",
    }
    return mapping.get(channel_eng, channel_eng)
