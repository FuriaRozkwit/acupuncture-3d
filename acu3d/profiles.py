"""Profile konfiguracji punktow - edycje uzytkownika obok fabrycznej, nie zamiast niej.

Fabryczna konfiguracja to `data/points.json` i `data/meridians.json`. Edytor nigdy
ich nie nadpisuje. Kazdy profil to NAKLADKA: tylko te punkty, ktore ktos przesunal,
z kotwica fabryczna z chwili edycji. Dzieki temu:

- poprawki fabryczne (np. z nowej wersji repo) nadal docieraja do punktow,
  ktorych nikt nie ruszal;
- gdy fabryka poprawi punkt, ktory ma tez wersje w profilu, widac to wprost
  ("fabryczna zmienila sie od edycji") zamiast cichego zamrozenia starego stanu;
- powrot do fabrycznej to wybor profilu, nie odtwarzanie pliku z kopii.

Aktywny profil (ten, ktory laduje sie domyslnie) zapisany jest w
`data/profiles/aktywny.txt`. Brak pliku albo nazwa "fabryczna" = sama fabryka.

    python -m acu3d.profile lista       zarzadzanie z wiersza polecen
"""
from __future__ import annotations

import copy
import hashlib
import json
import re
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
FACTORY_POINTS = DATA / "points.json"
PROFILES = DATA / "profiles"
ACTIVE = PROFILES / "aktywny.txt"
FACTORY = "fabryczna"
DEFAULT_NEW = "moja"               # profil zakladany przy pierwszym zapisie z fabrycznej


def slug(name: str) -> str:
    s = re.sub(r"[^a-z0-9_-]+", "-", name.strip().lower()).strip("-")
    if not s:
        raise ValueError("pusta nazwa profilu")
    return s


def path(name: str) -> Path:
    return PROFILES / f"{slug(name)}.json"


def names() -> list:
    return sorted(p.stem for p in PROFILES.glob("*.json")) if PROFILES.exists() else []


def exists(name: str) -> bool:
    return name == FACTORY or path(name).exists()


def active_name() -> str:
    try:
        n = ACTIVE.read_text(encoding="utf-8").strip()
    except OSError:
        return FACTORY
    return n if n and exists(n) else FACTORY


def set_active(name: str) -> None:
    if not exists(name):
        raise KeyError(f"nie ma profilu '{name}'")
    PROFILES.mkdir(parents=True, exist_ok=True)
    ACTIVE.write_text((FACTORY if name == FACTORY else slug(name)) + "\n", encoding="utf-8")


def load(name: str) -> dict:
    if name == FACTORY:
        return {"name": FACTORY, "points": {}}
    return json.loads(path(name).read_text(encoding="utf-8"))


def _write(doc: dict) -> None:
    PROFILES.mkdir(parents=True, exist_ok=True)
    doc["updated"] = time.strftime("%Y-%m-%d %H:%M")
    doc["points"] = dict(sorted(doc["points"].items(), key=lambda kv: _order(kv[0])))
    path(doc["name"]).write_text(json.dumps(doc, ensure_ascii=False, indent=1) + "\n",
                                 encoding="utf-8")


def _order(code: str):
    ch, _, n = code.rpartition("-")
    return (ch, int(n) if n.isdigit() else 0)


def create(name: str, copy_from: str | None = None) -> str:
    """Nowy profil - pusty albo kopia innego (takze kopia fabrycznej = pusty)."""
    n = slug(name)
    if n == FACTORY:
        raise ValueError("'fabryczna' to nazwa zastrzezona")
    if path(n).exists():
        raise FileExistsError(f"profil '{n}' juz istnieje")
    pts = {} if copy_from in (None, FACTORY) else copy.deepcopy(load(copy_from)["points"])
    _write({"name": n, "created": time.strftime("%Y-%m-%d %H:%M"),
            "copied_from": copy_from or FACTORY, "points": pts})
    return n


def delete(name: str) -> None:
    if name == FACTORY:
        raise ValueError("fabrycznej nie da sie usunac")
    path(name).unlink()
    if active_name() == FACTORY and ACTIVE.exists():
        ACTIVE.unlink()


def factory_points() -> list:
    return json.loads(FACTORY_POINTS.read_text(encoding="utf-8"))["points"]


def factory_point(code: str) -> dict:
    return copy.deepcopy(next(p for p in factory_points() if p["code"] == code))


def factory_anchor(code: str) -> dict:
    return factory_point(code)["anchor"]


def apply(points: list, name: str) -> dict:
    """Naklada profil na liste punktow (w miejscu). Zwraca, co sie zmienilo.

    `stale` to punkty, ktorych kotwica FABRYCZNA zmienila sie od czasu edycji -
    wtedy profil nadal wygrywa (to byla swiadoma decyzja), ale warto to zobaczyc.
    """
    info = {"name": name, "changed": [], "stale": [], "missing": []}
    if name == FACTORY:
        return info
    by = {p["code"]: p for p in points}
    for code, rec in load(name)["points"].items():
        p = by.get(code)
        if p is None:
            info["missing"].append(code)
            continue
        if rec.get("factory_anchor") is not None and rec["factory_anchor"] != p["anchor"]:
            info["stale"].append(code)
        p["anchor"] = copy.deepcopy(rec["anchor"])
        if rec.get("source"):
            p["az_source"] = rec["source"]
        p["_profile"] = name
        info["changed"].append(code)
    return info


def save_point(name: str, code: str, anchor: dict, source: str) -> str:
    """Zapisuje kotwice punktu do profilu. Gdy rowna fabrycznej - usuwa nadpisanie.

    Zwraca nazwe profilu, do ktorego trafil zapis. Z fabrycznej nie zapisujemy
    nigdy: zakladamy (albo bierzemy istniejacy) profil DEFAULT_NEW i robimy go
    domyslnym - wtedy edytowana wersja wlacza sie przy nastepnym uruchomieniu.
    """
    if name == FACTORY:
        name = DEFAULT_NEW if path(DEFAULT_NEW).exists() else create(DEFAULT_NEW)
        set_active(name)
    doc = load(name)
    fa = factory_anchor(code)
    if anchor == fa:
        doc["points"].pop(code, None)                 # powrot do fabrycznej
    else:
        doc["points"][code] = {"anchor": copy.deepcopy(anchor), "source": source,
                               "factory_anchor": fa,
                               "edited": time.strftime("%Y-%m-%d %H:%M")}
    _write(doc)
    return doc["name"]


def remove_point(name: str, code: str) -> bool:
    if name == FACTORY:
        return False
    doc = load(name)
    hit = doc["points"].pop(code, None) is not None
    if hit:
        _write(doc)
    return hit


def fingerprint(name: str) -> str:
    """Skrot tresci profilu - do kluczy cache (inny profil = inne przebiegi)."""
    if name == FACTORY:
        return FACTORY
    pts = load(name)["points"]
    return hashlib.md5(json.dumps({k: v["anchor"] for k, v in pts.items()},
                                  sort_keys=True).encode()).hexdigest()[:10]
