"""Import opisow punktow z repo `fiszki` (OCR Deadmana) do data/point_meta.json.

Normalizuje kody kanalow do standardu uzywanego w tej aplikacji:
    L.I. -> LI, KID -> KI, HE -> HT, P -> PC, LIV -> LR
Uruchomienie:
    python tools/import_meta.py [sciezka/do/points.json]
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_SRC = Path(r"C:\repos\fiszki\dane\points.json")

ALIAS = {"L.I.": "LI", "KID": "KI", "HE": "HT", "P": "PC", "LIV": "LR"}
KEEP = ("channel", "pinyin", "english", "categories", "location", "locationNote",
        "needling", "actions", "indications", "commentary", "combinations")


def norm_code(code: str) -> str:
    m = re.match(r"^(.*?)-(\d+)$", code.strip())
    if not m:
        return code.strip()
    ch, num = m.group(1), int(m.group(2))
    return f"{ALIAS.get(ch, ch)}-{num}"


def main() -> None:
    src = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_SRC
    if not src.exists():
        raise SystemExit(f"Brak pliku zrodlowego: {src}")
    raw = json.loads(src.read_text(encoding="utf-8-sig"))
    out = {}
    for p in raw:
        code = norm_code(p["code"])
        rec = {k: p[k] for k in KEEP if k in p and p[k]}
        rec["source"] = "Deadman, A Manual of Acupuncture (OCR)"
        out[code] = rec
    dst = ROOT / "data" / "point_meta.json"
    dst.write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"{len(out)} punktow -> {dst}")


if __name__ == "__main__":
    main()
