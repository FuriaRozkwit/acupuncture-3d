"""Czy do repozytorium nie wsiakl tekst z data/point_meta.json (material chroniony).

Sprawdza ZAWARTOSC, nie nazwy plikow: bierze wszystkie ciagi 8 kolejnych slow z
point_meta i szuka ich w plikach, ktore ida do repozytorium. Tak wlasnie wyszedl
na jaw wyciek w polu `categories` - kontrola nazw plikow go nie widziala.

    python tools/scan_leak.py              # pliki sledzone + nowe (bez .gitignore)
    python tools/scan_leak.py --history    # dodatkowo wszystkie bloby z historii
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
META = ROOT / "data" / "point_meta.json"
SHINGLE = 8
TEXT_SUFFIX = {".py", ".json", ".md", ".txt", ".yml", ".yaml", ".toml", ".cfg", ".csv", ".html"}


def words(text: str) -> list:
    return re.findall(r"[a-z]+", text.lower())


def shingles(text: str) -> set:
    w = words(text)
    return {" ".join(w[i:i + SHINGLE]) for i in range(len(w) - SHINGLE + 1)}


def meta_shingles() -> set:
    doc = json.loads(META.read_text(encoding="utf-8"))
    out = set()
    for rec in doc.values():
        for value in rec.values():
            if isinstance(value, str):
                out |= shingles(value)
            elif isinstance(value, list):
                out |= shingles(" ".join(str(v) for v in value))
    return out


def repo_files() -> list:
    cmd = ["git", "-C", str(ROOT), "ls-files", "--cached", "--others", "--exclude-standard"]
    names = subprocess.run(cmd, capture_output=True, text=True).stdout.split("\n")
    return [ROOT / n for n in names if n and Path(n).suffix.lower() in TEXT_SUFFIX]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--history", action="store_true", help="przeszukaj tez bloby z historii git")
    args = ap.parse_args()
    if not META.exists():
        print("brak data/point_meta.json - nie ma z czym porownywac", file=sys.stderr)
        return 2
    ref = meta_shingles()
    print(f"fragmentow odniesienia (po {SHINGLE} slow): {len(ref)}")

    bad = 0
    for f in repo_files():
        if f.resolve() == META.resolve():
            continue
        try:
            text = f.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            continue
        hit = shingles(text) & ref
        if hit:
            bad += 1
            print(f"WYCIEK: {f.relative_to(ROOT)} - {len(hit)} fragmentow, np. \"{sorted(hit)[0]}\"")
    print(f"plikow sprawdzonych: {len(repo_files())}, z tekstem z podrecznika: {bad}")

    if args.history:
        blobs = subprocess.run(["git", "-C", str(ROOT), "rev-list", "--objects", "--all"],
                               capture_output=True, text=True).stdout.split("\n")
        seen = 0
        for line in blobs:
            parts = line.split(" ", 1)
            if len(parts) != 2 or Path(parts[1]).suffix.lower() not in TEXT_SUFFIX:
                continue
            body = subprocess.run(["git", "-C", str(ROOT), "cat-file", "-p", parts[0]],
                                  capture_output=True, text=True).stdout
            seen += 1
            hit = shingles(body) & ref
            if hit:
                bad += 1
                print(f"WYCIEK W HISTORII: {parts[1]} ({parts[0][:8]}) - {len(hit)} fragmentow")
        print(f"blobow w historii sprawdzonych: {seen}")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
