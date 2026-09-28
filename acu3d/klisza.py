"""Moduł Kliszy Klinicznej (Klisza): Zapis, analiza i obróbka łączeń napięciowych.

Klisza to jednostka zapisu diagnostycznego w modelu 3D:
1. Przechowuje narysowane na ciele linie napięć (trajektorie 3D w metrach i cunach).
2. Wykrywa punkty akupunkturowe leżące w korytarzu narysowanego pasma, uporządkowane
   sekwencyjnie od początku do końca pociągnięcia (kierunek badania / napięcia).
3. Klasyfikuje ciągłość meridianową (zgodność z kanałami głównymi, parami Yuan-Luo,
   8 Naczyniami Nadzwyczajnymi oraz taśmami powięziowymi Jing Jin).
4. Przygotowuje dane w ustrukturyzowanym formacie grafowo-tensorowym do dalszej
   algorytmicznej obróbki w modelu (sieci neuronowe, wnioskowanie regułowe, macierze przejść).
"""
from __future__ import annotations

import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

from .translate_pl import get_point_pl, polish_point_name
from .clinical_kb import (
    get_extraordinary_vessel_info,
    get_xi_cleft_info,
    get_yuan_luo_info,
)

ROOT = Path(__file__).resolve().parent.parent
KLISZE_DIR = ROOT / "data" / "klisze"
INDEX_FILE = ROOT / "data" / "klisze.json"

DEFAULT_CORRIDOR_M = 0.032  # 32 mm promień korytarza poszukiwania punktów wokół linii
FALLBACK_CUN_M = 0.022      # domyślny przelicznik 1 cun = 22 mm gdy brak segmentu


def _safe_print(text: str) -> None:
    enc = getattr(sys.stdout, "encoding", None) or "utf-8"
    sys.stdout.write(text.encode(enc, errors="replace").decode(enc, errors="replace") + "\n")
    sys.stdout.flush()


def _point_to_segment_distance(p: np.ndarray, a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    """Zwraca (odleglosc_min, rzut_t_od_0_do_1) punktu p wzgledem odcinka ab."""
    ab = b - a
    ab_len2 = float(np.dot(ab, ab))
    if ab_len2 < 1e-12:
        return float(np.linalg.norm(p - a)), 0.0
    t = float(np.dot(p - a, ab) / ab_len2)
    t = max(0.0, min(1.0, t))
    proj = a + t * ab
    return float(np.linalg.norm(p - proj)), t


class Klisza:
    """Pojedyncza klisza diagnostyczna sesji klinicznej."""

    def __init__(
        self,
        klisza_id: str,
        name: str,
        timestamp: str | None = None,
        strokes: list[dict[str, Any]] | None = None,
        clinical_notes: str = "",
        model_height_m: float = 1.72,
    ):
        self.id = klisza_id
        self.name = name
        self.timestamp = timestamp or datetime.now(timezone.utc).isoformat(timespec="seconds")
        self.strokes = strokes or []
        self.clinical_notes = clinical_notes
        self.model_height_m = model_height_m

    @property
    def total_length_m(self) -> float:
        return sum(s.get("length_m", 0.0) for s in self.strokes)

    @property
    def total_length_cun(self) -> float:
        return sum(s.get("length_cun", 0.0) for s in self.strokes)

    @property
    def all_matched_points(self) -> list[dict[str, Any]]:
        pts = []
        for s in self.strokes:
            pts.extend(s.get("matched_points", []))
        return pts

    def to_dict(self) -> dict[str, Any]:
        return {
            "id": self.id,
            "name": self.name,
            "timestamp": self.timestamp,
            "model_height_m": round(self.model_height_m, 3),
            "total_strokes": len(self.strokes),
            "total_length_m": round(self.total_length_m, 4),
            "total_length_cun": round(self.total_length_cun, 2),
            "clinical_notes": self.clinical_notes,
            "strokes": self.strokes,
            "graph_representation": self.build_graph_representation(),
        }

    def build_graph_representation(self) -> dict[str, Any]:
        """Buduje graf do dalszej algorytmicznej obróbki w modelu."""
        nodes = []
        node_set = set()
        edges = []

        for s_idx, s in enumerate(self.strokes):
            m_pts = s.get("matched_points", [])
            for p in m_pts:
                k = p["key"]
                if k not in node_set:
                    node_set.add(k)
                    nodes.append({
                        "key": k,
                        "code": p["code"],
                        "channel": p["channel"],
                        "side": p["side"],
                        "pos": p.get("pos"),
                    })

            for i in range(len(m_pts) - 1):
                p1, p2 = m_pts[i], m_pts[i + 1]
                intra = (p1["channel"] == p2["channel"])
                d_cun = round(abs(p2["along_cun"] - p1["along_cun"]), 2)
                d_m = round(abs(p2["along_m"] - p1["along_m"]), 4)
                edges.append({
                    "stroke_idx": s_idx,
                    "source": p1["key"],
                    "target": p2["key"],
                    "dist_cun": d_cun,
                    "dist_m": d_m,
                    "intra_meridian": intra,
                    "source_channel": p1["channel"],
                    "target_channel": p2["channel"],
                })

        return {
            "nodes": nodes,
            "edges": edges,
            "num_nodes": len(nodes),
            "num_edges": len(edges),
        }

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> "Klisza":
        return cls(
            klisza_id=d.get("id", "klisza_unknown"),
            name=d.get("name", "Bez nazwy"),
            timestamp=d.get("timestamp"),
            strokes=d.get("strokes", []),
            clinical_notes=d.get("clinical_notes", ""),
            model_height_m=d.get("model_height_m", 1.72),
        )


class KliszaManager:
    """Zarządca bazy klisz diagnostycznych oraz analizy trajektorii 3D."""

    def __init__(self, atlas, wire=None, storage_dir: Path | None = None):
        self.atlas = atlas
        self.wire = wire
        self.storage_dir = Path(storage_dir) if storage_dir else KLISZE_DIR
        self.storage_dir.mkdir(parents=True, exist_ok=True)
        self.index_file = INDEX_FILE
        self._mpos = np.array([m["pos"] for m in atlas.markers], dtype=float)

    def analyze_stroke(
        self,
        points: np.ndarray | list,
        stroke_id: int = 1,
        corridor_m: float = DEFAULT_CORRIDOR_M,
    ) -> dict[str, Any]:
        """Dokładna analiza geometryczna i anatomiczna pojedynczego pociągnięcia."""
        pts = np.asarray(points, dtype=float)
        if len(pts) < 2:
            return {
                "stroke_id": stroke_id,
                "points_count": len(pts),
                "length_m": 0.0,
                "length_cun": 0.0,
                "points_3d": pts.tolist(),
                "matched_points": [],
                "matched_channels": {},
                "dominant_channel": None,
            }

        # 1. Obliczenie długości metrycznej i lokalnej skali cun
        dists = np.linalg.norm(np.diff(pts, axis=0), axis=1)
        cum_dist = np.insert(np.cumsum(dists), 0, 0.0)
        total_m = float(cum_dist[-1])

        # Skala cun wzdłuż trajektorii
        cum_cun = [0.0]
        for i in range(len(dists)):
            p_mid = 0.5 * (pts[i] + pts[i + 1])
            cun_scale = self._get_cun_scale(p_mid)
            step_cun = dists[i] / cun_scale
            cum_cun.append(cum_cun[-1] + step_cun)
        cum_cun = np.array(cum_cun, dtype=float)
        total_cun = float(cum_cun[-1])

        # 2. Kotwica początkowa i końcowa (wire describe)
        start_anchor = self._describe_anchor(pts[0])
        end_anchor = self._describe_anchor(pts[-1])

        # 3. Wyszukanie punktów w korytarzu wokół polilinii
        matched_points = []
        n_segments = len(pts) - 1

        for m_idx, marker in enumerate(self.atlas.markers):
            m_pos = self._mpos[m_idx]
            # Szybki filtr zgrubny AABB
            if (m_pos[0] < pts[:, 0].min() - corridor_m or m_pos[0] > pts[:, 0].max() + corridor_m or
                m_pos[1] < pts[:, 1].min() - corridor_m or m_pos[1] > pts[:, 1].max() + corridor_m or
                m_pos[2] < pts[:, 2].min() - corridor_m or m_pos[2] > pts[:, 2].max() + corridor_m):
                continue

            min_dist = float("inf")
            best_along_m = 0.0
            best_along_cun = 0.0

            for i in range(n_segments):
                d_seg, t = _point_to_segment_distance(m_pos, pts[i], pts[i + 1])
                if d_seg < min_dist:
                    min_dist = d_seg
                    best_along_m = cum_dist[i] + t * dists[i]
                    best_along_cun = cum_cun[i] + t * (cum_cun[i + 1] - cum_cun[i])

            if min_dist <= corridor_m:
                code = marker["code"]
                key = marker["key"]
                side = marker["side"]
                ch = marker["channel"]
                meta = self.atlas.meta.get(code, {})

                # Wzbogacenie o polską bazę Deadmana
                clean_code = code.replace("·L", "").replace("·R", "").replace(":L", "").replace(":R", "").strip()
                pl_rec = get_point_pl(clean_code)
                pl_name = (pl_rec.get("nazwa_pl") if pl_rec else None) or polish_point_name(code, meta.get("english", ""))
                pinyin = (pl_rec.get("pinyin") if pl_rec else None) or marker.get("pinyin", "")

                actions = (pl_rec.get("dzialania_pl") if pl_rec else None) or meta.get("actions", [])
                action_str = "; ".join(actions[:2]) if isinstance(actions, list) else str(actions)

                matched_points.append({
                    "key": key,
                    "code": code,
                    "side": side,
                    "channel": ch,
                    "nazwa_pl": pl_name,
                    "pinyin": pinyin,
                    "dist_mm": round(min_dist * 1000.0, 1),
                    "along_m": round(best_along_m, 4),
                    "along_cun": round(best_along_cun, 2),
                    "along_pct": round((best_along_m / max(total_m, 1e-6)) * 100.0, 1),
                    "pos": [round(float(c), 5) for c in m_pos],
                    "dzialania": action_str[:160] if action_str else "",
                    "naczynia_nadzwyczajne": get_extraordinary_vessel_info(clean_code) or [],
                    "yuan_luo": get_yuan_luo_info(clean_code) or [],
                    "xi_cleft": get_xi_cleft_info(clean_code) or "",
                })

        # Sortujemy punkty według ich pozycji na ścieżce od początku do końca!
        matched_points.sort(key=lambda p: p["along_m"])

        # 4. Tally kanałów
        chan_tally: dict[str, int] = {}
        for p in matched_points:
            c = p["channel"]
            chan_tally[c] = chan_tally.get(c, 0) + 1

        dominant_chan = max(chan_tally, key=chan_tally.get) if chan_tally else None

        return {
            "stroke_id": stroke_id,
            "points_count": len(pts),
            "length_m": round(total_m, 4),
            "length_cun": round(total_cun, 2),
            "start_anchor": start_anchor,
            "end_anchor": end_anchor,
            "points_3d": [[round(float(c), 5) for c in p] for p in pts],
            "matched_points": matched_points,
            "matched_channels": chan_tally,
            "dominant_channel": dominant_chan,
        }

    def _get_cun_scale(self, pos: np.ndarray) -> float:
        """Zwraca metryczną wielkość 1 cun (w metrach) w danym rejonie ciała."""
        if self.wire is not None:
            try:
                d = self.wire.describe(pos)
                seg = d["segment"]
                an = getattr(self.atlas, "anatomy", None)
                if an is not None and hasattr(an, "seg"):
                    sg = an.seg(seg)
                    canon_cuns = {
                        "forearm": 12.0, "upper_arm": 9.0, "hand": 7.5,
                        "shank": 16.0, "thigh": 19.0, "foot": 10.0,
                    }.get(seg)
                    if canon_cuns and sg and getattr(sg, "length", None):
                        return float(sg.length / canon_cuns)
            except Exception:
                pass
        return FALLBACK_CUN_M

    def _describe_anchor(self, pos: np.ndarray) -> dict[str, Any]:
        """Zwraca zwięzły opis anatomiczny pozycji (segment, cun, mięsień)."""
        if self.wire is None:
            return {"pos": [round(float(c), 4) for c in pos]}
        try:
            d = self.wire.describe(pos)
            mus = d.get("muscle")
            mus_name = mus.get("name") if mus else None
            return {
                "segment": d.get("segment", "ciało"),
                "cun": round(d.get("cun", 0.0), 2),
                "az": round(d.get("az", 0.0), 1),
                "flesh_mm": round(d.get("flesh_mm", 0.0), 1),
                "muscle": mus_name,
            }
        except Exception:
            return {"pos": [round(float(c), 4) for c in pos]}

    def create_klisza(
        self,
        strokes: list[np.ndarray | list],
        name: str | None = None,
        notes: str = "",
        corridor_m: float = DEFAULT_CORRIDOR_M,
    ) -> Klisza:
        """Tworzy pełną kliszę diagnostyczną ze wszystkich narysowanych pociągnięć."""
        now = datetime.now(timezone.utc)
        ts_str = now.strftime("%Y%m%d_%H%M%S")
        klisza_id = f"klisza_{ts_str}"

        analyzed_strokes = []
        for idx, s_pts in enumerate(strokes, start=1):
            analyzed_strokes.append(self.analyze_stroke(s_pts, stroke_id=idx, corridor_m=corridor_m))

        if not name:
            all_pts = []
            for s in analyzed_strokes:
                for p in s.get("matched_points", []):
                    all_pts.append(p["code"])
            if len(all_pts) >= 2:
                name = f"Klisza {ts_str[-6:]}: {all_pts[0]} -> {all_pts[-1]} ({len(all_pts)} pkt)"
            elif len(all_pts) == 1:
                name = f"Klisza {ts_str[-6:]}: rejon {all_pts[0]}"
            else:
                name = f"Klisza {ts_str[-6:]}: {len(strokes)} pasm(a)"

        return Klisza(
            klisza_id=klisza_id,
            name=name,
            timestamp=now.isoformat(timespec="seconds"),
            strokes=analyzed_strokes,
            clinical_notes=notes,
            model_height_m=round(float(getattr(self.atlas.rig, "H", 1.72)), 3) if hasattr(self.atlas, "rig") else 1.72,
        )

    def save_klisza(self, klisza: Klisza) -> Path:
        """Zapisuje kliszę do pliku JSON oraz aktualizuje indeks."""
        file_path = self.storage_dir / f"{klisza.id}.json"
        doc = klisza.to_dict()
        file_path.write_text(json.dumps(doc, ensure_ascii=False, indent=2), encoding="utf-8")

        # Aktualizacja indeksu data/klisze.json
        index_data = self.list_all_metadata()
        # Usun jesli juz istnieje ten sam id
        index_data = [item for item in index_data if item["id"] != klisza.id]
        index_data.insert(0, {
            "id": klisza.id,
            "name": klisza.name,
            "timestamp": klisza.timestamp,
            "file": f"klisze/{klisza.id}.json",
            "total_strokes": len(klisza.strokes),
            "total_length_m": round(klisza.total_length_m, 4),
            "total_length_cun": round(klisza.total_length_cun, 2),
            "points_count": len(klisza.all_matched_points),
            "points_summary": [p["code"] for p in klisza.all_matched_points[:8]],
        })
        self.index_file.write_text(json.dumps(index_data, ensure_ascii=False, indent=2), encoding="utf-8")
        return file_path

    def list_all_metadata(self) -> list[dict[str, Any]]:
        """Zwraca metadane wszystkich zapisanych klisz z indeksu lub dysku."""
        if self.index_file.exists():
            try:
                return json.loads(self.index_file.read_text(encoding="utf-8"))
            except Exception:
                pass

        # Jesli brak indeksu, skanujemy katalog
        items = []
        for p in sorted(self.storage_dir.glob("klisza_*.json"), reverse=True):
            try:
                d = json.loads(p.read_text(encoding="utf-8"))
                items.append({
                    "id": d["id"],
                    "name": d.get("name", p.stem),
                    "timestamp": d.get("timestamp", ""),
                    "file": f"klisze/{p.name}",
                    "total_strokes": d.get("total_strokes", len(d.get("strokes", []))),
                    "total_length_m": d.get("total_length_m", 0.0),
                    "total_length_cun": d.get("total_length_cun", 0.0),
                    "points_count": len(d.get("graph_representation", {}).get("nodes", [])),
                })
            except Exception:
                pass
        return items

    def load_klisza(self, klisza_id: str) -> Klisza | None:
        """Wczytuje kliszę o podanym ID."""
        target_path = self.storage_dir / f"{klisza_id}.json"
        if not target_path.exists():
            # sprobuj bez prefixu
            for p in self.storage_dir.glob(f"*{klisza_id}*.json"):
                target_path = p
                break
        if not target_path.exists():
            return None
        doc = json.loads(target_path.read_text(encoding="utf-8"))
        return Klisza.from_dict(doc)

    def print_klisza_report(self, klisza: Klisza) -> None:
        """Drukuje czytelny, pełny raport kliniczny w konsoli terminala."""
        out = [
            "",
            "╔" + "═" * 76 + "╗",
            f"║  KLISZA DIAGNOSTYCZNA: {klisza.name}".ljust(77) + "║",
            f"║  ID: {klisza.id:<24} Data: {klisza.timestamp}".ljust(77) + "║",
            f"║  Pasma napięć: {len(klisza.strokes)}   Długość: {klisza.total_length_m*100.0:.1f} cm ({klisza.total_length_cun:.1f} cun)".ljust(77) + "║",
            "╠" + "═" * 76 + "╣",
        ]

        for s in klisza.strokes:
            s_id = s.get("stroke_id", 1)
            l_m = s.get("length_m", 0.0)
            l_cun = s.get("length_cun", 0.0)
            dom_ch = s.get("dominant_channel", "-")
            start_anc = s.get("start_anchor", {})
            end_anc = s.get("end_anchor", {})

            out.append(f"║  PASMO #{s_id}: {l_m*100.0:.1f} cm | {l_cun:.1f} cun | Kanał wiodący: [{dom_ch}]".ljust(77) + "║")
            st_txt = f"Początek: {start_anc.get('segment','-')} {start_anc.get('cun',0):.1f}c ({start_anc.get('muscle','-')})"
            en_txt = f"Koniec: {end_anc.get('segment','-')} {end_anc.get('cun',0):.1f}c ({end_anc.get('muscle','-')})"
            out.append(f"║    {st_txt[:70]}".ljust(77) + "║")
            out.append(f"║    {en_txt[:70]}".ljust(77) + "║")

            m_pts = s.get("matched_points", [])
            if m_pts:
                out.append("╟" + "─" * 76 + "╢")
                out.append(f"║    PUNKTY PRZECIĘTE WZDŁUŻ TORU ({len(m_pts)} w korytarzu 32 mm):".ljust(77) + "║")
                for p in m_pts:
                    code = p["code"]
                    pl_name = p["nazwa_pl"]
                    pinyin = p["pinyin"]
                    dist_mm = p["dist_mm"]
                    al_cun = p["along_cun"]
                    al_pct = p["along_pct"]
                    disp = f"{code:<7} ({pinyin} • {pl_name})" if pinyin else f"{code:<7} ({pl_name})"
                    line_str = f"{disp:<36} +{al_cun:4.1f}c ({al_pct:4.1f}%) [odstęp {dist_mm:.0f}mm]"
                    out.append(f"║      * {line_str}".ljust(77) + "║")
                    if p.get("dzialania"):
                        act = p["dzialania"]
                        if len(act) > 64:
                            act = act[:64] + "..."
                        out.append(f"║        Działanie: {act}".ljust(77) + "║")
            else:
                out.append("║    (Brak punktów akupunkturowych w bezpośrednim korytarzu pasma)".ljust(77) + "║")
            out.append("╟" + "─" * 76 + "╢")

        # Podsumowanie powiązań grafowych do dalszej obróbki
        graph = klisza.build_graph_representation()
        out.append(f"║  STRUKTURA MODELU: {graph['num_nodes']} węzłów (punktów), {graph['num_edges']} krawędzi napięcia".ljust(77) + "║")
        out.append("║  Gotowe do dalszej obróbki algorytmicznej i inferencji modelu.".ljust(77) + "║")
        out.append("╚" + "═" * 76 + "╝")
        out.append("")
        _safe_print("\n".join(out))
