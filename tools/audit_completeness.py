import json
from pathlib import Path
from acu3d.model import Atlas

atlas = Atlas()
print('=== ATLAS MARKERS & POINTS ===')
print(f'Liczba wygenerowanych markerow 3D na ciele: {len(atlas.markers)}')
print(f'Liczba definicji punktow w points.json: {len(atlas.points)}')
print(f'Liczba kanalow w points.json: {len(atlas.channels)}')

canonical = {
    'LU': 11,
    'LI': 20,
    'ST': 45,
    'SP': 21,
    'HT': 9,
    'SI': 19,
    'BL': 67,
    'KI': 27,
    'PC': 9,
    'SJ': 23,
    'GB': 44,
    'LR': 14,
    'REN': 24,
    'DU': 28,
}

by_ch_points = {}
for p in atlas.points:
    ch = p.get('channel', '')
    by_ch_points[ch] = by_ch_points.get(ch, 0) + 1

by_ch_markers = {}
for m in atlas.markers:
    ch = m['channel']
    by_ch_markers[ch] = by_ch_markers.get(ch, 0) + 1

print('\nZestawienie kanalow w modelu:')
print(f"{'Kanal':<8} | {'W modelu':<10} | {'Kanonicznie':<12} | {'Kompletnosc':<12} | {'Markery 3D':<15}")
print('-' * 70)
total_model = 0
total_canon = sum(canonical.values())

for ch, canon_count in canonical.items():
    model_count = by_ch_points.get(ch, 0)
    total_model += model_count
    pct = (model_count / canon_count) * 100.0 if canon_count else 0
    mark_cnt = by_ch_markers.get(ch, 0)
    print(f"{ch:<8} | {model_count:<10} | {canon_count:<12} | {pct:5.1f}%       | {mark_cnt:<15}")

print('-' * 70)
print(f"SUMA     | {total_model:<10} | {total_canon:<12} | {total_model/total_canon*100:5.1f}%")

print('\n=== SCIEZKI MERIDIANOW (meridians.json) ===')
print(f'Liczba kanalow w meridians.json: {len(atlas.paths)}')
total_path_points = 0
for ch, branches in atlas.paths.items():
    pts_in_channel_paths = set()
    for branch in branches:
        pts_in_channel_paths.update(branch)
    canon = canonical.get(ch, 0)
    total_path_points += len(pts_in_channel_paths)
    status = "KOMPLETNY" if len(pts_in_channel_paths) == canon else f"brakuje {canon - len(pts_in_channel_paths)}"
    print(f'  Kanał {ch:<5}: {len(pts_in_channel_paths)}/{canon} pkt w ścieżce ({len(branches)} gałęzi) -> {status}')

print(f'\nŁącznie punktów w ścieżkach: {total_path_points}/{total_canon} ({total_path_points/total_canon*100:.1f}%)')

print('\n=== KOTWICE ANATOMICZNE I DOKŁADNOŚĆ POZYCJONOWANIA ===')
anchor_types = {}
accuracy_types = {}
for p in atlas.points:
    acc = p.get('accuracy', 'nieznana')
    accuracy_types[acc] = accuracy_types.get(acc, 0) + 1
    anc = p.get('anchor', {})
    if 'struct' in anc:
        t = 'strukturalna (kostna / Z-Anatomy)'
    elif 'seg' in anc:
        t = 'odcinkowa (proporcjonalna cunowa)'
    else:
        t = 'inna'
    anchor_types[t] = anchor_types.get(t, 0) + 1

print('Typy kotwiczenia punktów na ciele:')
for t, cnt in anchor_types.items():
    print(f'  - {t}: {cnt} punktów ({cnt/len(atlas.points)*100:.1f}%)')

print('Deklarowana dokładność (accuracy):')
for acc, cnt in accuracy_types.items():
    print(f'  - {acc}: {cnt} punktów ({cnt/len(atlas.points)*100:.1f}%)')

print('\n=== METADANE KLINICZNE (PL & Deadman) ===')
meta_pl_path = Path('data/point_meta_pl.json')
if meta_pl_path.exists():
    meta_pl = json.loads(meta_pl_path.read_text(encoding='utf-8'))
    print(f'Liczba punktów w bazie Deadmana PL: {len(meta_pl)}')
    with_loc = sum(1 for v in meta_pl.values() if v.get('lokalizacja_pl'))
    with_palp = sum(1 for v in meta_pl.values() if v.get('palpacja_pl'))
    with_nakl = sum(1 for v in meta_pl.values() if v.get('naklucie_pl'))
    with_acts = sum(1 for v in meta_pl.values() if v.get('dzialania_pl'))
    with_ind = sum(1 for v in meta_pl.values() if v.get('wskazania_pl'))
    print(f'  - lokalizacja anatomiczna PL: {with_loc}/{len(atlas.points)} ({with_loc/len(atlas.points)*100:.1f}%)')
    print(f'  - wskazówki palpacyjne PL: {with_palp}/{len(atlas.points)} ({with_palp/len(atlas.points)*100:.1f}%)')
    print(f'  - technika nakłucia PL: {with_nakl}/{len(atlas.points)} ({with_nakl/len(atlas.points)*100:.1f}%)')
    print(f'  - działania TCM PL: {with_acts}/{len(atlas.points)} ({with_acts/len(atlas.points)*100:.1f}%)')
    print(f'  - wskazania kliniczne PL: {with_ind}/{len(atlas.points)} ({with_ind/len(atlas.points)*100:.1f}%)')

# Sprawdzenie punktów nadzwyczajnych (Extra Points)
extra_pts = [p.get('code', '') for p in atlas.points if not any(p.get('code', '').startswith(c + '-') for c in canonical.keys())]
print(f'\nPunkty dodatkowe (Extra Points / poza 14 meridianami): {len(extra_pts)}')


