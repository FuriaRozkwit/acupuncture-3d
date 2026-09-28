import sys
from pathlib import Path
REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

import numpy as np
from acu3d.model import Atlas
from acu3d.wire import Wire
from acu3d.klisza import KliszaManager

atlas = Atlas()
wire = Wire()
mgr = KliszaManager(atlas, wire=wire)

# Symulujemy pociagniecie wzdluz ramienia: LU-1 -> LU-5 -> LU-7 -> LU-11
p_lu1 = atlas.by_key["LU-1·L"]["pos"]
p_lu5 = atlas.by_key["LU-5·L"]["pos"]
p_lu7 = atlas.by_key["LU-7·L"]["pos"]
p_lu11 = atlas.by_key["LU-11·L"]["pos"]

# 1. Odcinek od ramienia do lokcia
t1 = np.linspace(0, 1, 15)
s1 = [p_lu1 + (p_lu5 - p_lu1) * t for t in t1]
# 2. Od lokcia do nadgarstka i kciuka
t2 = np.linspace(0, 1, 20)[1:]
s2 = [p_lu5 + (p_lu7 - p_lu5) * t for t in t2]
t3 = np.linspace(0, 1, 10)[1:]
s3 = [p_lu7 + (p_lu11 - p_lu7) * t for t in t3]

stroke_pts = np.array(s1 + s2 + s3)
print(f"Stworzono trajektorie testowa: {len(stroke_pts)} punktow 3D")

klisza = mgr.create_klisza([stroke_pts], name="Testowa Klisza Taiyin Reki (LU)")
saved_path = mgr.save_klisza(klisza)
print(f"Zapisano klisze do: {saved_path}")

mgr.print_klisza_report(klisza)

# Sprawdzamy wczytanie z indeksu
all_meta = mgr.list_all_metadata()
print(f"Liczba klisz w indeksie: {len(all_meta)}")
assert len(all_meta) >= 1
assert all_meta[0]["id"] == klisza.id

loaded = mgr.load_klisza(klisza.id)
assert loaded is not None
assert len(loaded.strokes) == 1
print(f"Wczytano poprawnie klisze o {len(loaded.all_matched_points)} dopasowanych punktach!")
