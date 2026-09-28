import sys
from acu3d.model import Atlas
atlas = Atlas()
print("Atlas mode:", atlas.mode, flush=True)
print("Cache path:", atlas._cache_path(), flush=True)
print("Cache exists:", atlas._cache_path().exists(), flush=True)
