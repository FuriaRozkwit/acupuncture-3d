"""Eksport pojedynczych obiektow Z-Anatomy po nazwie - uruchamiany WEWNATRZ Blendera.

    blender --background --factory-startup --python tools/blender_export_named.py -- <blend> <outdir> <Nazwa1,Nazwa2,...>

Sluzy do wyciagniecia kosci referencyjnych, po ktorych dopasowujemy do siebie
uklady wspolrzednych Z-Anatomy i BodyParts3D.
"""
import os
import sys

import bpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from blender_export import mesh_of, slug, write_ply  # noqa: E402


def main():
    a = sys.argv[sys.argv.index("--") + 1:]
    blend, outdir, names = a[0], a[1], [n.strip() for n in a[2].split(",")]
    bpy.ops.wm.open_mainfile(filepath=blend)
    dg = bpy.context.evaluated_depsgraph_get()
    os.makedirs(outdir, exist_ok=True)
    for name in names:
        obj = bpy.data.objects.get(name)
        if obj is None or obj.type != "MESH":
            print(f"[named] brak obiektu: {name}")
            continue
        v, t = mesh_of(obj, dg)
        if not v:
            print(f"[named] pusta siatka: {name}")
            continue
        path = os.path.join(outdir, f"{slug(name)}.ply")
        write_ply(path, v, t)
        print(f"[named] {name}: {len(t)} trojkatow -> {path}")


main()
