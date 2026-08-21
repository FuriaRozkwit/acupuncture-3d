"""Spis tresci modelu Z-Anatomy - uruchamiany WEWNATRZ Blendera.

    blender --background --factory-startup --python tools/blender_inventory.py -- <blend> <out.json>

Wypisuje drzewo kolekcji i liste obiektow siatkowych (nazwa, sciezka kolekcji,
liczba trojkatow, bounding box), zeby dalo sie wybrac, co w ogole eksportowac.
"""
import json
import sys

import bpy


def arg(i, default=None):
    a = sys.argv[sys.argv.index("--") + 1:] if "--" in sys.argv else []
    return a[i] if len(a) > i else default


def collection_path(obj):
    """Sciezka kolekcji, w ktorej siedzi obiekt (od korzenia sceny)."""
    parent = {}
    for c in bpy.data.collections:
        for ch in c.children:
            parent[ch.name] = c.name
    if not obj.users_collection:
        return ""
    name = obj.users_collection[0].name
    path = [name]
    seen = set()
    while name in parent and name not in seen:
        seen.add(name)
        name = parent[name]
        path.append(name)
    return "/".join(reversed(path))


def main():
    blend = arg(0)
    out = arg(1, "inventory.json")
    bpy.ops.wm.open_mainfile(filepath=blend)

    objects = []
    for obj in bpy.data.objects:
        if obj.type != "MESH" or obj.data is None:
            continue
        me = obj.data
        bb = [list(v) for v in obj.bound_box]
        xs = [p[0] for p in bb]
        ys = [p[1] for p in bb]
        zs = [p[2] for p in bb]
        objects.append({
            "name": obj.name,
            "collection": collection_path(obj),
            "verts": len(me.vertices),
            "polys": len(me.polygons),
            "bounds": [round(min(xs), 4), round(max(xs), 4),
                       round(min(ys), 4), round(max(ys), 4),
                       round(min(zs), 4), round(max(zs), 4)],
            "hide": bool(obj.hide_get()) if obj.name in bpy.context.view_layer.objects else None,
        })

    tree = {}
    for c in bpy.data.collections:
        tree[c.name] = {"children": [ch.name for ch in c.children],
                        "objects": len(c.objects)}

    doc = {"blend": blend, "n_objects": len(objects), "collections": tree, "objects": objects}
    with open(out, "w", encoding="utf-8") as f:
        json.dump(doc, f, ensure_ascii=False, indent=1)

    total = sum(o["polys"] for o in objects)
    print(f"[inventory] obiektow siatkowych: {len(objects)}, wielokatow razem: {total}")
    print(f"[inventory] kolekcji: {len(tree)} -> {out}")


main()
