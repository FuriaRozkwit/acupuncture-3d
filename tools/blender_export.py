"""Eksport siatek Z-Anatomy do plikow PLY - uruchamiany WEWNATRZ Blendera.

    blender --background --factory-startup --python tools/blender_export.py -- <blend> <outdir> [grupy]

Grupy: muscles | skeleton | insertions | lines  (domyslnie: muscles,skeleton,lines)

Nie uzywa operatorow eksportu ani zaznaczania - czyta dane siatek wprost i sam
zapisuje binarne PLY. Dzieki temu nie przeszkadzaja ukryte obiekty ani wylaczone
kolekcje, a wynik jest powtarzalny. Decymacja odbywa sie pozniej po stronie
PyVisty (tools/prepare_meshes.py), bo VTK robi to lepiej i bez modyfikatorow.

Wspolrzedne zapisujemy w ukladzie swiata (macierz obiektu zaaplikowana).
"""
import json
import os
import re
import struct
import sys

import bpy

GROUPS = {
    "muscles": "4:",
    "skeleton": "1:",
    "insertions": "2:",
    "joints": "3:",
    "lines": "Bonus",
}


def args():
    return sys.argv[sys.argv.index("--") + 1:] if "--" in sys.argv else []


def slug(name):
    s = name.strip().lower()
    s = s.replace("(", "").replace(")", "")
    s = re.sub(r"[^a-z0-9]+", "-", s).strip("-")
    return s or "part"


def write_ply(path, verts, tris):
    """Binarne PLY: lista wierzcholkow float32 + trojkaty uint32."""
    with open(path, "wb") as f:
        f.write(b"ply\nformat binary_little_endian 1.0\n")
        f.write(b"comment Z-Anatomy (CC BY-SA 4.0) via acu3d exporter\n")
        f.write(f"element vertex {len(verts)}\n".encode())
        f.write(b"property float x\nproperty float y\nproperty float z\n")
        f.write(f"element face {len(tris)}\n".encode())
        f.write(b"property list uchar uint vertex_indices\nend_header\n")
        f.write(struct.pack(f"<{len(verts) * 3}f", *[c for v in verts for c in v]))
        for t in tris:
            f.write(struct.pack("<B3I", 3, t[0], t[1], t[2]))


def mesh_of(obj, depsgraph):
    """Trojkaty i wierzcholki obiektu w ukladzie swiata."""
    ev = obj.evaluated_get(depsgraph)
    try:
        me = ev.to_mesh()
    except RuntimeError:
        return None, None
    if me is None or len(me.vertices) == 0:
        ev.to_mesh_clear()
        return None, None
    mw = obj.matrix_world
    verts = [tuple(mw @ v.co) for v in me.vertices]
    me.calc_loop_triangles()
    tris = [tuple(t.vertices) for t in me.loop_triangles]
    ev.to_mesh_clear()
    return verts, tris


def main():
    a = args()
    blend, outdir = a[0], a[1]
    groups = (a[2] if len(a) > 2 else "muscles,skeleton,lines").split(",")
    bpy.ops.wm.open_mainfile(filepath=blend)
    depsgraph = bpy.context.evaluated_depsgraph_get()

    os.makedirs(outdir, exist_ok=True)
    catalog = {"source": "Z-Anatomy, CC BY-SA 4.0", "groups": {}, "parts": []}
    world = [1e9, -1e9, 1e9, -1e9, 1e9, -1e9]

    def track(verts):
        for v in verts:
            for i in range(3):
                world[i * 2] = min(world[i * 2], v[i])
                world[i * 2 + 1] = max(world[i * 2 + 1], v[i])

    for group in groups:
        prefix = GROUPS.get(group)
        if prefix is None:
            print(f"[export] nieznana grupa: {group}")
            continue
        objs = [o for o in bpy.data.objects
                if o.type == "MESH" and o.users_collection
                and o.users_collection[0].name.startswith(prefix)]
        gdir = os.path.join(outdir, group)
        os.makedirs(gdir, exist_ok=True)

        if group == "skeleton":                      # kosci lecza w jedna siatke
            allv, allt, names = [], [], []
            for o in objs:
                v, t = mesh_of(o, depsgraph)
                if not v:
                    continue
                off = len(allv)
                allv += v
                allt += [(x + off, y + off, z + off) for x, y, z in t]
                names.append(o.name)
            track(allv)
            path = os.path.join(gdir, "skeleton.ply")
            write_ply(path, allv, allt)
            catalog["groups"][group] = {"file": "skeleton/skeleton.ply",
                                        "parts": len(names), "polys": len(allt)}
            print(f"[export] {group}: {len(names)} obiektow -> jedna siatka, {len(allt)} trojkatow")
            continue

        n = 0
        for o in objs:
            v, t = mesh_of(o, depsgraph)
            if not v or not t:
                continue
            track(v)
            name = o.name
            side = "L" if name.endswith(".l") else "R" if name.endswith(".r") else "C"
            base = re.sub(r"\.(l|r|j|g)$", "", name)
            fn = f"{slug(name)}.ply"
            write_ply(os.path.join(gdir, fn), v, t)
            xs = [p[0] for p in v]
            ys = [p[1] for p in v]
            zs = [p[2] for p in v]
            catalog["parts"].append({
                "name": base, "raw": name, "group": group, "side": side,
                "file": f"{group}/{fn}", "polys": len(t), "verts": len(v),
                "bounds": [round(min(xs), 5), round(max(xs), 5), round(min(ys), 5),
                           round(max(ys), 5), round(min(zs), 5), round(max(zs), 5)],
                "centroid": [round(sum(xs) / len(xs), 5), round(sum(ys) / len(ys), 5),
                             round(sum(zs) / len(zs), 5)],
            })
            n += 1
        catalog["groups"][group] = {"parts": n}
        print(f"[export] {group}: {n} plikow")

    catalog["world_bounds"] = [round(w, 5) for w in world]
    with open(os.path.join(outdir, "catalog.json"), "w", encoding="utf-8") as f:
        json.dump(catalog, f, ensure_ascii=False, indent=1)
    print(f"[export] zakres swiata: {catalog['world_bounds']}")
    print(f"[export] katalog -> {os.path.join(outdir, 'catalog.json')}")


main()
