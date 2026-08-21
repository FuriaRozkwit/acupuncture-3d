"""Budowa siatek 3D warstw ciala z rigu proporcjonalnego.

Model jest stylizowany (proceduralny), nie skanem anatomicznym - jego zadaniem jest
poprawna *topografia* powierzchni, na ktorej siadaja punkty. Kazda warstwa powstaje
z tych samych segmentow, tylko z inna skala promieni:

    skora   1.00   - powierzchnia, na ktorej leza punkty
    miesnie 0.88   - warstwa podskorna (uproszczona)
    kosci   ~0.3   - osie kosci, stawy, kregoslup, zebra, miednica, czaszka
"""
from __future__ import annotations

import numpy as np
import pyvista as pv

from .rig import Rig, Segment

LIMB_SEGMENTS = ("upper_arm", "forearm", "hand", "thigh", "shank", "foot")
MIRRORED = LIMB_SEGMENTS
CENTRAL = ("trunk", "neck", "head")


def loft(rings: list[np.ndarray], cap_start: bool = True, cap_end: bool = True) -> pv.PolyData:
    """Zszywa kolejne przekroje w powloke czworokatna."""
    n = rings[0].shape[0]
    pts = np.vstack(rings)
    faces = []
    for i in range(len(rings) - 1):
        a, b = i * n, (i + 1) * n
        for j in range(n):
            k = (j + 1) % n
            faces.append([4, a + j, a + k, b + k, b + j])
    mesh = pv.PolyData(pts, np.hstack(faces))
    caps = []
    if cap_start:
        caps.append(_cap(rings[0], flip=True))
    if cap_end:
        caps.append(_cap(rings[-1], flip=False))
    for c in caps:
        mesh = mesh.merge(c)
    return mesh.clean()


def _cap(ring: np.ndarray, flip: bool) -> pv.PolyData:
    c = ring.mean(axis=0)
    pts = np.vstack([ring, c])
    n = ring.shape[0]
    faces = []
    for j in range(n):
        k = (j + 1) % n
        tri = [3, j, k, n] if not flip else [3, k, j, n]
        faces.append(tri)
    return pv.PolyData(pts, np.hstack(faces))


def segment_mesh(seg: Segment, scale: float = 1.0, steps: int = 26, res: int = 40) -> pv.PolyData:
    ts = np.linspace(0.0, 1.0, steps)
    rings = [seg.ring(float(t), n=res, scale=scale) for t in ts]
    return loft(rings)


def _joint(seg: Segment, t: float, scale: float, k: float = 1.0) -> pv.PolyData:
    rx, ry, dy = seg.radii(t)
    c = seg.axis_point(t) + dy * seg.e_ant
    r = 0.5 * (rx + ry) * scale * k
    return pv.Sphere(radius=r, center=c, theta_resolution=24, phi_resolution=24)


def _mirror(mesh: pv.PolyData) -> pv.PolyData:
    m = mesh.copy()
    m.points[:, 0] *= -1.0
    return m.flip_faces()


def soft_layer(rig: Rig, scale: float) -> pv.PolyData:
    """Warstwa 'objetosciowa': skora (scale=1) lub miesnie (scale<1)."""
    parts: list[pv.PolyData] = []
    for name in CENTRAL:
        parts.append(segment_mesh(rig.seg(name), scale))
    right: list[pv.PolyData] = []
    for name in MIRRORED:
        seg = rig.seg(name)
        right.append(segment_mesh(seg, scale))
        prox_k = {"thigh": 0.92, "upper_arm": 0.92, "hand": 0.8, "foot": 0.8}.get(name, 1.0)
        right.append(_joint(seg, 0.0, scale, prox_k))   # kula stawowa nie moze wystawac
        #                                                 poza sasiedni segment - inaczej
        #                                                 chowa punkty na faldzie nadgarstka
        dist_k = 0.85 if name in ("shank", "forearm") else 1.0
        right.append(_joint(seg, 1.0, scale, dist_k))
    parts += right + [_mirror(p) for p in right]
    merged = parts[0].merge(parts[1:]) if len(parts) > 1 else parts[0]
    return merged.clean().compute_normals(split_vertices=False)


def skeleton(rig: Rig) -> pv.PolyData:
    """Uproszczony uklad kostny: kosci dlugie, stawy, kregoslup, zebra, miednica, czaszka."""
    parts: list[pv.PolyData] = []
    right: list[pv.PolyData] = []
    for name in LIMB_SEGMENTS:
        seg = rig.seg(name)
        s = 0.34 if name in ("upper_arm", "forearm", "thigh", "shank") else 0.5
        right.append(segment_mesh(seg, s, steps=8, res=18))
        right.append(_joint(seg, 0.0, s, 1.5))
        right.append(_joint(seg, 1.0, s, 1.4))
    parts += right + [_mirror(p) for p in right]

    # kregoslup: kregi jako splaszczone kule na tylnej scianie tulowia
    trunk = rig.seg("trunk")
    for name, z in rig.vertebra_z.items():
        t = rig.trunk_t_for_z(z)
        rx, ry, dy = trunk.radii(t)
        c = trunk.axis_point(t) + (dy - ry * 0.55) * trunk.e_ant
        parts.append(pv.Sphere(radius=0.012 * rig.H, center=c, theta_resolution=16, phi_resolution=16))
    # szyjne C1-C6 w gore po szyi
    neck = rig.seg("neck")
    for i in range(6):
        t = 0.1 + i * 0.15
        rx, ry, dy = neck.radii(t)
        c = neck.axis_point(t) + (dy - ry * 0.5) * neck.e_ant
        parts.append(pv.Sphere(radius=0.010 * rig.H, center=c, theta_resolution=14, phi_resolution=14))

    # zebra: 10 par obreczy w klatce piersiowej
    for i, zf in enumerate(np.linspace(0.805, 0.665, 10)):
        t = rig.trunk_t_for_z(zf * rig.H)
        ring = trunk.ring(t, n=64, scale=0.93 - 0.01 * i)
        keep = ring[ring[:, 1] > trunk.axis_point(t)[1] - trunk.radii(t)[1] * 0.55]
        if len(keep) > 3:
            parts.append(pv.Spline(np.vstack([keep, keep[:1]]), 200).tube(radius=0.006 * rig.H))
    # mostek
    zs = np.linspace(0.80, 0.685, 8) * rig.H
    stern = []
    for z in zs:
        t = rig.trunk_t_for_z(z)
        stern.append(trunk.surface_by_lateral(t, 0.0, "anterior") - trunk.e_ant * 0.012 * rig.H)
    parts.append(pv.Spline(np.array(stern), 60).tube(radius=0.014 * rig.H))
    # miednica
    for zf in (0.520, 0.560):
        t = rig.trunk_t_for_z(zf * rig.H)
        ring = trunk.ring(t, n=64, scale=0.92)
        parts.append(pv.Spline(np.vstack([ring, ring[:1]]), 240).tube(radius=0.010 * rig.H))
    # obojczyki i lopatki (zarys)
    for sgn in (1, -1):
        a = np.array([sgn * 0.012 * rig.H, rig.D["chest"] * 0.75, rig.levels["sternal_notch"]])
        b = np.array([sgn * rig.W["shoulder_joint_x"], rig.D["chest"] * 0.15,
                      rig.levels["acromion"] - 0.005 * rig.H])
        parts.append(pv.Spline(np.array([a, (a + b) / 2 + np.array([0, 0.01 * rig.H, 0]), b]), 40)
                     .tube(radius=0.008 * rig.H))
    # czaszka
    head = rig.seg("head")
    parts.append(segment_mesh(head, 0.86, steps=20, res=32))

    merged = parts[0].merge(parts[1:])
    return merged.clean()


def build_layers(rig: Rig) -> dict:
    return {
        "skin": soft_layer(rig, 1.0),
        "muscle": soft_layer(rig, 0.88),
        "skeleton": skeleton(rig),
    }
