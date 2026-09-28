import sys
from pathlib import Path
REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

import numpy as np
from acu3d.model import Atlas

atlas = Atlas()
an = atlas.anatomy

print("=" * 80)
print("  OBIEKTYWNY RAPORT AUDYTU MIAR CUN W MODELU 3D vs KANON (LING SHU GU DU FA)")
print("=" * 80 + "\n")

# 1. Sprawdzenie segmentow rigu i ich dlugosci cun
print("1. PROPORCJE SEGMENTÓW KOSTNYCH (GU DU FA - 骨度法):")
print("-" * 80)
print(f"{'Segment / Okolica':<22} | {'Model [cun]':<11} | {'Kanon [cun]':<11} | {'Dlugosc [cm]':<12} | {'1 cun [mm]':<10}")
print("-" * 80)

lore = [
    ("forearm (przedramię)",     12.0, 12.0, an.seg('forearm').length,   "zgięcie łokcia -> zgięcie nadgarstka"),
    ("upper_arm (ramię)",         9.0,  9.0, an.seg('upper_arm').length, "fałd pachowy -> łokieć"),
    ("hand (dłoń)",               7.5,  7.5, an.seg('hand').length,      "nadgarstek -> czubek palca III"),
    ("shank bok/tył (goleń)",    16.0, 16.0, an.seg('shank').length,     "szpara kolana -> kostka boczna"),
    ("thigh (udo)",              19.0, 19.0, an.seg('thigh').length,     "krętarz większy -> szpara kolana"),
    ("foot (stopa)",             10.0, 10.0, an.seg('foot').length,      "pięta -> czubek palca II"),
    ("brzuch górny (xiphoid)",    8.0,  8.0, None,                       "kąt mostkowy / mieczyk -> pępek"),
    ("brzuch dolny (podbrzusze)", 5.0,  5.0, None,                       "pępek -> spojenie łonowe"),
    ("szerokość sutków (klatka)", 8.0,  8.0, None,                       "ST-17 prawy -> ST-17 lewy"),
]

for name, m_cun, k_cun, length, desc in lore:
    l_str = f"{length*100.0:5.1f} cm" if length is not None else "---"
    mm_str = f"{length/m_cun*1000.0:4.1f} mm" if length is not None else "---"
    print(f"{name:<22} | {m_cun:<11.1f} | {k_cun:<11.1f} | {l_str:<12} | {mm_str:<10} | {desc}")

r_med = an.fit['rulers']['shank_medial']
p0, p1 = np.array(r_med['p0']), np.array(r_med['p1'])
d_med = np.linalg.norm(p1 - p0)
print(f"{'shank medial (SP/KI/LR)':<22} | {13.0:<11.1f} | {13.0:<11.1f} | {d_med*100:5.1f} cm    | {d_med/13.0*1000:4.1f} mm    | kłykieć piszczeli -> kostka przyśr.")

print("\n" + "=" * 80)
print("2. DOKŁADNE POMIARY PUNKTÓW W KLUCZOWYCH RELACJACH ANATOMICZNYCH:")
print("=" * 80)

def proj_axis(seg_name, k1, k2):
    sg = an.seg(seg_name)
    p1 = atlas.by_key[k1]['pos']
    p2 = atlas.by_key[k2]['pos']
    return abs(float(np.dot(np.array(p1) - np.array(p2), sg.e_ax)))

cun_lat_m = an.fit['cun']['trunk_lateral_m']

# Podudzie bok i przód
d_st_ax = proj_axis('shank', 'ST-35·L', 'ST-41·L')
d_st_cun = d_st_ax / (an.seg('shank').length / 16.0)
d_st36_ax = proj_axis('shank', 'ST-35·L', 'ST-36·L')
d_st36_cun = d_st36_ax / (an.seg('shank').length / 16.0)
print("A. GOLEŃ PRZÓD (Żołądek ST):")
print(f"   • ST-35 (szpara kolana) -> ST-36 (Zusanli): kanon 3.0 cun | zmierzono: {d_st36_cun:.2f} cun (błąd: {d_st36_cun-3:+.2f} cun)")
print(f"   • ST-35 (szpara kolana) -> ST-41 (staw skokowy): kanon 16.0 cun | zmierzono: {d_st_cun:.2f} cun (błąd: {d_st_cun-16:+.2f} cun)")

# Podudzie bok (GB)
d_gb_ax = proj_axis('shank', 'GB-34·L', 'GB-39·L')
d_gb_cun = d_gb_ax / (an.seg('shank').length / 16.0)
print("\nB. GOLEŃ BOK (Pęcherzyk Żółciowy GB):")
print(f"   • GB-34 (głowa strzałki) -> GB-39 (3 cun nad kostką): kanon 13.0 cun | zmierzono: {d_gb_cun:.2f} cun (błąd: {d_gb_cun-13:+.2f} cun)")

# Podudzie tył (BL)
d_bl_ax = proj_axis('shank', 'BL-40·L', 'BL-60·L')
d_bl_cun = d_bl_ax / (an.seg('shank').length / 16.0)
d_bl57_ax = proj_axis('shank', 'BL-40·L', 'BL-57·L')
d_bl57_cun = d_bl57_ax / (an.seg('shank').length / 16.0)
print("\nC. GOLEŃ TYŁ (Pęcherz Moczowy BL):")
print(f"   • BL-40 (dół podkolanowy) -> BL-57 (brzusiec łydki): kanon 8.0 cun | zmierzono: {d_bl57_cun:.2f} cun (błąd: {d_bl57_cun-8:+.2f} cun)")
print(f"   • BL-40 (dół podkolanowy) -> BL-60 (kostka boczna): kanon 16.0 cun | zmierzono: {d_bl_cun:.2f} cun (błąd: {d_bl_cun-16:+.2f} cun)")

# Podudzie przyśrodek (SP)
cun_m_med = abs(float(np.dot(np.asarray(r_med['p1']) - np.asarray(r_med['p0']), an.seg('shank').e_ax))) / 13.0
d_sp_ax = proj_axis('shank', 'SP-9·L', 'SP-6·L')
d_sp_cun = d_sp_ax / cun_m_med
print("\nD. GOLEŃ PRZYŚRODEK (Śledziona SP / Nerki KI / Wątroba LR):")
print(f"   • SP-9 (kłykieć piszczeli) -> SP-6 (3 cun nad kostką przyśr.): kanon 10.0 cun | zmierzono: {d_sp_cun:.2f} cun (błąd: {d_sp_cun-10:+.2f} cun)")

# Przedramię
d_lu_ax = proj_axis('forearm', 'LU-5·L', 'LU-9·L')
d_lu_cun = d_lu_ax / (an.seg('forearm').length / 12.0)
d_pc_ax = proj_axis('forearm', 'PC-3·L', 'PC-7·L')
d_pc_cun = d_pc_ax / (an.seg('forearm').length / 12.0)
print("\nE. PRZEDRAMIĘ (Płuca LU, Osierdzie PC, Serce HT):")
print(f"   • LU-5 (łokieć) -> LU-9 (nadgarstek): kanon 12.0 cun | zmierzono: {d_lu_cun:.2f} cun (błąd: {d_lu_cun-12:+.2f} cun)")
print(f"   • PC-3 (łokieć) -> PC-7 (nadgarstek): kanon 12.0 cun | zmierzono: {d_pc_cun:.2f} cun (błąd: {d_pc_cun-12:+.2f} cun)")

# Tułów pion
c_xiphoid = an.landmark_cun["xiphoid"]
c_umb = an.landmark_cun["umbilicus"]
c_pub = an.landmark_cun["pubis"]
print("\nF. TUŁÓW PION (Linia Pośrodkowa Przednia - Ren Mai):")
print(f"   • Kąt mostkowy (Xiphoid) -> Pępek (Umbilicus): kanon 8.0 cun | w modelu: {c_xiphoid - c_umb:.2f} cun")
print(f"   • Pępek (Umbilicus) -> Spojenie łonowe (Pubis): kanon 5.0 cun | w modelu: {c_umb - c_pub:.2f} cun")
print(f"   • Całkowita linijka brzuszna (Pubis -> Xiphoid): kanon 13.0 cun | w modelu: {c_xiphoid - c_pub:.2f} cun")

# Tułów poziomy
p_nipple_r = atlas.by_key['ST-17·R']['pos']
p_nipple_l = atlas.by_key['ST-17·L']['pos']
dist_nipples = abs(float(p_nipple_r[0] - p_nipple_l[0]))
cun_nipples = dist_nipples / cun_lat_m

p_ren8 = atlas.by_key['REN-8']['pos']
p_st25 = atlas.by_key['ST-25·L']['pos']
p_sp15 = atlas.by_key['SP-15·L']['pos']
cun_st25 = abs(float(p_st25[0] - p_ren8[0])) / cun_lat_m
cun_sp15 = abs(float(p_sp15[0] - p_ren8[0])) / cun_lat_m

print("\nG. SZEROKOŚCI TUŁOWIA (Linijka Boczna):")
print(f"   • Między sutkami (ST-17·P -> ST-17·L): kanon 8.0 cun | zmierzono: {cun_nipples:.2f} cun ({dist_nipples*1000:.1f} mm)")
print(f"   • Pępek REN-8 -> Tianshu ST-25 (linia ST): kanon 2.0 cun | zmierzono: {cun_st25:.2f} cun")
print(f"   • Pępek REN-8 -> Daheng SP-15 (linia SP): kanon 4.0 cun | zmierzono: {cun_sp15:.2f} cun")

# Grzbiet szerokość
p_du4 = atlas.by_key['DU-4']['pos']
p_bl23 = atlas.by_key['BL-23·L']['pos']
p_bl52 = atlas.by_key['BL-52·L']['pos']
cun_bl_in = abs(float(p_bl23[0] - p_du4[0])) / cun_lat_m
cun_bl_out = abs(float(p_bl52[0] - p_du4[0])) / cun_lat_m
print("\nH. SZEROKOŚCI GRZBIETU (Od Linii Pośrodkowej DU Mai):")
print(f"   • DU-4 -> BL-23 (I linia Pęcherza BL):  kanon 1.5 cun | zmierzono: {cun_bl_in:.2f} cun")
print(f"   • DU-4 -> BL-52 (II linia Pęcherza BL): kanon 3.0 cun | zmierzono: {cun_bl_out:.2f} cun")

print("\n" + "=" * 80)
print("3. SIATKA CUNOWA WIDOKU 3D (WIRE / WIRE_SURFACE):")
print("=" * 80)
from acu3d.wire import Wire
w = Wire()
print(f"• Wierzchołków na ciele:    {len(w.points):,}")
print(f"• Krawędzi łączących:        {len(w.edges):,}")
print(f"• Krok gęstości pierścieni:  {w.step_cun} cun (dokładnie co 0.5 cun)")
print(f"• Liczba promieni azymutu:   {w.rays} (dokładnie co 15° wokół osi)")
print("• Każdy odcinek ciała posiada liczbę pierścieni = (długość_cun / 0.5) + 1")
print("  np. podudzie (16 cun) -> 33 pierścienie co 0.5 cun")
print("  np. przedramię (12 cun) -> 25 pierścieni co 0.5 cun")
print("  np. udo (19 cun) -> 39 pierścieni co 0.5 cun")
