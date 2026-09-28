from acu3d.extraordinary import EXTRAORDINARY_VESSELS, get_vessel, get_point_extraordinary_memberships
from acu3d.jing_jin import JING_JIN_CHANNELS, MYERS_LINES, analyze_stroke_fascia
from acu3d.clinical_kb import get_shu_mu_info, get_yuan_luo_info, get_xi_cleft_info

print("8 Naczynia count:", len(EXTRAORDINARY_VESSELS))
for k, v in EXTRAORDINARY_VESSELS.items():
    name = v["name_pl"]
    seq_len = len(v["sequence"])
    print(f"  {k}: {name} ({seq_len} pkt)")

print("\nMyers lines count:", len(MYERS_LINES))
for k, v in MYERS_LINES.items():
    name = v["name_pl"]
    jj_name = v["jing_jin_name"]
    print(f"  {k}: {name} -> {jj_name}")

print("\nShu-Mu test LU:", get_shu_mu_info("LU"))
print("Yuan-Luo test LU:", get_yuan_luo_info("LU-9"))
print("Xi-cleft test LU:", get_xi_cleft_info("LU-6"))
print("Point extraordinary test GB-41:", get_point_extraordinary_memberships("GB-41"))
print("Stroke fascia test BL:", analyze_stroke_fascia([], "BL"))
