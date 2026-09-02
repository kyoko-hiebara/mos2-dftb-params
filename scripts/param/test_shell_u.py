#!/usr/bin/env python3
"""skf_v3 に殻別 Hubbard U を入れ ShellResolvedSCC=Yes で バンド損失と V_S 準位を評価。"""
import json, os, shutil, subprocess, sys
SCRIPTS = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, SCRIPTS)
import optimize_multi as om
from compare_bands import band_loss
import extended_loss as ext
ROOT = om.ROOT
U = json.load(open(f"{ROOT}/local_opt/shell_hubbard.json"))
OUT = f"{ROOT}/local_opt/skf_v3_shU"
os.makedirs(OUT, exist_ok=True)
for f in os.listdir(f"{ROOT}/local_opt/skf_v3"):
    if f.endswith(".skf"): shutil.copy(f"{ROOT}/local_opt/skf_v3/{f}", f"{OUT}/{f}")
for el, shells in [("Mo", ["4d", "5p", "5s"]), ("S", ["3d", "3p", "3s"])]:
    p = f"{OUT}/{el}-{el}.skf"; L = open(p).readlines(); t = L[1].split()
    old = t[4:7]
    for i, nl in enumerate(shells):
        t[4 + i] = f"{U[el]['U_fd'][nl]:.6f}"
    L[1] = " ".join(t) + "\n"; open(p, "w").writelines(L)
    print(el, "U d/p/s:", old, "->", t[4:7])
om.DEFECT_HSD = om.DEFECT_HSD.replace("SCC = Yes\n", "SCC = Yes\n  ShellResolvedSCC = Yes\n")
wd = f"{ROOT}/local_opt/shU_eval"; os.makedirs(wd, exist_ok=True)
r = subprocess.run([sys.executable, f"{SCRIPTS}/dftb_bands.py", OUT, "3.16", f"{wd}/dftb", "--json", f"{wd}/dftb.json",
                    "--thickness", "3.127", "--s-lmax", "d", "--shell-resolved", "--extra-kpts", f"{SCRIPTS}/mesh_kpts_12x12.json"],
                   capture_output=True, text=True)
print(r.stdout.strip()[-120:])
res = band_loss(om.VASP_JSON, f"{wd}/dftb.json", verbose=False)
e = ext.all_terms(om.VASP_JSON, f"{wd}/dftb.json", target_midgap=-5.068)
depth = om.defect_depth(OUT, f"{wd}/defect")
print(f"shell-resolved U: band loss={res['loss']:.4f} (v3 0.2514) rmsVB={res['rms_vb']:.3f} rmsCB={res['rms_cb']:.3f} gapK={res['gap_K_dftb']:.3f} "
      f"mVB={e['m_VB_KM_dftb']:.2f} mCB={e['m_CB_KM_dftb']:.2f} dQ={e['dq_dftb']:.3f} mesh={e['rms_mesh']:.3f} midgap={e['midgap_abs']:.3f}")
print(f"V_S depth = {depth} (v3 0.690, target 0.554)")
