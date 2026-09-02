#!/usr/bin/env python3
"""Sb 反発込みセットの検証: バルク Sb E(scale) と Sb2S3 E(scale) を PBE 参照と比較。
使い方: python validate_sb_rep.py <skf_dir_with_splines>"""
import os, sys
import numpy as np
from ase.db import connect
SCRIPTS = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, SCRIPTS)
from dftb_energy_sb import run_one
skf = os.path.abspath(sys.argv[1])
db = connect(f"{os.path.dirname(SCRIPTS)}/local_opt/repfit_sb/sb_ref_dft.db")
wd = "/private/tmp/claude-501/-Users-crocus-uhuhu-MoS2-DFTB/6416faee-14e1-48c8-9838-35c06f3e7813/scratchpad/val_sb_rep"
for prefix, label in [("sb_s", "bulk Sb (u=0.2336)"), ("sb2s3_s", "Sb2S3")]:
    rows = [(r.name, r.toatoms(), r.energy) for r in db.select() if r.name.startswith(prefix) and ("u0.2336" in r.name or prefix == "sb2s3_s")]
    if not rows:
        continue
    sc, ed, et = [], [], []
    for name, atoms, e_dft in sorted(rows):
        e, f, err = run_one(atoms, skf, os.path.join(wd, name), True, 0.2)
        if e is None:
            print("FAILED", name, err); continue
        sc.append(float(name.split("_s")[1].split("_")[0])); ed.append(e_dft); et.append(e)
    sc, ed, et = np.array(sc), np.array(ed), np.array(et)
    ed -= ed.min(); et -= et.min()
    def smin(x, y):
        c = np.polyfit(x, y, 2); return -c[1] / (2 * c[0])
    print(f"{label}: scale_eq PBE {smin(sc, ed):.4f} vs DFTB {smin(sc, et):.4f}; RMS(E, min-aligned) = {np.sqrt(((ed-et)**2).mean())*1000:.1f} meV/cell")
    for x, y, z in zip(sc, ed, et):
        print(f"  s={x:.2f}  PBE {y:+.4f}  DFTB {z:+.4f}  diff {z-y:+.4f}")
