#!/usr/bin/env python3
"""spl を SKF に結合して <skf>rep を作り、E(a) を LAK と比較 (mac 版)。
使い方: python attach_validate.py <skf_dir> <spl_dir> <out_dir> [--pairs Mo-Mo,Mo-S,S-Mo,S-S] [--dft-db pairset3_dft.db]"""
import argparse
import os
import shutil
import subprocess
import sys

import numpy as np
from ase.db import connect
from ase.io import write as ase_write

SCRIPTS = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, SCRIPTS)
from dftb_energy_local import HSD, parse_results_tag, DFTB  # noqa: E402
HA = 27.211386245988

ap = argparse.ArgumentParser()
ap.add_argument("skf"); ap.add_argument("spl"); ap.add_argument("out")
ap.add_argument("--pairs", default="Mo-Mo,Mo-S,S-Mo,S-S")
ap.add_argument("--dft-db", default=f"{os.path.dirname(SCRIPTS)}/local_opt/repfit3/pairset3_dft.db")
ap.add_argument("--no-validate", action="store_true")
args = ap.parse_args()
os.makedirs(args.out, exist_ok=True)
for f in os.listdir(args.skf):
    if f.endswith(".skf"):
        shutil.copy(os.path.join(args.skf, f), os.path.join(args.out, f))
for pair in args.pairs.split(","):
    spl = os.path.join(args.spl, pair + ".spl")
    if not os.path.exists(spl):
        a, b = pair.split("-"); spl = os.path.join(args.spl, f"{b}-{a}.spl")
    base = open(os.path.join(args.skf, pair + ".skf")).read()
    open(os.path.join(args.out, pair + ".skf"), "w").write(base + open(spl).read())
    print("attached", pair, "<-", os.path.basename(spl))
if args.no_validate:
    sys.exit(0)
# --- E(a) validation ---
db = connect(args.dft_db)
rows = [(r.get("name"), r.toatoms(), r.energy) for r in db.select() if r.get("name", "").startswith("a_")]
wd = f"/private/tmp/claude-501/-Users-crocus-uhuhu-MoS2-DFTB/6416faee-14e1-48c8-9838-35c06f3e7813/scratchpad/attach_val"
res = []
for name, atoms, e_dft in rows:
    d = os.path.join(wd, name); os.makedirs(d, exist_ok=True)
    ase_write(os.path.join(d, "geo.gen"), atoms, format="gen")
    hsd = HSD.format(skf=os.path.abspath(args.out), kblock="  KPointsAndWeights = SupercellFolding {\n    12 0 0\n    0 12 0\n    0 0 1\n    0.5 0.5 0.0\n  }")
    hsd = hsd.replace("  PolynomialRepulsive = SetForAll { Yes }\n", "")
    open(os.path.join(d, "dftb_in.hsd"), "w").write(hsd)
    r = subprocess.run([DFTB], cwd=d, capture_output=True, text=True, env=dict(os.environ, OMP_NUM_THREADS="4"))
    e, _ = parse_results_tag(os.path.join(d, "results.tag"))
    res.append((float(name[2:]), e_dft, e * HA))
res.sort()
a = np.array([x[0] for x in res]); ed = np.array([x[1] for x in res]); et = np.array([x[2] for x in res])
ed -= ed.min(); et -= et.min()
msk = (a >= 3.04) & (a <= 3.32)
rms = np.sqrt(((ed - et)[msk] ** 2).mean())
def amin(a, e):
    m = (a >= 3.08) & (a <= 3.26); c = np.polyfit(a[m], e[m], 2); return -c[1] / (2 * c[0])
print(f"E(a) RMS (3.04-3.32, min-aligned) = {rms*1000:.1f} meV/cell; a_eq LAK {amin(a, ed):.4f} vs DFTB {amin(a, et):.4f} A "
      f"({(amin(a, et)/amin(a, ed)-1)*100:+.2f}%)")
for x, y, z in zip(a, ed, et):
    print(f"  a={x:.2f}  LAK {y:+.4f}  DFTB {z:+.4f}  diff {z-y:+.4f}")
