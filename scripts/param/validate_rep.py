#!/usr/bin/env python3
"""反発込み SKF (skf_v2rep) の検証: E(a) カーブを LAK と比較。"""
import os
import subprocess
import sys

import numpy as np
from ase.db import connect

sys.path.insert(0, "/workspace/MoS2_DFTB")
from dftb_energy import HSD, parse_results_tag  # noqa: E402

SKF = "/workspace/MoS2_DFTB/dftb/skf_v2rep"
DFTB = "/opt/dftbplus/bin/dftb+"
HA = 27.211386245988


def dftb_energy_rep(atoms, workdir):
    from ase.io import write as ase_write
    os.makedirs(workdir, exist_ok=True)
    ase_write(os.path.join(workdir, "geo.gen"), atoms, format="gen")
    kblock = ("  KPointsAndWeights = SupercellFolding {\n"
              "    12 0 0\n    0 12 0\n    0 0 1\n    0.5 0.5 0.0\n  }")
    hsd = HSD.format(skf=SKF, kblock=kblock)
    hsd = hsd.replace("  PolynomialRepulsive = SetForAll { Yes }\n", "")
    with open(os.path.join(workdir, "dftb_in.hsd"), "w") as f:
        f.write(hsd)
    env = dict(os.environ, OMP_NUM_THREADS="4", OPENBLAS_NUM_THREADS="1")
    r = subprocess.run([DFTB], cwd=workdir, capture_output=True, text=True,
                       timeout=900, env=env)
    tag = os.path.join(workdir, "results.tag")
    if r.returncode != 0 or not os.path.exists(tag):
        return None
    e, _ = parse_results_tag(tag)
    return e * HA


db = connect("/workspace/MoS2_DFTB/repfit/fit_dft.db")
rows = [(row.get("name"), row.toatoms(), row.energy) for row in db.select()
        if row.get("name", "").startswith("a_")]
rows.sort(key=lambda t: t[0])

print("  a      E_LAK-min   E_DFTB-min  (eV/cell)")
avals, e_lak, e_dftb = [], [], []
for name, atoms, edft in rows:
    a = float(name.split("_")[1])
    ed = dftb_energy_rep(atoms, f"/root/valrep/{name}")
    if ed is None:
        print(f"{name}: DFTB failed")
        continue
    avals.append(a)
    e_lak.append(edft)
    e_dftb.append(ed)

avals = np.array(avals)
e_lak = np.array(e_lak) - min(e_lak)
e_dftb = np.array(e_dftb) - min(e_dftb)
for a, el, ed in zip(avals, e_lak, e_dftb):
    print(f"{a:.2f}   {el:9.4f}   {ed:9.4f}")

# 平衡格子: 平衡近傍 (3.08-3.32) で 4 次フィット
msk = (avals >= 3.04) & (avals <= 3.32)
for label, e in [("LAK", e_lak), ("DFTB", e_dftb)]:
    c = np.polyfit(avals[msk], e[msk], 4)
    xs = np.linspace(3.05, 3.3, 1001)
    a0 = xs[np.polyval(c, xs).argmin()]
    print(f"{label}: a0 = {a0:.4f} A")
rms = np.sqrt(((e_lak - e_dftb) ** 2).mean())
print(f"E(a) curve RMS deviation: {rms*1000:.1f} meV/cell over {len(avals)} points")
