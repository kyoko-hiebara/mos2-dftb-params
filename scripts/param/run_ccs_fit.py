#!/usr/bin/env python3
"""ccs_fit で Mo-S / S-S / Mo-Mo の反発ポテンシャルをフィット。

venv_ccs で実行。前提: repfit/pairset_dft.db と repfit/pairset_dftb.db。
出力: structures.json, CCS_params.json (+ error analysis)
"""
import json
import os
import sys

PREFIX = sys.argv[1] if len(sys.argv) > 1 else "pairset"
os.chdir("/workspace/MoS2_DFTB/repfit")

# 分子はスピン状態が DFT(triplet)/DFTB(非偏極) で不整合なため除外
from ase.calculators.singlepoint import SinglePointCalculator  # noqa: E402
from ase.db import connect  # noqa: E402

for suffix in ["dft", "dftb"]:
    srcdb = connect(f"{PREFIX}_{suffix}.db")
    dstp = f"fit_{suffix}.db"
    if os.path.exists(dstp):
        os.remove(dstp)
    dst = connect(dstp)
    n = 0
    for row in srcdb.select():
        name = row.get("name", "")
        if name.startswith(("S2_", "SO_", "O2_")):
            continue
        atoms = row.toatoms()
        atoms.calc = SinglePointCalculator(
            atoms, energy=row.energy,
            forces=(row.forces if "forces" in row else None))
        dst.write(atoms, name=name)
        n += 1
    print(f"{dstp}: {n} structures (molecules excluded)")

from ccs_fit.scripts.ccs_fetch import ccs_fetch  # noqa: E402

print("fetching pair data ...", flush=True)
ccs_fetch(mode="DFTB", R_c=5.0, Ns="all",
          DFT_DB="fit_dft.db", DFTB_DB="fit_dftb.db")
print("structures.json written", flush=True)

inp = {
    "General": {"interface": "DFTB"},
    "Twobody": {
        "Mo-S": {"Rcut": 3.6, "Resolution": 0.05, "Swtype": "rep"},
        "S-S": {"Rcut": 4.2, "Resolution": 0.05, "Swtype": "rep"},
        "Mo-Mo": {"Rcut": 4.6, "Resolution": 0.05, "Swtype": "rep"},
    },
}
json.dump(inp, open("CCS_input.json", "w"), indent=2)

from ccs_fit import ccs_fit  # noqa: E402

print("fitting ...", flush=True)
ccs_fit("CCS_input.json")
print("CCS fit done", flush=True)
for f in ["CCS_params.json", "ccs_error.out", "error.out"]:
    if os.path.exists(f):
        print(f"-> {f} ({os.path.getsize(f)} bytes)")
# 簡単なエラーサマリ
try:
    import numpy as np
    err = np.loadtxt("error.out")
    de = err[:, 0] - err[:, 1]
    print(f"fit residual: RMS={np.sqrt((de**2).mean())*1000:.2f} meV, "
          f"max={abs(de).max()*1000:.2f} meV over {len(de)} structures")
except Exception as e:
    print("error summary skipped:", e)
print("CCS_FIT_DONE")
