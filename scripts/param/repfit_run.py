#!/usr/bin/env python3
"""CCS 反発フィット (venv_ccs_local で実行): fetch(力込み可) -> fit -> spl 出力。
使い方: python repfit_run.py <dir> <prefix> --pairs "Mo-S:3.6,S-S:4.2,Mo-Mo:4.6" [--forces] [--rc 5.0]"""
import argparse
import json
import os
import sys

ap = argparse.ArgumentParser()
ap.add_argument("dir"); ap.add_argument("prefix")
ap.add_argument("--pairs", required=True)
ap.add_argument("--forces", action="store_true")
ap.add_argument("--rc", type=float, default=5.0)
ap.add_argument("--res", type=float, default=0.05)
ap.add_argument("--swtype", default="rep")
args = ap.parse_args()
os.chdir(args.dir)
from ccs_fit.scripts.ccs_fetch import ccs_fetch
print("fetch ...", flush=True)
ccs_fetch(mode="DFTB", R_c=args.rc, Ns="all", DFT_DB=f"{args.prefix}_dft.db",
          DFTB_DB=f"{args.prefix}_dftb.db", include_forces=args.forces)
two = {}
for item in args.pairs.split(","):
    pair, rc = item.split(":")
    two[pair] = {"Rcut": float(rc), "Resolution": args.res, "Swtype": args.swtype}
inp = {"General": {"interface": "DFTB"}, "Twobody": two}
json.dump(inp, open("CCS_input.json", "w"), indent=2)
from ccs_fit import ccs_fit
print("fit ...", flush=True)
try:
    ccs_fit("CCS_input.json")
except Exception as ex:
    if not args.forces:
        raise
    print("force fit failed:", str(ex)[:200], "-> retry energies only", flush=True)
    ccs_fetch(mode="DFTB", R_c=args.rc, Ns="all", DFT_DB=f"{args.prefix}_dft.db",
              DFTB_DB=f"{args.prefix}_dftb.db", include_forces=False)
    ccs_fit("CCS_input.json")
os.get_terminal_size = lambda *a: os.terminal_size((100, 40))
sys.argv = ["x", "CCS_params.json"]
from ccs_fit.scripts.ccs_export_sktable import main as export_main
export_main()
print("spl files:", sorted(f for f in os.listdir(".") if f.endswith(".spl")))
try:
    import numpy as np
    err = np.loadtxt("error.out")
    de = err[:, 0] - err[:, 1]
    print(f"fit residual (energies): RMS={np.sqrt((de**2).mean())*1000:.1f} meV, max={abs(de).max()*1000:.1f} meV, n={len(de)}")
except Exception as e:
    print("error summary skipped:", e)
print("REPFIT_DONE")
