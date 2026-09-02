#!/usr/bin/env python3
"""DFTB+ で単層 MoS2 のバンド構造を計算し、固有値配列を返す。

使い方:
  python3 dftb_bands.py <skf_dir> <a_lattice> <workdir> [--json out.json]

手順: SCC 計算 (12x12x1) -> 電荷読込みでバンドパス 1 発計算 -> band.out パース。
出力 JSON: {"kpts": [[...]], "eigs": [[e1, e2, ...], ...] (eV), "efermi": float}
"""
import argparse
import json
import os
import subprocess
import sys

import numpy as np
from ase.build import mx2
from ase.io import write as ase_write

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bandpath_common import kpath_points  # noqa: E402

DFTB = "/Users/crocus/uhuhu/MoS2_DFTB/sw_local/dftbplus-install/bin/dftb+"

HSD_COMMON = """
Hamiltonian = DFTB {{
  SCC = Yes
{shell}  SCCTolerance = {scctol}
  MaxSCCIterations = {maxscc}
  Mixer = Broyden {{}}
  SlaterKosterFiles = Type2FileNames {{
    Prefix = "{skf}/"
    Separator = "-"
    Suffix = ".skf"
  }}
  MaxAngularMomentum {{
    Mo = "d"
    S = "{slmax}"
  }}
  PolynomialRepulsive = SetForAll {{ Yes }}
  Filling = Fermi {{ Temperature [K] = 40 }}
{kblock}
{extra}
}}
Options {{ WriteResultsTag = Yes }}
Analysis {{ WriteBandOut = Yes }}
ParserOptions {{ ParserVersion = 14 }}
"""


def SHELL(args):
    return "  ShellResolvedSCC = Yes\n" if getattr(args, "shell", False) else ""


def make_geometry(a, workdir, thickness=3.13):
    atoms = mx2(formula="MoS2", kind="2H", a=a, thickness=thickness, vacuum=10.0)
    ase_write(os.path.join(workdir, "geo.gen"), atoms, format="gen")
    return atoms


def run_dftb(workdir, hsd):
    with open(os.path.join(workdir, "dftb_in.hsd"), "w") as f:
        f.write(hsd)
    env = dict(os.environ, OMP_NUM_THREADS="4", OPENBLAS_NUM_THREADS="1")
    r = subprocess.run([DFTB], cwd=workdir, capture_output=True, text=True,
                       timeout=1200, env=env)
    ok = "ERROR" not in r.stdout and r.returncode == 0
    return ok, r.stdout + r.stderr


def parse_band_out(path):
    """band.out -> (nk, nbands) eV 配列と占有数配列"""
    eigs, occs = [], []
    cur_e, cur_o = None, None
    with open(path) as f:
        for line in f:
            t = line.split()
            if not t:
                continue
            if t[0] == "KPT":
                if cur_e:
                    eigs.append(cur_e)
                    occs.append(cur_o)
                cur_e, cur_o = [], []
            elif len(t) >= 3 and cur_e is not None:
                cur_e.append(float(t[1]))
                cur_o.append(float(t[2]))
    if cur_e:
        eigs.append(cur_e)
        occs.append(cur_o)
    return np.array(eigs), np.array(occs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("skf_dir")
    ap.add_argument("a", type=float)
    ap.add_argument("workdir")
    ap.add_argument("--json", default=None)
    ap.add_argument("--thickness", type=float, default=3.13)
    ap.add_argument("--s-lmax", choices=["p", "d"], default="p", dest="slmax")
    ap.add_argument("--shell-resolved", action="store_true", dest="shell")
    ap.add_argument("--extra-kpts", default=None, dest="extra_kpts",
                    help="JSON list of extra k-points appended after the path")
    args = ap.parse_args()

    os.makedirs(args.workdir, exist_ok=True)
    skf = os.path.abspath(args.skf_dir)
    make_geometry(args.a, args.workdir, args.thickness)

    geoblock = 'Geometry = GenFormat {\n  <<< "geo.gen"\n}\n'

    # 1) SCC
    kscc = ("  KPointsAndWeights = SupercellFolding {\n"
            "    12 0 0\n    0 12 0\n    0 0 1\n    0.5 0.5 0.0\n  }")
    hsd1 = geoblock + HSD_COMMON.format(skf=skf, kblock=kscc, extra="", maxscc=200, scctol="1e-7", slmax=args.slmax, shell=SHELL(args))
    ok, log = run_dftb(args.workdir, hsd1)
    if not ok:
        print("SCC FAILED", log[-1500:])
        sys.exit(1)

    # 2) band path (read charges, 1 iteration)
    if not os.path.exists(os.path.join(args.workdir, "charges.bin")):
        print("NO charges.bin AFTER SCC")
        sys.exit(1)
    pts = np.array(kpath_points())
    if args.extra_kpts:
        pts = np.vstack([pts, np.array(json.load(open(args.extra_kpts)))])
    klines = "  KPointsAndWeights = {\n" + "\n".join(
        f"    {p[0]:.10f} {p[1]:.10f} {p[2]:.10f} 1.0" for p in pts) + "\n  }"
    extra = "  ReadInitialCharges = Yes"
    hsd2 = geoblock + HSD_COMMON.format(skf=skf, kblock=klines, extra=extra, maxscc=1, scctol="1e8", slmax=args.slmax, shell=SHELL(args))
    ok, log = run_dftb(args.workdir, hsd2)
    if not ok:
        print("BAND RUN FAILED", log[-1500:])
        sys.exit(1)

    eigs, occs = parse_band_out(os.path.join(args.workdir, "band.out"))
    # VBM/CBM from occupations
    occupied = occs > 0.5
    vbm = eigs[occupied].max()
    cbm = eigs[~occupied].min()
    out = {
        "kpts": pts.tolist(),
        "eigs": eigs.tolist(),
        "occs": occs.tolist(),
        "vbm": float(vbm),
        "cbm": float(cbm),
        "gap": float(cbm - vbm),
    }
    if args.json:
        json.dump(out, open(args.json, "w"))
    print(f"nk={eigs.shape[0]} nb={eigs.shape[1]} VBM={vbm:.4f} "
          f"CBM={cbm:.4f} gap={cbm-vbm:.4f} eV")


if __name__ == "__main__":
    main()
