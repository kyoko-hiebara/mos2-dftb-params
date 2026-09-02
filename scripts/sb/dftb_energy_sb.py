#!/usr/bin/env python3
"""任意元素・3D 周期系対応の DFTB エネルギー+力計算 (反発ゼロ、または --keep-rep で SKF 内スプライン使用)。
使い方: python dftb_energy_sb.py <dft_db> <out_prefix> <skf_dir> [--keep-rep] [--kspacing 0.2]
出力: <out_prefix>_dft.db / <out_prefix>_dftb.db (行順一致)"""
import argparse
import os
import shutil
import subprocess
import sys

import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator
from ase.db import connect
from ase.io import write as ase_write

SCRIPTS = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, SCRIPTS)
from dftb_energy_local import parse_results_tag, DFTB  # noqa: E402
LMAX = {"Mo": "d", "S": "d", "O": "p", "Sb": "d", "H": "s"}
HA, BOHR = 27.211386245988, 0.529177210903

HSD = """Geometry = GenFormat {{
  <<< "geo.gen"
}}
Hamiltonian = DFTB {{
  SCC = Yes
  SCCTolerance = 1e-6
  MaxSCCIterations = 300
  Mixer = Broyden {{}}
  SlaterKosterFiles = Type2FileNames {{
    Prefix = "{skf}/"
    Separator = "-"
    Suffix = ".skf"
  }}
  MaxAngularMomentum {{
{lmax}
  }}
{rep}
  Filling = Fermi {{ Temperature [K] = 100 }}
  KPointsAndWeights = SupercellFolding {{
    {k0} 0 0
    0 {k1} 0
    0 0 {k2}
    0.5 0.5 0.5
  }}
}}
Options {{ WriteResultsTag = Yes }}
Analysis {{ PrintForces = Yes }}
ParserOptions {{ ParserVersion = 14 }}
"""
HSD = HSD.replace(
    "  PolynomialRepulsive", (os.environ.get("DFTB_EXTRA_HSD", "").replace("{", "{{").replace("}", "}}") + "\n"
                              if os.environ.get("DFTB_EXTRA_HSD") else "") + "  PolynomialRepulsive", 1)


def kmesh(atoms, spacing):
    rc = atoms.cell.reciprocal() * 2 * np.pi
    return [max(1, int(np.ceil(np.linalg.norm(b) / spacing))) for b in rc]


def run_one(atoms, skf, wd, keep_rep, spacing):
    os.makedirs(wd, exist_ok=True)
    ase_write(os.path.join(wd, "geo.gen"), atoms, format="gen")
    els = sorted(set(atoms.get_chemical_symbols()))
    lmax = "\n".join(f'    {e} = "{LMAX[e]}"' for e in els)
    k = kmesh(atoms, spacing)
    if keep_rep:
        lines = []
        for i, e1 in enumerate(els):
            for e2 in els[i:]:
                has = "Spline" in open(os.path.join(skf, f"{e1}-{e2}.skf")).read()
                lines.append(f"    {e1}-{e2} = {'No' if has else 'Yes'}")
                if e1 != e2:
                    lines.append(f"    {e2}-{e1} = {'No' if has else 'Yes'}")
        rep = "  PolynomialRepulsive = {\n" + "\n".join(lines) + "\n  }"
    else:
        rep = "  PolynomialRepulsive = SetForAll { Yes }"
    open(os.path.join(wd, "dftb_in.hsd"), "w").write(HSD.format(skf=skf, lmax=lmax, rep=rep, k0=k[0], k1=k[1], k2=k[2]))
    r = subprocess.run([DFTB], cwd=wd, capture_output=True, text=True, timeout=1800,
                       env=dict(os.environ, OMP_NUM_THREADS="6"))
    tag = os.path.join(wd, "results.tag")
    if r.returncode != 0 or not os.path.exists(tag):
        return None, None, r.stdout[-300:]
    e, f = parse_results_tag(tag)
    return e * HA, (f * HA / BOHR if f is not None else None), None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dft_db"); ap.add_argument("out_prefix"); ap.add_argument("skf_dir")
    ap.add_argument("--keep-rep", action="store_true")
    ap.add_argument("--kspacing", type=float, default=0.2)
    ap.add_argument("--workdir", default="/private/tmp/claude-501/-Users-crocus-uhuhu-MoS2-DFTB/6416faee-14e1-48c8-9838-35c06f3e7813/scratchpad/dftb_energy_sb")
    args = ap.parse_args()
    skf = os.path.abspath(args.skf_dir)
    for p in (args.out_prefix + "_dft.db", args.out_prefix + "_dftb.db"):
        if os.path.exists(p):
            os.remove(p)
    src = connect(args.dft_db)
    d1 = connect(args.out_prefix + "_dft.db"); d2 = connect(args.out_prefix + "_dftb.db")
    n = 0
    for row in src.select():
        atoms = row.toatoms(); name = row.get("name", f"row{row.id}")
        e, f, err = run_one(atoms, skf, os.path.join(args.workdir, name), args.keep_rep, args.kspacing)
        if e is None:
            print("FAILED", name, err, flush=True); continue
        a1 = atoms.copy(); a1.calc = SinglePointCalculator(a1, energy=row.energy, forces=row.forces)
        a2 = atoms.copy(); a2.calc = SinglePointCalculator(a2, energy=e, forces=f)
        d1.write(a1, name=name); d2.write(a2, name=name); n += 1
        print(f"{name}: E_dft={row.energy:.4f} E_dftb={e:.4f} (k={kmesh(atoms, args.kspacing)})", flush=True)
    print("DONE", n, "structures")


if __name__ == "__main__":
    main()
