#!/usr/bin/env python3
"""skf_v2 で参照構造群の DFTB 電子エネルギー+力を計算し ASE DB に格納。

- 反発ゼロ (PolynomialRepulsive=SetForAll Yes) の SCC 計算
- E_rep^target = E_LAK - E_elec を ccs_fit に渡すための DFTB 側データ

使い方:
  python3 dftb_energy.py <dft_db> <out_prefix> <skf_dir> [--workdir tmp_dftb]
成功した構造のみを <out_prefix>_dft.db / <out_prefix>_dftb.db のペアで出力
(ccs_fetch が id で突き合わせるため、両 DB の行順を一致させる)。
"""
import argparse
import os
import shutil
import subprocess
import sys

import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator
from ase.db import connect
from ase.io import write as ase_write

DFTB = "/Users/crocus/uhuhu/MoS2_DFTB/sw_local/dftbplus-install/bin/dftb+"

HSD = """Geometry = GenFormat {{
  <<< "geo.gen"
}}
Hamiltonian = DFTB {{
  SCC = Yes
  SCCTolerance = 1e-6
  MaxSCCIterations = 250
  Mixer = Broyden {{}}
  SlaterKosterFiles = Type2FileNames {{
    Prefix = "{skf}/"
    Separator = "-"
    Suffix = ".skf"
  }}
  MaxAngularMomentum {{
    Mo = "d"
    S = "d"
  }}
  PolynomialRepulsive = SetForAll {{ Yes }}
  Filling = Fermi {{ Temperature [K] = 100 }}
{kblock}
}}
Options {{ WriteResultsTag = Yes }}
Analysis {{ PrintForces = Yes }}
ParserOptions {{ ParserVersion = 14 }}
"""


def kblock_for(atoms):
    if not any(atoms.pbc):
        return "  # cluster"
    n = len(atoms)
    if n <= 4:
        k = 12
    elif n <= 30:
        k = 4
    else:
        k = 3
    return ("  KPointsAndWeights = SupercellFolding {\n"
            f"    {k} 0 0\n    0 {k} 0\n    0 0 1\n    0.5 0.5 0.0\n  }}")


def parse_results_tag(path):
    """results.tag から total_energy と forces を抽出"""
    energy, forces = None, None
    with open(path) as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        ln = lines[i]
        if ln.startswith("extrapolated0_energy"):
            energy = float(lines[i + 1].split()[0])
            i += 2
            continue
        if ln.startswith("forces"):
            # forces   :real:2:3,N
            shape = ln.split(":")[-1].strip().split(",")
            nat = int(shape[1])
            vals = []
            j = i + 1
            while len(vals) < 3 * nat:
                vals.extend(float(x) for x in lines[j].split())
                j += 1
            forces = np.array(vals).reshape(nat, 3)
            i = j
            continue
        i += 1
    return energy, forces


def run_one(atoms, skf, workdir):
    os.makedirs(workdir, exist_ok=True)
    # 分子はセルを保持しつつ pbc=True の Γ 点計算にする (gen 'S' 形式)
    a = atoms.copy()
    if not any(a.pbc):
        a.pbc = True
    ase_write(os.path.join(workdir, "geo.gen"), a, format="gen")
    kblock = kblock_for(atoms if any(atoms.pbc) else _gamma_only())
    hsd = HSD.format(skf=skf, kblock=kblock)
    if not any(atoms.pbc):
        hsd = hsd.replace("  # cluster",
                          "  KPointsAndWeights = SupercellFolding {\n"
                          "    1 0 0\n    0 1 0\n    0 0 1\n    0.0 0.0 0.0\n  }")
    with open(os.path.join(workdir, "dftb_in.hsd"), "w") as f:
        f.write(hsd)
    env = dict(os.environ, OMP_NUM_THREADS="8")
    r = subprocess.run([DFTB], cwd=workdir, capture_output=True, text=True,
                       timeout=1800, env=env)
    tag = os.path.join(workdir, "results.tag")
    if r.returncode != 0 or not os.path.exists(tag):
        return None, None, r.stdout[-300:]
    e_ha, f_ha = parse_results_tag(tag)
    HA = 27.211386245988
    BOHR = 0.529177210903
    return e_ha * HA, f_ha * HA / BOHR if f_ha is not None else None, None


class _gamma_only:
    pbc = [False]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dft_db")
    ap.add_argument("out_prefix")
    ap.add_argument("skf_dir")
    ap.add_argument("--workdir", default="/tmp/dftb_energy")
    args = ap.parse_args()

    skf = os.path.abspath(args.skf_dir)
    dft_out = args.out_prefix + "_dft.db"
    dftb_out = args.out_prefix + "_dftb.db"
    for p in (dft_out, dftb_out):
        if os.path.exists(p):
            os.remove(p)
    src = connect(args.dft_db)
    dst_dft = connect(dft_out)
    dst_dftb = connect(dftb_out)
    nok = nfail = 0
    for row in src.select():
        atoms = row.toatoms()
        name = row.get("name", str(row.id))
        wd = os.path.join(args.workdir, name)
        shutil.rmtree(wd, ignore_errors=True)
        e, f, err = run_one(atoms, skf, wd)
        if e is None:
            print(f"{name}: DFTB FAILED {err}", flush=True)
            nfail += 1
            continue
        dft_atoms = row.toatoms()
        dft_atoms.calc = SinglePointCalculator(
            dft_atoms, energy=row.energy,
            forces=row.forces if "forces" in row else None)
        dst_dft.write(dft_atoms, name=name)
        out = atoms.copy()
        out.calc = SinglePointCalculator(out, energy=e, forces=f)
        dst_dftb.write(out, name=name)
        print(f"{name}: E_elec = {e:.6f} eV", flush=True)
        nok += 1
        shutil.rmtree(wd, ignore_errors=True)
    print(f"done: {nok} ok, {nfail} failed -> {dft_out} + {dftb_out}")


if __name__ == "__main__":
    main()
