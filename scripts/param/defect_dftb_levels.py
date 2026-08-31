#!/usr/bin/env python3
"""欠陥スーパーセルの DFTB 固有値 (ギャップ周辺) を計算し LAK と比較。

使い方: python3 defect_dftb_levels.py <POSCAR> <skf_dir> <workdir> <nelec_dftb>
"""
import os
import subprocess
import sys

import numpy as np
from ase.io import read, write as ase_write

DFTB = "/opt/dftbplus/bin/dftb+"
HA = 27.211386245988

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
{angmom}
  }}
  PolynomialRepulsive = SetForAll {{ Yes }}
  Filling = Fermi {{ Temperature [K] = 50 }}
  KPointsAndWeights = SupercellFolding {{
    2 0 0
    0 2 0
    0 0 1
    0.0 0.0 0.0
  }}
}}
Options {{ WriteResultsTag = Yes }}
ParserOptions {{ ParserVersion = 14 }}
Parallel {{ UseOmpThreads = Yes }}
"""


def parse_tag_eigs(path):
    lines = open(path).readlines()
    eig, occ, shape_e = None, None, None
    i = 0
    while i < len(lines):
        ln = lines[i]
        if ln.startswith("eigenvalues"):
            shape = [int(x) for x in ln.split(":")[-1].strip().split(",")]
            n = np.prod(shape)
            vals = []
            j = i + 1
            while len(vals) < n:
                vals.extend(float(x) for x in lines[j].split())
                j += 1
            eig = np.array(vals).reshape(shape[::-1])  # (spin, nk, nb) 逆順
            shape_e = shape
            i = j
            continue
        if ln.startswith("occupations") or ln.startswith("filling"):
            shape = [int(x) for x in ln.split(":")[-1].strip().split(",")]
            n = np.prod(shape)
            vals = []
            j = i + 1
            while len(vals) < n:
                vals.extend(float(x) for x in lines[j].split())
                j += 1
            occ = np.array(vals).reshape(shape[::-1])
            i = j
            continue
        i += 1
    return eig, occ


def main():
    poscar, skf, workdir = sys.argv[1], os.path.abspath(sys.argv[2]), sys.argv[3]
    os.makedirs(workdir, exist_ok=True)
    atoms = read(poscar)
    ase_write(os.path.join(workdir, "geo.gen"), atoms, format="gen")
    lmax = {"Mo": "d", "S": "d", "O": "p", "H": "s"}
    angmom = "\n".join(f'    {el} = "{lmax[el]}"'
                        for el in sorted(set(atoms.get_chemical_symbols())))
    with open(os.path.join(workdir, "dftb_in.hsd"), "w") as f:
        f.write(HSD.format(skf=skf, angmom=angmom))
    env = dict(os.environ, OMP_NUM_THREADS="8", OPENBLAS_NUM_THREADS="1")
    r = subprocess.run([DFTB], cwd=workdir, capture_output=True, text=True,
                       timeout=3600, env=env)
    if r.returncode != 0:
        print("DFTB FAILED", r.stdout[-500:])
        sys.exit(1)
    eig, occ = parse_tag_eigs(os.path.join(workdir, "results.tag"))
    # eig shape: (1, nk, nb) -> (nk, nb)
    e = eig[0] if eig.ndim == 3 else eig
    o = occ[0] if occ is not None and occ.ndim == 3 else occ
    e = e * HA
    nk, nb = e.shape[-2], e.shape[-1]
    if e.shape[0] == nb:  # 転置チェック
        e = e.T
        o = o.T if o is not None else None
    vbm = e[o > 0.5].max()
    cbm = e[o < 0.5].min()
    print(f"nk={e.shape[0]} nb={e.shape[1]}")
    print(f"highest occ = {vbm:.4f}, lowest empty = {cbm:.4f}, gap = {cbm-vbm:.4f} eV")
    nocc_g = int((o[0] > 0.5).sum())
    print("Γ-point levels around gap:")
    for ib in range(nocc_g - 4, min(nocc_g + 5, e.shape[1])):
        print(f"  band {ib+1}: {e[0, ib]:8.4f}  occ={o[0, ib]:.1f}")
    print("dispersion (min..max over k):")
    for ib in range(nocc_g - 3, min(nocc_g + 4, e.shape[1])):
        print(f"  band {ib+1}: {e[:, ib].min():8.4f} .. {e[:, ib].max():8.4f}"
              f"  (width {e[:, ib].max()-e[:, ib].min():.3f})  occ~{o[:, ib].mean():.1f}")


if __name__ == "__main__":
    main()
