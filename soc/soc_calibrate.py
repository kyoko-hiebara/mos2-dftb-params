#!/usr/bin/env python3
"""DFTB+ SpinOrbit 定数の較正: K 点 VB 分裂を LAK+SOC 参照 (150.2 meV) に合わせる。

ξ_S(3p) は原子的な値に固定し、ξ_Mo(4d) を 2 点計算 + 線形補間で決定。
"""
import os
import subprocess
import sys

import numpy as np

SCRIPTS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPTS)
from optimize_multi import DFTB, HA, parse_tag_eigs  # noqa: E402

from ase.build import mx2  # noqa: E402
from ase.io import write as ase_write  # noqa: E402

ROOT = "/Users/crocus/uhuhu/MoS2_DFTB"
SKF = f"{ROOT}/local_opt/skf_v3"
TARGET_VB_SPLIT = 0.1502  # eV (LAK+SOC @K)
XI_MO_P = 0.036           # eV, 固定 (原子的初期値)
XI_S_P = 0.055            # eV, 固定

HSD = """Geometry = GenFormat {{
  <<< "geo.gen"
}}
Hamiltonian = DFTB {{
  SCC = Yes
  SCCTolerance = {tol}
  MaxSCCIterations = {mx}
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
  SpinOrbit = {{
    Dual = Yes
    Mo [eV] = {{0.0 {xi_mo_p} {xi_mo_d}}}
    S [eV] = {{0.0 {xi_s_p} 0.0}}
  }}
  Filling = Fermi {{ Temperature [K] = 40 }}
{kblock}
{extra}
}}
Options {{ WriteResultsTag = Yes }}
ParserOptions {{ ParserVersion = 14 }}
"""


def run_soc(xi_mo_d, wd):
    os.makedirs(wd, exist_ok=True)
    atoms = mx2(formula="MoS2", kind="2H", a=3.16, thickness=3.127, vacuum=10.0)
    ase_write(os.path.join(wd, "geo.gen"), atoms, format="gen")
    env = dict(os.environ, OMP_NUM_THREADS="8")
    # 1) SCC
    kscc = ("  KPointsAndWeights = SupercellFolding {\n"
            "    12 0 0\n    0 12 0\n    0 0 1\n    0.5 0.5 0.0\n  }")
    hsd = HSD.format(skf=SKF, kblock=kscc, tol="1e-6", mx="200", extra="",
                     xi_mo_p=XI_MO_P, xi_mo_d=xi_mo_d, xi_s_p=XI_S_P)
    open(os.path.join(wd, "dftb_in.hsd"), "w").write(hsd)
    r = subprocess.run([DFTB], cwd=wd, capture_output=True, text=True,
                       timeout=900, env=env)
    if r.returncode != 0:
        raise RuntimeError(r.stdout[-500:])
    # 2) K 点 1 発
    kK = ("  KPointsAndWeights = {\n"
          "    0.33333333 0.33333333 0.0 1.0\n  }")
    hsd2 = HSD.format(skf=SKF, kblock=kK, tol="1e6", mx="1",
                      extra="  ReadInitialCharges = Yes",
                      xi_mo_p=XI_MO_P, xi_mo_d=xi_mo_d, xi_s_p=XI_S_P)
    open(os.path.join(wd, "dftb_in.hsd"), "w").write(hsd2)
    r = subprocess.run([DFTB], cwd=wd, capture_output=True, text=True,
                       timeout=900, env=env)
    if r.returncode != 0:
        raise RuntimeError(r.stdout[-500:])
    eig, occ = parse_tag_eigs(os.path.join(wd, "results.tag"))
    e = np.sort(eig.flatten()) * HA
    o = occ.flatten()
    nelec = int(round(o.sum()))
    vb = e[:nelec]
    cb = e[nelec:]
    return vb[-1] - vb[-2], cb[1] - cb[0], cb[0] - vb[-1]


def main():
    # LAK 参照の CB 分裂も表示 (EIGENVAL から)
    lines = open(f"{ROOT}/results/ref_calc/soc_lak/EIGENVAL").readlines()
    nelec, nk, nb = (int(x) for x in lines[5].split())
    idx = 7
    for ik in range(nk):
        k = [float(x) for x in lines[idx].split()[:3]]
        if abs(k[0] - 1 / 3) < 1e-3 and abs(k[1] - 1 / 3) < 1e-3:
            eigs = sorted(float(lines[idx + 1 + ib].split()[1])
                          for ib in range(nb))
            vb, cb = eigs[:nelec], eigs[nelec:]
            print(f"LAK+SOC @K: VB split {1000*(vb[-1]-vb[-2]):.1f} meV, "
                  f"CB split {1000*(cb[1]-cb[0]):.1f} meV, "
                  f"gap {cb[0]-vb[-1]:.4f} eV")
            break
        idx += nb + 2

    xis = [0.070, 0.100]
    splits = []
    for xi in xis:
        vbs, cbs, gap = run_soc(xi, f"/tmp/soc_{int(xi*1000)}")
        splits.append(vbs)
        print(f"xi_Mo(4d)={xi:.3f} eV -> VB split {vbs*1000:.1f} meV, "
              f"CB split {cbs*1000:.1f} meV, K gap {gap:.4f} eV", flush=True)
    # 線形補間
    slope = (splits[1] - splits[0]) / (xis[1] - xis[0])
    xi_opt = xis[0] + (TARGET_VB_SPLIT - splits[0]) / slope
    vbs, cbs, gap = run_soc(xi_opt, "/tmp/soc_opt")
    print(f"\ncalibrated xi_Mo(4d) = {xi_opt:.4f} eV")
    print(f"  -> VB split {vbs*1000:.1f} meV (target {TARGET_VB_SPLIT*1000:.1f})")
    print(f"  -> CB split {cbs*1000:.1f} meV, K gap (SOC) {gap:.4f} eV")
    print("\nDFTB+ snippet:")
    print("  SpinOrbit = {")
    print("    Dual = Yes")
    print(f"    Mo [eV] = {{0.0 {XI_MO_P} {xi_opt:.4f}}}")
    print(f"    S [eV] = {{0.0 {XI_S_P} 0.0}}")
    print("  }")


if __name__ == "__main__":
    main()
