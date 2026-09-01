#!/usr/bin/env python3
"""Sb の SpinOrbit 定数 ξ_5p を GPAW+SOC 参照バンドに較正。

A7 Sb は反転対称なので分裂ではなく SOC 込みバンド全体の RMS でフィット。
"""
import json
import os
import subprocess
import sys

import numpy as np

SCRIPTS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPTS)
from optimize_multi import DFTB, HA, parse_tag_eigs  # noqa: E402

from ase.spacegroup import crystal  # noqa: E402
from ase.io import write as ase_write  # noqa: E402

ROOT = "/Users/crocus/uhuhu/MoS2_DFTB"
SKF = f"{ROOT}/local_opt/skf_v3sb"
REF = json.load(open(f"{ROOT}/local_opt/sb_ref_bands_soc.json"))
REF_E = np.array(REF["eigs"])
REF_K = np.array(REF["kpts"])

SB_ATOMS = crystal("Sb", [(0, 0, 0.2336)], spacegroup=166,
                   cellpar=[4.3084, 4.3084, 11.274, 90, 90, 120],
                   primitive_cell=True)

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
    Sb = "d"
  }}
  PolynomialRepulsive = SetForAll {{ Yes }}
  SpinOrbit = {{
    Dual = Yes
    Sb [eV] = {{0.0 {xi_p} {xi_d}}}
  }}
  Filling = Fermi {{ Temperature [K] = 300 }}
  {extra}
{kblock}
}}
Options {{ WriteResultsTag = Yes }}
ParserOptions {{ ParserVersion = 14 }}
"""


def sb_soc_bands(xi_p, xi_d, wd):
    os.makedirs(wd, exist_ok=True)
    ase_write(os.path.join(wd, "geo.gen"), SB_ATOMS, format="gen")
    env = dict(os.environ, OMP_NUM_THREADS="8")
    kscc = ("  KPointsAndWeights = SupercellFolding {\n"
            "    8 0 0\n    0 8 0\n    0 0 8\n    0.5 0.5 0.5\n  }")
    open(os.path.join(wd, "dftb_in.hsd"), "w").write(
        HSD.format(skf=SKF, kblock=kscc, tol="1e-6", mx="250", extra="",
                   xi_p=xi_p, xi_d=xi_d))
    r = subprocess.run([DFTB], cwd=wd, capture_output=True, text=True,
                       timeout=1800, env=env)
    if r.returncode != 0:
        raise RuntimeError("SCC: " + r.stdout[-300:])
    lines = open(os.path.join(wd, "results.tag")).readlines()
    ef = None
    for i, ln in enumerate(lines):
        if ln.startswith("fermi_level"):
            ef = float(lines[i + 1].split()[0]) * HA
            break
    kl = "  KPointsAndWeights = {\n" + "\n".join(
        f"    {k[0]:.8f} {k[1]:.8f} {k[2]:.8f} 1.0" for k in REF_K) + "\n  }"
    open(os.path.join(wd, "dftb_in.hsd"), "w").write(
        HSD.format(skf=SKF, kblock=kl, tol="1e6", mx="1",
                   extra="ReadInitialCharges = Yes", xi_p=xi_p, xi_d=xi_d))
    r = subprocess.run([DFTB], cwd=wd, capture_output=True, text=True,
                       timeout=1800, env=env)
    if r.returncode != 0:
        raise RuntimeError("bands: " + r.stdout[-300:])
    eig, _ = parse_tag_eigs(os.path.join(wd, "results.tag"))
    e = np.squeeze(eig) * HA
    if e.shape[0] != len(REF_K):
        e = e.T
    return np.sort(e, axis=1) - ef


def soc_loss(e_dftb):
    nk = min(len(e_dftb), len(REF_E))
    bmeans = REF_E[:nk].mean(axis=0)
    sel = [ib for ib, m in enumerate(bmeans) if -15.0 < m < 6.0]
    ref = REF_E[:nk][:, sel[:e_dftb.shape[1]]]
    n = min(ref.shape[1], e_dftb.shape[1])
    ref, dt = ref[:, :n], e_dftb[:nk, :n]
    w = np.where(np.abs(ref) < 2.0, 1.0, 0.25)
    rms_all = float(np.sqrt(((dt - ref) ** 2 * w).sum() / w.sum()))
    msk = np.abs(ref) < 2.0
    rms_ef = float(np.sqrt(((dt[msk] - ref[msk]) ** 2).mean()))
    return rms_all, rms_ef


def main():
    print(f"reference SOC bands: {REF_E.shape}", flush=True)
    xis = [0.30, 0.45, 0.60, 0.75, 0.90]
    losses = []
    for xi in xis:
        e = sb_soc_bands(xi, 0.0, f"/tmp/sbsoc_{int(xi*100)}")
        ra, re = soc_loss(e)
        losses.append(re)
        print(f"xi_5p = {xi:.2f} eV -> rms_EF {re:.4f}, rms_all {ra:.4f}",
              flush=True)
    c = np.polyfit(xis, losses, 2)
    xi_opt = float(np.clip(-c[1] / (2 * c[0]), 0.2, 1.0))
    e = sb_soc_bands(xi_opt, 0.0, "/tmp/sbsoc_opt")
    ra, re = soc_loss(e)
    print(f"\ncalibrated xi_Sb(5p) = {xi_opt:.3f} eV "
          f"-> rms_EF {re:.4f}, rms_all {ra:.4f}")
    np.save(f"{ROOT}/local_opt/sb_dftb_bands_soc.npy", e)
    json.dump({"xi_5p": xi_opt, "rms_ef": re, "rms_all": ra},
              open(f"{ROOT}/local_opt/sb_soc_calib.json", "w"), indent=1)
    print("\nDFTB+ snippet (full 4-element set):")
    print("  SpinOrbit = {")
    print("    Dual = Yes")
    print("    Mo [eV] = {0.0 0.036 0.0931}")
    print("    S [eV] = {0.0 0.055 0.0}")
    print("    O [eV] = {0.0 0.02 0.0}")
    print(f"    Sb [eV] = {{0.0 {xi_opt:.3f} 0.0}}")
    print("  }")


if __name__ == "__main__":
    main()
