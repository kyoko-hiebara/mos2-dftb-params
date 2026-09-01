#!/usr/bin/env python3
"""GPAW: バルク Sb の SOC 込み参照バンド (非自己無撞着 SOC 後処理)。"""
import json
import os

import numpy as np
from ase.spacegroup import crystal
from gpaw import GPAW, PW, FermiDirac
from gpaw.spinorbit import soc_eigenstates

ROOT = "/Users/crocus/uhuhu/MoS2_DFTB"
os.environ.setdefault(
    "GPAW_SETUP_PATH",
    f"{ROOT}/sw_local/gpaw-setups/gpaw-setups-24.11.0")

sb = crystal("Sb", [(0, 0, 0.2336)], spacegroup=166,
             cellpar=[4.3084, 4.3084, 11.274, 90, 90, 120],
             primitive_cell=True)

calc = GPAW(mode=PW(350), xc="PBE", kpts=(10, 10, 10),
            occupations=FermiDirac(0.02), txt="sb_scf2.txt", symmetry="off")
sb.calc = calc
sb.get_potential_energy()
print("SCF done", flush=True)

path = sb.cell.bandpath(npoints=80)
bcalc = calc.fixed_density(kpts=path.kpts, symmetry="off",
                           nbands=26, convergence={"bands": 22},
                           txt="sb_bands2.txt")

soc = soc_eigenstates(bcalc, n2=22)
e_soc = soc.eigenvalues()          # (nk, 2*n)
ef_soc = soc.fermi_level
e_soc = e_soc - ef_soc
print(f"SOC done: shape {e_soc.shape}, E_F(SOC) = {ef_soc:.4f} eV", flush=True)

xcoords, label_x, labels = path.get_linear_kpoint_axis()
json.dump({
    "kpts": path.kpts.tolist(),
    "eigs": e_soc.tolist(),
    "xcoords": xcoords.tolist(),
    "label_x": list(map(float, label_x)),
    "labels": list(labels),
    "path": path.path,
}, open(f"{ROOT}/local_opt/sb_ref_bands_soc.json", "w"))
print("wrote sb_ref_bands_soc.json;",
      "Gamma window:", np.round(e_soc[0][(np.abs(e_soc[0]) < 6)], 2))
