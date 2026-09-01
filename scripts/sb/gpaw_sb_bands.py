#!/usr/bin/env python3
"""GPAW で バルク Sb (A7, R-3m) の PBE 参照バンドを計算し JSON 出力。

半金属なのでフェルミレベル整列 (E_F = 0)。
出力: sb_ref_bands.json {kpts, eigs (E-E_F, eV), path_labels, nelec}
"""
import json
import os

import numpy as np
from ase.spacegroup import crystal
from gpaw import GPAW, PW, FermiDirac

ROOT = "/Users/crocus/uhuhu/MoS2_DFTB"
os.environ.setdefault("GPAW_SETUP_PATH",
                      f"{ROOT}/sw_local/gpaw-setups/gpaw-setups-24.11.0")

# 実験構造 (hex setting): a=4.3084, c=11.274, z(Sb)=0.2336
sb = crystal("Sb", [(0, 0, 0.2336)], spacegroup=166,
             cellpar=[4.3084, 4.3084, 11.274, 90, 90, 120],
             primitive_cell=True)
print("atoms:", len(sb), "cell:", np.round(sb.cell.lengths(), 4))

calc = GPAW(mode=PW(350), xc="PBE", kpts=(10, 10, 10),
            occupations=FermiDirac(0.02), txt="sb_scf.txt", symmetry="off")
sb.calc = calc
sb.get_potential_energy()
ef = calc.get_fermi_level()
print(f"SCF done, E_F = {ef:.4f} eV")

path = sb.cell.bandpath(npoints=80)
print("path:", path.path)
bcalc = calc.fixed_density(kpts=path.kpts, symmetry="off",
                           nbands=26, convergence={"bands": 22},
                           txt="sb_bands.txt")
e_kn = np.array([bcalc.get_eigenvalues(kpt=k) for k in range(len(path.kpts))])
e_kn -= ef

xcoords, label_x, labels = path.get_linear_kpoint_axis()
out = {
    "kpts": path.kpts.tolist(),
    "eigs": e_kn.tolist(),
    "xcoords": xcoords.tolist(),
    "label_x": list(map(float, label_x)),
    "labels": list(labels),
    "nelec": 30,
    "path": path.path,
}
json.dump(out, open(f"{ROOT}/local_opt/sb_ref_bands.json", "w"))
print("wrote sb_ref_bands.json;",
      f"bands at Gamma (E-EF): {np.round(e_kn[0][8:22], 2)}")
