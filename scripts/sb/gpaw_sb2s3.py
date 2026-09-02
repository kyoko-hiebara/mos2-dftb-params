#!/usr/bin/env python3
"""Sb2S3 (輝安鉱, Pnma) の PBE 参照: 実験セルで内部座標を緩和 → Γ-X-S-Y-Γ-Z バンド。
S-Sb ペアの転移性検証用。出力: sb2s3_relaxed.traj, sb2s3_ref_bands.json"""
import json
import os

import numpy as np
from ase.io import read, write
from ase.optimize import BFGS
from gpaw import GPAW, PW, FermiDirac

ROOT = "/Users/crocus/uhuhu/MoS2_DFTB"
os.environ.setdefault("GPAW_SETUP_PATH", f"{ROOT}/sw_local/gpaw-setups/gpaw-setups-24.11.0")
atoms = read(f"{ROOT}/local_opt/sb2s3_expt.traj")
calc = GPAW(mode=PW(350), xc="PBE", kpts=(2, 6, 2), occupations=FermiDirac(0.02),
            txt=f"{ROOT}/local_opt/sb2s3_relax.txt")
atoms.calc = calc
opt = BFGS(atoms, logfile=f"{ROOT}/local_opt/sb2s3_bfgs.log")
opt.run(fmax=0.05, steps=40)
write(f"{ROOT}/local_opt/sb2s3_relaxed.traj", atoms)
e0 = atoms.get_potential_energy()
print("relaxed; E =", e0, "max|F| =", np.abs(atoms.get_forces()).max(), flush=True)
path = atoms.cell.bandpath("GXSYGZ", npoints=60)
nelec = calc.get_number_of_electrons()
nocc = int(round(nelec)) // 2
bcalc = calc.fixed_density(kpts=path.kpts, symmetry="off", nbands=nocc + 10,
                           convergence={"bands": nocc + 6}, txt=f"{ROOT}/local_opt/sb2s3_bands.txt")
eigs = np.array([bcalc.get_eigenvalues(kpt=i) for i in range(len(path.kpts))])
vbm = eigs[:, :nocc].max(); cbm = eigs[:, nocc:].min()
out = dict(kpts=path.kpts.tolist(), eigs=eigs.tolist(), nocc=nocc, vbm=float(vbm), cbm=float(cbm),
           gap=float(cbm - vbm), labels=path.path, xcoords=path.get_linear_kpoint_axis()[0].tolist(),
           label_x=path.get_linear_kpoint_axis()[1].tolist(), label_names=path.get_linear_kpoint_axis()[2])
json.dump(out, open(f"{ROOT}/local_opt/sb2s3_ref_bands.json", "w"))
print(f"Sb2S3 PBE: gap = {cbm - vbm:.3f} eV (expt ~1.7, PBE lit ~1.2-1.4); nocc={nocc}", flush=True)
