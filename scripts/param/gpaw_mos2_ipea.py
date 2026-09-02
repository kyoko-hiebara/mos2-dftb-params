#!/usr/bin/env python3
"""単層 MoS2 (a=3.16, LAK 緩和形状) の PBE 真空準位基準バンド端 (GPAW, mac)。
出力: local_opt/mos2_pbe_ipea.json {vacuum, vbm, cbm, IP, EA, midgap_abs, gap_K, ...}"""
import json
import os

import numpy as np
from ase.io import read
from gpaw import GPAW, PW, FermiDirac
from gpaw.mpi import world

ROOT = "/Users/crocus/uhuhu/MoS2_DFTB"
os.environ.setdefault("GPAW_SETUP_PATH",
                      f"{ROOT}/sw_local/gpaw-setups/gpaw-setups-24.11.0")

atoms = read(f"{ROOT}/results/ref_calc/bands_lak_a316/POSCAR")
if world.rank == 0: print("geometry:", atoms.cell.lengths(), atoms.get_chemical_formula(), flush=True)
calc = GPAW(mode=PW(450), xc="PBE", kpts={"size": (12, 12, 1), "gamma": True},
            occupations=FermiDirac(0.01), nbands=24,
            convergence={"bands": 16},
            txt=f"{ROOT}/local_opt/mos2_pbe_ipea.txt")
atoms.calc = calc
e = atoms.get_potential_energy()
ef = calc.get_fermi_level()
v = calc.get_electrostatic_potential()
vz = v.mean(axis=(0, 1))
nz = len(vz)
vac = float(np.mean(vz[[0, 1, 2, nz - 3, nz - 2, nz - 1]]))
nk = len(calc.get_ibz_k_points())
eigs = np.array([calc.get_eigenvalues(kpt=k) for k in range(nk)])
wk = np.array(calc.get_k_point_weights())
occ = np.array([calc.get_occupation_numbers(kpt=k) for k in range(nk)])
occ = occ / wk[:, None]          # k 重みを除いた占有数 (0..2)
omax = occ.max()
vbm = float(eigs[occ > 0.5 * omax].max())
cbm = float(eigs[occ < 0.5 * omax].min())
kpts = calc.get_ibz_k_points()
vb_k = np.array([eigs[k][occ[k] > 0.5 * omax].max() for k in range(nk)])
cb_k = np.array([eigs[k][occ[k] < 0.5 * omax].min() for k in range(nk)])
iK = int(np.argmin(cb_k - vb_k))
# 厳密な Γ / K / M 点は fixed-density の非 SCF 計算で取得
nelec = calc.get_number_of_electrons()
nocc = int(round(nelec)) // 2
bcalc = calc.fixed_density(kpts=[[0, 0, 0], [1 / 3, 1 / 3, 0], [0.5, 0, 0]],
                           symmetry="off", nbands=20, convergence={"bands": 16},
                           txt=f"{ROOT}/local_opt/mos2_pbe_ipea_bands.txt")
eG, eK, eM = [np.array(bcalc.get_eigenvalues(kpt=i)) for i in range(3)]
vbmK, cbmK = float(eK[nocc - 1]), float(eK[nocc])
vbmG, cbmG = float(eG[nocc - 1]), float(eG[nocc])
out = dict(vacuum=vac, efermi=ef, vbm=vbm, cbm=cbm, IP=vac - vbm, EA=vac - cbm,
           midgap_abs=(vbmK + cbmK) / 2 - vac, midgap_mesh_abs=(vbm + cbm) / 2 - vac, gap=cbm - vbm,
           gap_K=cbmK - vbmK, vbm_K_rel_vac=vbmK - vac, cbm_K_rel_vac=cbmK - vac,
           vbm_G_rel_vac=vbmG - vac, cbm_G_rel_vac=cbmG - vac, nelec=float(nelec),
           midgap_K_abs=(vbmK + cbmK) / 2 - vac,
           kK=kpts[iK].tolist(), energy=e, ibz_kpts=np.array(kpts).tolist(),
           ibz_weights=wk.tolist(), vb_k=(vb_k - vac).tolist(), cb_k=(cb_k - vac).tolist(),
           vbm_k_frac=kpts[int(vb_k.argmax())].tolist())
if world.rank == 0:
    json.dump(out, open(f"{ROOT}/local_opt/mos2_pbe_ipea.json", "w"), indent=1)
if world.rank == 0: print(f"MoS2 PBE @a=3.16: V_vac={vac:.3f} VBM={vbm:.3f} CBM={cbm:.3f} -> IP={vac-vbm:.3f} EA={vac-cbm:.3f} "
      f"gap={cbm-vbm:.3f} (K gap {cbmK-vbmK:.3f}); midgap rel. vac = {(vbm+cbm)/2-vac:.3f} eV", flush=True)
