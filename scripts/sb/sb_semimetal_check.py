#!/usr/bin/env python3
"""Sb 半金属性の決定的検証: BZ 全体 (16^3) で状態カウントベースの
バンドオーバーラップ (max E_30th - min E_31st) を GPAW+SOC で計算。"""
import os

import numpy as np
from ase.spacegroup import crystal
from gpaw import GPAW, PW, FermiDirac
from gpaw.spinorbit import soc_eigenstates

ROOT = "/Users/crocus/uhuhu/MoS2_DFTB"
os.environ.setdefault(
    "GPAW_SETUP_PATH", f"{ROOT}/sw_local/gpaw-setups/gpaw-setups-24.11.0")

sb = crystal("Sb", [(0, 0, 0.2336)], spacegroup=166,
             cellpar=[4.3084, 4.3084, 11.274, 90, 90, 120],
             primitive_cell=True)

calc = GPAW(mode=PW(350), xc="PBE", kpts=(10, 10, 10),
            occupations=FermiDirac(0.02), txt="sb_scf3.txt", symmetry="off")
sb.calc = calc
sb.get_potential_energy()

# 16^3 全 BZ メッシュ (uniform, no symmetry)
n = 16
kpts = np.array([[i / n, j / n, l / n]
                 for i in range(n) for j in range(n) for l in range(n)])
bcalc = calc.fixed_density(kpts=kpts, symmetry="off", nbands=22,
                           convergence={"bands": 18}, txt="sb_mesh.txt")
soc = soc_eigenstates(bcalc, n2=18)
e = soc.eigenvalues()          # (nk, 36)
ef = soc.fermi_level
e = np.sort(e, axis=1) - ef

# 30 電子/セル -> SOC では 30 状態占有 (index 29 が最高占有)
E_ho = e[:, 29]
E_lu = e[:, 30]
overlap = E_ho.max() - E_lu.min()
k_h = kpts[E_ho.argmax()]
k_e = kpts[E_lu.argmin()]
print(f"max E_30th = {E_ho.max():+.4f} eV at k = {np.round(k_h, 3)}")
print(f"min E_31st = {E_lu.min():+.4f} eV at k = {np.round(k_e, 3)}")
print(f"band overlap (state counting) = {overlap*1000:+.1f} meV "
      f"(semimetal if > 0; expt ~ +180 meV)")
print(f"direct gap min over BZ = {(E_lu - E_ho).min()*1000:.1f} meV")
nh = (E_ho > 0).sum()
ne = (E_lu < 0).sum()
print(f"k-points with holes: {nh}/{len(kpts)}, with electrons: {ne}/{len(kpts)}")
np.savez(f"{ROOT}/local_opt/sb_mesh_soc.npz", kpts=kpts, e=e)
