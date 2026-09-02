#!/usr/bin/env python3
"""Sb(111) 6 バイレイヤースラブの PBE 仕事関数 (GPAW, mac)。
出力: local_opt/sb_slab_wf.json {W, efermi, vacuum, ...} と sb_slab.traj"""
import json
import os

import numpy as np
from ase.io import write
from ase.spacegroup import crystal
from gpaw import GPAW, PW, FermiDirac
from gpaw.mpi import world

ROOT = "/Users/crocus/uhuhu/MoS2_DFTB"
os.environ.setdefault("GPAW_SETUP_PATH",
                      f"{ROOT}/sw_local/gpaw-setups/gpaw-setups-24.11.0")

conv = crystal("Sb", [(0, 0, 0.2336)], spacegroup=166,
               cellpar=[4.3084, 4.3084, 11.274, 90, 90, 120],
               primitive_cell=False)           # 6 原子 = 3 バイレイヤー
slab = conv.repeat((1, 1, 2))                 # 12 原子 = 6 バイレイヤー
zs = slab.positions[:, 2]
thick = zs.max() - zs.min()
VAC = 15.0
cell = slab.cell.copy()
cell[2] = [0, 0, thick + 2 * VAC]
slab.set_cell(cell)
slab.positions[:, 2] += VAC - zs.min()
slab.pbc = (True, True, True)
write(f"{ROOT}/local_opt/sb_slab.traj", slab)   # 全ランクで呼ぶ (ASE 側で master のみ書く)
if world.rank == 0: print("slab:", len(slab), "atoms, thickness %.2f A, cell c %.2f A" % (thick, cell[2, 2]), flush=True)

calc = GPAW(mode=PW(350), xc="PBE", kpts=(8, 8, 1),
            occupations=FermiDirac(0.05),
            txt=f"{ROOT}/local_opt/sb_slab_wf.txt")
slab.calc = calc
e = slab.get_potential_energy()
ef = calc.get_fermi_level()
v = calc.get_electrostatic_potential()       # (nx, ny, nz) eV
vz = v.mean(axis=(0, 1))
nz = len(vz)
# 真空中央 (z=0 と z=c の境界付近) の平均
imid = [0, 1, 2, nz - 3, nz - 2, nz - 1]
vac = float(np.mean(vz[imid]))
W = vac - ef
out = dict(W=W, efermi=ef, vacuum=vac, energy=e, natoms=len(slab),
           thickness=thick, vz=vz.tolist())
if world.rank == 0:
    json.dump(out, open(f"{ROOT}/local_opt/sb_slab_wf.json", "w"), indent=1)
if world.rank == 0: print(f"Sb(111) 6BL PBE: E_F = {ef:.3f} eV, V_vac = {vac:.3f} eV, W = {W:.3f} eV "
      f"(expt 4.55-4.7)", flush=True)
