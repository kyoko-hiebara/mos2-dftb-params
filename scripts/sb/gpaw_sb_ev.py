#!/usr/bin/env python3
"""Sb 系反発フィット用 PBE 参照 (GPAW, mac): バルク Sb の等方スケール x 内部座標 u、
Sb2S3 (緩和構造) の等方スケール + ランダム変位 (力込み)。出力: repfit_sb/sb_ref_dft.db"""
import os

import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator
from ase.db import connect
from ase.io import read
from ase.spacegroup import crystal
from gpaw import GPAW, PW, FermiDirac

ROOT = "/Users/crocus/uhuhu/MoS2_DFTB"
os.environ.setdefault("GPAW_SETUP_PATH", f"{ROOT}/sw_local/gpaw-setups/gpaw-setups-24.11.0")
OUT = f"{ROOT}/local_opt/repfit_sb"
os.makedirs(OUT, exist_ok=True)
dbp = f"{OUT}/sb_ref_dft.db"
db = connect(dbp)
done = {r.name for r in db.select()}


def store(name, atoms, calc_kw, txt):
    if name in done:
        print("skip", name, flush=True); return
    atoms = atoms.copy()
    atoms.calc = GPAW(mode=PW(350), xc="PBE", occupations=FermiDirac(0.05), txt=txt, **calc_kw)
    e = atoms.get_potential_energy(); f = atoms.get_forces()
    atoms.calc = SinglePointCalculator(atoms, energy=e, forces=f)
    db.write(atoms, name=name)
    print(f"{name}: E={e:.4f} eV maxF={np.abs(f).max():.3f}", flush=True)


# --- bulk Sb (A7) ---
def sb_bulk(scale, u):
    return crystal("Sb", [(0, 0, u)], spacegroup=166,
                   cellpar=[4.3084 * scale, 4.3084 * scale, 11.274 * scale, 90, 90, 120], primitive_cell=True)


for s in [0.94, 0.96, 0.98, 1.00, 1.02, 1.04, 1.07, 1.10]:
    store(f"sb_s{s:.2f}_u0.2336", sb_bulk(s, 0.2336), dict(kpts=(10, 10, 10)), f"{OUT}/gpaw_bulk.txt")
for u in [0.222, 0.228, 0.240, 0.246]:
    store(f"sb_s1.00_u{u:.3f}", sb_bulk(1.0, u), dict(kpts=(10, 10, 10)), f"{OUT}/gpaw_bulk.txt")
for s, u in [(0.97, 0.228), (1.03, 0.240)]:
    store(f"sb_s{s:.2f}_u{u:.3f}", sb_bulk(s, u), dict(kpts=(10, 10, 10)), f"{OUT}/gpaw_bulk.txt")

# --- Sb2S3 ---
traj = f"{ROOT}/local_opt/sb2s3_relaxed.traj"
if os.path.exists(traj):
    base = read(traj)
    for s in [0.96, 0.98, 1.00, 1.02, 1.05]:
        a = base.copy(); a.set_cell(base.cell * s, scale_atoms=True)
        store(f"sb2s3_s{s:.2f}", a, dict(kpts=(2, 6, 2)), f"{OUT}/gpaw_sb2s3.txt")
    rng = np.random.default_rng(7)
    for i in range(3):
        a = base.copy(); a.positions += rng.normal(0, 0.10, a.positions.shape)
        store(f"sb2s3_rattle{i}", a, dict(kpts=(2, 6, 2)), f"{OUT}/gpaw_sb2s3.txt")
else:
    print("sb2s3_relaxed.traj missing -> Sb2S3 part skipped")
print("SB_EV_DONE", len(list(db.select())), "structures")
