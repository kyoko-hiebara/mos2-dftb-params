#!/usr/bin/env python3
"""Sb 反発フィット用参照の追加: 異方的ひずみ (a のみ / c のみ ±3%) x u。既存 DB に追記。"""
import os
import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator
from ase.db import connect
from ase.spacegroup import crystal
from gpaw import GPAW, PW, FermiDirac
ROOT = "/Users/crocus/uhuhu/MoS2_DFTB"
os.environ.setdefault("GPAW_SETUP_PATH", f"{ROOT}/sw_local/gpaw-setups/gpaw-setups-24.11.0")
db = connect(f"{ROOT}/local_opt/repfit_sb/sb_ref_dft.db")
done = {r.name for r in db.select()}
def sb(sa, sc, u):
    return crystal("Sb", [(0, 0, u)], spacegroup=166, cellpar=[4.3084 * sa, 4.3084 * sa, 11.274 * sc, 90, 90, 120], primitive_cell=True)
for sa, sc, u in [(0.97, 1.0, 0.2336), (1.03, 1.0, 0.2336), (1.0, 0.97, 0.2336), (1.0, 1.03, 0.2336),
                  (0.97, 1.03, 0.240), (1.03, 0.97, 0.228), (1.02, 1.05, 0.2336), (0.98, 0.95, 0.2336)]:
    name = f"sb_a{sa:.2f}_c{sc:.2f}_u{u:.4f}"
    if name in done:
        continue
    a = sb(sa, sc, u)
    a.calc = GPAW(mode=PW(350), xc="PBE", kpts=(10, 10, 10), occupations=FermiDirac(0.05), txt=f"{ROOT}/local_opt/repfit_sb/gpaw_bulk2.txt")
    e = a.get_potential_energy(); f = a.get_forces()
    a.calc = SinglePointCalculator(a, energy=e, forces=f); db.write(a, name=name)
    print(f"{name}: E={e:.4f} maxF={np.abs(f).max():.3f}", flush=True)
print("SB_EV2_DONE", len(list(db.select())))
