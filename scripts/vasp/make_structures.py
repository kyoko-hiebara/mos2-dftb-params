#!/usr/bin/env python3
"""MoS2 参照構造の生成 (DFTB パラメータ化用 VASP 参照計算の入力)。

生成物:
  ref_calc/mono/POSCAR      1H-MoS2 単層 (真空 20 A)
  ref_calc/bulk/POSCAR      2H-MoS2 バルク
  ref_calc/molecules/*/POSCAR  二原子/小分子 (ボックス 15 A)
"""
import os
from ase.build import mx2
from ase.spacegroup import crystal
from ase import Atoms
from ase.io import write

ROOT = "/workspace/MoS2_DFTB/ref_calc"

A0 = 3.160          # 面内格子定数の初期値 (実験値)
C0 = 12.295         # 2H バルク c 軸 (実験値)
THICK = 3.127       # S-S 層厚の初期値


def save(atoms, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    write(path, atoms, format="vasp", direct=True, sort=True)
    print(f"wrote {path}: {atoms.get_chemical_formula()}")


# --- 1H-MoS2 単層 ---
mono = mx2(formula="MoS2", kind="2H", a=A0, thickness=THICK, vacuum=10.0)
save(mono, f"{ROOT}/mono/POSCAR")

# --- 2H-MoS2 バルク (P6_3/mmc, #194; Mo 2c, S 4f z=0.621) ---
bulk = crystal(
    ("Mo", "S"),
    basis=[(1 / 3, 2 / 3, 1 / 4), (1 / 3, 2 / 3, 0.621)],
    spacegroup=194,
    cellpar=[A0, A0, C0, 90, 90, 120],
)
save(bulk, f"{ROOT}/bulk/POSCAR")

# --- 小分子 (スピン多重度は INCAR 側で設定) ---
mol_specs = {
    "S2": Atoms("S2", positions=[(0, 0, 0), (0, 0, 1.89)]),
    "O2": Atoms("O2", positions=[(0, 0, 0), (0, 0, 1.21)]),
    "SO": Atoms("SO", positions=[(0, 0, 0), (0, 0, 1.48)]),
    "H2": Atoms("H2", positions=[(0, 0, 0), (0, 0, 0.74)]),
    "H2S": Atoms(
        "SH2",
        positions=[(0, 0, 0), (0.96, 0, 0.93), (-0.96, 0, 0.93)],
    ),
    "H2O": Atoms(
        "OH2",
        positions=[(0, 0, 0), (0.76, 0, 0.59), (-0.76, 0, 0.59)],
    ),
}
for name, atoms in mol_specs.items():
    atoms.set_cell([15.0, 15.5, 16.0])  # 対称性を落とした直方体
    atoms.center()
    save(atoms, f"{ROOT}/molecules/{name}/POSCAR")

print("done")
