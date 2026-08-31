#!/usr/bin/env python3
"""LAK 参照計算 (OUTCAR) を集めて ASE DB (dft.db) を構築。

使い方:
  python3 build_dft_db.py <out_db> <dir1> [dir2 ...]
各 dir 内の */OUTCAR を走査 (最終イオンステップの構造+E+F)。
"""
import glob
import os
import sys

from ase.db import connect
from ase.io import read


def main():
    out_db = sys.argv[1]
    roots = sys.argv[2:]
    if os.path.exists(out_db):
        os.remove(out_db)
    db = connect(out_db)
    nok = nfail = 0
    for root in roots:
        for outcar in sorted(glob.glob(os.path.join(root, "*", "OUTCAR"))):
            name = os.path.basename(os.path.dirname(outcar))
            try:
                atoms = read(outcar, index=-1)
                e = atoms.get_potential_energy()
                atoms.get_forces()
                db.write(atoms, name=name)
                print(f"{name}: E = {e:.6f} eV  nat={len(atoms)}")
                nok += 1
            except Exception as ex:
                print(f"{name}: SKIP ({ex})")
                nfail += 1
    print(f"done: {nok} entries, {nfail} skipped -> {out_db}")


if __name__ == "__main__":
    main()
