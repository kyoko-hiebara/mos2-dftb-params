#!/usr/bin/env python3
"""LAK 参照バンド構造計算 (0-weight k 点方式) と JSON 抽出。

使い方:
  python3 run_bands_vasp.py <a_lattice> <workdir> [--zS <z_S offset>] [--np 24]

手順:
  1) 構造生成 (a, 内部座標は E(a) 緩和済みの値を CONTCAR から流用するのが理想だが、
     ここでは thickness を引数で受ける)
  2) NELM=1 のダミー実行で IBZKPT を取得
  3) IBZKPT + バンドパス (重み 0) を連結した KPOINTS で LAK SCF
  4) EIGENVAL から重み 0 の k 点の固有値を抽出 -> bands_lak.json
"""
import argparse
import json
import os
import shutil
import subprocess
import sys

import numpy as np
from ase.build import mx2
from ase.io import write as ase_write

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bandpath_common import kpath_points  # noqa: E402

VASP = "/root/vasp.6.4.2/bin/vasp_std"
POTCAR = "/workspace/MoS2_DFTB/ref_calc/POTCAR_Mo_S"

INCAR_SCF = """SYSTEM = mono MoS2 LAK bands
PREC = Accurate
ENCUT = 900
EDIFF = 1E-7
NELM = 200
ISMEAR = 0
SIGMA = 0.03
METAGGA = LIBXC
LIBXC1 = MGGA_X_LAK
LIBXC2 = MGGA_C_LAK
LASPH = .TRUE.
LREAL = .FALSE.
LWAVE = .FALSE.
LCHARG = .FALSE.
LMAXMIX = 6
NBANDS = 28
NCORE = 4
"""


def run_vasp(workdir, np_ranks):
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    with open(os.path.join(workdir, "stdout.log"), "w") as f:
        subprocess.run(
            ["mpirun", "--allow-run-as-root", "-np", str(np_ranks), VASP],
            cwd=workdir, stdout=f, stderr=subprocess.STDOUT, timeout=7200,
            stdin=subprocess.DEVNULL, env=env)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("a", type=float)
    ap.add_argument("workdir")
    ap.add_argument("--thickness", type=float, default=3.13)
    ap.add_argument("--np", type=int, default=24, dest="np_ranks")
    args = ap.parse_args()

    wd = args.workdir
    os.makedirs(wd, exist_ok=True)
    atoms = mx2(formula="MoS2", kind="2H", a=args.a,
                thickness=args.thickness, vacuum=10.0)
    ase_write(os.path.join(wd, "POSCAR"), atoms, format="vasp", direct=True,
              sort=True)
    shutil.copy(POTCAR, os.path.join(wd, "POTCAR"))

    # --- step 1: get IBZKPT with a cheap PBE 1-step run
    with open(os.path.join(wd, "INCAR"), "w") as f:
        f.write("SYSTEM = ibz\nPREC = Accurate\nENCUT = 400\nNELM = 1\n"
                "EDIFF = 1E-4\nISMEAR = 0\nSIGMA = 0.03\nGGA = PE\n"
                "LWAVE = .FALSE.\nLCHARG = .FALSE.\nNCORE = 4\n")
    with open(os.path.join(wd, "KPOINTS"), "w") as f:
        f.write("mono\n0\nGamma\n12 12 1\n0 0 0\n")
    run_vasp(wd, args.np_ranks)
    ibz = open(os.path.join(wd, "IBZKPT")).readlines()
    nibz = int(ibz[1].split()[0])
    ibz_pts = [ibz[3 + i] for i in range(nibz)]

    # --- step 2: combined KPOINTS (IBZ weighted + path with weight 0)
    path = kpath_points()
    with open(os.path.join(wd, "KPOINTS"), "w") as f:
        f.write("scf + zero-weight band path\n")
        f.write(f"{nibz + len(path)}\n")
        f.write("Reciprocal\n")
        for line in ibz_pts:
            f.write(line if line.endswith("\n") else line + "\n")
        for p in path:
            f.write(f"{p[0]:.10f} {p[1]:.10f} {p[2]:.10f} 0.000000\n")

    with open(os.path.join(wd, "INCAR"), "w") as f:
        f.write(INCAR_SCF)
    run_vasp(wd, args.np_ranks)

    # --- step 3: parse EIGENVAL, keep only zero-weight kpts
    with open(os.path.join(wd, "EIGENVAL")) as f:
        lines = f.readlines()
    nelec, nk, nb = (int(x) for x in lines[5].split())
    eigs, kw, kpts = [], [], []
    idx = 7
    for ik in range(nk):
        kline = lines[idx].split()
        kpts.append([float(kline[0]), float(kline[1]), float(kline[2])])
        kw.append(float(kline[3]))
        band_e = []
        for ib in range(nb):
            band_e.append(float(lines[idx + 1 + ib].split()[1]))
        eigs.append(band_e)
        idx += nb + 2
    eigs = np.array(eigs)
    kw = np.array(kw)
    sel = kw < 1e-10
    eigs_path = eigs[sel]
    kpts_path = np.array(kpts)[sel]

    nocc = nelec // 2
    vbm = float(eigs_path[:, :nocc].max())
    cbm = float(eigs_path[:, nocc:].min())
    out = {
        "a": args.a,
        "kpts": kpts_path.tolist(),
        "eigs": eigs_path.tolist(),
        "nelec": nelec,
        "nocc": nocc,
        "vbm": vbm,
        "cbm": cbm,
        "gap": cbm - vbm,
    }
    jpath = os.path.join(wd, "bands_lak.json")
    json.dump(out, open(jpath, "w"))
    print(f"nk_path={len(eigs_path)} nb={nb} nelec={nelec} "
          f"VBM={vbm:.4f} CBM={cbm:.4f} gap={cbm - vbm:.4f} eV")
    print(f"wrote {jpath}")


if __name__ == "__main__":
    main()
