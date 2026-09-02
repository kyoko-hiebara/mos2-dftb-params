#!/usr/bin/env python3
"""Sb(111) 2BL / 4x4 MoS2 ヘテロ構造の DFTB+ チェック: E_F と MoS2 バンド端 (Mo d PDOS) の相対位置、Mulliken 電荷移動。
使い方: python dftb_interface_check.py <skf_dir> [--soc] [--dz 3.0]"""
import argparse
import os
import subprocess
import sys

import numpy as np
from ase.build import mx2
from ase.io import write as ase_write
from ase.spacegroup import crystal

SCRIPTS = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, SCRIPTS)
from optimize_multi import DFTB, HA, parse_tag_eigs  # noqa: E402


def build(dz):
    ml = mx2("MoS2", kind="2H", a=3.16, thickness=3.127, vacuum=0).repeat((4, 4, 1))
    conv = crystal("Sb", [(0, 0, 0.2336)], spacegroup=166,
                   cellpar=[4.3084, 4.3084, 11.274, 90, 90, 120], primitive_cell=False)
    sb = conv.repeat((3, 3, 1))
    # 2 バイレイヤー (z 昇順 4 層) を取り出す
    z = sb.positions[:, 2]; zs = np.unique(np.round(z, 3))
    keep = np.isin(np.round(z, 3), zs[:4]); sb = sb[keep]
    # MoS2 の面内格子に合わせて Sb を歪ませる (2.2%)
    cell = ml.cell.copy()
    sb.set_cell([cell[0], cell[1], sb.cell[2]], scale_atoms=True)
    top_s = ml.positions[:, 2].max()
    sb.positions[:, 2] += top_s + dz - sb.positions[:, 2].min()
    het = ml + sb
    zmin, zmax = het.positions[:, 2].min(), het.positions[:, 2].max()
    cell[2] = [0, 0, zmax - zmin + 20.0]
    het.set_cell(cell); het.positions[:, 2] += 10.0 - zmin; het.pbc = True
    return het, len(ml)


HSD = """Geometry = GenFormat {{
  <<< "geo.gen"
}}
Hamiltonian = DFTB {{
  SCC = Yes
  SCCTolerance = 1e-5
  MaxSCCIterations = 400
  Mixer = Broyden {{ MixingParameter = 0.1 }}
  SlaterKosterFiles = Type2FileNames {{
    Prefix = "{skf}/"
    Separator = "-"
    Suffix = ".skf"
  }}
  MaxAngularMomentum {{
    Mo = "d"
    S = "d"
    Sb = "d"
  }}
  PolynomialRepulsive = SetForAll {{ Yes }}
  Filling = Fermi {{ Temperature [K] = 300 }}
{soc}
  KPointsAndWeights = SupercellFolding {{
    4 0 0
    0 4 0
    0 0 1
    0.5 0.5 0.0
  }}
}}
Options {{ WriteResultsTag = Yes }}
Analysis {{
  ProjectStates {{
    Region {{
      Atoms = 1:{nml}
      Label = "mos2"
    }}
    Region {{
      Atoms = {nsb0}:{ntot}
      Label = "sb"
    }}
  }}
}}
ParserOptions {{ ParserVersion = 14 }}
"""
SOC = """  SpinOrbit = {
    Dual = Yes
    Mo [eV] = {0.0 0.036 0.0931}
    S  [eV] = {0.0 0.055 0.0}
    Sb [eV] = {0.0 0.57 0.0}
  }"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("skf"); ap.add_argument("--soc", action="store_true"); ap.add_argument("--dz", type=float, default=3.0)
    ap.add_argument("--wd", default=None)
    args = ap.parse_args()
    het, nml = build(args.dz)
    wd = args.wd or f"/private/tmp/claude-501/-Users-crocus-uhuhu-MoS2-DFTB/6416faee-14e1-48c8-9838-35c06f3e7813/scratchpad/iface_{os.path.basename(args.skf)}"
    os.makedirs(wd, exist_ok=True)
    ase_write(os.path.join(wd, "geo.gen"), het, format="gen")
    ase_write(os.path.join(wd, "het.traj"), het)
    open(os.path.join(wd, "dftb_in.hsd"), "w").write(HSD.format(skf=os.path.abspath(args.skf), nml=nml, nsb0=nml + 1, ntot=len(het),
                                                                 soc=SOC if args.soc else ""))
    r = subprocess.run([DFTB], cwd=wd, capture_output=True, text=True, timeout=3600, env=dict(os.environ, OMP_NUM_THREADS="8"))
    if r.returncode != 0:
        print("DFTB FAILED", r.stdout[-500:]); sys.exit(1)
    lines = open(os.path.join(wd, "results.tag")).readlines()
    ef = None
    for i, ln in enumerate(lines):
        if ln.startswith("fermi_level"):
            ef = float(lines[i + 1].split()[0]) * HA
    # Mulliken: detailed.out から原子電荷
    q = []
    det = open(os.path.join(wd, "detailed.out")).read().splitlines()
    for i, ln in enumerate(det):
        if ln.strip().startswith("Atom") and "Net charge" in ln or "Net atomic charges" in ln:
            j = i + 2
            while j < len(det) and det[j].strip():
                t = det[j].split()
                if len(t) >= 2:
                    try: q.append(float(t[-1]))
                    except ValueError: pass
                j += 1
            break
    q = np.array(q)
    def parse_proj(fn):
        E, Wt = [], []
        for ln in open(fn):
            t = ln.split()
            if len(t) == 2:
                try:
                    E.append(float(t[0])); Wt.append(float(t[1]))
                except ValueError:
                    pass
        return np.array(E), np.array(Wt)
    print(f"E_F = {ef:.3f} eV  (natoms {len(het)}: MoS2 {nml}, Sb {len(het)-nml}, dz={args.dz})")
    E, Wt = parse_proj(os.path.join(wd, "mos2.out"))
    e = np.sort(E[Wt > 0.6])
    win = e[(e > ef - 3) & (e < ef + 3)]
    gaps = np.diff(win); i = int(np.argmax(gaps))
    print(f"MoS2-projected gap: VBM {win[i]-ef:+.3f} -> CBM {win[i+1]-ef:+.3f} eV rel. E_F (gap {gaps[i]:.3f} eV); "
          f"E_F - CBM(MoS2) = {ef-win[i+1]:+.3f} eV")
    det = open(os.path.join(wd, "detailed.out")).read().splitlines()
    idx = [k for k, ln in enumerate(det) if "Atomic gross charges" in ln]
    if idx:
        j = idx[-1] + 2; q = []
        while j < len(det) and det[j].strip():
            try: q.append(float(det[j].split()[-1]))
            except ValueError: pass
            j += 1
        q = np.array(q)
        if len(q) == len(het):
            print(f"Mulliken net charge: MoS2 {q[:nml].sum():+.3f} e, Sb {q[nml:].sum():+.3f} e per cell")


if __name__ == "__main__":
    main()
