#!/usr/bin/env python3
"""Sb 反発込み SKF の緩和テスト: バルク Sb (格子+内部 u)、Sb2S3 (固定セル内部座標) を PBE と比較。
使い方: python relax_sb_test.py <skf_dir>"""
import os, subprocess, sys
import numpy as np
from ase.io import read, write as ase_write
from ase.spacegroup import crystal
SCRIPTS = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, SCRIPTS)
from optimize_multi import DFTB
HSD = """Geometry = GenFormat {{
  <<< "geo.gen"
}}
Driver = GeometryOptimization {{
  Optimiser = Rational {{}}
{latblock}  Convergence {{ GradElem [eV/A] = 3e-3 }}
  MaxSteps = 300
}}
Hamiltonian = DFTB {{
  SCC = Yes
  SCCTolerance = 1e-6
  MaxSCCIterations = 300
  Mixer = Broyden {{}}
  SlaterKosterFiles = Type2FileNames {{
    Prefix = "{skf}/"
    Separator = "-"
    Suffix = ".skf"
  }}
  MaxAngularMomentum {{
{lmax}
  }}
  Filling = Fermi {{ Temperature [K] = 100 }}
  KPointsAndWeights = SupercellFolding {{
    {k0} 0 0
    0 {k1} 0
    0 0 {k2}
    0.5 0.5 0.5
  }}
}}
Options {{ WriteResultsTag = Yes }}
ParserOptions {{ ParserVersion = 14 }}
"""
W = "/private/tmp/claude-501/-Users-crocus-uhuhu-MoS2-DFTB/6416faee-14e1-48c8-9838-35c06f3e7813/scratchpad/relax_sb"
skf = os.path.abspath(sys.argv[1])
def run(name, atoms, latopt, k):
    wd = os.path.join(W, name); os.makedirs(wd, exist_ok=True)
    ase_write(os.path.join(wd, "geo.gen"), atoms, format="gen")
    lmax = "\n".join(f'    {e} = "d"' for e in sorted(set(atoms.get_chemical_symbols())))
    latblock = "  LatticeOpt = Yes\n  FixAngles = Yes\n" if latopt == "Yes" else ""
    open(os.path.join(wd, "dftb_in.hsd"), "w").write(HSD.format(skf=skf, latblock=latblock, lmax=lmax, k0=k[0], k1=k[1], k2=k[2]))
    r = subprocess.run([DFTB], cwd=wd, capture_output=True, text=True, timeout=3600, env=dict(os.environ, OMP_NUM_THREADS="4"))
    if not os.path.exists(os.path.join(wd, "geo_end.gen")):
        print(name, "FAILED", r.stdout[-200:].replace("\n", " | ")); return None
    return read(os.path.join(wd, "geo_end.gen"))
# --- bulk Sb: 実験構造から出発 (PBE の等方平衡スケール 1.023) ---
sb = crystal("Sb", [(0, 0, 0.2336)], spacegroup=166, cellpar=[4.3084, 4.3084, 11.274, 90, 90, 120], primitive_cell=True)
out = run("sb_bulk", sb, "Yes", (10, 10, 10))
if out is not None:
    v0, v1 = sb.get_volume(), out.get_volume()
    from ase.neighborlist import neighbor_list
    i, j, d = neighbor_list("ijd", out, 3.8); d1 = np.sort(d[i == 0])
    cp = out.cell.cellpar()
    print(f"bulk Sb: V/V0 = {v1/v0:.4f} (PBE isotropic scan -> {1.023**3:.4f}); a_rh={cp[0]:.3f} alpha={cp[3]:.2f} (expt 4.507, 57.1); "
          f"Sb-Sb shells {d1[:3].round(3)} / {d1[3:6].round(3)} (expt 2.908 x3 / 3.355 x3)")
# --- Sb2S3: PBE 緩和構造から出発、固定セル ---
s23 = read(f"{os.path.dirname(SCRIPTS)}/local_opt/sb2s3_relaxed.traj")
out = run("sb2s3", s23, "No", (2, 6, 2))
if out is not None:
    disp = np.linalg.norm(out.positions - s23.positions, axis=1)
    print(f"Sb2S3 (fixed cell): RMS displacement from PBE positions = {np.sqrt((disp**2).mean()):.3f} A, max {disp.max():.3f} A")
