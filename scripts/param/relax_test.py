#!/usr/bin/env python3
"""反発込み SKF で単層 MoS2 を格子+内部座標ごと緩和し、LAK 平衡 (a=3.2199, 厚さ 3.115) と比較。
使い方: python relax_test.py <skf_dir> [<skf_dir> ...]"""
import os, subprocess, sys
import numpy as np
from ase.build import mx2
from ase.io import read, write as ase_write
SCRIPTS = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, SCRIPTS)
from optimize_multi import DFTB
HSD = """Geometry = GenFormat {{
  <<< "geo.gen"
}}
Driver = GeometryOptimization {{
  Optimiser = Rational {{}}
  LatticeOpt = Yes
  FixAngles = Yes
  FixLengths = {{ No No Yes }}
  Convergence {{ GradElem [eV/A] = 2e-3 }}
  MaxSteps = 200
}}
Hamiltonian = DFTB {{
  SCC = Yes
  SCCTolerance = 1e-7
  MaxSCCIterations = 250
  Mixer = Broyden {{}}
  SlaterKosterFiles = Type2FileNames {{
    Prefix = "{skf}/"
    Separator = "-"
    Suffix = ".skf"
  }}
  MaxAngularMomentum {{
    Mo = "d"
    S = "d"
  }}
  Filling = Fermi {{ Temperature [K] = 100 }}
  KPointsAndWeights = SupercellFolding {{
    12 0 0
    0 12 0
    0 0 1
    0.5 0.5 0.0
  }}
}}
Options {{ WriteResultsTag = Yes }}
ParserOptions {{ ParserVersion = 14 }}
"""
HSD = HSD.replace(
    "  PolynomialRepulsive", (os.environ.get("DFTB_EXTRA_HSD", "").replace("{", "{{").replace("}", "}}") + "\n"
                              if os.environ.get("DFTB_EXTRA_HSD") else "") + "  PolynomialRepulsive", 1)
W = "/private/tmp/claude-501/-Users-crocus-uhuhu-MoS2-DFTB/6416faee-14e1-48c8-9838-35c06f3e7813/scratchpad/relax"
for skf in sys.argv[1:]:
    wd = os.path.join(W, os.path.basename(skf.rstrip("/").replace("/", "_")))
    os.makedirs(wd, exist_ok=True)
    atoms = mx2("MoS2", kind="2H", a=3.16, thickness=3.127, vacuum=10.0)
    ase_write(os.path.join(wd, "geo.gen"), atoms, format="gen")
    open(os.path.join(wd, "dftb_in.hsd"), "w").write(HSD.format(skf=os.path.abspath(skf)))
    r = subprocess.run([DFTB], cwd=wd, capture_output=True, text=True, timeout=1800, env=dict(os.environ, OMP_NUM_THREADS="4"))
    if r.returncode != 0 or not os.path.exists(os.path.join(wd, "geo_end.gen")):
        print(skf, "FAILED:", r.stdout[-300:].replace("\n", " | ")); continue
    out = read(os.path.join(wd, "geo_end.gen"))
    a = out.cell.lengths()[0]; z = out.positions[:, 2]; t = z.max() - z.min()
    print(f"{skf}: a = {a:.4f} A (LAK 3.2199, {100*(a/3.2199-1):+.2f}%), S-S thickness = {t:.4f} A (LAK 3.115, {100*(t/3.115-1):+.2f}%)")
