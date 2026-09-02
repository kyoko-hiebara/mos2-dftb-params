#!/usr/bin/env python3
"""反発ターゲット D(a) = E_LAK(a) − E_elec(a) の形状指標 (単調減少・凸性) を電子部 SKF について計算。
使い方: python erep_shape.py <skf_dir> [<skf_dir> ...]   (モジュールとしても利用: shape_terms(skf))"""
import os, subprocess, sys
import numpy as np
from ase.db import connect
from ase.io import write as ase_write
SCRIPTS = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, SCRIPTS)
from dftb_energy_local import HSD, parse_results_tag, DFTB
HA = 27.211386245988
DB = f"{os.path.dirname(SCRIPTS)}/local_opt/repfit3/pairset3_dft.db"
_ROWS = None
def rows():
    global _ROWS
    if _ROWS is None:
        _ROWS = sorted([(float(r.name[2:]), r.toatoms(), r.energy) for r in connect(DB).select() if r.name.startswith("a_")])
    return _ROWS
def e_elec(atoms, skf, wd, extra=""):
    os.makedirs(wd, exist_ok=True)
    ase_write(os.path.join(wd, "geo.gen"), atoms, format="gen")
    hsd = HSD.format(skf=skf, kblock="  KPointsAndWeights = SupercellFolding {\n    12 0 0\n    0 12 0\n    0 0 1\n    0.5 0.5 0.0\n  }")
    if extra:
        hsd = hsd.replace("  PolynomialRepulsive", extra + "\n  PolynomialRepulsive")
    open(os.path.join(wd, "dftb_in.hsd"), "w").write(hsd)
    r = subprocess.run([DFTB], cwd=wd, capture_output=True, text=True, env=dict(os.environ, OMP_NUM_THREADS="4"))
    tag = os.path.join(wd, "results.tag")
    if not os.path.exists(tag):
        return None
    e, _ = parse_results_tag(tag); return e * HA
def shape_terms(skf, wd, extra="", sel=(2.95, 3.04, 3.10, 3.16, 3.22, 3.26, 3.32, 3.40)):
    a, D = [], []
    for av, atoms, e_dft in rows():
        if sel and not any(abs(av - s) < 1e-3 for s in sel):
            continue
        ee = e_elec(atoms, os.path.abspath(skf), os.path.join(wd, f"a_{av:.2f}"), extra)
        if ee is None:
            return None
        a.append(av); D.append(e_dft - ee)
    a, D = np.array(a), np.array(D)
    dD = np.diff(D)                      # a 増加方向: 単調減少なら全て負
    mono = float((np.maximum(0.0, dD) ** 2).sum())
    d2 = D[2:] - 2 * D[1:-1] + D[:-2]    # 凸なら正
    conv = float((np.maximum(0.0, -d2) ** 2).sum())
    return dict(a=a.tolist(), D=(D - D[-1]).tolist(), mono=mono, conv=conv, span=float(D.max() - D.min()))
if __name__ == "__main__":
    W = "/private/tmp/claude-501/-Users-crocus-uhuhu-MoS2-DFTB/6416faee-14e1-48c8-9838-35c06f3e7813/scratchpad/erep_shape"
    for skf in sys.argv[1:]:
        t = shape_terms(skf, os.path.join(W, os.path.basename(skf.rstrip("/"))), sel=())
        print(f"{skf}: mono={t['mono']:.4f} conv={t['conv']:.4f} span={t['span']:.2f} eV")
        print("   a:", " ".join(f"{x:.2f}" for x in t["a"]))
        print("   D(a)-D(end):", " ".join(f"{x:+.3f}" for x in t["D"]))
