#!/usr/bin/env python3
"""自由原子 (弱 confinement) の殻別 Hubbard U を hotcent で計算 (解析的カーネル法と有限差分)。
出力: local_opt/shell_hubbard.json"""
import json
import os

from hotcent.atomic_dft import AtomicDFT
from hotcent.confinement import PowerConfinement

ROOT = "/Users/crocus/uhuhu/MoS2_DFTB"
XC = "GGA_X_PBE+GGA_C_PBE"
SPECS = {
    "Mo": dict(configuration="[Kr] 4d5 5s1 5p0", valence=["4d", "5s", "5p"], scalarrel=True, kw={}),
    "S": dict(configuration="[Ne] 3s2 3p4 3d0", valence=["3s", "3p", "3d"], scalarrel=False, kw={}),
    "O": dict(configuration="[He] 2s2 2p4", valence=["2s", "2p"], scalarrel=False, kw={}),
    "Sb": dict(configuration="[Kr] 4d10 5s2 5p3 5d0", valence=["5s", "5p", "5d"], scalarrel=True,
               kw=dict(mix=0.12, nodegpts=180)),
}
out = {}
for el, sp in SPECS.items():
    atom = AtomicDFT(el, xc=XC, configuration=sp["configuration"], valence=sp["valence"],
                     scalarrel=sp["scalarrel"], confinement=PowerConfinement(r0=40.0, s=4),
                     perturbative_confinement=False, txt=None, **sp["kw"])
    atom.run()
    res = {"eig": {nl: float(atom.get_eigenvalue(nl)) for nl in sp["valence"]}, "U_analytic": {}, "U_fd": {}}
    for nl in sp["valence"]:
        try:
            res["U_analytic"][nl] = float(atom.get_analytical_hubbard_value(nl))
        except Exception as ex:
            res["U_analytic"][nl] = None
            print(el, nl, "analytic failed:", str(ex)[:80], flush=True)
        try:
            res["U_fd"][nl] = float(atom.get_hubbard_value(nl, scheme="central", maxstep=1))
        except Exception as ex:
            res["U_fd"][nl] = None
            print(el, nl, "fd failed:", str(ex)[:80], flush=True)
    out[el] = res
    print(el, {k: (round(v, 4) if v else None) for k, v in res["U_analytic"].items()},
          "fd:", {k: (round(v, 4) if v else None) for k, v in res["U_fd"].items()}, flush=True)
json.dump(out, open(f"{ROOT}/local_opt/shell_hubbard.json", "w"), indent=1)
print("DONE")
