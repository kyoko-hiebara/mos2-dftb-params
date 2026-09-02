#!/usr/bin/env python3
"""study の最良試行の Mo_UJ から DFTB+ の OrbitalPotential ブロックを出力 (無ければ空)。使い方: python u_block_from_study.py optm4,optm4c"""
import os, sys
import optuna
optuna.logging.set_verbosity(optuna.logging.WARNING)
ROOT = "/Users/crocus/uhuhu/MoS2_DFTB"
best = None
for name in sys.argv[1].split(","):
    p = f"{ROOT}/local_opt/{name}/study.db"
    if not os.path.exists(p):
        continue
    st = optuna.load_study(study_name=name, storage=f"sqlite:///{p}")
    try:
        t = st.best_trial
    except ValueError:
        continue
    if best is None or t.value < best.value:
        best = t
uj = best.params.get("Mo_UJ", 0.0) if best else 0.0
if uj > 1e-4:
    print(f"  OrbitalPotential = {{\n    Functional = FLL\n    Mo = {{\n      Shells = {{3}}\n      UJ = {uj:.5f}\n    }}\n  }}")
