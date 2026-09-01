#!/usr/bin/env python3
"""skf_v3sb: Mo/S(v3) + O(v4o) + Sb(optsb2) の統合セットを密グリッドで生成し、
バルク Sb バンドを最終検証 (GPAW vs DFTB vs PTBP)。"""
import json
import os
import shutil
import subprocess
import sys

import numpy as np

SCRIPTS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPTS)

import optimize_confinement as oc  # noqa: E402
import optimize_sb2 as sb2  # noqa: E402  (SPECS["Sb"] spd, make_atom tuned, FREE["Sb"])

ROOT = "/Users/crocus/uhuhu/MoS2_DFTB"
OUT = f"{ROOT}/local_opt/skf_v3sb"

oc.SPECS["O"] = dict(configuration="[He] 2s2 2p4", valence=["2s", "2p"],
                     scalarrel=False, occupations={"2s": 2, "2p": 4})

# 密グリッド
oc.RMIN, oc.DR, oc.NPTS = 0.4, 0.02, 980
oc.NTHETA, oc.NR = 150, 50


def main():
    import optuna
    optuna.logging.set_verbosity(optuna.logging.WARNING)

    # --- パラメータ収集 ---
    st1 = optuna.load_study(study_name="optm1",
                            storage=f"sqlite:///{ROOT}/local_opt/optm1/study.db")
    fp = dict(st1.best_trial.params)
    pMo = {"r_dens": fp["Mo_rd"], "s_wf": fp["Mo_swf"],
           "r_wf": {"4d": fp["Mo_r4d"], "5s": fp["Mo_r5s"], "5p": fp["Mo_r5p"]}}
    pS = {"r_dens": fp["S_rd"], "s_wf": fp["S_swf"],
          "r_wf": {"3s": fp["S_r3s"], "3p": fp["S_r3p"], "3d": fp["S_r3d"]}}
    e3d = fp["S_e3d"]
    sto = optuna.load_study(study_name="opto2",
                            storage=f"sqlite:///{ROOT}/local_opt/opto2/study.db")
    fo = dict(sto.best_trial.params)
    pO = {"r_dens": fo["O_rd"], "s_wf": fo["O_swf"],
          "r_wf": {"2s": fo["O_r2s"], "2p": fo["O_r2p"]}}
    ssb = optuna.load_study(study_name="optsb2",
                            storage=f"sqlite:///{ROOT}/local_opt/optsb2/study.db")
    fs = dict(ssb.best_trial.params)
    pSb = {"r_dens": fs["Sb_rd"], "s_wf": fs["Sb_swf"],
           "r_wf": {"5s": fs["Sb_r5s"], "5p": fs["Sb_r5p"],
                    "5d": fs["Sb_r5d"]}}
    shifts = {"Mo": {"4d": fp["Mo_sh4d"], "5s": fp["Mo_sh5s"],
                     "5p": fp["Mo_sh5p"]},
              "S": {"3s": fp["S_sh3s"], "3p": fp["S_sh3p"]},
              "O": {"2s": fo["O_sh2s"], "2p": fo["O_sh2p"]},
              "Sb": {"5s": fs["Sb_sh5s"], "5p": fs["Sb_sh5p"]}}
    oc.FREE["Sb"]["eig"]["5d"] = fs["Sb_e5d"]
    oc.FREE["Sb"]["U"]["5d"] = oc.FREE["Sb"]["U"]["5p"]

    # --- 既存 (v3 Mo/S + v4o O 系) をコピー ---
    os.makedirs(OUT, exist_ok=True)
    for f in os.listdir(f"{ROOT}/local_opt/skf_v4o"):
        if f.endswith(".skf"):
            shutil.copy(f"{ROOT}/local_opt/skf_v4o/{f}", f"{OUT}/{f}")

    # --- Sb ペアを生成 ---
    for task in [("Sb", "Sb", pSb, pSb, OUT, e3d, shifts),
                 ("Mo", "Sb", pMo, pSb, OUT, e3d, shifts),
                 ("S", "Sb", pS, pSb, OUT, e3d, shifts),
                 ("O", "Sb", pO, pSb, OUT, e3d, shifts)]:
        print("pair", oc.gen_pair(task), "done", flush=True)
    print("skf_v3sb:", sorted(os.listdir(OUT)))

    # --- バルク Sb 検証 (密グリッド) ---
    e = sb2.sb_bands(OUT, "/tmp/v3sb_check")
    rms_all, rms_ef = sb2.band_loss_sb(e)
    print(f"\nfinal Sb bulk: rms_EF = {rms_ef:.3f} eV, rms_all = {rms_all:.3f} eV")
    json.dump({"params_sb": fs, "rms_ef": rms_ef, "rms_all": rms_all},
              open(f"{OUT}/sb_metrics.json", "w"), indent=1)
    np.save(f"{ROOT}/local_opt/sb_dftb_bands.npy", e)


if __name__ == "__main__":
    main()
