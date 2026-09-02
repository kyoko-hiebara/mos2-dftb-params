#!/usr/bin/env python3
"""O 系 confinement/onsite の最適化 (mac ローカル)。

ターゲット (sub_O 5x5 の LAK 参照):
  - ギャップ 1.826 eV + ギャップ内無準位の維持
  - O 2s 孤立準位: VBM - 20.122 eV
  - VB スペクトル (VBM..-8 eV, 全 k ソート) の RMS
変数 6: O r_dens, r_wf(2s), r_wf(2p), s_wf, onsite shift(2s), shift(2p)
Mo/S は skf_v3 (optm1 ベスト) に固定。

使い方: python3 optimize_o.py <n_trials> [--tag opto1]
"""
import argparse
import json
import os
import shutil
import subprocess
import sys

import numpy as np

SCRIPTS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPTS)

import optimize_confinement as oc  # noqa: E402
from optimize_multi import DFTB, HA, parse_tag_eigs  # noqa: E402

from ase.io import read  # noqa: E402
from ase.io import write as ase_write  # noqa: E402

ROOT = "/Users/crocus/uhuhu/MoS2_DFTB"
SUBO_POSCAR = f"{ROOT}/results/ref_calc/defects/sub_O_gpu/CONTCAR"

oc.SPECS["O"] = dict(configuration="[He] 2s2 2p4", valence=["2s", "2p"],
                     scalarrel=False, occupations={"2s": 2, "2p": 4})

# ---- LAK 参照 (sub_O EIGENVAL) ----
def lak_reference():
    lines = open(f"{ROOT}/results/ref_calc/defects/sub_O_lak/EIGENVAL").readlines()
    nelec, nk, nb = (int(x) for x in lines[5].split())
    idx = 7
    eigs = []
    for ik in range(nk):
        eigs.append([float(lines[idx + 1 + ib].split()[1]) for ib in range(nb)])
        idx += nb + 2
    e = np.array(eigs)
    nocc = nelec // 2
    vbm = e[:, :nocc].max()
    cbm = e[:, nocc:].min()
    o2s = e[:, 100].mean() - vbm          # band 101 (0-based 100)
    # VB 窓: バンドごとの k 平均 (Mo セミコア除外)、-8 eV より浅いもの
    bm = e[:, 100:nocc].mean(axis=0) - vbm
    win = np.sort(bm[bm > -8.0])
    return dict(gap=float(cbm - vbm), o2s=float(o2s), vb=win)


REF = lak_reference()
print(f"LAK ref: gap={REF['gap']:.4f}, O2s={REF['o2s']:.3f}, "
      f"VB window n={len(REF['vb'])}", flush=True)

HSD = """Geometry = GenFormat {{
  <<< "geo.gen"
}}
Hamiltonian = DFTB {{
  SCC = Yes
  SCCTolerance = 1e-5
  MaxSCCIterations = 200
  Mixer = Broyden {{}}
  SlaterKosterFiles = Type2FileNames {{
    Prefix = "{skf}/"
    Separator = "-"
    Suffix = ".skf"
  }}
  MaxAngularMomentum {{
    Mo = "d"
    S = "d"
    O = "p"
  }}
  PolynomialRepulsive = SetForAll {{ Yes }}
  Filling = Fermi {{ Temperature [K] = 50 }}
  KPointsAndWeights = SupercellFolding {{
    2 0 0
    0 2 0
    0 0 1
    0.0 0.0 0.0
  }}
}}
Options {{ WriteResultsTag = Yes }}
ParserOptions {{ ParserVersion = 14 }}
"""
HSD = HSD.replace(
    "  PolynomialRepulsive", (os.environ.get("DFTB_EXTRA_HSD", "").replace("{", "{{").replace("}", "}}") + "\n"
                              if os.environ.get("DFTB_EXTRA_HSD") else "") + "  PolynomialRepulsive", 1)


def v3_params():
    import optuna
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    base = os.environ.get("O_BASE_STUDY", "optm1")
    st = optuna.load_study(study_name=base,
                           storage=f"sqlite:///{ROOT}/local_opt/{base}/study.db")
    fp = dict(st.best_trial.params)
    print(f"O base: {base} best trial {st.best_trial.number}", flush=True)
    pMo = {"r_dens": fp["Mo_rd"], "s_wf": fp["Mo_swf"],
           "r_wf": {"4d": fp["Mo_r4d"], "5s": fp["Mo_r5s"],
                    "5p": fp["Mo_r5p"]}}
    if "Mo_s4d" in fp:
        pMo["s_wf_shell"] = {"4d": fp["Mo_s4d"]}
    pS = {"r_dens": fp["S_rd"], "s_wf": fp["S_swf"],
          "r_wf": {"3s": fp["S_r3s"], "3p": fp["S_r3p"], "3d": fp["S_r3d"]}}
    return pMo, pS, fp["S_e3d"]


PMO, PS, E3D = v3_params()


def evaluate(po, oshifts, wd):
    skf = os.path.join(wd, "skf")
    os.makedirs(skf, exist_ok=True)
    for f in ["Mo-Mo.skf", "Mo-S.skf", "S-Mo.skf", "S-S.skf"]:
        shutil.copy(f"{ROOT}/local_opt/{os.environ.get('O_BASE_SKF', 'skf_v3')}/{f}", f"{skf}/{f}")
    shifts = {"O": oshifts}
    for task in [("Mo", "O", PMO, po, skf, E3D, shifts),
                 ("S", "O", PS, po, skf, E3D, shifts),
                 ("O", "O", po, po, skf, E3D, shifts)]:
        oc.gen_pair(task)
    # sub_O DFTB
    dwd = os.path.join(wd, "defect")
    os.makedirs(dwd, exist_ok=True)
    atoms = read(SUBO_POSCAR)
    ase_write(os.path.join(dwd, "geo.gen"), atoms, format="gen")
    open(os.path.join(dwd, "dftb_in.hsd"), "w").write(HSD.format(skf=skf))
    env = dict(os.environ, OMP_NUM_THREADS="8")
    r = subprocess.run([DFTB], cwd=dwd, capture_output=True, text=True,
                       timeout=1200, env=env)
    tag = os.path.join(dwd, "results.tag")
    if r.returncode != 0 or not os.path.exists(tag):
        return None
    eig, occ = parse_tag_eigs(tag)
    e = (eig[0] if eig.ndim == 3 else eig) * HA
    o = (occ[0] if occ.ndim == 3 else occ)
    nocc = int((o[0] > 0.5).sum())
    vbm = e[o > 0.5].max()
    cbm = e[o < 0.5].min()
    gap = cbm - vbm
    # ギャップ内フラット準位チェック: 最低空バンド群がフラットで、その上に
    # >0.25 eV のギャップがある場合のみ (V_S 型の孤立準位) を検出
    ingap = 0.0
    emeta = [(e[:, ib].mean(), e[:, ib].max() - e[:, ib].min())
             for ib in range(nocc, min(nocc + 8, e.shape[1]))]
    emeta.sort()
    for i in range(len(emeta) - 1):
        if emeta[i + 1][0] - emeta[i][0] > 0.25:
            if all(w < 0.08 for _, w in emeta[:i + 1]):
                ingap = 1.0
            break
    # O2s: 最深占有バンドが孤立しているか
    occ_means = np.array([e[:, ib].mean() for ib in range(nocc)])
    order = np.argsort(occ_means)
    deepest, second = occ_means[order[0]], occ_means[order[1]]
    if second - deepest > 1.0:
        o2s = deepest - vbm
    else:
        o2s = None
    # VB 窓スペクトル (バンド平均ベース)
    bm = occ_means - vbm
    win = np.sort(bm[bm > -8.0])
    n = min(len(win), len(REF["vb"]))
    vb_rms = float(np.sqrt(((win[-n:] - REF["vb"][-n:]) ** 2).mean()))
    return dict(gap=float(gap), o2s=o2s, vb_rms=vb_rms, ingap=ingap)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("n_trials", type=int)
    ap.add_argument("--tag", default="opto1")
    args = ap.parse_args()

    import optuna
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    base = os.path.join(ROOT, "local_opt", args.tag)
    os.makedirs(base, exist_ok=True)
    study = optuna.create_study(study_name=args.tag,
                                storage=f"sqlite:///{base}/study.db",
                                load_if_exists=True, direction="minimize")
    if len(study.trials) == 0:
        study.enqueue_trial({"O_rd": 3.74, "O_r2s": 2.49, "O_r2p": 2.49,
                             "O_swf": 2.0, "O_sh2s": 0.0, "O_sh2p": 0.0})

    def objective(trial):
        po = {"r_dens": trial.suggest_float("O_rd", 2.2, 9.0),
              "s_wf": trial.suggest_float("O_swf", 1.5, 4.5),
              "r_wf": {"2s": trial.suggest_float("O_r2s", 1.6, 5.0),
                       "2p": trial.suggest_float("O_r2p", 1.6, 5.0)}}
        oshifts = {"2s": trial.suggest_float("O_sh2s", -0.10, 0.10),
                   "2p": trial.suggest_float("O_sh2p", -0.12, 0.12)}
        wd = os.path.join(base, f"trial_{trial.number}")
        try:
            res = evaluate(po, oshifts, wd)
        except Exception as ex:
            print(f"trial {trial.number}: FAILED {str(ex)[:150]}", flush=True)
            res = None
        finally:
            shutil.rmtree(wd, ignore_errors=True)
        if res is None:
            return 1e3
        loss = 5.0 * (res["gap"] - REF["gap"]) ** 2 + 2.0 * res["ingap"]
        if res["o2s"] is None:
            loss += 3.0
            o2s_txt = "lost"
        else:
            loss += 2.0 * (res["o2s"] - REF["o2s"]) ** 2
            o2s_txt = f"{res['o2s']:.2f}"
        loss += 1.0 * res["vb_rms"] ** 2
        for k, v in res.items():
            if v is not None:
                trial.set_user_attr(k, float(v))
        print(f"trial {trial.number}: loss={loss:.4f} gap={res['gap']:.3f} "
              f"O2s={o2s_txt} (ref {REF['o2s']:.2f}) vbRMS={res['vb_rms']:.3f}",
              flush=True)
        return loss

    study.optimize(objective, n_trials=args.n_trials)
    t = study.best_trial
    print("BEST:", round(study.best_value, 4),
          {k: round(v, 4) for k, v in t.user_attrs.items()})
    print({k: round(v, 3) for k, v in t.params.items()})


if __name__ == "__main__":
    main()
