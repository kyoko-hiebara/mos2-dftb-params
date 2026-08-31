#!/usr/bin/env python3
"""マルチターゲット最適化 (ローカル Mac 版):
単層バンド構造 (LAK) + V_S 欠陥準位深さを同時フィット。

2 段階: バンド損失 <= BAND_GATE の trial のみ 5x5 欠陥スーパーセルを評価。
使い方:
  python3 optimize_multi.py <n_trials> [--tag optm] [--seed-db path/to/study.db]
複数プロセス並列可 (同じ sqlite storage)。
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
from multiprocessing import Pool

ROOT = "/Users/crocus/uhuhu/MoS2_DFTB"
SCRIPTS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPTS)

from compare_bands import band_loss  # noqa: E402
from optimize_confinement import gen_pair  # noqa: E402  (SKF 生成 16 変数版)

import numpy as np  # noqa: E402
from ase.io import read  # noqa: E402
from ase.io import write as ase_write  # noqa: E402

DFTB = f"{ROOT}/sw_local/dftbplus-install/bin/dftb+"
VASP_JSON = f"{ROOT}/results/ref_calc/bands_lak_a316/bands_lak.json"
DEFECT_POSCAR = f"{ROOT}/results/ref_calc/defects/vac_S_gpu/CONTCAR"
TARGET = json.load(open(f"{ROOT}/results/vac_s_target.json"))
TARGET_DEPTH = TARGET["depth_below_cbm"]          # 0.554 eV
A_LAT, THICKNESS = 3.16, 3.127

BAND_GATE = 0.32      # これ以下のバンド損失で欠陥評価に進む
SKIP_PENALTY = 0.46   # 未評価時の代替ペナルティ (現 v2 の欠陥誤差相当)
W_DEF = 8.0
HA = 27.211386245988

DEFECT_HSD = """Geometry = GenFormat {{
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


def parse_tag_eigs(path):
    lines = open(path).readlines()
    eig = occ = None
    i = 0
    while i < len(lines):
        ln = lines[i]
        if ln.startswith("eigenvalues") or ln.startswith("filling"):
            shape = [int(x) for x in ln.split(":")[-1].strip().split(",")]
            n = int(np.prod(shape))
            vals = []
            j = i + 1
            while len(vals) < n:
                vals.extend(float(x) for x in lines[j].split())
                j += 1
            arr = np.array(vals).reshape(shape[::-1])
            if ln.startswith("eigenvalues"):
                eig = arr
            else:
                occ = arr
            i = j
            continue
        i += 1
    return eig, occ


def defect_depth(skf_dir, workdir):
    """5x5 V_S スーパーセルの DFTB 計算 -> ホスト CBM(min) - 欠陥準位(mean)"""
    os.makedirs(workdir, exist_ok=True)
    atoms = read(DEFECT_POSCAR)
    ase_write(os.path.join(workdir, "geo.gen"), atoms, format="gen")
    with open(os.path.join(workdir, "dftb_in.hsd"), "w") as f:
        f.write(DEFECT_HSD.format(skf=skf_dir))
    env = dict(os.environ, OMP_NUM_THREADS="6")
    r = subprocess.run([DFTB], cwd=workdir, capture_output=True, text=True,
                       timeout=2400, env=env)
    tag = os.path.join(workdir, "results.tag")
    if r.returncode != 0 or not os.path.exists(tag):
        return None
    eig, occ = parse_tag_eigs(tag)
    e = (eig[0] if eig.ndim == 3 else eig) * HA
    o = occ[0] if occ.ndim == 3 else occ
    if e.shape[0] != o.shape[0] or e.shape[0] < e.shape[1]:
        pass
    nocc = int((o[0] > 0.5).sum())
    nb = e.shape[1]
    means = sorted(e[:, ib].mean() for ib in range(nocc, min(nocc + 10, nb)))
    split = None
    for i in range(len(means) - 1):
        if means[i + 1] - means[i] > 0.25:
            split = i
            break
    if split is None:
        return None  # 分離不能 (欠陥準位が CB に埋まった等) -> ペナルティ側で処理
    defect_means = means[:split + 1]
    host_cbm = min(e[:, ib].min() for ib in range(nocc, min(nocc + 10, nb))
                   if e[:, ib].mean() > means[split] + 0.25)
    return float(host_cbm - np.mean(defect_means))


def suggest(trial):
    p = {"Mo": {"r_dens": trial.suggest_float("Mo_rd", 6.0, 16.0),
                "r_wf": {"4d": trial.suggest_float("Mo_r4d", 3.5, 9.0),
                         "5s": trial.suggest_float("Mo_r5s", 3.5, 9.0),
                         "5p": trial.suggest_float("Mo_r5p", 3.0, 9.0)}},
         "S": {"r_dens": trial.suggest_float("S_rd", 4.0, 12.0),
               "r_wf": {"3s": trial.suggest_float("S_r3s", 2.5, 7.0),
                        "3p": trial.suggest_float("S_r3p", 2.5, 7.0),
                        "3d": trial.suggest_float("S_r3d", 2.5, 7.0)}}}
    p["Mo"]["s_wf"] = trial.suggest_float("Mo_swf", 1.5, 5.0)
    p["S"]["s_wf"] = trial.suggest_float("S_swf", 1.5, 5.0)
    e3d = trial.suggest_float("S_e3d", 0.05, 0.80)
    shifts = {"Mo": {"4d": trial.suggest_float("Mo_sh4d", -0.06, 0.06),
                     "5s": trial.suggest_float("Mo_sh5s", -0.06, 0.06),
                     "5p": trial.suggest_float("Mo_sh5p", -0.06, 0.06)},
              "S": {"3s": trial.suggest_float("S_sh3s", -0.06, 0.06),
                    "3p": trial.suggest_float("S_sh3p", -0.06, 0.06)}}
    return p, e3d, shifts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("n_trials", type=int)
    ap.add_argument("--tag", default="optm")
    ap.add_argument("--seed-db", default=None)
    args = ap.parse_args()

    import optuna
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    base = os.path.join(ROOT, "local_opt", args.tag)
    os.makedirs(base, exist_ok=True)
    storage = f"sqlite:///{base}/study.db"
    study = optuna.create_study(study_name=args.tag, storage=storage,
                                load_if_exists=True, direction="minimize")

    if args.seed_db and len(study.trials) == 0:
        old = optuna.load_study(study_name="opt_r2",
                                storage=f"sqlite:///{args.seed_db}")
        done = sorted([t for t in old.trials
                       if t.state.name == "COMPLETE" and t.value < 900],
                      key=lambda t: t.value)
        for t in done[:6]:
            study.enqueue_trial(t.params)
        print(f"seeded {min(6, len(done))} from opt_r2 "
              f"(best {done[0].value:.4f})", flush=True)

    def objective(trial):
        p, e3d, shifts = suggest(trial)
        reg = sum(v ** 2 for d in shifts.values() for v in d.values())
        wd = os.path.join(base, f"trial_{trial.number}")
        skf = os.path.join(wd, "skf")
        os.makedirs(skf, exist_ok=True)
        try:
            tasks = [("Mo", "Mo", p["Mo"], p["Mo"], skf, e3d, shifts),
                     ("Mo", "S", p["Mo"], p["S"], skf, e3d, shifts),
                     ("S", "S", p["S"], p["S"], skf, e3d, shifts)]
            with Pool(3) as pool:
                pool.map(gen_pair, tasks)
            dftb_json = os.path.join(wd, "dftb.json")
            r = subprocess.run(
                [sys.executable, os.path.join(SCRIPTS, "dftb_bands.py"),
                 skf, str(A_LAT), os.path.join(wd, "dftb"),
                 "--json", dftb_json, "--thickness", str(THICKNESS),
                 "--s-lmax", "d"],
                capture_output=True, text=True, timeout=900)
            if r.returncode != 0 or not os.path.exists(dftb_json):
                raise RuntimeError("band step failed: " + r.stdout[-200:])
            res = band_loss(VASP_JSON, dftb_json, verbose=False)
            core = res["loss"] + 1.0 * reg
            if core <= BAND_GATE:
                depth = defect_depth(skf, os.path.join(wd, "defect"))
                if depth is None:
                    loss = core + W_DEF * 0.35 ** 2
                    res["depth"] = -1.0
                else:
                    loss = core + W_DEF * (depth - TARGET_DEPTH) ** 2
                    res["depth"] = depth
                res["defect_evaluated"] = 1.0
            else:
                loss = core + SKIP_PENALTY
                res["defect_evaluated"] = 0.0
            res["core"] = core
            for k, v in res.items():
                trial.set_user_attr(k, float(v))
            dtxt = (f"depth={res.get('depth', -9):.3f}"
                    if res["defect_evaluated"] else "skip")
            print(f"trial {trial.number}: loss={loss:.4f} core={core:.4f} "
                  f"gapK={res['gap_K_dftb']:.3f} {dtxt}", flush=True)
            return loss
        except Exception as e:
            print(f"trial {trial.number}: FAILED {str(e)[:200]}", flush=True)
            return 1e3
        finally:
            shutil.rmtree(wd, ignore_errors=True)

    study.optimize(objective, n_trials=args.n_trials)
    t = study.best_trial
    print("BEST:", round(study.best_value, 4),
          {k: round(v, 4) for k, v in t.user_attrs.items()})
    print({k: round(v, 3) for k, v in t.params.items()})


if __name__ == "__main__":
    main()
