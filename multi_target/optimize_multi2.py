#!/usr/bin/env python3
"""マルチターゲット最適化 第 2 世代 (mac):
  LAK バンド (パス + 12x12 メッシュ) + バンド端曲率 (有効質量) + Q 谷 + 絶対整列 + V_S 準位。
optimize_multi.py との差分:
  - 損失に extended_loss (質量/Q/メッシュ/整列) を追加
  - onsite シフトと S 3d onsite の探索範囲を拡張 (v3 で境界張り付き)
  - シフト正則化を廃止 (共通成分は整列項で決まる)
使い方: python3 optimize_multi2.py <n_trials> [--tag optm2] [--seed-from optm1]
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
from multiprocessing import Pool

SCRIPTS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPTS)

from compare_bands import band_loss  # noqa: E402
from optimize_confinement import gen_pair  # noqa: E402
from optimize_multi import (ROOT, VASP_JSON, TARGET_DEPTH, A_LAT, THICKNESS,  # noqa: E402
                            BAND_GATE, SKIP_PENALTY, W_DEF, defect_depth)
import extended_loss as ext  # noqa: E402

MESH_KPTS = os.path.join(SCRIPTS, "mesh_kpts_12x12.json")
MO_S4D = os.environ.get("MO_S4D", "0") == "1"
ext.write_mesh_kpts(MESH_KPTS)

ZOOM = {}


def set_zoom(best, frac):
    ZOOM.clear()
    ZOOM.update({"_best": best, "_frac": frac})


def sf(trial, name, lo, hi):
    if ZOOM and name in ZOOM["_best"]:
        b, w = ZOOM["_best"][name], ZOOM["_frac"] * (hi - lo)
        lo, hi = max(lo, b - w), min(hi, b + w)
    return trial.suggest_float(name, lo, hi)


def suggest(trial):
    p = {"Mo": {"r_dens": sf(trial, "Mo_rd", 6.0, 16.0),
                "r_wf": {"4d": sf(trial, "Mo_r4d", 3.0, 9.0),
                         "5s": sf(trial, "Mo_r5s", 3.0, 9.0),
                         "5p": sf(trial, "Mo_r5p", 3.0, 9.0)}},
         "S": {"r_dens": sf(trial, "S_rd", 4.0, 12.0),
               "r_wf": {"3s": sf(trial, "S_r3s", 2.2, 7.0),
                        "3p": sf(trial, "S_r3p", 2.5, 7.0),
                        "3d": sf(trial, "S_r3d", 2.5, 8.0)}}}
    p["Mo"]["s_wf"] = sf(trial, "Mo_swf", 1.5, 5.0)
    p["S"]["s_wf"] = sf(trial, "S_swf", 1.5, 5.0)
    if MO_S4D:
        p["Mo"]["s_wf_shell"] = {"4d": sf(trial, "Mo_s4d", 1.2, 6.0)}
    e3d = sf(trial, "S_e3d", 0.05, 1.50)
    shifts = {"Mo": {"4d": sf(trial, "Mo_sh4d", -0.06, 0.16),
                     "5s": sf(trial, "Mo_sh5s", -0.10, 0.12),
                     "5p": sf(trial, "Mo_sh5p", -0.10, 0.12)},
              "S": {"3s": sf(trial, "S_sh3s", -0.10, 0.10),
                    "3p": sf(trial, "S_sh3p", -0.10, 0.10)}}
    return p, e3d, shifts


def evaluate(p, e3d, shifts, wd):
    """SKF 生成 -> バンド (パス+メッシュ) -> 損失辞書 (欠陥は呼び出し側)"""
    skf = os.path.join(wd, "skf")
    os.makedirs(skf, exist_ok=True)
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
         "--s-lmax", "d", "--extra-kpts", MESH_KPTS],
        capture_output=True, text=True, timeout=900)
    if r.returncode != 0 or not os.path.exists(dftb_json):
        raise RuntimeError("band step failed: " + r.stdout[-200:])
    res = band_loss(VASP_JSON, dftb_json, verbose=False)
    res.update(ext.all_terms(VASP_JSON, dftb_json))
    return res, skf


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("n_trials", type=int)
    ap.add_argument("--tag", default="optm2")
    ap.add_argument("--seed-from", default=None)
    ap.add_argument("--nseed", type=int, default=6)
    ap.add_argument("--zoom", type=float, default=None,
                    help="seed-from の最良点周りに各変数の範囲を ±zoom*range に絞る")
    args = ap.parse_args()

    import optuna
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    base = os.path.join(ROOT, "local_opt", args.tag)
    os.makedirs(base, exist_ok=True)
    if args.zoom:
        src0 = args.seed_from.split(",")[0]
        src = optuna.load_study(study_name=src0, storage=f"sqlite:///{ROOT}/local_opt/{src0}/study.db")
        set_zoom(dict(src.best_trial.params), args.zoom)
        print(f"zoom {args.zoom} around {src0} best ({src.best_value:.4f})", flush=True)
    study = optuna.create_study(study_name=args.tag,
                                storage=f"sqlite:///{base}/study.db",
                                load_if_exists=True, direction="minimize")
    if args.seed_from and len(study.trials) == 0:
        old = optuna.load_study(study_name=args.seed_from,
                                storage=f"sqlite:///{ROOT}/local_opt/{args.seed_from}/study.db")
        done = sorted([t for t in old.trials
                       if t.state.name == "COMPLETE" and t.value < 900],
                      key=lambda t: t.value)
        for t in done[:args.nseed]:
            study.enqueue_trial(t.params)
        print(f"seeded {min(args.nseed, len(done))} from {args.seed_from}", flush=True)
    print(f"alignment target (K midgap rel. vacuum) = {ext.align_target():.3f} eV; weights W_M={ext.W_M} W_Q={ext.W_Q} "
          f"W_MESH={ext.W_MESH} W_AL={ext.W_AL}; MO_S4D={MO_S4D}", flush=True)

    def objective(trial):
        p, e3d, shifts = suggest(trial)
        wd = os.path.join(base, f"trial_{trial.number}")
        try:
            res, skf = evaluate(p, e3d, shifts, wd)
            core = res["loss"]                      # 従来のバンド損失 (ゲート用)
            total = core + res["loss_ext"]
            if core <= BAND_GATE:
                depth = defect_depth(skf, os.path.join(wd, "defect"))
                if depth is None:
                    total += W_DEF * 0.35 ** 2
                    res["depth"] = -1.0
                else:
                    total += W_DEF * (depth - TARGET_DEPTH) ** 2
                    res["depth"] = depth
                res["defect_evaluated"] = 1.0
            else:
                total += SKIP_PENALTY
                res["defect_evaluated"] = 0.0
            res["core"] = core
            for k, v in res.items():
                trial.set_user_attr(k, float(v))
            dtxt = (f"depth={res.get('depth', -9):.3f}" if res["defect_evaluated"] else "skip")
            print(f"trial {trial.number}: loss={total:.4f} core={core:.4f} ext={res['loss_ext']:.3f} "
                  f"gapK={res['gap_K_dftb']:.3f} mVB={res['m_VB_KM_dftb']:.2f} mCB={res['m_CB_KM_dftb']:.2f} "
                  f"dQ={res['dq_dftb']:.3f} mesh={res['rms_mesh']:.3f} al={res['align_err']:+.2f} {dtxt}",
                  flush=True)
            return total
        except Exception as e:
            print(f"trial {trial.number}: FAILED {str(e)[:200]}", flush=True)
            return 1e3
        finally:
            shutil.rmtree(wd, ignore_errors=True)

    study.optimize(objective, n_trials=args.n_trials)
    t = study.best_trial
    print("BEST:", round(study.best_value, 4),
          {k: round(v, 4) for k, v in t.user_attrs.items()})
    print({k: round(v, 4) for k, v in t.params.items()})


if __name__ == "__main__":
    main()
