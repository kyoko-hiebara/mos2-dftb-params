#!/usr/bin/env python3
"""confinement パラメータの最適化 (optuna TPE)。S は spd 基底。

9 変数: Mo r_dens, r_wf(4d,5s,5p); S r_dens, r_wf(3s,3p,3d), onsite(3d)
ターゲット: LAK バンド構造 (a=3.16, 実験格子)。

使い方:
  python3 optimize_confinement.py <a_lat> <thickness> <vasp_bands.json> <n_trials> [--tag opt_spd]
複数プロセス同時起動可 (同じ sqlite storage を共有)。
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

from hotcent.atomic_dft import AtomicDFT  # noqa: E402
from hotcent.confinement import PowerConfinement  # noqa: E402
from hotcent.offsite_twocenter import Offsite2cTable  # noqa: E402

XC = "GGA_X_PBE+GGA_C_PBE"
RMIN, DR, NPTS = 0.4, 0.025, 780   # 粗グリッド (最適化用)
NTHETA, NR = 100, 40

SPECS = {
    "Mo": dict(configuration="[Kr] 4d5 5s1 5p0", valence=["4d", "5s", "5p"],
               scalarrel=True, occupations={"4d": 5, "5s": 1, "5p": 0}),
    "S": dict(configuration="[Ne] 3s2 3p4 3d0", valence=["3s", "3p", "3d"],
              scalarrel=False, occupations={"3s": 2, "3p": 4, "3d": 0}),
}

FREE = json.load(open(os.path.join(SCRIPTS, "free_atom_props.json")))


def make_atom(el, p):
    spec = SPECS[el]
    s_wf = p.get("s_wf", 2.0)
    wf_conf = {nl: PowerConfinement(r0=p["r_wf"][nl], s=s_wf)
               for nl in spec["valence"]}
    atom = AtomicDFT(el, xc=XC,
                     configuration=spec["configuration"],
                     valence=spec["valence"],
                     scalarrel=spec["scalarrel"],
                     confinement=PowerConfinement(r0=p["r_dens"], s=p.get("s_dens", 2.0)),
                     wf_confinement=wf_conf,
                     perturbative_confinement=False,
                     txt=None)
    atom.run()
    return atom


def gen_pair(task):
    el1, el2, p1, p2, outdir, e3d, shifts = task
    atom1 = make_atom(el1, p1)
    atom2 = atom1 if el2 == el1 else make_atom(el2, p2)
    cwd = os.getcwd()
    os.chdir(outdir)
    try:
        off2c = Offsite2cTable(atom1, atom2)
        off2c.run(RMIN, DR, NPTS, ntheta=NTHETA, nr=NR,
                  superposition="density", xc=XC)
        kw = {}
        if el1 == el2:
            eig = dict(FREE[el1]["eig"])
            hub = dict(FREE[el1]["U"])
            for nl, dv in shifts.get(el1, {}).items():
                eig[nl] = eig[nl] + dv
            if el1 == "S":
                eig["3d"] = e3d
                hub["3d"] = hub["3p"]
            kw = dict(eigenvalues=eig, hubbardvalues=hub,
                      occupations=SPECS[el1]["occupations"], spe=0.0)
        off2c.write(filename_template="{el1}-{el2}.skf", **kw)
    finally:
        os.chdir(cwd)
    return f"{el1}-{el2}"


def evaluate(params, e3d, shifts, workdir, a_lat, thickness, vasp_json):
    skf = os.path.join(workdir, "skf")
    os.makedirs(skf, exist_ok=True)
    tasks = [("Mo", "Mo", params["Mo"], params["Mo"], skf, e3d, shifts),
             ("Mo", "S", params["Mo"], params["S"], skf, e3d, shifts),
             ("S", "S", params["S"], params["S"], skf, e3d, shifts)]
    with Pool(3) as pool:
        pool.map(gen_pair, tasks)

    dftb_json = os.path.join(workdir, "dftb.json")
    r = subprocess.run(
        [sys.executable, os.path.join(SCRIPTS, "dftb_bands.py"),
         skf, str(a_lat), os.path.join(workdir, "dftb"),
         "--json", dftb_json, "--thickness", str(thickness),
         "--s-lmax", "d"],
        capture_output=True, text=True, timeout=1200)
    if r.returncode != 0 or not os.path.exists(dftb_json):
        return None, r.stdout[-500:] + r.stderr[-300:]
    res = band_loss(vasp_json, dftb_json, verbose=False)
    return res, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("a_lat", type=float)
    ap.add_argument("thickness", type=float)
    ap.add_argument("vasp_json")
    ap.add_argument("n_trials", type=int)
    ap.add_argument("--tag", default="opt_spd")
    args = ap.parse_args()

    import optuna
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    base = f"/workspace/MoS2_DFTB/dftb/{args.tag}"
    os.makedirs(base, exist_ok=True)
    storage = f"sqlite:///{base}/study.db"
    study = optuna.create_study(study_name=args.tag, storage=storage,
                                load_if_exists=True, direction="minimize")

    def objective(trial):
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
        reg = sum(v ** 2 for d in shifts.values() for v in d.values())
        wd = os.path.join(base, f"trial_{trial.number}")
        os.makedirs(wd, exist_ok=True)
        try:
            res, err = evaluate(p, e3d, shifts, wd, args.a_lat,
                                args.thickness, args.vasp_json)
        except Exception as e:
            res, err = None, str(e)
        if res is None:
            print(f"trial {trial.number}: FAILED {err[:200]}", flush=True)
            shutil.rmtree(wd, ignore_errors=True)
            return 1e3
        res["loss"] = res["loss"] + 1.0 * reg
        for k, v in res.items():
            trial.set_user_attr(k, v)
        print(f"trial {trial.number}: loss={res['loss']:.4f} "
              f"rmsVB={res['rms_vb']:.3f} rmsCB={res['rms_cb']:.3f} "
              f"gapK_D={res['gap_K_dftb']:.3f}/min={res['gap_dftb']:.3f} "
              f"(V={res['gap_K_vasp']:.3f})", flush=True)
        shutil.rmtree(wd, ignore_errors=True)
        return res["loss"]

    study.optimize(objective, n_trials=args.n_trials)
    print("BEST:", study.best_params, study.best_value)


if __name__ == "__main__":
    main()
