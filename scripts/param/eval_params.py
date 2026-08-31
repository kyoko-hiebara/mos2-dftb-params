#!/usr/bin/env python3
"""与えた confinement パラメータで SKF 生成 + DFTB バンド + LAK 比較を単発実行。

使い方:
  python3 eval_params.py <params.json|optuna:TAG> <a> <thickness> <vasp_json> <outdir>
     [--fine] [--s-lmax p|d] [--plot out.png]

params.json は optuna フラット形式 {"Mo_rd":.., "Mo_r4d":.., ..., "S_r3d":(任意)}
optuna:TAG 指定で study のベストパラメータを使用。
--fine で本番グリッド (dr=0.02, N=980, ntheta=150, nr=50)。
"""
import argparse
import json
import os
import subprocess
import sys
from multiprocessing import Pool

SCRIPTS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPTS)
from compare_bands import band_loss, plot  # noqa: E402

from hotcent.atomic_dft import AtomicDFT  # noqa: E402
from hotcent.confinement import PowerConfinement  # noqa: E402
from hotcent.offsite_twocenter import Offsite2cTable  # noqa: E402

XC = "GGA_X_PBE+GGA_C_PBE"
FREE = json.load(open(os.path.join(SCRIPTS, "free_atom_props.json")))

GRID = {}  # set in main


def make_spec(el, with_s_d):
    if el == "Mo":
        return dict(configuration="[Kr] 4d5 5s1 5p0", valence=["4d", "5s", "5p"],
                    scalarrel=True, occupations={"4d": 5, "5s": 1, "5p": 0})
    if el == "S" and with_s_d:
        return dict(configuration="[Ne] 3s2 3p4 3d0", valence=["3s", "3p", "3d"],
                    scalarrel=False, occupations={"3s": 2, "3p": 4, "3d": 0})
    if el == "S":
        return dict(configuration="[Ne] 3s2 3p4", valence=["3s", "3p"],
                    scalarrel=False, occupations={"3s": 2, "3p": 4})
    raise ValueError(el)


def unflatten(p):
    prm = {"Mo": {"r_dens": p["Mo_rd"],
                  "r_wf": {"4d": p["Mo_r4d"], "5s": p["Mo_r5s"], "5p": p["Mo_r5p"]}},
           "S": {"r_dens": p["S_rd"],
                 "r_wf": {"3s": p["S_r3s"], "3p": p["S_r3p"]}}}
    if "S_r3d" in p:
        prm["S"]["r_wf"]["3d"] = p["S_r3d"]
    if "Mo_swf" in p:
        prm["Mo"]["s_wf"] = p["Mo_swf"]
    if "S_swf" in p:
        prm["S"]["s_wf"] = p["S_swf"]
    return prm


def get_shifts(p):
    sh = {"Mo": {}, "S": {}}
    for key, el, nl in [("Mo_sh4d", "Mo", "4d"), ("Mo_sh5s", "Mo", "5s"),
                        ("Mo_sh5p", "Mo", "5p"), ("S_sh3s", "S", "3s"),
                        ("S_sh3p", "S", "3p")]:
        if key in p:
            sh[el][nl] = p[key]
    return sh


def make_atom(el, p, with_s_d):
    spec = make_spec(el, with_s_d)
    s_wf = p.get("s_wf", 2.0)
    wf_conf = {nl: PowerConfinement(r0=p["r_wf"][nl], s=s_wf)
               for nl in spec["valence"]}
    atom = AtomicDFT(el, xc=XC, configuration=spec["configuration"],
                     valence=spec["valence"], scalarrel=spec["scalarrel"],
                     confinement=PowerConfinement(r0=p["r_dens"], s=2.0),
                     wf_confinement=wf_conf,
                     perturbative_confinement=False, txt=None)
    atom.run()
    return atom, spec


def free_props_with_d(el, spec):
    """S に 3d を含む場合は自由原子を再計算して 3d 固有値を得る"""
    atom = AtomicDFT(el, xc=XC, configuration=spec["configuration"],
                     valence=spec["valence"], scalarrel=spec["scalarrel"],
                     confinement=PowerConfinement(r0=40.0, s=4),
                     perturbative_confinement=False, txt=None)
    atom.run()
    eig = {nl: atom.get_eigenvalue(nl) for nl in spec["valence"]}
    U3p = FREE["S"]["U"]["3p"]
    hub = {nl: U3p for nl in spec["valence"]}
    return eig, hub


def gen_pair(task):
    el1, el2, prm, outdir, with_s_d, e3d, shifts = task
    atom1, spec1 = make_atom(el1, prm[el1], with_s_d)
    atom2 = atom1 if el2 == el1 else make_atom(el2, prm[el2], with_s_d)[0]
    cwd = os.getcwd()
    os.chdir(outdir)
    try:
        off2c = Offsite2cTable(atom1, atom2)
        off2c.run(GRID["rmin"], GRID["dr"], GRID["N"],
                  ntheta=GRID["ntheta"], nr=GRID["nr"],
                  superposition="density", xc=XC)
        kw = {}
        if el1 == el2:
            eig = dict(FREE[el1]["eig"])
            hub = dict(FREE[el1]["U"])
            for nl, dv in shifts.get(el1, {}).items():
                eig[nl] = eig[nl] + dv
            if el1 == "S" and with_s_d:
                eig["3d"] = e3d if e3d is not None else                     free_props_with_d(el1, spec1)[0]["3d"]
                hub["3d"] = hub["3p"]
            kw = dict(eigenvalues=eig, hubbardvalues=hub,
                      occupations=spec1["occupations"], spe=0.0)
        off2c.write(filename_template="{el1}-{el2}.skf", **kw)
    finally:
        os.chdir(cwd)
    return f"{el1}-{el2}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("params")
    ap.add_argument("a", type=float)
    ap.add_argument("thickness", type=float)
    ap.add_argument("vasp_json")
    ap.add_argument("outdir")
    ap.add_argument("--fine", action="store_true")
    ap.add_argument("--grid", default=None,
                    help="dr,N,ntheta,nr (e.g. 0.025,780,150,50)")
    ap.add_argument("--s-lmax", choices=["p", "d"], default="p")
    ap.add_argument("--plot", default=None)
    args = ap.parse_args()

    GRID.update(dict(rmin=0.4, dr=0.02, N=980, ntheta=150, nr=50) if args.fine
                else dict(rmin=0.4, dr=0.025, N=780, ntheta=100, nr=40))
    if args.grid:
        dr, N, ntheta, nr = args.grid.split(",")
        GRID.update(dict(rmin=0.4, dr=float(dr), N=int(N),
                         ntheta=int(ntheta), nr=int(nr)))

    if args.params.startswith("optuna:"):
        import optuna
        optuna.logging.set_verbosity(optuna.logging.WARNING)
        tag = args.params.split(":", 1)[1]
        study = optuna.load_study(
            study_name=tag,
            storage=f"sqlite:////workspace/MoS2_DFTB/dftb/{tag}/study.db")
        p = dict(study.best_trial.params)
        print("using optuna best:", {k: round(v, 3) for k, v in p.items()},
              "loss:", round(study.best_value, 4))
    else:
        p = json.load(open(args.params))

    with_s_d = args.s_lmax == "d"
    if with_s_d and "S_r3d" not in p:
        p["S_r3d"] = p["S_r3p"]

    prm = unflatten(p)
    skf = os.path.join(args.outdir, "skf")
    os.makedirs(skf, exist_ok=True)
    e3d = p.get("S_e3d")
    shifts = get_shifts(p)
    tasks = [("Mo", "Mo", prm, skf, with_s_d, e3d, shifts),
             ("Mo", "S", prm, skf, with_s_d, e3d, shifts),
             ("S", "S", prm, skf, with_s_d, e3d, shifts)]
    with Pool(3) as pool:
        pool.map(gen_pair, tasks)

    dftb_json = os.path.join(args.outdir, "dftb.json")
    cmd = [sys.executable, os.path.join(SCRIPTS, "dftb_bands.py"), skf,
           str(args.a), os.path.join(args.outdir, "dftb"),
           "--json", dftb_json, "--thickness", str(args.thickness)]
    if with_s_d:
        cmd += ["--s-lmax", "d"]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=1200)
    print(r.stdout[-400:])
    if r.returncode != 0:
        print("DFTB FAILED", r.stderr[-400:])
        sys.exit(1)

    band_loss(args.vasp_json, dftb_json, verbose=True)

    # CBM の k 位置診断
    import numpy as np
    d = json.load(open(dftb_json))
    e = np.array(d["eigs"])
    occ = np.array(d["occs"]) > 0.5
    nvb = int(occ[0].sum())
    cbm_k = e[:, nvb:].min(axis=1)
    vbm_k = e[:, :nvb].max(axis=1)
    print(f"DFTB CBM at k-idx {cbm_k.argmin()} (G=0, M=20, K=32, G=56); "
          f"VBM at k-idx {vbm_k.argmax()}")

    if args.plot:
        plot(args.vasp_json, dftb_json, args.plot)


if __name__ == "__main__":
    main()
