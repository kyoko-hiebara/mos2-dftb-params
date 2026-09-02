#!/usr/bin/env python3
"""Sb 第 3 世代最適化 (mac): SOC をループ内に入れ、
  (1) PBE+SOC パスバンド (E_F±2 eV 重視)
  (2) 16^3 全 BZ メッシュ (SOC) の E_F 近傍 RMS と状態カウント法オーバーラップ
  (3) Sb(111) 6BL スラブの E_F (= -仕事関数) の PBE 整列
を同時フィット。変数 8 (optsb2 と同じ、シフト範囲拡張)。
使い方: python3 optimize_sb3.py <n_trials> [--tag optsb3]
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

import optimize_sb2 as sb2  # noqa: E402  (SPECS/FREE/make_atom tuned, SB_ATOMS)
from optimize_multi import DFTB, HA, parse_tag_eigs  # noqa: E402

from ase.io import write as ase_write  # noqa: E402
from ase.spacegroup import crystal  # noqa: E402

ROOT = "/Users/crocus/uhuhu/MoS2_DFTB"
oc = sb2.oc
XI_5P = 0.57   # eV (2026-09-01 較正値、ループ内では固定)

REF_SOC = json.load(open(f"{ROOT}/local_opt/sb_ref_bands_soc.json"))
REF_SOC_E = np.array(REF_SOC["eigs"])
REF_SOC_K = np.array(REF_SOC["kpts"])
MESH = np.load(f"{ROOT}/local_opt/sb_mesh_soc.npz")
MESH_K, MESH_E = MESH["kpts"], MESH["e"]           # (4096, 3), (4096, 36) E-EF
assert MESH_E.shape[1] == 36, MESH_E.shape
MESH_VAL = MESH_E[:, 20:]        # 4d セミコア 20 状態を除いた 5s5p 由来 16 状態
N_OCC_VAL = 10                   # 2 Sb x 5 価電子 (SOC 2 成分)
OV_REF = float(MESH_VAL[:, N_OCC_VAL - 1].max() - MESH_VAL[:, N_OCC_VAL].min())
MESH_WIN = 1.5

SLAB_JSON = f"{ROOT}/local_opt/sb_slab_wf.json"


def slab_target():
    if os.environ.get("SB_WF_TARGET"):
        return -float(os.environ["SB_WF_TARGET"])
    return -json.load(open(SLAB_JSON))["W"]


def build_slab():
    conv = crystal("Sb", [(0, 0, 0.2336)], spacegroup=166,
                   cellpar=[4.3084, 4.3084, 11.274, 90, 90, 120],
                   primitive_cell=False)
    slab = conv.repeat((1, 1, 2))
    zs = slab.positions[:, 2]
    thick = zs.max() - zs.min()
    cell = slab.cell.copy()
    cell[2] = [0, 0, thick + 30.0]
    slab.set_cell(cell)
    slab.positions[:, 2] += 15.0 - zs.min()
    slab.pbc = True
    return slab


SLAB = build_slab()

HSD = """Geometry = GenFormat {{
  <<< "geo.gen"
}}
Hamiltonian = DFTB {{
  SCC = Yes
  SCCTolerance = {tol}
  MaxSCCIterations = {mx}
  Mixer = Broyden {{}}
  SlaterKosterFiles = Type2FileNames {{
    Prefix = "{skf}/"
    Separator = "-"
    Suffix = ".skf"
  }}
  MaxAngularMomentum {{
    Sb = "d"
  }}
  PolynomialRepulsive = SetForAll {{ Yes }}
{soc}
  Filling = Fermi {{ Temperature [K] = 300 }}
  {extra}
{kblock}
}}
Options {{ WriteResultsTag = Yes }}
ParserOptions {{ ParserVersion = 14 }}
"""
SOC_BLOCK = """  SpinOrbit = {{
    Dual = Yes
    Sb [eV] = {{0.0 {xi} 0.0}}
  }}"""


def _run(wd, hsd, nthreads="4"):
    open(os.path.join(wd, "dftb_in.hsd"), "w").write(hsd)
    env = dict(os.environ, OMP_NUM_THREADS=nthreads)
    r = subprocess.run([DFTB], cwd=wd, capture_output=True, text=True,
                       timeout=900, env=env)
    if r.returncode != 0:
        raise RuntimeError("dftb+: " + r.stdout[-300:])


def _ef(wd):
    lines = open(os.path.join(wd, "results.tag")).readlines()
    for i, ln in enumerate(lines):
        if ln.startswith("fermi_level"):
            return float(lines[i + 1].split()[0]) * HA
    raise RuntimeError("no fermi_level")


def _eigs(wd, nk):
    eig, _ = parse_tag_eigs(os.path.join(wd, "results.tag"))
    e = np.squeeze(eig) * HA
    if e.ndim == 1:
        e = e[None, :]
    if e.shape[0] != nk:
        e = e.T
    return np.sort(e, axis=1)


def _kblock(kpts):
    return "  KPointsAndWeights = {\n" + "\n".join(
        f"    {k[0]:.8f} {k[1]:.8f} {k[2]:.8f} 1.0" for k in kpts) + "\n  }"


def bulk_soc(skf, wd):
    """SCC(SOC, 8^3) -> E_F; パス固有値; 16^3 メッシュ固有値 (いずれも E-E_F)"""
    os.makedirs(wd, exist_ok=True)
    ase_write(os.path.join(wd, "geo.gen"), sb2.SB_ATOMS, format="gen")
    soc = SOC_BLOCK.format(xi=XI_5P)
    kscc = ("  KPointsAndWeights = SupercellFolding {\n"
            "    8 0 0\n    0 8 0\n    0 0 8\n    0.5 0.5 0.5\n  }")
    _run(wd, HSD.format(skf=skf, kblock=kscc, tol="1e-6", mx="250", extra="", soc=soc))
    ef = _ef(wd)
    _run(wd, HSD.format(skf=skf, kblock=_kblock(REF_SOC_K), tol="1e6", mx="1",
                        extra="ReadInitialCharges = Yes", soc=soc))
    e_path = _eigs(wd, len(REF_SOC_K)) - ef
    _run(wd, HSD.format(skf=skf, kblock=_kblock(MESH_K), tol="1e6", mx="1",
                        extra="ReadInitialCharges = Yes", soc=soc))
    e_mesh = _eigs(wd, len(MESH_K)) - ef
    return ef, e_path, e_mesh


def slab_ef(skf, wd):
    os.makedirs(wd, exist_ok=True)
    ase_write(os.path.join(wd, "geo.gen"), SLAB, format="gen")
    kscc = ("  KPointsAndWeights = SupercellFolding {\n"
            "    8 0 0\n    0 8 0\n    0 0 1\n    0.5 0.5 0.0\n  }")
    _run(wd, HSD.format(skf=skf, kblock=kscc, tol="1e-5", mx="300", extra="", soc=""))
    return _ef(wd)


def path_loss(e_dftb):
    nk = min(len(e_dftb), len(REF_SOC_E))
    bmeans = REF_SOC_E[:nk].mean(axis=0)
    sel = [ib for ib, m in enumerate(bmeans) if -15.0 < m < 6.0]
    ref = REF_SOC_E[:nk][:, sel[:e_dftb.shape[1]]]
    n = min(ref.shape[1], e_dftb.shape[1])
    ref, dt = ref[:, :n], e_dftb[:nk, :n]
    w = np.where(np.abs(ref) < 2.0, 1.0, 0.25)
    rms_all = float(np.sqrt(((dt - ref) ** 2 * w).sum() / w.sum()))
    msk = np.abs(ref) < 2.0
    rms_ef = float(np.sqrt(((dt[msk] - ref[msk]) ** 2).mean()))
    return rms_all, rms_ef


def mesh_loss(e_mesh):
    """16^3 メッシュ: 参照の 4d 抜き 36 状態 (5s5p5d x2) と DFTB 36 状態を対応"""
    n = MESH_VAL.shape[1]
    ref, dt = MESH_VAL, e_mesh[:, :n]
    msk = np.abs(ref) < MESH_WIN
    rms = float(np.sqrt(((dt[msk] - ref[msk]) ** 2).mean()))
    ov = float(dt[:, N_OCC_VAL - 1].max() - dt[:, N_OCC_VAL].min())
    ih, ie = int(dt[:, N_OCC_VAL - 1].argmax()), int(dt[:, N_OCC_VAL].argmin())
    return rms, ov, MESH_K[ih], MESH_K[ie]


W_PATH_ALL, W_MESH, W_OV, W_AL = 0.3, 1.0, 5.0, 2.0


def evaluate(psb, shifts, e5d, wd):
    skf = os.path.join(wd, "skf")
    os.makedirs(skf, exist_ok=True)
    oc.FREE["Sb"]["eig"]["5d"] = e5d
    oc.FREE["Sb"]["U"]["5d"] = oc.FREE["Sb"]["U"]["5p"]
    oc.gen_pair(("Sb", "Sb", psb, psb, skf, 0.0, {"Sb": shifts}))
    ef, e_path, e_mesh = bulk_soc(skf, os.path.join(wd, "bulk"))
    rms_all, rms_ef = path_loss(e_path)
    rms_mesh, ov, kh, ke = mesh_loss(e_mesh)
    try:
        efs = slab_ef(skf, os.path.join(wd, "slab"))
        al = efs - slab_target()
    except Exception as ex:
        efs, al = -99.0, 1.0
    loss = (rms_ef ** 2 + W_PATH_ALL * rms_all ** 2 + W_MESH * rms_mesh ** 2
            + W_OV * (ov - OV_REF) ** 2 + W_AL * al ** 2)
    return dict(rms_all=rms_all, rms_ef=rms_ef, rms_mesh=rms_mesh, overlap=ov,
                ef_bulk=ef, ef_slab=efs, align_err=al, loss=loss)

ZOOM = {}


def set_zoom(best, frac):
    ZOOM.clear()
    ZOOM.update({"_best": best, "_frac": frac})


def sf(trial, name, lo, hi):
    if ZOOM and name in ZOOM["_best"]:
        b, w = ZOOM["_best"][name], ZOOM["_frac"] * (hi - lo)
        lo, hi = max(lo, b - w), min(hi, b + w)
    return trial.suggest_float(name, lo, hi)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("n_trials", type=int)
    ap.add_argument("--tag", default="optsb3")
    ap.add_argument("--seed-from", default="optsb2")
    ap.add_argument("--nseed", type=int, default=3)
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
    print(f"targets: overlap_ref={OV_REF*1000:+.1f} meV, slab E_F target={slab_target():.3f} eV",
          flush=True)
    if len(study.trials) == 0:
        for src in args.seed_from.split(","):
            old = optuna.load_study(study_name=src, storage=f"sqlite:///{ROOT}/local_opt/{src}/study.db")
            done = [t for t in old.trials if t.state.name == "COMPLETE" and t.value is not None and t.value < 900]
            for t in sorted(done, key=lambda t: t.value)[:args.nseed]:
                study.enqueue_trial(t.params)
            print(f"seeded {min(args.nseed, len(done))} from {src}", flush=True)

    def objective(trial):
        psb = {"r_dens": sf(trial, "Sb_rd", 4.0, 14.0),
               "s_wf": sf(trial, "Sb_swf", 1.5, 5.0),
               "r_wf": {"5s": sf(trial, "Sb_r5s", 3.0, 9.0),
                        "5p": sf(trial, "Sb_r5p", 3.0, 9.0),
                        "5d": sf(trial, "Sb_r5d", 3.0, 9.0)}}
        shifts = {"5s": sf(trial, "Sb_sh5s", -0.15, 0.10),
                  "5p": sf(trial, "Sb_sh5p", -0.15, 0.10)}
        e5d = sf(trial, "Sb_e5d", 0.05, 0.90)
        wd = os.path.join(base, f"trial_{trial.number}")
        try:
            res = evaluate(psb, shifts, e5d, wd)
        except Exception as ex:
            print(f"trial {trial.number}: FAILED {str(ex)[:150]}", flush=True)
            res = None
        finally:
            shutil.rmtree(wd, ignore_errors=True)
        if res is None:
            return 1e3
        for k, v in res.items():
            trial.set_user_attr(k, float(v))
        print(f"trial {trial.number}: loss={res['loss']:.4f} rmsEF={res['rms_ef']:.3f} "
              f"mesh={res['rms_mesh']:.3f} ov={res['overlap']*1000:+.0f}meV "
              f"EFslab={res['ef_slab']:.2f} al={res['align_err']:+.2f}", flush=True)
        return res["loss"]

    study.optimize(objective, n_trials=args.n_trials)
    t = study.best_trial
    print("BEST:", round(study.best_value, 4), {k: round(v, 4) for k, v in t.user_attrs.items()})
    print({k: round(v, 4) for k, v in t.params.items()})


if __name__ == "__main__":
    main()
