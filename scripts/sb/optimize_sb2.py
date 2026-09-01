#!/usr/bin/env python3
"""Sb (A7 バルク) の confinement/onsite 最適化 — GPAW PBE 参照 (mac 完結)。

損失: フェルミ整列したバンドの重み付き RMS (E_F±2 eV を重視)。
変数 6: Sb r_dens, r_wf(5s), r_wf(5p), s_wf, shift(5s), shift(5p)

使い方: python3 optimize_sb.py <n_trials> [--tag optsb1]
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

from ase.spacegroup import crystal  # noqa: E402
from ase.io import write as ase_write  # noqa: E402
from hotcent.atomic_dft import AtomicDFT  # noqa: E402
from hotcent.confinement import PowerConfinement  # noqa: E402

ROOT = "/Users/crocus/uhuhu/MoS2_DFTB"
REF = json.load(open(f"{ROOT}/local_opt/sb_ref_bands.json"))
REF_E = np.array(REF["eigs"])          # (nk, nb), E-EF
REF_K = np.array(REF["kpts"])

oc.SPECS["Sb"] = dict(configuration="[Kr] 4d10 5s2 5p3 5d0",
                      valence=["5s", "5p", "5d"], scalarrel=True,
                      occupations={"5s": 2, "5p": 3, "5d": 0})

SB_ATOMS = crystal("Sb", [(0, 0, 0.2336)], spacegroup=166,
                   cellpar=[4.3084, 4.3084, 11.274, 90, 90, 120],
                   primitive_cell=True)

# ---- Sb 自由原子 (キャッシュ) ----
FREE_SB_FILE = f"{ROOT}/local_opt/free_sb.json"
if os.path.exists(FREE_SB_FILE):
    FREE_SB = json.load(open(FREE_SB_FILE))
else:
    atom = AtomicDFT("Sb", xc="GGA_X_PBE+GGA_C_PBE",
                     configuration="[Kr] 4d10 5s2 5p3",
                     valence=["5s", "5p"], scalarrel=True,
                     confinement=PowerConfinement(r0=40.0, s=4),
                     perturbative_confinement=False, mix=0.12,
                     nodegpts=180, txt=None)
    atom.run()
    eig = {nl: atom.get_eigenvalue(nl) for nl in ["5s", "5p"]}
    try:
        U = atom.get_hubbard_value("5p", scheme="forward", maxstep=1)
    except Exception as ex:
        print("Hubbard central failed, fallback:", str(ex)[:80])
        U = 0.21  # Ha, 代表的な 5p 値 (後で refine 可)
    FREE_SB = {"eig": eig, "U": {"5s": U, "5p": U}}
    json.dump(FREE_SB, open(FREE_SB_FILE, "w"))
oc.FREE["Sb"] = FREE_SB
print("Sb free atom:", {k: round(v, 4) for k, v in FREE_SB["eig"].items()},
      flush=True)


def _make_atom_tuned(el, p):
    spec = oc.SPECS[el]
    s_wf = p.get("s_wf", 2.0)
    wf_conf = {nl: PowerConfinement(r0=p["r_wf"][nl], s=s_wf)
               for nl in spec["valence"]}
    atom = AtomicDFT(el, xc="GGA_X_PBE+GGA_C_PBE",
                     configuration=spec["configuration"],
                     valence=spec["valence"],
                     scalarrel=spec["scalarrel"],
                     confinement=PowerConfinement(
                         r0=p["r_dens"], s=p.get("s_dens", 2.0)),
                     wf_confinement=wf_conf,
                     perturbative_confinement=False,
                     mix=0.12, nodegpts=180, txt=None)
    atom.run()
    return atom


oc.make_atom = _make_atom_tuned

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
  Filling = Fermi {{ Temperature [K] = 300 }}
  {extra}
{kblock}
}}
Options {{ WriteResultsTag = Yes }}
ParserOptions {{ ParserVersion = 14 }}
"""


def sb_bands(skf_dir, wd):
    os.makedirs(wd, exist_ok=True)
    ase_write(os.path.join(wd, "geo.gen"), SB_ATOMS, format="gen")
    env = dict(os.environ, OMP_NUM_THREADS="6")
    kscc = ("  KPointsAndWeights = SupercellFolding {\n"
            "    8 0 0\n    0 8 0\n    0 0 8\n    0.5 0.5 0.5\n  }")
    open(os.path.join(wd, "dftb_in.hsd"), "w").write(
        HSD.format(skf=skf_dir, kblock=kscc, tol="1e-6", mx="200", extra=""))
    r = subprocess.run([DFTB], cwd=wd, capture_output=True, text=True,
                       timeout=600, env=env)
    if r.returncode != 0:
        raise RuntimeError("SCC: " + r.stdout[-300:])
    # E_F from SCC
    ef = None
    for ln in open(os.path.join(wd, "results.tag")):
        if ln.startswith("fermi_level"):
            ef = None
        if ef is None and ln.strip() and not ln.startswith(("fermi_level",)):
            pass
    lines = open(os.path.join(wd, "results.tag")).readlines()
    for i, ln in enumerate(lines):
        if ln.startswith("fermi_level"):
            ef = float(lines[i + 1].split()[0]) * HA
            break
    # band pass
    kl = "  KPointsAndWeights = {\n" + "\n".join(
        f"    {k[0]:.8f} {k[1]:.8f} {k[2]:.8f} 1.0" for k in REF_K) + "\n  }"
    open(os.path.join(wd, "dftb_in.hsd"), "w").write(
        HSD.format(skf=skf_dir, kblock=kl, tol="1e6", mx="1",
                   extra="ReadInitialCharges = Yes"))
    r = subprocess.run([DFTB], cwd=wd, capture_output=True, text=True,
                       timeout=600, env=env)
    if r.returncode != 0:
        raise RuntimeError("bands: " + r.stdout[-300:])
    eig, _ = parse_tag_eigs(os.path.join(wd, "results.tag"))
    e = (eig[0] if eig.ndim == 3 else eig) * HA
    return e - ef


def band_loss_sb(e_dftb):
    """重み付き RMS: 参照側の各 (k,n) 状態にエネルギー依存の重み"""
    nk = min(len(e_dftb), len(REF_E))
    # 参照: E_F-15..+6 窓のバンド (バンド平均で選択)
    bmeans = REF_E[:nk].mean(axis=0)
    sel = [ib for ib, m in enumerate(bmeans) if -15.0 < m < 6.0]
    ref = REF_E[:nk][:, sel[:e_dftb.shape[1]]]
    n = min(ref.shape[1], e_dftb.shape[1])
    ref, dt = ref[:, :n], e_dftb[:nk, :n]
    w = np.where(np.abs(ref) < 2.0, 1.0, 0.25)
    rms_all = float(np.sqrt(((dt - ref) ** 2 * w).sum() / w.sum()))
    msk = np.abs(ref) < 2.0
    rms_ef = float(np.sqrt(((dt[msk] - ref[msk]) ** 2).mean()))
    return rms_all, rms_ef


def evaluate(psb, shifts, e5d, wd):
    skf = os.path.join(wd, "skf")
    os.makedirs(skf, exist_ok=True)
    oc.FREE["Sb"]["eig"]["5d"] = e5d
    oc.FREE["Sb"]["U"]["5d"] = oc.FREE["Sb"]["U"]["5p"]
    oc.gen_pair(("Sb", "Sb", psb, psb, skf, 0.0, {"Sb": shifts}))
    e = sb_bands(skf, os.path.join(wd, "dftb"))
    rms_all, rms_ef = band_loss_sb(e)
    return dict(rms_all=rms_all, rms_ef=rms_ef,
                loss=rms_ef ** 2 + 0.3 * rms_all ** 2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("n_trials", type=int)
    ap.add_argument("--tag", default="optsb1")
    args = ap.parse_args()

    import optuna
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    base = os.path.join(ROOT, "local_opt", args.tag)
    os.makedirs(base, exist_ok=True)
    study = optuna.create_study(study_name=args.tag,
                                storage=f"sqlite:///{base}/study.db",
                                load_if_exists=True, direction="minimize")
    if len(study.trials) == 0:
        # rule-of-thumb: rcov(Sb)=1.39 A = 2.63 Bohr
        study.enqueue_trial({"Sb_rd": 6.604, "Sb_r5s": 4.95, "Sb_r5p": 6.045,
                             "Sb_r5d": 6.0, "Sb_swf": 1.914,
                             "Sb_sh5s": 0.02, "Sb_sh5p": 0.013, "Sb_e5d": 0.30})

    def objective(trial):
        psb = {"r_dens": trial.suggest_float("Sb_rd", 4.0, 14.0),
               "s_wf": trial.suggest_float("Sb_swf", 1.5, 5.0),
               "r_wf": {"5s": trial.suggest_float("Sb_r5s", 3.0, 9.0),
                        "5p": trial.suggest_float("Sb_r5p", 3.0, 9.0),
                        "5d": trial.suggest_float("Sb_r5d", 3.0, 9.0)}}
        shifts = {"5s": trial.suggest_float("Sb_sh5s", -0.08, 0.08),
                  "5p": trial.suggest_float("Sb_sh5p", -0.08, 0.08)}
        e5d = trial.suggest_float("Sb_e5d", 0.05, 0.70)
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
        print(f"trial {trial.number}: loss={res['loss']:.4f} "
              f"rmsEF={res['rms_ef']:.3f} rmsAll={res['rms_all']:.3f}",
              flush=True)
        return res["loss"]

    study.optimize(objective, n_trials=args.n_trials)
    t = study.best_trial
    print("BEST:", round(study.best_value, 4),
          {k: round(v, 4) for k, v in t.user_attrs.items()})
    print({k: round(v, 3) for k, v in t.params.items()})


if __name__ == "__main__":
    main()
