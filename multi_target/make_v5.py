#!/usr/bin/env python3
"""skf_v5: optm2 (Mo/S, 整列+質量+メッシュ+V_S) + opto2 (O, 共通シフト再アンカー) + optsb3 (Sb, SOC ループ内・整列)
を密グリッドで 4 元素 16 ペア生成し、全ターゲットを検証する。
使い方: python3 make_v5.py [--out skf_v5] [--mos2-study optm2] [--sb-study optsb3] [--skip-gen]"""
import argparse
import json
import os
import shutil
import subprocess
import sys
from multiprocessing import Pool

import numpy as np

SCRIPTS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPTS)

import optimize_confinement as oc  # noqa: E402
_default_make_atom = oc.make_atom
import optimize_sb2 as sb2  # noqa: E402  (SPECS["Sb"], FREE["Sb"], make_atom -> tuned)
_tuned_make_atom = sb2._make_atom_tuned
import optimize_multi as om  # noqa: E402
import optimize_sb3 as sb3  # noqa: E402
from compare_bands import band_loss  # noqa: E402
import extended_loss as ext  # noqa: E402
from ase.io import read, write as ase_write  # noqa: E402

ROOT = "/Users/crocus/uhuhu/MoS2_DFTB"
oc.SPECS["O"] = dict(configuration="[He] 2s2 2p4", valence=["2s", "2p"],
                     scalarrel=False, occupations={"2s": 2, "2p": 4})
oc.RMIN, oc.DR, oc.NPTS = 0.4, 0.02, 980
oc.NTHETA, oc.NR = 150, 50


def _dispatch_make_atom(el, p):
    return _tuned_make_atom(el, p) if el == "Sb" else _default_make_atom(el, p)


oc.make_atom = _dispatch_make_atom
V3_SHIFT_MEAN = np.mean([0.05988774, 0.03579509, 0.03870145, 0.02528376, 0.02723118])


def best_params(study, tag):
    """study はカンマ区切り可: 最良値が最小の study のベスト試行を採用"""
    import optuna
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    best = None
    for name in study.split(","):
        path = f"{ROOT}/local_opt/{name}/study.db"
        if not os.path.exists(path):
            continue
        st = optuna.load_study(study_name=name, storage=f"sqlite:///{path}")
        try:
            t = st.best_trial
        except ValueError:
            continue
        if best is None or t.value < best[1].value:
            best = (name, t)
    name, t = best
    print(f"{tag}: best trial {t.number} of {name}, value {t.value:.4f}", flush=True)
    return dict(t.params), dict(t.user_attrs)


def collect(mos2_study, sb_study, o_study="opto2"):
    fp, ap = best_params(mos2_study, "Mo/S")
    fo, _ = best_params(o_study, "O")
    fs, asb = best_params(sb_study, "Sb")
    pMo = {"r_dens": fp["Mo_rd"], "s_wf": fp["Mo_swf"],
           "r_wf": {"4d": fp["Mo_r4d"], "5s": fp["Mo_r5s"], "5p": fp["Mo_r5p"]}}
    if "Mo_s4d" in fp:
        pMo["s_wf_shell"] = {"4d": fp["Mo_s4d"]}
    pS = {"r_dens": fp["S_rd"], "s_wf": fp["S_swf"],
          "r_wf": {"3s": fp["S_r3s"], "3p": fp["S_r3p"], "3d": fp["S_r3d"]}}
    pO = {"r_dens": fo["O_rd"], "s_wf": fo["O_swf"], "r_wf": {"2s": fo["O_r2s"], "2p": fo["O_r2p"]}}
    pSb = {"r_dens": fs["Sb_rd"], "s_wf": fs["Sb_swf"],
           "r_wf": {"5s": fs["Sb_r5s"], "5p": fs["Sb_r5p"], "5d": fs["Sb_r5d"]}}
    new_mean = np.mean([fp["Mo_sh4d"], fp["Mo_sh5s"], fp["Mo_sh5p"], fp["S_sh3s"], fp["S_sh3p"]])
    d_o = float(new_mean - V3_SHIFT_MEAN) if o_study.startswith("opto2") else 0.0
    shifts = {"Mo": {"4d": fp["Mo_sh4d"], "5s": fp["Mo_sh5s"], "5p": fp["Mo_sh5p"]},
              "S": {"3s": fp["S_sh3s"], "3p": fp["S_sh3p"]},
              "O": {"2s": fo["O_sh2s"] + d_o, "2p": fo["O_sh2p"] + d_o},
              "Sb": {"5s": fs["Sb_sh5s"], "5p": fs["Sb_sh5p"]}}
    print(f"O shifts re-anchored by {d_o:+.4f} Ha (Mo/S mean shift {V3_SHIFT_MEAN:.4f} -> {new_mean:.4f})")
    return dict(pMo=pMo, pS=pS, pO=pO, pSb=pSb, e3d=fp["S_e3d"], e5d=fs["Sb_e5d"], shifts=shifts,
                raw=dict(mos2=fp, o=fo, sb=fs), attrs=dict(mos2=ap, sb=asb))


def gen_all(P, out):
    os.makedirs(out, exist_ok=True)
    oc.FREE["Sb"]["eig"]["5d"] = P["e5d"]
    oc.FREE["Sb"]["U"]["5d"] = oc.FREE["Sb"]["U"]["5p"]
    sh, e3d = P["shifts"], P["e3d"]
    tasks = [("Mo", "Mo", P["pMo"], P["pMo"]), ("Mo", "S", P["pMo"], P["pS"]), ("S", "S", P["pS"], P["pS"]),
             ("Mo", "O", P["pMo"], P["pO"]), ("S", "O", P["pS"], P["pO"]), ("O", "O", P["pO"], P["pO"]),
             ("Sb", "Sb", P["pSb"], P["pSb"]), ("Mo", "Sb", P["pMo"], P["pSb"]),
             ("S", "Sb", P["pS"], P["pSb"]), ("O", "Sb", P["pO"], P["pSb"])]
    tasks = [(a, b, pa, pb, out, e3d, sh) for a, b, pa, pb in tasks]
    pool_tasks = [t for t in tasks if "Sb" not in (t[0], t[1])]
    sb_tasks = [t for t in tasks if "Sb" in (t[0], t[1])]
    with Pool(5) as pool:
        for name in pool.imap_unordered(oc.gen_pair, pool_tasks):
            print("  pair", name, "done", flush=True)
    for t in sb_tasks:      # spawn ワーカーでは oc.FREE["Sb"]["eig"]["5d"] が失われるため直列
        print("  pair", oc.gen_pair(t), "done (serial)", flush=True)
    l2 = open(os.path.join(out, "Sb-Sb.skf")).readlines()[1].split()
    assert abs(float(l2[0]) - P["e5d"]) < 1e-5, f"Sb 5d onsite not written: {l2[:3]}"
    print("SKF files:", sorted(f for f in os.listdir(out) if f.endswith(".skf")))


SUBO_HSD = """Geometry = GenFormat {{
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


def validate_subo(skf, wd):
    import optimize_o as oo   # LAK 参照 (gap, O2s, VB 窓)
    os.makedirs(wd, exist_ok=True)
    ase_write(os.path.join(wd, "geo.gen"), read(oo.SUBO_POSCAR), format="gen")
    open(os.path.join(wd, "dftb_in.hsd"), "w").write(SUBO_HSD.format(skf=skf))
    r = subprocess.run([om.DFTB], cwd=wd, capture_output=True, text=True, timeout=1800,
                       env=dict(os.environ, OMP_NUM_THREADS="8"))
    if r.returncode != 0:
        return {"error": r.stdout[-200:]}
    eig, occ = om.parse_tag_eigs(os.path.join(wd, "results.tag"))
    e = (eig[0] if eig.ndim == 3 else eig) * om.HA
    o = (occ[0] if occ.ndim == 3 else occ)
    nocc = int((o[0] > 0.5).sum())
    vbm, cbm = e[o > 0.5].max(), e[o < 0.5].min()
    emeta = sorted((e[:, ib].mean(), e[:, ib].max() - e[:, ib].min()) for ib in range(nocc, min(nocc + 8, e.shape[1])))
    ingap = 0.0
    for i in range(len(emeta) - 1):
        if emeta[i + 1][0] - emeta[i][0] > 0.25:
            if all(w < 0.08 for _, w in emeta[:i + 1]):
                ingap = 1.0
            break
    occ_means = np.array([e[:, ib].mean() for ib in range(nocc)])
    order = np.argsort(occ_means)
    o2s = float(occ_means[order[0]] - vbm) if occ_means[order[1]] - occ_means[order[0]] > 1.0 else None
    bm = occ_means - vbm
    win = np.sort(bm[bm > -8.0])
    n = min(len(win), len(oo.REF["vb"]))
    vb_rms = float(np.sqrt(((win[-n:] - oo.REF["vb"][-n:]) ** 2).mean()))
    return dict(gap=float(cbm - vbm), gap_ref=oo.REF["gap"], o2s=o2s, o2s_ref=oo.REF["o2s"],
                vb_rms=vb_rms, ingap=ingap)


def validate_sb2s3(skf, wd):
    ref_path = f"{ROOT}/local_opt/sb2s3_ref_bands.json"
    traj = f"{ROOT}/local_opt/sb2s3_relaxed.traj"
    if not (os.path.exists(ref_path) and os.path.exists(traj)):
        return {"skipped": "no GPAW reference yet"}
    ref = json.load(open(ref_path))
    kpts = np.array(ref["kpts"]); er = np.array(ref["eigs"]); nocc = ref["nocc"]
    os.makedirs(wd, exist_ok=True)
    ase_write(os.path.join(wd, "geo.gen"), read(traj), format="gen")
    hsd = """Geometry = GenFormat {{
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
    S = "d"
    Sb = "d"
  }}
  PolynomialRepulsive = SetForAll {{ Yes }}
  Filling = Fermi {{ Temperature [K] = 100 }}
  {extra}
{kblock}
}}
Options {{ WriteResultsTag = Yes }}
ParserOptions {{ ParserVersion = 14 }}
"""
    kscc = "  KPointsAndWeights = SupercellFolding {\n    2 0 0\n    0 6 0\n    0 0 2\n    0.5 0.5 0.5\n  }"
    env = dict(os.environ, OMP_NUM_THREADS="8")
    open(os.path.join(wd, "dftb_in.hsd"), "w").write(hsd.format(skf=skf, kblock=kscc, tol="1e-5", mx="300", extra=""))
    r = subprocess.run([om.DFTB], cwd=wd, capture_output=True, text=True, timeout=1800, env=env)
    if r.returncode != 0:
        return {"error": "SCC " + r.stdout[-200:]}
    kl = "  KPointsAndWeights = {\n" + "\n".join(f"    {k[0]:.8f} {k[1]:.8f} {k[2]:.8f} 1.0" for k in kpts) + "\n  }"
    open(os.path.join(wd, "dftb_in.hsd"), "w").write(hsd.format(skf=skf, kblock=kl, tol="1e6", mx="1",
                                                                  extra="ReadInitialCharges = Yes"))
    r = subprocess.run([om.DFTB], cwd=wd, capture_output=True, text=True, timeout=1800, env=env)
    if r.returncode != 0:
        return {"error": "bands " + r.stdout[-200:]}
    eig, _ = om.parse_tag_eigs(os.path.join(wd, "results.tag"))
    ed = np.sort(np.squeeze(eig) * om.HA, axis=1)
    if ed.shape[0] != len(kpts):
        ed = ed.T
    # DFTB 占有数: S 6 + Sb 5 (4d はセミコアで DFTB に無し) -> 12*6+8*5 = 112 e -> 56 バンド
    nocc_d = (12 * 6 + 8 * 5) // 2
    er_al = er - er[:, :nocc].max(); ed_al = ed - ed[:, :nocc_d].max()
    nvb, ncb = 8, 4
    dv = er_al[:, nocc - nvb:nocc] - ed_al[:, nocc_d - nvb:nocc_d]
    dc = er_al[:, nocc:nocc + ncb] - ed_al[:, nocc_d:nocc_d + ncb]
    gap_r = float(er_al[:, nocc:].min() - er_al[:, :nocc].max())
    gap_d = float(ed_al[:, nocc_d:].min() - ed_al[:, :nocc_d].max())
    np.save(os.path.join(wd, "dftb_bands.npy"), ed)
    return dict(gap_ref=gap_r, gap_dftb=gap_d, rms_vb8=float(np.sqrt((dv ** 2).mean())),
                rms_cb4=float(np.sqrt((dc ** 2).mean())))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="skf_v5")
    ap.add_argument("--mos2-study", default="optm2")
    ap.add_argument("--sb-study", default="optsb3")
    ap.add_argument("--o-study", default="opto2")
    ap.add_argument("--skip-gen", action="store_true")
    args = ap.parse_args()
    OUT = f"{ROOT}/local_opt/{args.out}"
    P = collect(args.mos2_study, args.sb_study, args.o_study)
    if not args.skip_gen:
        gen_all(P, OUT)
    wd = f"{ROOT}/local_opt/{args.out}_eval"
    os.makedirs(wd, exist_ok=True)
    M = {"params": P["raw"], "shifts": P["shifts"], "e3d": P["e3d"], "e5d": P["e5d"]}

    # --- MoS2 bands / masses / mesh / alignment ---
    r = subprocess.run([sys.executable, f"{SCRIPTS}/dftb_bands.py", OUT, "3.16", f"{wd}/dftb", "--json",
                        f"{wd}/dftb.json", "--thickness", "3.127", "--s-lmax", "d",
                        "--extra-kpts", f"{SCRIPTS}/mesh_kpts_12x12.json"], capture_output=True, text=True)
    res = band_loss(om.VASP_JSON, f"{wd}/dftb.json", verbose=False)
    res.update(ext.all_terms(om.VASP_JSON, f"{wd}/dftb.json"))
    M["mos2"] = res
    print(f"\n[MoS2] band loss {res['loss']:.4f} (v3 0.2514) rmsVB {res['rms_vb']:.3f} rmsCB {res['rms_cb']:.3f} "
          f"gapK {res['gap_K_dftb']:.3f}/{res['gap_K_vasp']:.3f} dGK {res['dGK']:+.3f}")
    print(f"       m*VB K->M {res['m_VB_KM_dftb']:.2f}/{res['m_VB_KM_lak']:.2f}  K->G {res['m_VB_KG_dftb']:.2f}/{res['m_VB_KG_lak']:.2f}  "
          f"m*CB K->M {res['m_CB_KM_dftb']:.2f}/{res['m_CB_KM_lak']:.2f}  K->G {res['m_CB_KG_dftb']:.2f}/{res['m_CB_KG_lak']:.2f}  "
          f"m*VB(G) {res['m_VBG_dftb']:.2f}/{res['m_VBG_lak']:.2f}")
    print(f"       Q-K {res['dq_dftb']:.3f}/{res['dq_lak']:.3f}  mesh RMS {res['rms_mesh']:.3f}  "
          f"midgap {res['midgap_abs']:.3f} (target {ext.align_target():.3f}, err {res['align_err']:+.3f})")
    depth = om.defect_depth(OUT, f"{wd}/defect")
    M["vs_depth"] = depth
    print(f"[V_S]  depth {depth if depth is None else round(depth, 4)} (target {om.TARGET_DEPTH:.3f}, v3 0.690)")
    # --- sub_O ---
    so = validate_subo(OUT, f"{wd}/subo")
    M["sub_o"] = so
    print(f"[O_S]  {so}")
    # --- Sb bulk (SOC 込み, 整列) ---
    ef, e_path, e_mesh = sb3.bulk_soc(OUT, f"{wd}/sb_bulk")
    rms_all, rms_ef = sb3.path_loss(e_path)
    rms_mesh, ov, kh, ke = sb3.mesh_loss(e_mesh)
    efs = sb3.slab_ef(OUT, f"{wd}/sb_slab")
    e_noso = sb2.sb_bands(OUT, f"{wd}/sb_noso")
    rms_all0, rms_ef0 = sb2.band_loss_sb(e_noso)
    tgt = sb3.slab_target() if os.path.exists(sb3.SLAB_JSON) else None
    M["sb"] = dict(rms_ef_soc=rms_ef, rms_all_soc=rms_all, rms_mesh_soc=rms_mesh, overlap=ov, overlap_ref=sb3.OV_REF,
                   k_hole=kh.tolist(), k_elec=ke.tolist(), ef_bulk=ef, ef_slab=efs, ef_slab_target=tgt,
                   rms_ef_noso=rms_ef0, rms_all_noso=rms_all0)
    print(f"[Sb]   SOC path rmsEF {rms_ef:.3f} (v3sb 0.327)  mesh {rms_mesh:.3f}  overlap {ov*1000:+.0f}/{sb3.OV_REF*1000:+.0f} meV  "
          f"hole k {np.round(kh,3)} elec k {np.round(ke,3)}")
    print(f"       no-SOC path rmsEF {rms_ef0:.3f} (v3sb 0.165)  E_F bulk {ef:.3f} slab {efs:.3f} (target {tgt})")
    # --- Sb2S3 ---
    s23 = validate_sb2s3(OUT, f"{wd}/sb2s3")
    M["sb2s3"] = s23
    print(f"[Sb2S3] {s23}")
    json.dump(M, open(f"{OUT}/v5_metrics.json", "w"), indent=1, default=float)
    print("metrics ->", f"{OUT}/v5_metrics.json")


if __name__ == "__main__":
    main()
