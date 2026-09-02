#!/usr/bin/env python3
"""バンド損失の拡張項 (2026-09-02 診断に基づく):
  - K 点 VB/CB と Γ 点 VB のバンド端曲率 (有効質量) の相対残差
  - Q 谷 (K-Γ 上の CB 極小) と K CBM の間隔
  - 12x12 SCF メッシュ (IBZ 19 点、EIGENVAL に同梱) の VB 上位 3 本 + CB 下位 2 本
  - 絶対エネルギー整列: K 点ミッドギャップの真空準位基準 (PBE/GPAW ターゲット)
"""
import json
import os

import numpy as np

ROOT = "/Users/crocus/uhuhu/MoS2_DFTB"
EIGENVAL = f"{ROOT}/results/ref_calc/bands_lak_a316/EIGENVAL"
ALIGN_JSON = f"{ROOT}/local_opt/mos2_pbe_ipea.json"
NSEMI, NVB, IDX_K = 4, 9, 32
NMESH = 19
HB = 3.80998  # eV A^2 = hbar^2/2m_e
DK_KM = DK_KG = 0.05516   # 1/A (bandpath_common, a=3.16)
DK_GM = 0.05740
W_M = float(os.environ.get("EXT_W_M", 0.5))
W_Q = float(os.environ.get("EXT_W_Q", 10.0))
W_MESH = float(os.environ.get("EXT_W_MESH", 1.0))
W_AL = float(os.environ.get("EXT_W_AL", 2.0))
Q_RANGE = (40, 49)


def load_vasp_mesh(path=EIGENVAL, nmesh=NMESH):
    lines = open(path).readlines()
    nelec, nk, nb = (int(x) for x in lines[5].split())
    idx = 7
    kpts, wts, eigs = [], [], []
    for ik in range(nk):
        t = lines[idx].split()
        kpts.append([float(t[0]), float(t[1]), float(t[2])])
        wts.append(float(t[3]))
        eigs.append([float(lines[idx + 1 + ib].split()[1]) for ib in range(nb)])
        idx += nb + 2
    return (np.array(kpts)[:nmesh], np.array(wts)[:nmesh],
            np.array(eigs)[:nmesh, NSEMI:])


MESH_K, MESH_W, MESH_E = load_vasp_mesh()


def write_mesh_kpts(path):
    json.dump(MESH_K.tolist(), open(path, "w"))
    return path


def align_target():
    if os.environ.get("MOS2_ALIGN_TARGET"):
        return float(os.environ["MOS2_ALIGN_TARGET"])
    if not os.path.exists(ALIGN_JSON):
        raise FileNotFoundError(f"alignment target missing: {ALIGN_JSON}")
    d = json.load(open(ALIGN_JSON))
    return float(d["midgap_abs"])


def _curv(e, ik, band, side, i):
    return e[ik + side * i, band] - e[ik, band]


def edge_terms(ev, ed):
    """ev, ed: パス固有値 (57, nb)、各自の VBM 整列済み。相対残差と質量を返す。"""
    res, masses = [], {}
    for band, lab in [(NVB - 1, "VB"), (NVB, "CB")]:
        for side, dlab in [(-1, "KM"), (1, "KG")]:
            for i in (1, 2):
                dv, dd = _curv(ev, IDX_K, band, side, i), _curv(ed, IDX_K, band, side, i)
                res.append((dd - dv) / abs(dv))
            k2 = (2 * DK_KM) ** 2
            masses[f"m_{lab}_{dlab}_lak"] = HB * k2 / _curv(ev, IDX_K, band, side, 2)
            masses[f"m_{lab}_{dlab}_dftb"] = HB * k2 / _curv(ed, IDX_K, band, side, 2)
    for i in (1, 2):
        dv, dd = _curv(ev, 0, NVB - 1, 1, i), _curv(ed, 0, NVB - 1, 1, i)
        res.append(0.7 * (dd - dv) / abs(dv))
    masses["m_VBG_lak"] = HB * (2 * DK_GM) ** 2 / _curv(ev, 0, NVB - 1, 1, 2)
    masses["m_VBG_dftb"] = HB * (2 * DK_GM) ** 2 / _curv(ed, 0, NVB - 1, 1, 2)
    res = np.array(res)
    return float((res ** 2).mean()), masses


def q_terms(ev, ed):
    a, b = Q_RANGE
    dqv = float(ev[a:b, NVB].min() - ev[IDX_K, NVB])
    dqd = float(ed[a:b, NVB].min() - ed[IDX_K, NVB])
    return (dqd - dqv) ** 2, dqv, dqd


def mesh_terms(ed_mesh, vbm_d):
    ev = MESH_E - MESH_E[:, :NVB].max()
    ed = ed_mesh - vbm_d
    sel = [NVB - 3, NVB - 2, NVB - 1, NVB, NVB + 1]
    bw = np.array([0.5, 1.0, 2.0, 2.0, 1.0])
    d = ev[:, sel] - ed[:, sel]
    w = MESH_W[:, None] * bw[None, :]
    return float(np.sqrt((w * d ** 2).sum() / w.sum()))


def all_terms(vasp_json, dftb_json, target_midgap=None):
    ev = np.array(json.load(open(vasp_json))["eigs"])[:, NSEMI:]
    dd = json.load(open(dftb_json))
    ed_all = np.array(dd["eigs"])
    npath = len(ev)
    ed = ed_all[:npath]
    vbm_v, vbm_d = ev[:, :NVB].max(), ed[:, :NVB].max()
    evr, edr = ev - vbm_v, ed - vbm_d
    loss_m, masses = edge_terms(evr, edr)
    loss_q, dqv, dqd = q_terms(evr, edr)
    out = dict(loss_m=loss_m, loss_q=loss_q, dq_lak=dqv, dq_dftb=dqd, **masses)
    if ed_all.shape[0] >= npath + NMESH:
        out["rms_mesh"] = mesh_terms(ed_all[npath:npath + NMESH], vbm_d)
    else:
        out["rms_mesh"] = 0.0
    midgap = 0.5 * (ed[IDX_K, NVB - 1] + ed[IDX_K, NVB])
    out["midgap_abs"] = float(midgap)
    if target_midgap is None:
        target_midgap = align_target()
    out["align_err"] = float(midgap - target_midgap)
    out["loss_ext"] = (W_M * loss_m + W_Q * loss_q + W_MESH * out["rms_mesh"] ** 2
                       + W_AL * out["align_err"] ** 2)
    return out


if __name__ == "__main__":
    import sys
    tgt = float(sys.argv[3]) if len(sys.argv) > 3 else None
    r = all_terms(sys.argv[1], sys.argv[2], tgt)
    for k, v in r.items():
        print(f"{k:16s} {v:+.4f}")
