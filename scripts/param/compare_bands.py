#!/usr/bin/env python3
"""LAK 参照バンド (VASP) と DFTB バンドの比較・損失評価。

- VASP: Mo_sv のセミコア (4s x1 + 4p x3 = 4 バンド) を除外
- 両者とも VBM = 0 に整列
- 損失 = VB 9 本 + CB 下位数本の重み付き RMS + ギャップ誤差ペナルティ

使い方:
  python3 compare_bands.py vasp_bands.json dftb_bands.json [--plot out.png]
"""
import argparse
import json
import sys

import numpy as np

NSEMICORE = 4     # VASP 側で除外する下位バンド数
NVB = 9           # 価電子バンド本数 (DFTB: Mo d5+s1 + 2S p4s2 -> 18 el / 2)
NCB_FIT = 4       # フィットに含める伝導バンド本数
W_VB = 1.0
W_CB = 1.5
W_GAP = 10.0
W_DIR = 20.0
W_GK = 20.0   # Γ-K VBM 順序ペナルティ
# バンド別重み (輸送に効く VBM/CBM 近傍を重視)
WB_VB = [0.2, 0.2, 0.2, 0.5, 1.0, 1.0, 1.0, 1.5, 2.0]
WB_CB = [2.0, 1.5, 1.0, 0.5]

def _kw(nk):
    import numpy as _np
    w = _np.ones(nk)
    w[max(0, 30 - 0):35] = 3.0     # K 谷近傍
    w[42:47] = 2.0                 # Q 谷
    w[0:3] = 2.0                   # Γ VBM
    w[54:nk] = 2.0
    return w
# K 点インデックス (bandpath_common: N_GM=20, N_MK=12 -> K = 32)
IDX_K = 32


def load_vasp(path):
    d = json.load(open(path))
    e = np.array(d["eigs"])[:, NSEMICORE:]
    return e, d


def load_dftb(path):
    d = json.load(open(path))
    return np.array(d["eigs"]), d


def band_loss(vasp_json, dftb_json, verbose=True):
    ev, dv = load_vasp(vasp_json)
    ed, dd = load_dftb(dftb_json)
    nk = min(len(ev), len(ed))
    ev, ed = ev[:nk], ed[:nk]

    vbm_v = ev[:, :NVB].max()
    vbm_d = ed[:, :NVB].max()
    ev = ev - vbm_v
    ed = ed - vbm_d

    ncb = min(NCB_FIT, ed.shape[1] - NVB, ev.shape[1] - NVB)
    dvb = ev[:, :NVB] - ed[:, :NVB]
    dcb = ev[:, NVB:NVB + ncb] - ed[:, NVB:NVB + ncb]

    wvb = np.array(WB_VB[:NVB])
    wcb = np.array(WB_CB[:ncb])
    wk = _kw(nk)
    kav_vb = (wk[:, None] * dvb ** 2).sum(axis=0) / wk.sum()
    kav_cb = (wk[:, None] * dcb ** 2).sum(axis=0) / wk.sum()
    rms_vb = float(np.sqrt((wvb * kav_vb).sum() / wvb.sum()))
    rms_cb = float(np.sqrt((wcb * kav_cb).sum() / wcb.sum()))

    # 直接ギャップ @ K
    gap_v = float(ev[IDX_K, NVB:].min() - ev[IDX_K, :NVB].max())
    gap_d = float(ed[IDX_K, NVB:].min() - ed[IDX_K, :NVB].max())
    # 全体ギャップ
    gg_v = float(ev[:, NVB:].min() - ev[:, :NVB].max())
    gg_d = float(ed[:, NVB:].min() - ed[:, :NVB].max())

    # 直接性ペナルティ: 参照が直接 (gap==gapK) のとき DFTB も直接に誘導
    indirectness = max(0.0, gap_d - gg_d) - max(0.0, gap_v - gg_v)
    # Γ-K の VBM 順序 (K が上 = 直接性の要) を参照に合わせる
    dGK_v = float(ev[IDX_K, :NVB].max() - ev[0, :NVB].max())
    dGK_d = float(ed[IDX_K, :NVB].max() - ed[0, :NVB].max())
    loss = (W_VB * rms_vb ** 2 + W_CB * rms_cb ** 2
            + W_GAP * (gap_v - gap_d) ** 2
            + W_DIR * indirectness ** 2
            + W_GK * (dGK_v - dGK_d) ** 2)
    if verbose:
        print(f"RMS(VB {NVB})   = {rms_vb:.4f} eV")
        print(f"RMS(CB {ncb})    = {rms_cb:.4f} eV")
        print(f"K-gap: VASP {gap_v:.4f} / DFTB {gap_d:.4f} eV")
        print(f"min gap: VASP {gg_v:.4f} / DFTB {gg_d:.4f} eV")
        print(f"VBM(K)-VBM(G): VASP {dGK_v:+.4f} / DFTB {dGK_d:+.4f} eV")
        print(f"LOSS = {loss:.6f}")
    return dict(loss=loss, rms_vb=rms_vb, rms_cb=rms_cb, dGK=dGK_d,
                gap_K_vasp=gap_v, gap_K_dftb=gap_d,
                gap_vasp=gg_v, gap_dftb=gg_d)


def plot(vasp_json, dftb_json, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ev, _ = load_vasp(vasp_json)
    ed, _ = load_dftb(dftb_json)
    nk = min(len(ev), len(ed))
    ev = ev[:nk] - ev[:nk, :NVB].max()
    ed = ed[:nk] - ed[:nk, :NVB].max()
    x = np.arange(nk)
    fig, ax = plt.subplots(figsize=(6, 5))
    for ib in range(ev.shape[1]):
        ax.plot(x, ev[:, ib], "-", color="k", lw=1.2,
                label="LAK (VASP)" if ib == 0 else None)
    for ib in range(ed.shape[1]):
        ax.plot(x, ed[:, ib], "--", color="tab:red", lw=1.2,
                label="DFTB" if ib == 0 else None)
    for xx, lab in [(0, "$\\Gamma$"), (20, "M"), (32, "K"), (56, "$\\Gamma$")]:
        ax.axvline(xx, color="gray", lw=0.5)
        ax.text(xx, ax.get_ylim()[0], lab, ha="center", va="top")
    ax.set_ylim(-8, 6)
    ax.set_ylabel("E - E$_{VBM}$ (eV)")
    ax.legend(loc="upper right", fontsize=8)
    ax.set_xticks([])
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("vasp_json")
    ap.add_argument("dftb_json")
    ap.add_argument("--plot", default=None)
    args = ap.parse_args()
    band_loss(args.vasp_json, args.dftb_json)
    if args.plot:
        plot(args.vasp_json, args.dftb_json, args.plot)
