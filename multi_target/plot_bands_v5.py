#!/usr/bin/env python3
"""LAK / v3 / v5 バンド比較図 (notes/bands_v5_vs_lak.png)。使い方: python plot_bands_v5.py <v5_dftb.json> [out.png]"""
import json, sys, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
NS, NVB = 4, 9
ev = np.array(json.load(open("../results/ref_calc/bands_lak_a316/bands_lak.json"))["eigs"])[:, NS:]
e3 = np.array(json.load(open("v3_eval/dftb.json"))["eigs"])[:57]
e5 = np.array(json.load(open(sys.argv[1]))["eigs"])[:57]
out = sys.argv[2] if len(sys.argv) > 2 else "../notes/bands_v5_vs_lak.png"
al = lambda e: e - e[:, :NVB].max()
ev, e3, e5 = al(ev), al(e3), al(e5)
HB, dk = 3.80998, 0.05516
m = lambda e, b, s: HB * (2 * dk) ** 2 / (e[32 + 2 * s, b] - e[32, b])
x = np.arange(57)
fig, axs = plt.subplots(1, 2, figsize=(11, 5), gridspec_kw={"width_ratios": [1.3, 1]})
ax = axs[0]
for ib in range(ev.shape[1]): ax.plot(x, ev[:, ib], "-", color="k", lw=1.3, label="LAK" if ib == 0 else None)
for ib in range(e3.shape[1]): ax.plot(x, e3[:, ib], "--", color="tab:blue", lw=1.0, label="skf_v3" if ib == 0 else None)
for ib in range(e5.shape[1]): ax.plot(x, e5[:, ib], "-", color="tab:red", lw=1.0, alpha=0.85, label="skf_v5" if ib == 0 else None)
for xx, lab in [(0, "Γ"), (20, "M"), (32, "K"), (56, "Γ")]:
    ax.axvline(xx, color="gray", lw=0.5); ax.text(xx, -7.9, lab, ha="center", va="bottom")
ax.set_ylim(-8, 5); ax.set_xticks([]); ax.set_ylabel("E − E$_{VBM}$ (eV)"); ax.legend(loc="upper right", fontsize=9)
ax.set_title("Monolayer MoS$_2$ @ a = 3.16 Å (LAK reference)")
ax = axs[1]; sl = slice(24, 50)
for e, c, ls, lab in [(ev, "k", "-", "LAK"), (e3, "tab:blue", "--", "v3"), (e5, "tab:red", "-", "v5")]:
    for ib in [NVB - 1, NVB]:
        ax.plot(x[sl], e[sl, ib], ls, color=c, lw=1.4, label=lab if ib == NVB - 1 else None)
ax.axvline(32, color="gray", lw=0.5)
for xx, lab in [(32, "K"), (44, "Q"), (24, "←M")]: ax.text(xx, -0.97, lab, ha="center", fontsize=9)
ax.set_ylim(-1.0, 2.6); ax.set_xticks([]); ax.set_title("band edges (M–K–Γ)"); ax.legend(fontsize=9, loc="center right")
dq = lambda e: e[40:49, NVB].min() - e[32, NVB]
txt = (f"m*(K, K→M) VB / CB  [m_e]\n  LAK {m(ev,8,-1):+.2f} / {m(ev,9,-1):+.2f}\n  v3  {m(e3,8,-1):+.2f} / {m(e3,9,-1):+.2f}\n  v5  {m(e5,8,-1):+.2f} / {m(e5,9,-1):+.2f}\n"
       f"Q−K [eV]: LAK {dq(ev):.2f}, v3 {dq(e3):.2f}, v5 {dq(e5):.2f}")
ax.text(0.03, 0.62, txt, transform=ax.transAxes, va="top", fontsize=8.5, family="monospace",
        bbox=dict(boxstyle="round", fc="white", ec="gray", alpha=0.9))
fig.tight_layout(); fig.savefig(out, dpi=150); print("wrote", out)
