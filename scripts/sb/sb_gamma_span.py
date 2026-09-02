#!/usr/bin/env python3
"""Sb SOC 定数の Γ スパン較正: Γ 点の価電子 5p 多重項 (最高占有 6 状態) のスパンを PBE+SOC 参照に合わせる。
使い方: python sb_gamma_span.py <skf_dir>"""
import json, os, sys
import numpy as np
SCRIPTS = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, SCRIPTS)
import optimize_sb3 as sb3
skf = os.path.abspath(sys.argv[1])
K = sb3.REF_SOC_K; E = sb3.REF_SOC_E
ig = int(np.argmin(np.linalg.norm(K, axis=1)))
ref = np.sort(E[ig])
refv = ref[20:36]   # 4d 20 状態を除いた価電子 16 状態 (占有 10 + 空 6)
span_ref = refv[9] - refv[4]
print(f"reference Γ (k={K[ig]}): top-6 occupied span = {span_ref:.4f} eV; levels {np.round(refv[4:10],3)}")
wd = "/private/tmp/claude-501/-Users-crocus-uhuhu-MoS2-DFTB/6416faee-14e1-48c8-9838-35c06f3e7813/scratchpad/sbspan"
res = []
for xi in [0.50, 0.57, 0.65]:
    sb3.XI_5P = xi
    ef, e_path, e_mesh = sb3.bulk_soc(skf, f"{wd}/xi{int(xi*100)}")
    d = np.sort(e_path[ig]); dv = d[:10]  # DFTB は 36 状態、下位 10 が占有
    span = dv[9] - dv[4]
    ra, re = sb3.path_loss(e_path)
    res.append((xi, span, re))
    print(f"xi={xi:.2f}: Γ span {span:.4f} eV (ref {span_ref:.4f}), levels {np.round(dv[4:10],3)}, rms_EF {re:.3f}")
xs = np.array([r[0] for r in res]); sp = np.array([r[1] for r in res])
xi_opt = float(np.interp(span_ref, sp, xs)) if np.all(np.diff(sp) > 0) else float(xs[np.argmin(abs(sp - span_ref))])
print(f"Γ-span calibrated xi_5p = {xi_opt:.3f} eV")
json.dump({"xi_5p": xi_opt, "span_ref": span_ref, "scan": res}, open(f"{skf}/sb_soc_gamma_span.json", "w"), indent=1)
