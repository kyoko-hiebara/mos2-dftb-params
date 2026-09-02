#!/usr/bin/env python3
"""DFTB+U (FLL) を Mo 4d に入れたときの V_S 準位・バンドの応答スキャン。使い方: python scan_u_mo.py <skf_dir> [UJ list in Ha]"""
import json, os, subprocess, sys
import numpy as np
SCRIPTS = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, SCRIPTS)
import optimize_multi as om
from compare_bands import band_loss
import extended_loss as ext
skf = os.path.abspath(sys.argv[1])
ujs = [float(x) for x in sys.argv[2:]] or [0.0, 0.02, 0.04, 0.06, 0.08, 0.10]
W = "/private/tmp/claude-501/-Users-crocus-uhuhu-MoS2-DFTB/6416faee-14e1-48c8-9838-35c06f3e7813/scratchpad/uscan_" + os.path.basename(skf)
base_defect = om.DEFECT_HSD
rows = []
for uj in ujs:
    blk = f"  OrbitalPotential = {{\n    Functional = FLL\n    Mo = {{\n      Shells = {{3}}\n      UJ = {uj:.4f}\n    }}\n  }}" if uj > 0 else ""
    wd = f"{W}/u{int(round(uj*1000)):03d}"; os.makedirs(wd, exist_ok=True)
    env = dict(os.environ, DFTB_EXTRA_HSD=blk)
    r = subprocess.run([sys.executable, f"{SCRIPTS}/dftb_bands.py", skf, "3.16", f"{wd}/dftb", "--json", f"{wd}/dftb.json",
                        "--thickness", "3.127", "--s-lmax", "d", "--extra-kpts", f"{SCRIPTS}/mesh_kpts_12x12.json"],
                       capture_output=True, text=True, env=env)
    if not os.path.exists(f"{wd}/dftb.json"):
        print(f"UJ={uj}: band step failed: {r.stdout[-200:]}"); continue
    res = band_loss(om.VASP_JSON, f"{wd}/dftb.json", verbose=False); res.update(ext.all_terms(om.VASP_JSON, f"{wd}/dftb.json"))
    blk_esc = blk.replace("{", "{{").replace("}", "}}")   # DEFECT_HSD は後で .format されるので波括弧をエスケープ
    om.DEFECT_HSD = base_defect.replace("  PolynomialRepulsive", (blk_esc + "\n" if blk else "") + "  PolynomialRepulsive")
    depth = om.defect_depth(skf, f"{wd}/defect")
    rows.append((uj, depth, res))
    print(f"UJ={uj:.3f} Ha ({uj*27.2114:.2f} eV): V_S depth={depth if depth is None else round(depth,3)} (target 0.554) | "
          f"band loss {res['loss']:.3f} gapK {res['gap_K_dftb']:.3f} dGK {res['dGK']:+.3f} dQ {res['dq_dftb']:.3f} "
          f"mVB {res['m_VB_KM_dftb']:+.2f} mCB {res['m_CB_KM_dftb']:+.2f} midgap {res['midgap_abs']:.3f} rmsVB/CB {res['rms_vb']:.3f}/{res['rms_cb']:.3f}", flush=True)
json.dump([(u, d, {k: v for k, v in r.items()}) for u, d, r in rows], open(f"{W}/scan.json", "w"), indent=1, default=float)
