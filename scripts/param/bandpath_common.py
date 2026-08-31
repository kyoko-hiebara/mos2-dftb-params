#!/usr/bin/env python3
"""単層 MoS2 の共通バンドパス (Γ-M-K-Γ) 定義。

VASP (KPOINTS_OPT) と DFTB+ (KPointsAndWeights) の両方に
同一の明示 k 点列を与えるためのユーティリティ。
"""
import numpy as np

# 分割数 (セグメントごと)
N_GM, N_MK, N_KG = 20, 12, 24


def kpath_points():
    """Γ-M-K-Γ の分数座標列 (重複端点なし、始点 Γ・終点 Γ 含む)"""
    G = np.array([0.0, 0.0, 0.0])
    M = np.array([0.5, 0.0, 0.0])
    K = np.array([1.0 / 3.0, 1.0 / 3.0, 0.0])
    pts = []
    for start, end, n in [(G, M, N_GM), (M, K, N_MK), (K, G, N_KG)]:
        for i in range(n):
            pts.append(start + (end - start) * i / n)
    pts.append(G)
    return np.array(pts)


def label_indices():
    """高対称点のインデックス {label: index}"""
    return {"G": 0, "M": N_GM, "K": N_GM + N_MK, "G2": N_GM + N_MK + N_KG}


def write_vasp_kpoints_opt(path):
    pts = kpath_points()
    with open(path, "w") as f:
        f.write("bandpath G-M-K-G explicit\n")
        f.write(f"{len(pts)}\n")
        f.write("Reciprocal\n")
        for p in pts:
            f.write(f"{p[0]:.10f} {p[1]:.10f} {p[2]:.10f} 1.0\n")


def dftb_klines_block():
    """DFTB+ の KPointsAndWeights 明示リスト (重み 1.0)"""
    pts = kpath_points()
    lines = ["  KPointsAndWeights = {"]
    for p in pts:
        lines.append(f"    {p[0]:.10f} {p[1]:.10f} {p[2]:.10f} 1.0")
    lines.append("  }")
    return "\n".join(lines)


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "vasp":
        write_vasp_kpoints_opt(sys.argv[2])
        print(f"wrote {sys.argv[2]} ({len(kpath_points())} kpts)")
    else:
        print(dftb_klines_block())
