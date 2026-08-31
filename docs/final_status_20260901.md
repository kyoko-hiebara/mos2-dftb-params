# MoS2 DFTB パラメータ化 — セッション最終状態 (2026-09-01)

成果物はすべてローカル `results/` に取得済み。リモート `/workspace/MoS2_DFTB` (RunPod 永続ボリューム) にも全データ+ビルド済みバイナリのバックアップ (`sw/bin_backup/`) が残っている。

## SKF セット (results/dftb/)

| セット | 内容 | 用途 |
|---|---|---|
| `skf_v2` | **電子部分の最終版** (Mo/S、spd 基底、16 変数最適化済み) | バンド計算・NEGF (反発不要の固定構造計算) |
| `skf_v2rep` | v2 + CCS 反発スプライン第 1 弾 | 構造緩和・エネルギー計算 (平衡格子 +1.1% 誤差あり、第 2 弾で改善予定) |
| `skf_v2o` | v2 + O/H 系ペア (暫定 confinement) | O 置換欠陥入り計算 (O_S 無準位性は検証済み) |
| `skf_v0` | 全 10 ペア (Mo,S,O,H) 初期版 | H 端終端など (未検証) |
| `skf_ptbp` | PTBP ベースライン (比較用) | — |

DFTB+ 設定: MaxAngularMomentum は Mo="d", S="d", O="p", H="s"。反発無しで使う場合は
`PolynomialRepulsive = SetForAll { Yes }`。

## 達成した精度 (LAK meta-GGA 参照、a=3.16 実験格子)

- バンド: 重み付き VB RMS 0.28 eV / CB RMS 0.29 eV、**K-K 直接ギャップ 1.933 vs 1.914 eV** (誤差 0.02)
- VBM 順序 (K-Γ): +0.019 vs +0.015 eV ✓、Q-K 谷間隔: 221 vs 242 meV (誤差 21 meV)
- PTBP 比較: 総合損失 12 分の 1 (4.52 → 0.37 → v2 で 0.21)、ギャップ誤差 30 分の 1
- **V_S 欠陥**: 二重縮退ギャップ内準位を再現 (分散 8 vs 20 meV)。位置は CBM−0.80 vs 参照 CBM−0.56 eV (0.24 eV 深い)
- **O_S 欠陥**: ギャップ内無準位 (電子的良性) を完全再現
- **SOC 参照**: LAK+SOC K 点 VB 分裂 = 150.2 meV (実験 145–150 meV) — DFTB SOC 定数の較正ターゲット確定

## 参照データ (results/ref_calc/)

- `bands_lak_a316/bands_lak.json` — フィットターゲットの単層バンド (Γ-M-K-Γ 57 点)
- `ascan/` E(a) 16 点 (a=2.90–3.40, 内部緩和済み) / `tscan/` 層厚 4 点 / `snapshots/` 変位構造 3 点 (力データ)
- `molecules_lak/` S2 (6 点), SO (5), O2 (4) 解離カーブ (スピン偏極)
- `defects/{vac_S,sub_O}_lak/` 5x5 スーパーセル LAK 静的 (EIGENVAL に準位)
- `soc_{lak,pbe}/` SOC 計算
- `repfit/` CCS フィット一式 (CCS_params.json, structures.json, ASE DB)

## 環境の再構築 (次回セッション)

1. 接続: `ssh -p <PORT> -i ~/.ssh/id_runpod_deploy root@<IP>` (鍵は id_runpod_deploy!)
2. `/workspace/MoS2_DFTB/sw/bin_backup/` に vasp_std_libxc / vasp_ncl / vasp_std_gpu / libxc / dftbplus のバイナリ・tgz あり。/opt に展開すれば数分で復元
3. **重要: RunPod A100 ポッドは cgroup CPU quota ≈ 26 コア** (nproc=252 は見かけ)。VASP は 24 ランク単独ジョブで運用すること。並列ジョブ多重投入は厳禁
4. LAK は libxc 経由のみ (GPU ポート不可)。PBE 系は GPU 版が 18 倍速

## 残タスク (優先順)

1. **マルチターゲット第 3 ラウンド**: V_S 準位 (CBM−0.56 eV) を損失に追加して 0.24 eV のずれを解消。
   2 段階スクリーニング + 4x4 セルで ~150 trials、26 コアで 4–6 時間 (64 vCPU なら ~2h)
2. **反発第 2 弾**: snapshots (力データ) を加えた ccs_fit 再実行 → 平衡格子 +1.1% を改善
   (`repfit_round2.sh` が雛形、pairset2 プレフィックス)
3. **O 系 confinement 最適化**: O_S の VB 構造 + SO/O2 カーブをターゲットに
4. **SOC 定数較正**: DFTB+ `SpinOrbit` ブロックの ξ_Mo(4d) 等を K 分裂 150 meV に合わせる
5. **NEGF スモークテスト**: 過去プロジェクト `PAW_NEGF/MoS2_comparison/DFTB-NEGF/runs/` の
   libNEGF 2 端子構成 (contact/transport hsd) を skf_v2rep で再実行 → PTBP 時代の T(E) と比較
6. H 系の検証 (エッジ終端リボン)
