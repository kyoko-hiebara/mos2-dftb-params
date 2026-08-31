# MoS2 DFTB パラメータ化プロジェクト — 実行計画と進捗

目的: pristine / S空孔 (Vs) / O置換 (O_S) の MoS2 の NEGF 輸送計算 (DFTB+ transport) 用の
高精度 Slater-Koster セット (Mo-S-O-H) を新規作成する。

## 参照レベル
- **LAK meta-GGA** (Lebeda–Aschebrock–Kümmel, PRL 133, 136402 (2024))
  - libxc 7.1.2 の `MGGA_X_LAK` (ID 342) + `MGGA_C_LAK` (ID 345)
  - VASP 6.4.2 + libxc 7.1.2 (`-DUSELIBXC`, libxc は `DISABLE_FHC=ON` でビルド)
  - INCAR: `METAGGA = LIBXC; LIBXC1 = MGGA_X_LAK; LIBXC2 = MGGA_C_LAK; LASPH = .TRUE.`
  - 検証済み: r2SCAN の native vs libxc 経由でエネルギー 3 μeV・力 1e-6 eV/Å・応力一致
  - 注: VASP 6.5+ なら `METAGGA = LAK` がネイティブ実装 (公式検証済み経路)。6.4.2+libxc は
    公式未検証の組み合わせなので r2SCAN 照合で自己検証した
- ギャップ精度: 半導体 (<4 eV) で HSE06 級、コストは半局所 (ユーザー要件: HSE06 は NG)
- POTCAR: Mo_sv / S_h / O_h / H_h (全て kinetic energy-density 対応 = meta-GGA OK)
- 収束: ENCUT = 900 eV (0.3 meV/cell)、k = 12x12x1 (0.3 µeV/cell) @ 単層プリミティブ

## 戦略 (文献調査の結論)
meta-GGA で SK 積分を直接生成できるツールは存在しない (hotcent/skprogs とも GGA まで)。
2026 年の「anti-transfer」論文 (arXiv:2608.14875) も直接継承に赤信号。
→ **実証済みの正攻法** (鉛ペロブスカイト DFTB, arXiv:2412.07016 が HSE06 超えを達成):
1. 電子部分: PBE + hotcent v2.0.1 (スカラー相対論) で SK 積分生成
2. confinement (r_dens, r_wf per shell) を **LAK バンド構造にフィット** (optuna TPE)
3. onsite 固有値・Hubbard U: PBE 自由原子から
4. 反発ポテンシャル: LAK のエネルギー・力データに CCS (ccs_fit) でフィット
5. (後段) SOC 定数: QUASINANO Part III (Jha & Heine, JCTC 2022) 流用または原子計算

## パイプライン (スクリプトは scripts/ に、リモートは /workspace/MoS2_DFTB)
1. `make_structures.py` — 単層/バルク/分子の POSCAR 生成 ✅
2. `run_conv.sh` — ENCUT/k 収束テスト (LAK) ✅
3. `run_ascan.sh` — E(a) スキャン 3.08–3.26 Å、内部緩和付き 🔄
4. `run_bands_vasp.py` — LAK 参照バンド (0-weight k 点方式、Γ-M-K-Γ 57点) ⏳
5. `gen_skf.py` — hotcent で全 10 ペア SKF 生成 (第0版: 経験則 confinement) 🔄
6. `dftb_bands.py` — DFTB+ 25.1 でバンド計算 + band.out パース
7. `compare_bands.py` — LAK↔DFTB バンド損失 (VASP 側は Mo 4s4p セミコア 4 本除外、
   VBM 整列、VB9本 + CB4本 + K点ギャップ x10 重み)
8. `optimize_confinement.py` — optuna で 7 パラメータ最適化 (粗グリッド SK、3ペア並列)
9. 反発フィット: 圧縮側 E(a) 追加 (a=2.90–3.06)、変位スナップショット、分子カーブ
   → ccs_fit (venv_ccs, numpy<1.23 制約のため専用 venv)
10. 欠陥参照: 5x5 スーパーセル Vs / O_S (PBE 緩和 → LAK 緩和 → 欠陥準位)
11. 検証: バンド/欠陥準位/フォノン? → DFTB+ NEGF スモークテスト (透過関数)

## 既知の注意点
- mpirun は stdin を吸う → リモートスクリプトはファイル化して実行 (`< /dev/null` 併用)
- /workspace は root squash (chown 不可)。tar は `--no-same-owner --no-same-permissions`
- DFTB+ 25.1 は gfortran-12 必須 (OMPI_FC=gfortran-12)
- pbe_libxc (GGA=LIBXC) は VASP 6.4.2 で動かず (原因未調査、meta-GGA 経路は問題なし)
- ccs_fit 0.22.5 は numpy<1.23 固定 → venv_ccs で分離
- hotcent SKF 出力名は `El1-El2_offsite2c.skf` → `filename_template="{el1}-{el2}.skf"` 指定
- VASP バンドは KPOINTS_OPT でなく 0-weight 方式を採用 (meta-GGA で確実)

## ベースライン比較
- PTBP (JCTC 2024, 全元素 SCC セット, Zenodo 14289468) を取得し、同じ損失関数で
  自作セットと比較する
