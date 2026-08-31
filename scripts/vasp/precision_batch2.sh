#!/bin/bash
# 精度バッチ v2: 現実的なコスト設定で再投入
# - snapshots: ENCUT 600, k 3x3, 8 構造 x 24 ランク
# - vac_S / sub_O: LAK 静的 1 点 (ENCUT 700, 96 ランク, 直列)
set -u
cd /workspace/MoS2_DFTB/ref_calc
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
VASP=/root/vasp.6.4.2/bin/vasp_std

# ---------- snapshots (8 構造) ----------
for i in 00 01 02 03 04 05 06 07; do
  d=snapshots/snap_$i
  sed -i "s/^ENCUT.*/ENCUT = 600/" $d/INCAR
  printf "s\n0\nGamma\n3 3 1\n0 0 0\n" > $d/KPOINTS
  rm -f $d/WAVECAR $d/CHGCAR
  ( cd $d && timeout 10800 mpirun --allow-run-as-root --bind-to none -np 24 $VASP > stdout.log 2>&1 < /dev/null ) &
done
echo "8 snapshots relaunched (ENCUT 600, 3x3, 24 ranks)"

# ---------- defect LAK static (sequential, 96 ranks) ----------
(
for name in vac_S sub_O; do
  d=defects/${name}_lak
  mkdir -p $d
  cp defects/${name}_gpu/CONTCAR $d/POSCAR
  cp defects/${name}_gpu/POTCAR $d/POTCAR
  printf "sc\n0\nGamma\n3 3 1\n0 0 0\n" > $d/KPOINTS
  cat > $d/INCAR << 'EOF'
SYSTEM = defect LAK static
PREC = Accurate
ENCUT = 700
EDIFF = 1E-6
NSW = 0
ISMEAR = 0
SIGMA = 0.03
METAGGA = LIBXC
LIBXC1 = MGGA_X_LAK
LIBXC2 = MGGA_C_LAK
LASPH = .TRUE.
LREAL = Auto
LWAVE = .FALSE.
LCHARG = .FALSE.
LMAXMIX = 6
NBANDS = 320
NCORE = 8
EOF
  ( cd $d && timeout 10800 mpirun --allow-run-as-root --bind-to none -np 96 $VASP > stdout.log 2>&1 < /dev/null )
  echo "DEFECT_LAK_STATIC_DONE_$name"
done
) &
echo "defect LAK static chain launched"
wait
echo PRECISION2_ALL_DONE
