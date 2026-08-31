#!/bin/bash
# 反発フィット用 LAK 参照データ第1弾:
#  (1) 圧縮側+遠方の E(a) (内部緩和付き)
#  (2) 層厚スキャン (a=3.16 固定、静的、力出力)
#  (3) 分子解離カーブ S2 / SO / O2 (スピン偏極)
set -u
source /workspace/MoS2_DFTB/venv/bin/activate
cd /workspace/MoS2_DFTB/ref_calc
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
VASP=/root/vasp.6.4.2/bin/vasp_std

lak_incar () {  # $1=NSW $2=extra
  cat << EOF
PREC = Accurate
ENCUT = 900
EDIFF = 1E-7
ISMEAR = 0
SIGMA = 0.03
METAGGA = LIBXC
LIBXC1 = MGGA_X_LAK
LIBXC2 = MGGA_C_LAK
LASPH = .TRUE.
LREAL = .FALSE.
LWAVE = .FALSE.
LCHARG = .FALSE.
LMAXMIX = 6
NCORE = 4
NSW = $1
$2
EOF
}

# ---------- (1) 圧縮側 + 遠方 E(a) ----------
python3 << 'PYEOF'
from ase.build import mx2
from ase.io import write
import os
for a in [2.90, 2.95, 3.00, 3.04, 3.32, 3.40]:
    d = f"/workspace/MoS2_DFTB/ref_calc/ascan/a_{a:.2f}"
    os.makedirs(d, exist_ok=True)
    atoms = mx2(formula="MoS2", kind="2H", a=a, thickness=3.127, vacuum=10.0)
    write(f"{d}/POSCAR", atoms, format="vasp", direct=True, sort=True)
print("iso-scan POSCARs written")
PYEOF
for a in 2.90 2.95 3.00 3.04 3.32 3.40; do
  d=ascan/a_$a
  cp POTCAR_Mo_S $d/POTCAR
  printf "m\n0\nGamma\n12 12 1\n0 0 0\n" > $d/KPOINTS
  { echo "SYSTEM = iso a=$a"; lak_incar 30 "IBRION = 2
ISIF = 2
EDIFFG = -5E-3"; } > $d/INCAR
  ( cd $d && timeout 14400 mpirun --allow-run-as-root -np 12 $VASP > stdout.log 2>&1 < /dev/null ) &
done

# ---------- (2) 層厚スキャン (静的、力を記録) ----------
python3 << 'PYEOF'
from ase.build import mx2
from ase.io import write
import os
for t in [2.80, 2.95, 3.30, 3.45]:
    d = f"/workspace/MoS2_DFTB/ref_calc/tscan/t_{t:.2f}"
    os.makedirs(d, exist_ok=True)
    atoms = mx2(formula="MoS2", kind="2H", a=3.16, thickness=t, vacuum=10.0)
    write(f"{d}/POSCAR", atoms, format="vasp", direct=True, sort=True)
print("t-scan POSCARs written")
PYEOF
for t in 2.80 2.95 3.30 3.45; do
  d=tscan/t_$t
  cp POTCAR_Mo_S $d/POTCAR
  printf "m\n0\nGamma\n12 12 1\n0 0 0\n" > $d/KPOINTS
  { echo "SYSTEM = tscan t=$t"; lak_incar 0 ""; } > $d/INCAR
  ( cd $d && timeout 7200 mpirun --allow-run-as-root -np 8 $VASP > stdout.log 2>&1 < /dev/null ) &
done

# ---------- (3) 分子解離カーブ ----------
python3 << 'PYEOF'
from ase import Atoms
from ase.io import write
import os
curves = {
    "S2": ("S2", [1.70, 1.80, 1.90, 2.00, 2.20, 2.50], 2),
    "SO": ("SO", [1.35, 1.48, 1.60, 1.80, 2.10], 2),
    "O2": ("O2", [1.10, 1.21, 1.35, 1.60], 2),
}
for name, (formula, dists, mag) in curves.items():
    for r in dists:
        d = f"/workspace/MoS2_DFTB/ref_calc/molecules_lak/{name}_r{r:.2f}"
        os.makedirs(d, exist_ok=True)
        atoms = Atoms(formula, positions=[(0, 0, 0), (0, 0, r)])
        atoms.set_cell([14.0, 14.5, 15.0])
        atoms.center()
        write(f"{d}/POSCAR", atoms, format="vasp", direct=True, sort=True)
print("molecule POSCARs written")
PYEOF
mol_potcar () {
  case $1 in
    S2) cat /workspace/MoS2_DFTB/sw/S_h/POTCAR ;;
    SO) cat /workspace/MoS2_DFTB/sw/O_h/POTCAR /workspace/MoS2_DFTB/sw/S_h/POTCAR ;;
    O2) cat /workspace/MoS2_DFTB/sw/O_h/POTCAR ;;
  esac
}
njobs=0
for m in S2 SO O2; do
  for d in molecules_lak/${m}_r*; do
    mol_potcar $m > $d/POTCAR
    printf "m\n0\nGamma\n1 1 1\n0 0 0\n" > $d/KPOINTS
    { echo "SYSTEM = $d"; lak_incar 0 "ISPIN = 2
NUPDOWN = 2
ISYM = 0"; } > $d/INCAR
    ( cd $d && timeout 7200 mpirun --allow-run-as-root -np 8 $VASP > stdout.log 2>&1 < /dev/null ) &
    njobs=$((njobs+1))
    # 分子は 6 個ずつ波状投入
    if [ $((njobs % 6)) -eq 0 ]; then wait; fi
  done
done
wait
echo REP_DATA1_DONE
