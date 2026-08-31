#!/bin/bash
# 5x5 単層スーパーセルの Vs (S空孔) / O_S (O置換) を PBE で予備緩和
set -u
source /workspace/MoS2_DFTB/venv/bin/activate
cd /workspace/MoS2_DFTB/ref_calc
mkdir -p defects
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

python3 << 'PYEOF'
import os
from ase.build import mx2
from ase.io import write

A_PBE = 3.183  # PBE 文献値 (予備緩和用; LAK 格子確定後に再構築する)
base = mx2(formula="MoS2", kind="2H", a=A_PBE, thickness=3.13, vacuum=10.0)
sc = base.repeat((5, 5, 1))

# 上面 S のインデックスを1つ選ぶ (最初の上面 S)
top_s = [i for i, at in enumerate(sc) if at.symbol == "S"
         and at.position[2] > sc.positions[:, 2].mean()]
idx = top_s[0]

vac = sc.copy()
del vac[idx]

osub = sc.copy()
osub[idx].symbol = "O"

for name, atoms in [("vac_S", vac), ("sub_O", osub)]:
    d = f"/workspace/MoS2_DFTB/ref_calc/defects/{name}_pbe"
    os.makedirs(d, exist_ok=True)
    write(f"{d}/POSCAR", atoms, format="vasp", direct=True, sort=True)
    print(name, atoms.get_chemical_formula(), "->", d)
PYEOF

# POTCAR: 元素順は POSCAR ヘッダに合わせる
head -1 defects/vac_S_pbe/POSCAR
head -1 defects/sub_O_pbe/POSCAR
cat /workspace/MoS2_DFTB/sw/Mo_sv/POTCAR /workspace/MoS2_DFTB/sw/S_h/POTCAR > defects/vac_S_pbe/POTCAR
cat /workspace/MoS2_DFTB/sw/Mo_sv/POTCAR /workspace/MoS2_DFTB/sw/O_h/POTCAR /workspace/MoS2_DFTB/sw/S_h/POTCAR > defects/sub_O_pbe/POTCAR

for d in defects/vac_S_pbe defects/sub_O_pbe; do
  printf "sc\n0\nGamma\n3 3 1\n0 0 0\n" > $d/KPOINTS
  cat > $d/INCAR << 'EOF'
SYSTEM = defect PBE pre-relax
PREC = Accurate
ENCUT = 900
EDIFF = 1E-6
EDIFFG = -0.02
IBRION = 2
ISIF = 2
NSW = 60
ISMEAR = 0
SIGMA = 0.03
GGA = PE
LASPH = .TRUE.
LREAL = Auto
LWAVE = .FALSE.
LCHARG = .FALSE.
NCORE = 6
EOF
  ( cd $d && timeout 21600 mpirun --allow-run-as-root -np 24 /root/vasp.6.4.2/bin/vasp_std > stdout.log 2>&1 < /dev/null ) &
done
wait
echo DEFECT_PBE_DONE
for d in defects/vac_S_pbe defects/sub_O_pbe; do
  echo "--- $d ---"; grep "F=" $d/stdout.log | tail -1
done
