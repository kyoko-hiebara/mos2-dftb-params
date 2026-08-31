#!/bin/bash
# 単層 MoS2 の E(a) スキャン (LAK)。各点で内部座標緩和 (ISIF=2)。
set -u
source /workspace/MoS2_DFTB/venv/bin/activate
cd /workspace/MoS2_DFTB/ref_calc
mkdir -p ascan
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

python3 << 'PYEOF'
from ase.build import mx2
from ase.io import write
import os
for a in [3.08, 3.10, 3.12, 3.14, 3.16, 3.18, 3.20, 3.22, 3.24, 3.26]:
    tag = f"{a:.2f}"
    d = f"/workspace/MoS2_DFTB/ref_calc/ascan/a_{tag}"
    os.makedirs(d, exist_ok=True)
    atoms = mx2(formula="MoS2", kind="2H", a=a, thickness=3.127, vacuum=10.0)
    write(f"{d}/POSCAR", atoms, format="vasp", direct=True, sort=True)
print("POSCARs written")
PYEOF

for a in 3.08 3.10 3.12 3.14 3.16 3.18 3.20 3.22 3.24 3.26; do
  d=ascan/a_$a
  cp POTCAR_Mo_S $d/POTCAR
  printf "mono\n0\nGamma\n12 12 1\n0 0 0\n" > $d/KPOINTS
  cat > $d/INCAR << EOF
SYSTEM = ascan a=$a LAK
PREC = Accurate
ENCUT = 900
EDIFF = 1E-7
EDIFFG = -5E-3
IBRION = 2
ISIF = 2
NSW = 30
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
EOF
  ( cd $d && timeout 5400 mpirun --allow-run-as-root -np 12 /root/vasp.6.4.2/bin/vasp_std > stdout.log 2>&1 < /dev/null ) &
done
wait

echo "=== E(a) results (LAK) ==="
for a in 3.08 3.10 3.12 3.14 3.16 3.18 3.20 3.22 3.24 3.26; do
  printf "a=%s  " $a
  grep "free  energy" ascan/a_$a/OUTCAR 2>/dev/null | tail -1 || echo MISSING
done
echo ASCAN_DONE
