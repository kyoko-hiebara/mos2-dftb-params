#!/bin/bash
# vasp_ncl ビルド完了を待って LAK+SOC / PBE+SOC の単層計算を実行し K 点分裂を抽出
set -u
until grep -qE "NCL_OK|NCL_FAILED" /workspace/MoS2_DFTB/logs/vasp_ncl_build.log 2>/dev/null; do sleep 60; done
if grep -q NCL_FAILED /workspace/MoS2_DFTB/logs/vasp_ncl_build.log; then
  echo "ncl build failed"; tail -20 /workspace/MoS2_DFTB/logs/vasp_ncl_build.log; exit 1
fi
echo "ncl ready"
cd /workspace/MoS2_DFTB/ref_calc
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

soc_incar () {  # $1 = xc block
  cat << EOF
PREC = Accurate
ENCUT = 900
EDIFF = 1E-6
NELM = 200
ISMEAR = 0
SIGMA = 0.03
LSORBIT = .TRUE.
ISYM = -1
GGA_COMPAT = .FALSE.
MAGMOM = 9*0.0
NBANDS = 60
$1
LASPH = .TRUE.
LREAL = .FALSE.
LWAVE = .FALSE.
LCHARG = .FALSE.
LMAXMIX = 6
NCORE = 4
EOF
}

for tag in lak pbe; do
  d=soc_$tag
  mkdir -p $d
  cp bands_lak_a316/POSCAR $d/POSCAR
  cp POTCAR_Mo_S $d/POTCAR
  printf "m\n0\nGamma\n12 12 1\n0 0 0\n" > $d/KPOINTS
  if [ $tag = lak ]; then
    soc_incar "METAGGA = LIBXC
LIBXC1 = MGGA_X_LAK
LIBXC2 = MGGA_C_LAK" > $d/INCAR
  else
    soc_incar "GGA = PE" > $d/INCAR
  fi
  ( cd $d && timeout 10000 mpirun --allow-run-as-root -np 24 /root/vasp.6.4.2/bin/vasp_ncl > stdout.log 2>&1 < /dev/null ) &
done
wait

source /workspace/MoS2_DFTB/venv/bin/activate
python3 << 'PYEOF'
import numpy as np
for tag in ["lak", "pbe"]:
    try:
        lines = open(f"/workspace/MoS2_DFTB/ref_calc/soc_{tag}/EIGENVAL").readlines()
        nelec, nk, nb = (int(x) for x in lines[5].split())
        idx = 7
        found = False
        for ik in range(nk):
            k = [float(x) for x in lines[idx].split()[:3]]
            eigs = [float(lines[idx + 1 + ib].split()[1]) for ib in range(nb)]
            if abs(k[0] - 1/3) < 1e-3 and abs(k[1] - 1/3) < 1e-3:
                # noncollinear: nelec bands occupied (1 e per band)
                vb = sorted(eigs)[:nelec]
                cb = sorted(eigs)[nelec:]
                print(f"{tag.upper()}+SOC @K: VB split = {(vb[-1]-vb[-3])*1000:.1f} meV "
                      f"(top pair vs next), CB split = {(cb[1]-cb[0])*1000:.1f} meV, "
                      f"gap = {cb[0]-vb[-1]:.4f} eV")
                found = True
                break
            idx += nb + 2
        if not found:
            print(f"{tag}: K point not found in EIGENVAL")
    except Exception as e:
        print(f"{tag}: FAILED ({e})")
PYEOF
echo SOC_CHAIN_DONE
