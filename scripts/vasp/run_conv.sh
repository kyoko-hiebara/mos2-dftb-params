#!/bin/bash
# LAK での ENCUT / k-mesh 収束テスト (単層 MoS2、ジョブ並列)
set -u
cd /workspace/MoS2_DFTB/ref_calc
mkdir -p conv
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

write_incar () {  # $1 = ENCUT
  cat << EOF
SYSTEM = conv test LAK
PREC = Accurate
ENCUT = $1
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
NCORE = 4
EOF
}

launch () {  # $1 = dir, $2 = encut, $3 = kmesh
  d=conv/$1
  mkdir -p $d
  cp mono/POSCAR $d/POSCAR
  cp POTCAR_Mo_S $d/POTCAR
  write_incar $2 > $d/INCAR
  printf "mono\n0\nGamma\n%s %s 1\n0 0 0\n" $3 $3 > $d/KPOINTS
  ( cd $d && timeout 1800 mpirun --allow-run-as-root -np 8 /root/vasp.6.4.2/bin/vasp_std > stdout.log 2>&1 < /dev/null ) &
}

# ENCUT scan @ 12x12x1
for E in 600 700 800 900 1000; do launch encut_$E $E 12; done
# k scan @ ENCUT=900
for K in 6 9 15 18; do launch kmesh_$K 900 $K; done
wait

echo "=== ENCUT convergence (12x12x1) ==="
for E in 600 700 800 900 1000; do
  printf "ENCUT=%s  " $E
  grep "free  energy" conv/encut_$E/OUTCAR 2>/dev/null | tail -1 || echo MISSING
done
echo "=== k convergence (ENCUT=900) ==="
for K in 6 9 12 15 18; do
  d=conv/kmesh_$K; [ $K = 12 ] && d=conv/encut_900
  printf "k=%sx%s  " $K $K
  grep "free  energy" $d/OUTCAR 2>/dev/null | tail -1 || echo MISSING
done
echo CONV_DONE
