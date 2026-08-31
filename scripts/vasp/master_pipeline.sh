#!/bin/bash
# CPU quota 26 コアに合わせた直列マスターパイプライン
# Stage1: 欠陥 LAK 静的 x2 -> Stage2: SOC -> Stage3: snapshots (時間の許す限り)
set -u
cd /workspace/MoS2_DFTB/ref_calc
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
VASP=/root/vasp.6.4.2/bin/vasp_std
LOG=/workspace/MoS2_DFTB/logs/master.log
echo "=== master pipeline start $(date) ===" >> $LOG

# ---------- Stage 1: defect LAK static (ENCUT 600, 24 ranks, sequential) ----------
for name in vac_S sub_O; do
  d=defects/${name}_lak
  mkdir -p $d
  cp defects/${name}_gpu/CONTCAR $d/POSCAR
  cp defects/${name}_gpu/POTCAR $d/POTCAR
  printf "sc\n0\nGamma\n2 2 1\n0 0 0\n" > $d/KPOINTS
  cat > $d/INCAR << 'EOF'
SYSTEM = defect LAK static
PREC = Accurate
ENCUT = 600
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
LCHGCAR = .FALSE.
LCHARG = .FALSE.
LMAXMIX = 6
NEDOS = 2000
NCORE = 4
EOF
  rm -f $d/stdout.log
  ( cd $d && timeout 7200 mpirun --allow-run-as-root -np 24 $VASP > stdout.log 2>&1 < /dev/null )
  echo "STAGE1_${name}_done $(grep -c 'F=' $d/stdout.log) $(date +%H:%M)" >> $LOG
done

# ---------- Stage 2: SOC (needs vasp_ncl from /root/vasp-ncl) ----------
if grep -q NCL2_OK /workspace/MoS2_DFTB/logs/vasp_ncl2_build.log 2>/dev/null; then
  NCL=/root/vasp-ncl/bin/vasp_ncl
  for tag in lak pbe; do
    d=soc_$tag
    mkdir -p $d
    cp bands_lak_a316/POSCAR $d/POSCAR
    cp POTCAR_Mo_S $d/POTCAR
    printf "m\n0\nGamma\n12 12 1\n0 0 0\n" > $d/KPOINTS
    { echo "PREC = Accurate"; echo "ENCUT = 900"; echo "EDIFF = 1E-6"; echo "NELM = 200"
      echo "ISMEAR = 0"; echo "SIGMA = 0.03"
      echo "LSORBIT = .TRUE."; echo "ISYM = -1"; echo "GGA_COMPAT = .FALSE."
      echo "MAGMOM = 9*0.0"; echo "NBANDS = 60"
      if [ $tag = lak ]; then
        echo "METAGGA = LIBXC"; echo "LIBXC1 = MGGA_X_LAK"; echo "LIBXC2 = MGGA_C_LAK"
      else
        echo "GGA = PE"
      fi
      echo "LASPH = .TRUE."; echo "LREAL = .FALSE."
      echo "LWAVE = .FALSE."; echo "LCHARG = .FALSE."; echo "LMAXMIX = 6"; echo "NCORE = 4"
    } > $d/INCAR
    rm -f $d/stdout.log
    ( cd $d && timeout 7200 mpirun --allow-run-as-root -np 24 $NCL > stdout.log 2>&1 < /dev/null )
    echo "STAGE2_soc_${tag}_done $(date +%H:%M)" >> $LOG
  done
else
  echo "STAGE2 skipped (ncl not ready)" >> $LOG
fi

# ---------- Stage 3: snapshots sequential ----------
for i in 00 01 02 03 04 05 06 07; do
  d=snapshots/snap_$i
  rm -f $d/stdout.log $d/OUTCAR $d/OSZICAR
  ( cd $d && timeout 3600 mpirun --allow-run-as-root -np 24 $VASP > stdout.log 2>&1 < /dev/null )
  echo "STAGE3_snap_${i}_done $(grep -c 'F=' $d/stdout.log) $(date +%H:%M)" >> $LOG
done
echo MASTER_PIPELINE_DONE >> $LOG
