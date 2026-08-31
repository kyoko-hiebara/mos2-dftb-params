#!/bin/bash
# 欠陥スーパーセル PBE 緩和の GPU 実行 (usage: run_defect_gpu.sh <name>)
set -u
name=$1
NVVER=$(ls /opt/nvidia/hpc_sdk/Linux_x86_64/ | grep -E "^[0-9]+\.[0-9]+$" | sort -V | tail -1)
NVBASE=/opt/nvidia/hpc_sdk/Linux_x86_64/$NVVER
export PATH=$NVBASE/compilers/bin:$NVBASE/comm_libs/mpi/bin:$PATH
export LD_LIBRARY_PATH=$NVBASE/compilers/lib:$NVBASE/compilers/extras/qd/lib:$NVBASE/comm_libs/mpi/lib:$NVBASE/math_libs/lib64:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=1
export OMPI_MCA_pml=ob1
export OMPI_MCA_btl=self,vader
export UCX_TLS=self,sm,tcp

src=/workspace/MoS2_DFTB/ref_calc/defects/${name}_pbe
dst=/workspace/MoS2_DFTB/ref_calc/defects/${name}_gpu
rm -rf $dst && mkdir -p $dst
cp $src/POSCAR $src/POTCAR $src/KPOINTS $dst/
# NCORE は GPU 版では 1 に (GPU ポートの制約)
sed 's/^NCORE.*/NCORE = 1/; s/^LREAL.*/LREAL = Auto/' $src/INCAR > $dst/INCAR
cd $dst
t0=$(date +%s)
timeout 43200 mpirun --allow-run-as-root --bind-to none -np 1 /root/vasp.6.4.2-gpu/bin/vasp_std > stdout.log 2>&1 < /dev/null
rc=$?
t1=$(date +%s)
echo "GPU defect $name finished rc=$rc wall=$((t1-t0))s"
grep "F=" stdout.log | tail -1
echo DEFECT_GPU_DONE_$name
