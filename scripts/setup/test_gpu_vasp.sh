#!/bin/bash
# GPU 版 VASP の検証: CPU 版 test_pbe と同一条件で TOTEN 一致 + 速度確認
set -u
NVVER=$(ls /opt/nvidia/hpc_sdk/Linux_x86_64/ | grep -E "^[0-9]+\.[0-9]+$" | sort -V | tail -1)
NVBASE=/opt/nvidia/hpc_sdk/Linux_x86_64/$NVVER
export PATH=$NVBASE/compilers/bin:$NVBASE/comm_libs/mpi/bin:$PATH
export LD_LIBRARY_PATH=$NVBASE/compilers/lib:$NVBASE/compilers/extras/qd/lib:$NVBASE/comm_libs/mpi/lib:$NVBASE/math_libs/lib64:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=1

cd /workspace/MoS2_DFTB/ref_calc
rm -rf test_gpu && mkdir test_gpu
cp test_pbe/POSCAR test_pbe/POTCAR test_pbe/KPOINTS test_pbe/INCAR test_gpu/
cd test_gpu
t0=$(date +%s)
timeout 900 mpirun --allow-run-as-root -np 1 /root/vasp.6.4.2-gpu/bin/vasp_std > stdout.log 2>&1 < /dev/null
t1=$(date +%s)
echo "GPU run wall time: $((t1-t0)) s"
grep "free  energy" OUTCAR | tail -1
echo "CPU reference:"
grep "free  energy" ../test_pbe/OUTCAR | tail -1
grep -i "GPU" stdout.log | head -3
echo GPU_TEST_DONE
