#!/bin/bash
# nvhpc インストール完了を待って GPU 版 VASP (OpenACC, cc80) をビルド
# 注: GPU ポートは libxc 非対応のため、このビルドは PBE/r2SCAN 等ネイティブ汎函数専用
set -u
until grep -q "NVHPC_INSTALL_OK" /workspace/MoS2_DFTB/logs/nvhpc_install.log 2>/dev/null; do
  if grep -q "NVHPC_INSTALL_FAILED" /workspace/MoS2_DFTB/logs/nvhpc_install.log 2>/dev/null; then
    echo "nvhpc install failed, aborting"; exit 1
  fi
  sleep 30
done
echo "nvhpc ready"

NVVER=$(ls /opt/nvidia/hpc_sdk/Linux_x86_64/ | grep -E "^[0-9]+\.[0-9]+$" | sort -V | tail -1)
NVBASE=/opt/nvidia/hpc_sdk/Linux_x86_64/$NVVER
export PATH=$NVBASE/compilers/bin:$NVBASE/comm_libs/mpi/bin:$PATH
echo "using nvhpc $NVVER"; which nvfortran mpif90

rsync -a --exclude 'build/*' /root/vasp.6.4.2/ /root/vasp.6.4.2-gpu/
cat > /root/vasp.6.4.2-gpu/makefile.include << 'MEOF'
# GPU (OpenACC) build for A100 — native functionals only (no libxc!)
CPP_OPTIONS = -DHOST=\"LinuxNV\" \
              -DMPI -DMPI_INPLACE -DMPI_BLOCK=8000 -Duse_collective \
              -DscaLAPACK \
              -DCACHE_SIZE=4000 \
              -Davoidalloc \
              -Dvasp6 \
              -Duse_bse_te \
              -Dtbdyn \
              -Dqd_emulate \
              -Dfock_dblbuf \
              -D_OPENMP \
              -D_OPENACC \
              -DUSENCCL -DUSENCCLP2P

CPP         = nvfortran -Mpreprocess -Mfree -Mextend -E $(CPP_OPTIONS) $*$(FUFFIX)  > $*$(SUFFIX)

FC          = mpif90 -acc -gpu=cc80 -mp
FCL         = mpif90 -acc -gpu=cc80 -mp -c++libs

FREE        = -Mfree

FFLAGS      = -Mbackslash -Mlarge_arrays

OFLAG       = -fast

DEBUG       = -Mfree -O0 -traceback

OBJECTS     = fftmpiw.o fftmpi_map.o fftw3d.o fft3dlib.o

LLIBS       = -cudalib=cublas,cusolver,cufft,nccl -cuda

SOURCE_O1  := pade_fit.o minimax_dependence.o
SOURCE_O2  := pead.o

CPP_LIB     = $(CPP)
FC_LIB      = nvfortran
CC_LIB      = nvc -w
CFLAGS_LIB  = -O
FFLAGS_LIB  = -O1 -Mfixed
FREE_LIB    = $(FREE)

OBJECTS_LIB = linpack_double.o

CXX_PARS    = nvc++ --no_warnings

VASP_TARGET_CPU ?= -tp host
FFLAGS     += $(VASP_TARGET_CPU)

NVROOT      =$(shell which nvfortran | awk -F /compilers/bin/nvfortran '{ print $$1 }')

QD         ?= $(NVROOT)/compilers/extras/qd
LLIBS      += -L$(QD)/lib -lqdmod -lqd
INCS       += -I$(QD)/include/qd

BLAS        = -lblas
LAPACK      = -llapack
SCALAPACK   = -Mscalapack
LLIBS      += $(SCALAPACK) $(LAPACK) $(BLAS)

FFTW_ROOT  ?= /usr
LLIBS      += -L$(FFTW_ROOT)/lib/x86_64-linux-gnu -lfftw3 -lfftw3_omp
INCS       += -I$(FFTW_ROOT)/include
MEOF

cd /root/vasp.6.4.2-gpu
make DEPS=1 -j32 std > /workspace/MoS2_DFTB/logs/vasp_gpu_build.log 2>&1 \
  && echo VASP_GPU_OK >> /workspace/MoS2_DFTB/logs/vasp_gpu_build.log \
  || echo VASP_GPU_FAILED >> /workspace/MoS2_DFTB/logs/vasp_gpu_build.log
tail -2 /workspace/MoS2_DFTB/logs/vasp_gpu_build.log
mkdir -p /workspace/MoS2_DFTB/sw/bin
cp bin/vasp_std /workspace/MoS2_DFTB/sw/bin/vasp_std_gpu 2>/dev/null || true
cp makefile.include /workspace/MoS2_DFTB/sw/makefile.include.gpu 2>/dev/null || true
echo GPU_CHAIN_DONE
