#!/bin/bash
# libxc 7.1.2 shared build + pylibxc + hotcent を venv に導入
set -e
cd /workspace/MoS2_DFTB/sw/libxc-7.1.2
ls setup.py pylibxc 2>/dev/null || echo "no setup.py/pylibxc at root"

cmake -B build-shared -DCMAKE_INSTALL_PREFIX=/opt/libxc-7.1.2-shared \
  -DENABLE_FORTRAN=OFF -DBUILD_SHARED_LIBS=ON -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_TESTING=OFF -DDISABLE_FHC=ON > /workspace/MoS2_DFTB/logs/libxc_shared_cmake.log 2>&1
cmake --build build-shared -j 48 > /workspace/MoS2_DFTB/logs/libxc_shared_build.log 2>&1
cmake --install build-shared > /dev/null
echo "shared libxc installed:"
ls /opt/libxc-7.1.2-shared/lib/ | head -5

source /workspace/MoS2_DFTB/venv/bin/activate
# pylibxc: libxc ソースルートの setup.py を利用
pip install . 2>&1 | tail -2 || echo "pip install of pylibxc from source root failed"
export LD_LIBRARY_PATH=/opt/libxc-7.1.2-shared/lib:${LD_LIBRARY_PATH:-}
python3 -c "import pylibxc; print('pylibxc OK, libxc version:', pylibxc.util.xc_version_string()); f=pylibxc.LibXCFunctional('MGGA_X_LAK','unpolarized'); print('LAK via pylibxc OK')" || echo PYLIBXC_CHECK_FAILED

# hotcent (C 拡張の再生成には cython)
pip install -q cython
cd /workspace/MoS2_DFTB/sw/hotcent
pip install . 2>&1 | tail -2
python3 -c "import hotcent; from hotcent.atomic_dft import AtomicDFT; print('hotcent OK')" || echo HOTCENT_CHECK_FAILED
echo SETUP_PYLIBXC_DONE
