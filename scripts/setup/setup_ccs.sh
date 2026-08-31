#!/bin/bash
# ccs_fit 専用 venv (numpy<1.23) + PTBP ベースラインセット取得
set -e
if [ ! -f /workspace/MoS2_DFTB/venv_ccs/bin/pip ]; then
  python3 -m venv /workspace/MoS2_DFTB/venv_ccs
fi
source /workspace/MoS2_DFTB/venv_ccs/bin/activate
pip install -q --upgrade pip 2>&1 | tail -1
pip install -q "numpy<1.23" ccs_fit "ase==3.22.1" 2>&1 | tail -1
python3 -c "import ccs_fit, numpy; print('ccs venv OK, numpy', numpy.__version__)"
deactivate

cd /workspace/MoS2_DFTB/sw
if [ ! -d ptbp ]; then
  curl -sL -o ParameterSets.zip "https://zenodo.org/api/records/14289468/files/ParameterSets.zip/content"
  unzip -q -o ParameterSets.zip -d ptbp
fi
echo "--- Mo-related PTBP files ---"
find ptbp -iname "*Mo*" | head -8
echo CCS_PTBP_SETUP_DONE
