#!/bin/bash
# 反発スプラインを SKF に結合し、E(a) カーブで検証
set -e
cd /workspace/MoS2_DFTB/repfit
source /workspace/MoS2_DFTB/venv_ccs/bin/activate
# ccs_export_sktable は端末サイズを要求するので script で疑似端末を与える
script -qec "ccs_export_sktable CCS_params.json" /dev/null > /dev/null 2>&1 || \
  COLUMNS=80 LINES=24 python3 -c "
import sys, os
sys.argv = ['x', 'CCS_params.json']
os.get_terminal_size = lambda *a: os.terminal_size((80, 24))
from ccs_fit.scripts.ccs_export_sktable import main
main()"
ls *.spl
deactivate

mkdir -p /workspace/MoS2_DFTB/dftb/skf_v2rep
for pair in Mo-Mo Mo-S S-Mo S-S; do
  spl=$pair.spl
  # S-Mo.spl が無ければ Mo-S.spl を使う (対称)
  if [ ! -f "$spl" ]; then
    alt=$(echo $pair | awk -F- '{print $2"-"$1}').spl
    spl=$alt
  fi
  cat /workspace/MoS2_DFTB/dftb/skf_v2/$pair.skf "$spl" > /workspace/MoS2_DFTB/dftb/skf_v2rep/$pair.skf
done
echo "combined SKF written:"
grep -c Spline /workspace/MoS2_DFTB/dftb/skf_v2rep/*.skf

# --- E(a) validation ---
source /workspace/MoS2_DFTB/venv/bin/activate
python3 /workspace/MoS2_DFTB/validate_rep.py
echo ATTACH_VALIDATE_DONE
