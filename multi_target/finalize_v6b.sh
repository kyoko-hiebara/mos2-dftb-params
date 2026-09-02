#!/bin/zsh
# v6 第 2 段: Mo/S rep+forces 結合 → Sb 反発 → SOC → 界面 → 緩和テスト
cd /Users/crocus/uhuhu/MoS2_DFTB/local_opt
source ../venv_local/bin/activate
export DYLD_LIBRARY_PATH=/Users/crocus/uhuhu/MoS2_DFTB/sw_local/libxc-install/lib
OUT=${1:-skf_v6}; MOSREP=${2:-repfit6/rep_forces}; MOS2=${3:-optm4,optm4c}
export DFTB_EXTRA_HSD="$(python u_block_from_study.py $MOS2)"
echo "[1] attach $(date)"; rm -rf ${OUT}_rep ${OUT}_rep_full
python attach_validate.py $OUT $MOSREP ${OUT}_rep --no-validate > /dev/null 2>&1 || echo "attach FAILED"
echo "[2] Sb repulsion $(date)"; mkdir -p repfit_sb_$OUT
python dftb_energy_sb.py repfit_sb/sb_ref_dft.db repfit_sb_$OUT/sbset ${OUT}_rep --keep-rep > repfit_sb_$OUT/dftb_energy.log 2>&1 || echo "dftb_energy_sb FAILED"
source ../venv_ccs_local/bin/activate
python repfit_run.py repfit_sb_$OUT sbset --pairs "S-Sb:3.6,Sb-Sb:3.9" --swtype sw > repfit_sb_$OUT/repfit.log 2>&1 || echo "repfit_sb FAILED"
grep -E "MSE" repfit_sb_$OUT/CCS_error.out
source ../venv_local/bin/activate
python attach_validate.py ${OUT}_rep repfit_sb_$OUT ${OUT}_rep_full --pairs Sb-Sb,Sb-S,S-Sb --no-validate > /dev/null 2>&1 || echo "attach_sb FAILED"
python validate_sb_rep.py ${OUT}_rep_full 2>&1 | grep scale_eq
python relax_sb_test.py ${OUT}_rep_full 2>&1 | tail -2
python relax_test.py ${OUT}_rep_full 2>&1 | tail -1
echo "[3] SOC $(date)"
python soc_calibrate.py /Users/crocus/uhuhu/MoS2_DFTB/local_opt/$OUT > soc_${OUT}_mos.log 2>&1 || echo "soc FAILED"
python sb_gamma_span.py $OUT > soc_${OUT}_sb_span.log 2>&1 || echo "span FAILED"
grep calibrated soc_${OUT}_mos.log soc_${OUT}_sb_span.log
echo "[4] interface $(date)"; python dftb_interface_check.py $OUT > iface_$OUT.log 2>&1 || echo "iface FAILED"
cat iface_$OUT.log
python plot_bands_v5.py ${OUT}_eval/dftb.json ../notes/bands_${OUT}_vs_lak.png | tail -1
echo "FINALIZE_${OUT}_DONE $(date)"
