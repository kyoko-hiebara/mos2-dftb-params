#!/bin/zsh
# skf_v5m: Mo/S 反発再フィット → Sb 反発再フィット → SOC → 界面 → 緩和テスト
cd /Users/crocus/uhuhu/MoS2_DFTB/local_opt
source ../venv_local/bin/activate
export DYLD_LIBRARY_PATH=/Users/crocus/uhuhu/MoS2_DFTB/sw_local/libxc-install/lib
OUT=skf_v5m
echo "[1] Mo/S repulsion $(date)"; mkdir -p repfit5m
python dftb_energy_local.py repfit3/pairset3_dft.db repfit5m/pairset5m $OUT --workdir /private/tmp/claude-501/-Users-crocus-uhuhu-MoS2-DFTB/6416faee-14e1-48c8-9838-35c06f3e7813/scratchpad/dftb_energy_v5m > repfit5m/dftb_energy.log 2>&1 || echo "dftb_energy FAILED"
source ../venv_ccs_local/bin/activate
python repfit_run.py repfit5m pairset5m --pairs "Mo-S:3.6,S-S:4.2,Mo-Mo:4.6" --forces > repfit5m/repfit.log 2>&1 || echo "repfit FAILED"
grep -E "MSE|Maxerror" repfit5m/CCS_error.out
source ../venv_local/bin/activate
python attach_validate.py $OUT repfit5m ${OUT}_rep > repfit5m/attach.log 2>&1 || echo "attach FAILED"
grep "^E(a)" repfit5m/attach.log
python relax_test.py ${OUT}_rep 2>&1 | tail -1
echo "[2] Sb repulsion $(date)"; mkdir -p repfit_sb_v5m
python dftb_energy_sb.py repfit_sb/sb_ref_dft.db repfit_sb_v5m/sbset ${OUT}_rep --keep-rep > repfit_sb_v5m/dftb_energy.log 2>&1 || echo "dftb_energy_sb FAILED"
source ../venv_ccs_local/bin/activate
python repfit_run.py repfit_sb_v5m sbset --pairs "S-Sb:3.6,Sb-Sb:3.9" --swtype sw > repfit_sb_v5m/repfit.log 2>&1 || echo "repfit_sb FAILED"
grep -E "MSE|Maxerror" repfit_sb_v5m/CCS_error.out
source ../venv_local/bin/activate
python attach_validate.py ${OUT}_rep repfit_sb_v5m ${OUT}_rep_full --pairs Sb-Sb,Sb-S,S-Sb --no-validate > repfit_sb_v5m/attach.log 2>&1 || echo "attach_sb FAILED"
python validate_sb_rep.py ${OUT}_rep_full > repfit_sb_v5m/validate.log 2>&1 || echo "validate FAILED"
grep scale_eq repfit_sb_v5m/validate.log
python relax_sb_test.py ${OUT}_rep_full 2>&1 | tail -2
echo "[3] SOC $(date)"
python soc_calibrate.py /Users/crocus/uhuhu/MoS2_DFTB/local_opt/$OUT > soc_${OUT}_mos.log 2>&1 || echo "soc FAILED"
python sb_gamma_span.py $OUT > soc_${OUT}_sb_span.log 2>&1 || echo "span FAILED"
grep calibrated soc_${OUT}_mos.log soc_${OUT}_sb_span.log
echo "[4] interface $(date)"; python dftb_interface_check.py $OUT > iface_$OUT.log 2>&1 || echo "iface FAILED"
cat iface_$OUT.log
echo "FINALIZE_V5M_DONE $(date)"
