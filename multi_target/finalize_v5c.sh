#!/bin/zsh
# 最終組み立て: make_v5 (O=opto3, Sb=best of optsb4/optsb4c) → Mo/S スプライン結合 → Sb 反発 (sw, energies) → SOC Γ スパン → 界面
cd /Users/crocus/uhuhu/MoS2_DFTB/local_opt
source ../venv_local/bin/activate
export DYLD_LIBRARY_PATH=/Users/crocus/uhuhu/MoS2_DFTB/sw_local/libxc-install/lib
MOS2=${1:-optm2}; SB=${2:-optsb4,optsb4c}; O=${3:-opto3}; OUT=${4:-skf_v5}; MOSREP=${5:-repfit5}
echo "[1] make_v5 $(date) ( $MOS2 / $SB / $O -> $OUT ; Mo/S splines from $MOSREP )"
python make_v5.py --out $OUT --mos2-study $MOS2 --sb-study $SB --o-study $O > make_$OUT.log 2>&1 || echo "make_v5 FAILED"
grep -E "^\[|re-anchored|best trial|Assert" make_$OUT.log
echo "[2] attach Mo/S splines $(date)"
python attach_validate.py $OUT $MOSREP ${OUT}rep --no-validate > ${OUT}_attach.log 2>&1 || echo "attach FAILED"
echo "[3] Sb repulsion (sw, energies) $(date)"; mkdir -p repfit_sb_final
python dftb_energy_sb.py repfit_sb/sb_ref_dft.db repfit_sb_final/sbset ${OUT}rep --keep-rep > repfit_sb_final/dftb_energy.log 2>&1 || echo "dftb_energy_sb FAILED"
source ../venv_ccs_local/bin/activate
python repfit_run.py repfit_sb_final sbset --pairs "S-Sb:3.6,Sb-Sb:3.9" --swtype sw > repfit_sb_final/repfit.log 2>&1 || echo "repfit_sb FAILED"
grep -E "MSE|Maxerror" repfit_sb_final/CCS_error.out
source ../venv_local/bin/activate
python attach_validate.py ${OUT}rep repfit_sb_final ${OUT}rep_full --pairs Sb-Sb,Sb-S,S-Sb --no-validate > repfit_sb_final/attach.log 2>&1 || echo "attach_sb FAILED"
python validate_sb_rep.py ${OUT}rep_full > repfit_sb_final/validate.log 2>&1 || echo "validate_sb_rep FAILED"
grep -E "scale_eq" repfit_sb_final/validate.log
echo "[4] SOC $(date)"
python soc_calibrate.py /Users/crocus/uhuhu/MoS2_DFTB/local_opt/$OUT > soc_${OUT}_mos.log 2>&1 || echo "soc_calibrate FAILED"
python sb_gamma_span.py $OUT > soc_${OUT}_sb_span.log 2>&1 || echo "sb_gamma_span FAILED"
grep -E "calibrated" soc_${OUT}_mos.log soc_${OUT}_sb_span.log
echo "[5] interface $(date)"; python dftb_interface_check.py $OUT > iface_$OUT.log 2>&1 || echo "iface FAILED"
grep -E "E_F =" iface_$OUT.log
echo "FINALIZE_DONE $(date)"
