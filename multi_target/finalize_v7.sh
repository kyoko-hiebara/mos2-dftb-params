#!/bin/zsh
# v6: optm4c (+U, 単調性) → 暫定 skf_v6 → O 再フィット (opto6) → 最終 skf_v6 → Mo/S 反発 (rep+力) → Sb 反発 → SOC → 界面 → 緩和
cd /Users/crocus/uhuhu/MoS2_DFTB/local_opt
source ../venv_local/bin/activate
export DYLD_LIBRARY_PATH=/Users/crocus/uhuhu/MoS2_DFTB/sw_local/libxc-install/lib
MOS2=optm6,optm6c; SB=optsb4c; OUT=skf_v7
export DFTB_EXTRA_HSD="$(python u_block_from_study.py $MOS2)"
echo "[0] U block:"; echo "$DFTB_EXTRA_HSD"
echo "[1] provisional $OUT $(date)"
python make_v5.py --out $OUT --mos2-study $MOS2 --sb-study $SB --o-study opto4 > make_${OUT}_prov.log 2>&1 || echo "make prov FAILED"
grep -E "^\[MoS2\]|^\[V_S\]|best trial" make_${OUT}_prov.log | cut -c1-200
echo "[2] opto6 $(date)"; rm -rf opto6
O_BASE_STUDY=$(python - <<'PY'
import optuna, os; optuna.logging.set_verbosity(optuna.logging.WARNING)
best=None
for n in ["optm6","optm6c"]:
    p=f"/Users/crocus/uhuhu/MoS2_DFTB/local_opt/{n}/study.db"
    if not os.path.exists(p): continue
    st=optuna.load_study(study_name=n, storage=f"sqlite:///{p}")
    try: t=st.best_trial
    except ValueError: continue
    if best is None or t.value<best[1].value: best=(n,t)
print(best[0])
PY
)
echo "O base study: $O_BASE_STUDY"
O_BASE_STUDY=$O_BASE_STUDY O_BASE_SKF=$OUT nohup python optimize_o.py 80 --tag opto6 > opto6_w1.log 2>&1 &
sleep 20
O_BASE_STUDY=$O_BASE_STUDY O_BASE_SKF=$OUT nohup python optimize_o.py 80 --tag opto6 > opto6_w2.log 2>&1 &
sleep 90
while [ "$(pgrep -f 'tag opto6' | wc -l)" -gt 0 ]; do sleep 60; done
echo "opto6: $(grep -h '^BEST' opto6_w1.log opto6_w2.log | head -1 | cut -c1-120)"
echo "[3] final $OUT $(date)"
python make_v5.py --out $OUT --mos2-study $MOS2 --sb-study $SB --o-study opto6 > make_$OUT.log 2>&1 || echo "make_v5 FAILED"
grep -E "^\[|best trial|re-anchored|Assert" make_$OUT.log | cut -c1-220; sed -n '/m\*VB/p;/Q-K/p' make_$OUT.log | cut -c1-200
echo "[4] Mo/S repulsion $(date)"; mkdir -p repfit7
python dftb_energy_local.py repfit3/pairset3_dft.db repfit7/pairset6 $OUT --workdir /private/tmp/claude-501/-Users-crocus-uhuhu-MoS2-DFTB/6416faee-14e1-48c8-9838-35c06f3e7813/scratchpad/dftb_energy_v7 > repfit7/dftb_energy.log 2>&1 || echo "dftb_energy FAILED"
source ../venv_ccs_local/bin/activate
for mode in rep sw; do for fo in "" "--forces"; do d=repfit7/${mode}${fo/--/_}; mkdir -p $d; cp repfit7/pairset6_dft.db repfit7/pairset6_dftb.db $d/
  python repfit_run.py $d pairset6 --pairs "Mo-S:3.6,S-S:4.2,Mo-Mo:4.6" --swtype $mode $fo > $d/repfit.log 2>&1; echo "  $d: $(grep -E 'MSE' $d/CCS_error.out)"; done; done
source ../venv_local/bin/activate
for d in repfit7/rep repfit7/rep_forces repfit7/sw_forces; do python attach_validate.py $OUT $d $d/skf > $d/attach.log 2>&1; echo "  $d: $(grep '^E(a)' $d/attach.log)"; python relax_test.py $d/skf 2>&1 | tail -1; done
echo "MILESTONE V7_ELECTRONIC_AND_REP_DONE $(date)"
