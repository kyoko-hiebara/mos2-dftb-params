# Spin-orbit coupling constants (calibrated against LAK+SOC)

Add to the `Hamiltonian = DFTB { }` block (works with skf_v3 / skf_v3o):

```
SpinOrbit = {
  Dual = Yes
  Mo [eV] = {0.0 0.036 0.0931}   # shells in s, p, d order
  S  [eV] = {0.0 0.055 0.0}
  O  [eV] = {0.0 0.02 0.0}       # atomic estimate (negligible for O_S physics)
  Sb [eV] = {0.0 0.57 0.0}       # calibrated vs GPAW PBE+SOC bulk bands
}
```

**Sb**: calibrated against GPAW PBE+SOC bulk Sb bands (non-self-consistent SOC).
Since A7 Sb is inversion-symmetric (Kramers-degenerate bands), the calibration
target is the Gamma-point valence multiplet span: 0.952 eV (DFTB) vs 0.950 eV
(reference); the SOC-split level pattern (2+2+2) is reproduced. See
`docs/sb_soc_comparison.png` and `scripts/sb/calibrate_sb_soc.py`.

Validation at K (monolayer, a = 3.16 Å): VB splitting 150.2 meV (= LAK+SOC
reference, fitted), CB splitting 2.8 meV (reference 2.8 meV — independent check,
not fitted). xi_Mo(4d) = 0.093 eV is consistent with the atomic value.
Recalibrate with `soc_calibrate.py` if the electronic set changes.
