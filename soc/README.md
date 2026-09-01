# Spin-orbit coupling constants (calibrated against LAK+SOC)

Add to the `Hamiltonian = DFTB { }` block (works with skf_v3 / skf_v3o):

```
SpinOrbit = {
  Dual = Yes
  Mo [eV] = {0.0 0.036 0.0931}   # shells in s, p, d order
  S  [eV] = {0.0 0.055 0.0}
}
```

Validation at K (monolayer, a = 3.16 Å): VB splitting 150.2 meV (= LAK+SOC
reference, fitted), CB splitting 2.8 meV (reference 2.8 meV — independent check,
not fitted). xi_Mo(4d) = 0.093 eV is consistent with the atomic value.
Recalibrate with `soc_calibrate.py` if the electronic set changes.
