#!/usr/bin/env python3
"""Mo-S-O-H の Slater-Koster ファイル生成 (hotcent v2, PBE)。

使い方:
  python3 gen_skf.py params.json outdir

params.json 形式 (Bohr 単位):
  {"Mo": {"r_dens": 8.7, "s_dens": 2.0,
          "r_wf": {"4d": 5.8, "5s": 5.8, "5p": 5.8}, "s_wf": 2.0},
   "S": {...}, "O": {...}, "H": {...}}

- 自由原子 (弱 confinement) から onsite 固有値・Hubbard U・スピン分極誤差を取得
- confined 原子で 2 中心 SK 積分 (density superposition)
- 全ペア (homo + hetero) の .skf を outdir に出力
"""
import json
import os
import sys
import time
from ase.units import Ha
from hotcent.atomic_dft import AtomicDFT
from hotcent.confinement import PowerConfinement
from hotcent.offsite_twocenter import Offsite2cTable

XC = "GGA_X_PBE+GGA_C_PBE"

SPECS = {
    "Mo": dict(configuration="[Kr] 4d5 5s1 5p0", valence=["4d", "5s", "5p"],
               scalarrel=True, occupations={"4d": 5, "5s": 1, "5p": 0},
               u_shell="4d"),
    "S": dict(configuration="[Ne] 3s2 3p4", valence=["3s", "3p"],
              scalarrel=False, occupations={"3s": 2, "3p": 4},
              u_shell="3p"),
    "O": dict(configuration="[He] 2s2 2p4", valence=["2s", "2p"],
              scalarrel=False, occupations={"2s": 2, "2p": 4},
              u_shell="2p"),
    "H": dict(configuration="1s1", valence=["1s"],
              scalarrel=False, occupations={"1s": 1},
              u_shell="1s"),
}

# SK テーブル範囲 (Bohr)
RMIN, DR, NPTS = 0.4, 0.02, 980


def free_atom_properties(el):
    """弱 confinement の自由原子: onsite 固有値と Hubbard U"""
    spec = SPECS[el]
    atom = AtomicDFT(el, xc=XC,
                     configuration=spec["configuration"],
                     valence=spec["valence"],
                     scalarrel=spec["scalarrel"],
                     confinement=PowerConfinement(r0=40.0, s=4),
                     perturbative_confinement=False,
                     txt=None)
    atom.run()
    eigs = {nl: atom.get_eigenvalue(nl) for nl in spec["valence"]}
    U = atom.get_hubbard_value(spec["u_shell"], scheme="central", maxstep=1)
    hubbards = {nl: U for nl in spec["valence"]}
    return eigs, hubbards


def confined_atom(el, p):
    spec = SPECS[el]
    wf_conf = {nl: PowerConfinement(r0=p["r_wf"][nl], s=p.get("s_wf", 2.0))
               for nl in spec["valence"]}
    atom = AtomicDFT(el, xc=XC,
                     configuration=spec["configuration"],
                     valence=spec["valence"],
                     scalarrel=spec["scalarrel"],
                     confinement=PowerConfinement(r0=p["r_dens"], s=p.get("s_dens", 2.0)),
                     wf_confinement=wf_conf,
                     perturbative_confinement=False,
                     txt=None)
    atom.run()
    return atom


def main():
    params = json.load(open(sys.argv[1]))
    outdir = sys.argv[2]
    os.makedirs(outdir, exist_ok=True)
    elements = [el for el in ["Mo", "S", "O", "H"] if el in params]

    t0 = time.time()
    free_props = {}
    for el in elements:
        eigs, hubbards = free_atom_properties(el)
        free_props[el] = (eigs, hubbards)
        print(f"[{time.time()-t0:6.1f}s] {el}: eig={ {k: round(v,5) for k,v in eigs.items()} } "
              f"U={ {k: round(v,5) for k,v in hubbards.items()} }", flush=True)

    atoms = {}
    for el in elements:
        atoms[el] = confined_atom(el, params[el])
        print(f"[{time.time()-t0:6.1f}s] confined {el} done", flush=True)

    cwd = os.getcwd()
    os.chdir(outdir)
    try:
        for i, e1 in enumerate(elements):
            for e2 in elements[i:]:
                off2c = Offsite2cTable(atoms[e1], atoms[e2])
                off2c.run(RMIN, DR, NPTS, superposition="density", xc=XC)
                eigs1, hub1 = free_props[e1]
                off2c.write(eigenvalues=eigs1,
                            hubbardvalues=hub1,
                            occupations=SPECS[e1]["occupations"],
                            spe=0.0)
                print(f"[{time.time()-t0:6.1f}s] SK {e1}-{e2} written", flush=True)
    finally:
        os.chdir(cwd)
    print(f"ALL DONE in {time.time()-t0:.1f}s; files:", sorted(os.listdir(outdir)))


if __name__ == "__main__":
    main()
