# reaxff-crosslinking

reactive molecular dynamics simulations of thermoset polymer crosslinking using ReaxFF in LAMMPS. an end-to-end open-source workflow for generating realistic crosslinked network structures from first principles, built as a demonstration of computational materials science applied to industrially relevant coating and structural polymer systems.

## overview

thermoset polymers form irreversible covalent networks during cure. standard molecular dynamics uses fixed bond topologies and cannot capture the bond formation and breaking that defines crosslinking. this project uses the ReaxFF reactive force field, which describes bonding through continuous bond order functions rather than fixed connectivity, allowing crosslinking to emerge naturally from the interatomic potential without artificial bond insertion algorithms.

two systems are implemented:

| system | resin | hardener/crosslinker | application |
|--------|-------|----------------------|-------------|
| DGEBF/DETDA | bisphenol F diglycidyl ether | diethyltoluenediamine | structural epoxy |
| IPA-NPG polyester / HMMA | isophthalic acid / neopentyl glycol polyester | hexamethoxymethylmelamine | coil coating |

the DGEBF/DETDA system directly replicates the Vashisth et al. (2018) benchmark, providing a validated starting point. the polyester/HMMA system is a novel extension relevant to industrial coating applications where melamine crosslinkers are standard.

## scientific background

### the epoxide-amine reaction

in the DGEBF/DETDA system, the primary amine nitrogen in DETDA attacks the terminal carbon of the epoxide ring in DGEBF. the C-O bond in the strained three-membered ring breaks, a new C-N bond forms, and the oxygen accepts a proton from the N-H group to yield a hydroxyl group. each DETDA molecule carries two amine groups each with two N-H bonds, giving a theoretical functionality of 4 crosslinks per hardener molecule.

in ReaxFF this process emerges from the bond order potential — as N approaches the epoxide C, the C-N bond order rises continuously from 0 while the C-O bond order in the ring drops, with no algorithmic intervention required.

### the polyester-melamine reaction

in the polyester/HMMA system, hydroxyl end groups on the polyester react with the N-CH2-OCH3 groups on the HMMA melamine crosslinker via acid-catalysed transetherification. the methoxy group leaves as methanol (CH3OH) and a new C-O ether bond forms between the HMMA methylene carbon and the polyester oxygen. formation of methanol as a byproduct is tracked in the analysis as mechanistic confirmation that the correct reaction pathway is occurring.

### why ReaxFF

classical force fields (OPLS-AA, COMPASS, PCFF) assign fixed bond topologies at the start of a simulation. simulating crosslinking with these fields requires artificial bond insertion algorithms that periodically search for reactive pairs within a cutoff radius and manually create bonds. ReaxFF eliminates this approximation at the cost of roughly 10-100x greater computational expense, and is most appropriate when reaction pathways, transition states, or the detailed energetics of network formation are of interest.

### force field

the Vashisth et al. (2018) re-optimised CHNO parameter set is used throughout. this field was specifically validated against epoxide ring-opening reaction barriers and thermo-mechanical properties of DGEBF/DETDA, making it the most appropriate available parameter set for these systems. the field covers C, H, O, N and is directly applicable to both systems without reparameterisation.

## results

### DGEBF/DETDA epoxy-amine

50% crosslink conversion achieved in 55 ps during the cure ramp, plateauing as reactive sites become separated in the densifying network.

![epoxy crosslinks](epoxy_amine/epoxy_crosslinks.png)

| property | simulation | experiment | reference |
|----------|-----------|------------|-----------|
| crosslink conversion | 50% | 85-95% | Vashisth 2018 |
| density | ~1.1 g/cm³ | 1.16 g/cm³ | Vashisth 2018 |
| time to 50% conversion | 55 ps | - | - |

### polyester/HMMA coil coating

2.7% crosslink conversion observed, with 2 methanol molecules formed matching 2 crosslink events — mass balance confirmed. the low conversion reflects the high activation barrier of uncatalysed transetherification and the nanosecond simulation timescale.

![polyester crosslinks](polyester_system/polyester_crosslinks.png)

| property | simulation | note |
|----------|-----------|------|
| crosslink conversion | 2.7% | uncatalysed, 105 ps |
| methanol formed | 2 | matches crosslink count |
| mass balance | 73/73 | confirmed |

## limitations

the principal limitation of this approach is the timescale gap between MD simulation (nanoseconds) and real cure cycles (seconds to minutes). conversion plateaus as reactive sites become separated in the densifying network and thermal diffusion is insufficient to bring them together within the simulation window.

Vashisth et al. address this using an accelerated ReaxFF method that applies biased forces to reactive site pairs within a defined distance window, achieving 85-95% conversion. this represents the most significant methodological extension of the current work.

for the polyester/HMMA system, the low conversion additionally reflects the absence of an acid catalyst. industrial cure of melamine-crosslinked coatings uses p-toluenesulfonic acid or dodecylbenzenesulfonic acid to lower the reaction barrier significantly.

## repository structure

```
reaxff-crosslinking/
├── README.md
├── epoxy_amine/
│   ├── dgebf.xyz
│   ├── detda.xyz
│   ├── packbox.inp
│   ├── in.reaxff.epoxy
│   ├── convert_to_LAMMPS.py
│   └── epoxy_analysis.py
├── polyester_system/
│   ├── build_polyester_molecules.py
│   ├── polyester_ipa_npg.xyz
│   ├── hmma.xyz
│   ├── packbox_polyester.inp
│   ├── in.reaxff_polyester
│   ├── convert_polyester_to_lammps.py
│   └── polyester_analysis.py
└── forcefield/
    └── ffield.reax.epoxy
```

## workflow

### epoxy-amine system

```bash
# pack simulation box
packmol < packbox.inp

# convert to lammps format
python convert_to_LAMMPS.py

# run crosslinking simulation
mpirun -np 4 lmp -in in.reaxff.epoxy

# analyse results
python epoxy_analysis.py bonds.reaxff
```

### polyester/HMMA system

```bash
# generate molecule xyz files
python build_polyester_molecules.py

# pack simulation box
packmol < packbox_polyester.inp

# convert to lammps format
python convert_polyester_to_lammps.py

# run crosslinking simulation
mpirun -np 4 lmp -in in.reaxff_polyester

# analyse results
python polyester_analysis.py bonds.reaxff
```

## requirements

- LAMMPS compiled with REAXFF, MOLECULE, MANYBODY, RIGID packages
- Python 3.8+
- RDKit: `conda install -c conda-forge rdkit`
- Packmol: `conda install -c conda-forge packmol`
- ASE: `pip install ase`

GPU acceleration via the LAMMPS Kokkos package with CUDA is strongly recommended for ReaxFF simulations.

## future work

**increase crosslink conversion** — implement the Vashisth et al. accelerated ReaxFF method, applying biased forces to reactive site pairs to overcome the diffusion limitation and achieve >80% conversion.

**property validation** — use the cured network structures for downstream calculations: glass transition temperature (Tg) via NPT density-temperature ramp, Young's modulus via nonequilibrium MD box deformation, thermal conductivity via the Müller-Plathe reverse NEMD method.

**explicit acid catalyst** — include p-toluenesulfonic acid in the polyester/HMMA simulation to lower the transetherification barrier and achieve meaningful conversion.

**multiscale CG→ReaxFF** — use MARTINI coarse-grained MD to equilibrate the system at microsecond timescales, then back-map to atomistic resolution before applying ReaxFF crosslinking. this gives realistic chain configurations and reactive site distributions without the artificial spatial homogeneity of Packmol packing. a natural extension given the author's PhD experience with OPLS-AA/MARTINI back-mapping for polyester coil coating systems.

## references

- Vashisth, A. et al. *Polymer* 2018, 158, 354-363
- Provenzano, M. et al. *ACS Appl. Polym. Mater.* 2025, 7, 4876-4884
- van Duin, A.C.T. et al. *J. Phys. Chem. A* 2001, 105, 9396-9409
- Thompson, A.P. et al. *Comput. Phys. Commun.* 2022, 271, 108171

## author

computational materials scientist with PhD in materials science (University of Surrey, 2018-2022). research focus: polymer blend miscibility, coatings chemistry, multiscale molecular simulation. published work on polyester/melamine surface segregation in coil coating systems using OPLS-AA/MARTINI multiscale workflow.
