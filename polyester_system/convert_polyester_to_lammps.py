# convert polyester to lammps using ase 

from ase.io import read
from ase.io.lammpsdata import write_lammps_data

#specify object as xyz coordinates
atoms = read("system_polyester.xyz")
atoms.set_cell([50, 50, 50])
atoms.set_pbc(True)
atoms.set_initial_charges([0.0] * len(atoms))

write_lammps_data("system_polyester.data", atoms,
                  specorder=["C", "H", "O", "N"],
                  atom_style="charge")

