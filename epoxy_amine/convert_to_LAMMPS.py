#taking packmol xyz into lammps format using ase

from ase.io import read, write
from ase.io.lammpsdata import write_lammps_data

atoms = read("system.xyz")
atoms.set_cell([40, 40, 40])
atoms.set_pbc(True)

# ReaxFF needs charges — initialise to 0 (QEq handles them at runtime)
atoms.set_initial_charges([0.0] * len(atoms))

write_lammps_data("system.data", atoms, specorder=["C","H","O","N"],
                  atom_style="charge")
