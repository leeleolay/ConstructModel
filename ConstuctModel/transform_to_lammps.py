import os
import numpy as np
from collections import Counter

def distance(coord1, coord2):
    coord1, coord2 = np.array(coord1), np.array(coord2)
    return np.linalg.norm(coord1 - coord2)

def triangle_area(a, b, c):
    a, b, c = np.array(a), np.array(b), np.array(c)
    return 0.5 * np.linalg.norm(np.cross(b - a, c - a))

def load_file(filename, data_type):
    with open(filename, 'r') as file:
        return {str(i+1): list(map(data_type, line.strip().split())) for i, line in enumerate(file)}

# Load data
box_ranges = load_file("box.txt", float)
atoms = load_file("coord.txt", float)
bonds = load_file("bond.txt", int)
angles = load_file("angle.txt", int)

xlo, xhi = box_ranges['1']
ylo, yhi = box_ranges['2']
zlo, zhi = box_ranges['3']

bond_lengths = {bond_id: distance(atoms[str(bond[0])], atoms[str(bond[1])]) for bond_id, bond in bonds.items()}

angle_areas = {angle_id: triangle_area(atoms[str(angle[0])], atoms[str(angle[1])], atoms[str(angle[2])]) for angle_id, angle in angles.items()}

# Building connectedatoms and DihedralID dictionaries
connectedatoms = {str(atom_id): [bond[1 - bond.index(atom_id)] for bond in bonds.values() if atom_id in bond] for atom_id in range(1, len(atoms) + 1)}

dihedral_IDs = {}
for bond_id, bond in bonds.items():
    common_atoms = [atom for atom in connectedatoms[str(bond[0])] + connectedatoms[str(bond[1])] if atom != bond[0] and atom != bond[1]]
    counts = Counter(common_atoms)
    dihedrallisttmp = [bond[0], bond[1]] + [atom for atom, count in counts.items() if count > 1]
    dihedral_IDs[bond_id] = dihedrallisttmp

# Generate LAMMPS data file
current_dir_name = os.path.basename(os.getcwd())
with open(f"{current_dir_name}.data", "w") as lammps_file:
    lammps_file.write(
        f"LAMMPS data file\n\n{len(atoms)} atoms\n1 atom types\n{len(bonds)} bonds\n1 bond types\n{len(angles)} angles\n1 angle types\n{len(dihedral_IDs)} dihedrals\n1 dihedral types\n\n"
        f"{xlo} {xhi} xlo xhi\n{ylo} {yhi} ylo yhi\n{zlo} {zhi} zlo zhi\n\nMasses\n\n1 1\n\n"
    )

    lammps_file.write("\nAtoms\n\n")
    for atom_id, atom in atoms.items():
        lammps_file.write(f"{atom_id} 1001 1 {atom[0]} {atom[1]} {atom[2]}\n")

    lammps_file.write("\nBonds\n\n")
    for bond_id, bond in bonds.items():
        lammps_file.write(f"{bond_id} 1 {bond[0]} {bond[1]} {bond_lengths[bond_id]}\n")

    lammps_file.write("\nAngles\n\n")
    for angle_id, angle in angles.items():
        lammps_file.write(f"{angle_id} 1 {angle[0]} {angle[1]} {angle[2]} {angle_areas[angle_id]}\n")

    lammps_file.write("\nDihedrals\n\n")
    for dihedral_id, dihedral in dihedral_IDs.items():
        lammps_file.write(f"{dihedral_id} 1 {dihedral[0]} {dihedral[1]} {dihedral[2]} {dihedral[3]}\n")
