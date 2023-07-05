from collections import Counter
import networkx as nx
import re
import math
import scipy.spatial as ss
import numpy as np

AtomID = {}
Atombond = {}
BondID = {}
Bondlenth = {}
AngleID = {}
Anglearea = {}
DihedralID = {}

def distance(coord1, coord2):
    # Compute Euclidean distance between two coordinates
    return math.sqrt(sum([(i - j)**2 for i, j in zip(coord1, coord2)]))

def triangle_area(a, b, c):
    # Compute area of triangle given coordinates of three points
    ab = np.array(b) - np.array(a)
    ac = np.array(c) - np.array(a)
    return np.linalg.norm(np.cross(ab, ac))/2.0

# Read box.txt
with open("box.txt", 'r') as box_file:
    box_ranges = [list(map(float, line.strip().split(','))) for line in box_file.readlines()]
    xlo, xhi = box_ranges[0]
    ylo, yhi = box_ranges[1]
    zlo, zhi = box_ranges[2]

# Read coord.txt
with open('coord.txt', 'r') as coord_file:
    for idx, line in enumerate(coord_file.readlines()):
        x, y, z = map(float, line.strip().split())
        AtomID[str(idx+1)] = [x, y, z]

# Read bond.txt
with open('bond.txt', 'r') as bond_file:
    for idx, line in enumerate(bond_file.readlines()):
        atom1, atom2 = map(int, line.strip().split())
        BondID[str(idx+1)] = [atom1, atom2]
        Bondlenth[str(idx+1)] = distance(AtomID[str(atom1)], AtomID[str(atom2)])

# Read angle.txt
with open('angle.txt', 'r') as angle_file:
    for idx, line in enumerate(angle_file.readlines()):
        atom1, atom2, atom3 = map(int, line.strip().split())
        AngleID[str(idx+1)] = [atom1, atom2, atom3]
        # Calculating area of triangle formed by the three atoms
        Anglearea[str(idx+1)] = triangle_area(AtomID[str(atom1)], AtomID[str(atom2)], AtomID[str(atom3)])

# Create dihedral
for i in range(1,len(BondID)+1):
    tmp = []
    for j in range(len(Atombond[str(BondID[str(i)][0])])):
        tmp.append(Atombond[str(BondID[str(i)][0])][j])
    for j in range(len(Atombond[str(BondID[str(i)][1])])):
        tmp.append(Atombond[str(BondID[str(i)][1])][j])
    for j in tmp:
        if j == BondID[str(i)][0] or j == BondID[str(i)][1]:
            tmp.remove(j)
    d1 = Counter(tmp)
    for key in list(d1.keys()):
        if d1[key] == 1:
            del d1[key]
            continue
    dihedrallisttmp = [BondID[str(i)][0],BondID[str(i)][1]]
    for key in list(d1.keys()):
        dihedrallisttmp.append(key)
    DihedralID[str(i)] = dihedrallisttmp


# Write to LAMMPS data file
with open("capsule_step1.data", "w") as lammps_file:
    lammps_file.write(
        f"LAMMPS data file\n\n{len(AtomID)} atoms\n1 atom types\n{len(BondID)} bonds\n1 bond types\n{len(AngleID)} angles\n1 angle types\n{len(DihedralID)} dihedrals\n1 dihedral types\n"
        f"{xlo} {xhi} xlo xhi\n{ylo} {yhi} ylo yhi\n{zlo} {zhi} zlo zhi\n\nMasses\n\n1 1\n\n"
        f"Atoms # full\n"
    )

    for key, value in AtomID.items():
        lammps_file.write(f"{key} 1 1 0.0 {value[0]} {value[1]} {value[2]}\n")
    
    lammps_file.write("\nBonds\n\n")
    for key, value in BondID.items():
        lammps_file.write(f"{key} 1 {value[0]} {value[1]} {Bondlenth[key]}\n")

    lammps_file.write("\nAngles\n\n")
    for key, value in AngleID.items():
        lammps_file.write(f"{key} 1 {value[0]} {value[1]} {value[2]} {Anglearea[key]}\n")

    lammps_file.write("\nDihedrals\n\n")
    for key, value in DihedralID.items():
        lammps_file.write(f"{key} 1 {value[2]} {value[1]} {value[0]} {value[2]}\n")
