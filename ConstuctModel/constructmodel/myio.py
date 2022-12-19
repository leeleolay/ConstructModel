import re
import argparse
import json

from .system import System
from .box import Box
from .data import *

class MyIO:
    def __init__(self, input_file_name:str, output_file_name:str):
        self.input_file_name = input_file_name
        self.output_file_name = output_file_name
        self._label_tag_atom = []
        self._label_tag_bond = []
        self._label_tag_angle = []
        self._label_tag_dihedral = []

    def read_file_LAMMPS(self, system:System):
        with open(self.input_file_name, 'r') as f:
            ncount = 0
            for line in f.readlines():
                ncount = ncount + 1
                if re.search('.*xlo xhi',line):
                    line = line.strip().split()
                    system.box.xlo = float(line[0])
                    system.box.xhi = float(line[1])
                elif re.search('.*ylo yhi',line):
                    line = line.strip().split()
                    system.box.xlo= float(line[0])
                    system.box.yhi = float(line[1])
                elif re.search('.*zlo zhi',line):
                    line = line.strip().split()
                    system.box.zlo = float(line[0])
                    system.box.zhi = float(line[1])
                elif re.search('.*atoms',line):
                    line = line.strip().split()
                    system.natoms = int(line[0])
                elif re.search('.*bonds',line):
                    line = line.strip().split()
                    system.nbonds = int(line[0])
                elif re.search('.*angles',line):
                    line = line.strip().split()
                    system.nangles = int(line[0])
                elif re.search('.*dihedrals',line):
                    line = line.strip().split()
                    system.ndihedrals = int(line[0])
                elif re.search('Atoms',line):
                    self._label_tag_atom.append(ncount + 2)
                    self._label_tag_atom.append(ncount + system.natoms + 1)
                elif re.search('Bonds',line):
                    self._label_tag_bond.append(ncount + 2)
                    self._label_tag_bond.append(ncount + system.nbonds + 1)
                elif re.search('Angles',line):
                    self._label_tag_angle.append(ncount + 2)
                    self._label_tag_angle.append(ncount + system.nangles + 1)
                elif re.search('Dihedrals',line):
                    self._label_tag_dihedral.append(ncount + 2)
                    self._label_tag_dihedral.append(ncount + system.ndihedrals + 1)
        with open(self.input_file_name,'r') as f:
            ncount = 0
            atom = Atom()
            bond = Bond()
            angle = Angle()
            dihedral = Dihedral()
            for line in f.readlines():
                ncount = ncount + 1
                if ncount >= self._label_tag_atom[0] and ncount <= self._label_tag_atom[1]:
                    line = line.strip().split()
                    atom.idx = line[0]
                    atom.x = line[3]
                    atom.y = line[4]
                    atom.z = line[5]
                    system.atoms.append(atom)
                if ncount >= self._label_tag_bond[0] and ncount <= self._label_tag_bond[1]:
                    line = line.strip().split()
                    bond.idx = line[0]
                    bond.atom1 = line[2]
                    bond.atom1 = line[3]
                    bond.length = line[4]
                    system.angles.append(bond)
                if ncount >= self._label_tag_angle[0] and ncount <= self._label_tag_angle[1]:
                    line = line.strip().split()
                    angle.idx = line[0]
                    angle.atom1 = line[2]
                    angle.atom2 = line[3]
                    angle.atom3 = line[4]
                    angle.area = line[5]
                    system.angles.append(angle)
                if ncount >= self._label_tag_angle[0] and ncount <= self._label_tag_angle[1]:
                    line = line.strip().split()
                    dihedral.idx = line[0]
                    dihedral.atom1 = line[2]
                    dihedral.atom2 = line[3]
                    dihedral.atom3 = line[4]
                    dihedral.atom4 = line[5]
                    system.dihedrals.append(dihedral)
    
    def output_LAMMPS(self, system:System):
        with open(self.output_file_name, 'w') as f:
            f.write("LAMMPS Description\n\n")
            f.write("%d atoms\n" % system.natoms)
            f.write("%d bonds\n" % system.nbonds)
            f.write("%d angles\n" % system.nangles)
            f.write("%d dihedrals\n" % system.ndihedrals)
            f.write("\n")
            f.write("%d atom types\n" % system.atomtypes)
            f.write("%d bond types\n" % system.bondtypes)
            f.write("%d angle types\n" % system.angletypes)
            f.write("%d dihedral types\n" % system.dihedraltypes)
            f.write("\n")
            f.write("%f %f xlo xhi\n" % (system.box.xlo, system.box.xhi))
            f.write("%f %f ylo yhi\n" % (system.box.ylo, system.box.yhi))
            f.write("%f %f zlo zhi\n" % (system.box.zlo, system.box.zhi))
            f.write("\n")
            f.write("Masses\n\n")
            for i in range(1, system.atomtypes):
                f.write("%d %f\n" % (i, system.masstypes[i-1]))
            f.write("\n")
            f.write("Atoms\n\n")
            for i in range(system.natoms):
                f.write("%d %d %d %f %f %f\n" % (system.atoms[i].idx, system.atoms[i].molidx, system.atoms[i].type, system.atoms[i].x, system.atoms[i].y, system.atoms[i].z))
            f.write("\n")
            f.write("Bonds\n\n")
            for i in range(system.nbonds):
                f.write("%d %d %d %d %f\n" % (system.bonds[i].idx, system.bonds[i].type, system.bonds[i].atom1, system.bonds[i].atom2, system.bonds[i].length))
            f.write("\n")
            f.write("Angles\n\n")
            for i in range(system.angles):
                f.write("%d %d %d %d %d %f\n" % (system.angles[i].idx, system.angles[i].type, system.angles[i].atom1, system.angles[i].atom2, system.angles[i].atom3, system.angles[i].area))
            f.write("\n")
            f.write("Dihedrals\n\n")
            for i in range(system.ndihedrals):
                f.write("%d %d %d %d %d %d\n" % (system.dihedrals[i].idx, system.dihedrals[i].type, system.dihedrals[i].atom1, system.dihedrals[i].atom2, system.dihedrals[i].atom3, system.dihedrals[i].atom4))
            f.write("\n")

    def load_json(json_file):
        with open(json_file, "r") as f:
            data = json.load(f)
        return data

    def parse_args():
        parser = argparse.ArgumentParser()
        parser.add_argument("json_file", type=str, help="a JSON file")
        parser.add_argument("-i","--input",required=True, help="input file")
        parser.add_argument("--output", type=str, default="output.txt", help="output file")
        return parser.parse_args()