import scipy.spatial as ss
import numpy as np
import re

class Box:
    def __init__(self, xlo=0, xhi=0, ylo=0, yhi=0, zlo=0, zhi=0):
        self.xlo = xlo
        self.xhi = xhi
        self.ylo = ylo
        self.yhi = yhi
        self.zlo = zlo
        self.zhi = zhi

    @property
    def size(self):
        return (self.xhi - self.xlo, self.yhi - self.ylo, self.zhi - self.zlo)

class System:
    def __init__(self, box:Box):
        # Initialize system properties
        self.box = box
        self.natoms = 0
        self.nbonds = 0
        self.nangles = 0
        self.ndihedrals = 0
        # Initialize atom, bond, angle and dihedral properties
        self.atoms = []
        self.bonds = []
        self.angles = []
        self.dihedrals = []      

class Atom:
    def __init__(self, idx, x, y, z, mass, type):
        self.idx = idx
        self.x = x
        self.y = y
        self.z = z
        self.mass = mass
        self.type = type

class Bond:
    def __init__(self, idx, atom1, atom2, bond_length):
        self.idx = idx
        self.atom1 = atom1
        self.atom2 = atom2
        self.bond_length = bond_length

class Angle:
    def __init__(self, idx, atom1, atom2, atom3, area):
        self.idx = idx
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.area = area

class Dihedral:
    def __init__(self, idx, atom1, atom2, atom3, atom4):
        self.idx = idx
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4

class IO:
    def __init__(self, input_file_name:str, output_file_name:str, system=System()):
        self.input_file_name = input_file_name
        self.output_file_name = output_file_name
        self.system = system
        self._label_tag_atom = []
        self._label_tag_bond = []
        self._label_tag_angle = []
        self._label_tag_dihedral = []

    def read_file(self,system):
        with open(self.input_file_name, 'r') as f:
            ncount = 0
            for line in f.readlines():
                ncount = ncount + 1
                if re.search('.*xlo xhi',line):
                    line = line.strip().split()
                    self.system.box.xlo = float(line[0])
                    self.system.box.xhi = float(line[1])
                elif re.search('.*ylo yhi',line):
                    line = line.strip().split()
                    self.system.box.xlo= float(line[0])
                    self.system.box.yhi = float(line[1])
                elif re.search('.*zlo zhi',line):
                    line = line.strip().split()
                    self.system.box.zlo = float(line[0])
                    self.system.box.zhi = float(line[1])
                elif re.search('.*atoms',line):
                    line = line.strip().split()
                    self.system.natoms = int(line[0])
                elif re.search('.*bonds',line):
                    line = line.strip().split()
                    self.system.nbonds = int(line[0])
                elif re.search('.*angles',line):
                    line = line.strip().split()
                    self.system.nangles = int(line[0])
                elif re.search('.*dihedrals',line):
                    line = line.strip().split()
                    self.system.ndihedrals = int(line[0])
                elif re.search('Atoms',line):
                    self._label_tag_atom.append(ncount + 2)
                    self._label_tag_atom.append(ncount + self.system.natoms + 1)
                elif re.search('Bonds',line):
                    self._label_tag_bond.append(ncount + 2)
                    self._label_tag_bond.append(ncount + self.system.nbonds + 1)
                elif re.search('Angles',line):
                    self._label_tag_angle.append(ncount + 2)
                    self._label_tag_angle.append(ncount + self.system.nangles + 1)
                elif re.search('Dihedrals',line):
                    self._label_tag_dihedral.append(ncount + 2)
                    self._label_tag_dihedral.append(ncount + self.system.ndihedrals + 1)
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
                    self.system.atoms.append(atom)
                if ncount >= self._label_tag_bond[0] and ncount <= self._label_tag_bond[1]:
                    line = line.strip().split()
                    bond.idx = line[0]
                    bond.atom1 = line[2]
                    bond.atom1 = line[3]
                    bond.bond_length = line[4]
                    self.system.angles.append(bond)
                if ncount >= self._label_tag_angle[0] and ncount <= self._label_tag_angle[1]:
                    line = line.strip().split()
                    angle.idx = line[0]
                    angle.atom1 = line[2]
                    angle.atom2 = line[3]
                    angle.atom3 = line[4]
                    angle.area = line[5]
                    self.system.angles.append(angle)
                if ncount >= self._label_tag_angle[0] and ncount <= self._label_tag_angle[1]:
                    line = line.strip().split()
                    dihedral.idx = line[0]
                    dihedral.atom1 = line[2]
                    dihedral.atom2 = line[3]
                    dihedral.atom3 = line[4]
                    dihedral.atom4 = line[5]
                    self.system.dihedrals.append(dihedral)
    
    def output(self,system):
        with open(self.output_file_name, 'w') as f:
            f.write("LAMMPS Description\n\n")
            f.write("%d atoms\n" % cnatom)
            f.write("%d bonds\n" % cnbond)
            f.write("%d angles\n" % cnangle)
            f.write("%d dihedrals\n" % cndihedral)
            f.write("\n")
            f.write("%d atom types\n" % cnparticle)
            f.write("%d bond types\n" % cnparticle)
            f.write("%d angle types\n" % cnparticle)
            f.write("%d dihedral types\n" % cnparticle)
            f.write("\n")
            f.write("%f %f xlo xhi\n" % (xlo, xhi))
            f.write("%f %f ylo yhi\n" % (ylo, yhi))
            f.write("%f %f zlo zhi\n" % (zlo, zhi))
            f.write("\n")
            f.write("Masses\n\n")
            for i in range(1, self.nparticle):
                f.write("%d %f\n" % (i, self.mass[i]))
            f.write("\n")
            f.write("Atoms\n\n")
            for i in range(1, self.natoms+1):
                f.write("%d %d %d %d %f %f %f\n" % (i, self.molecules[i], self.type[i], self.particle[i], self.x[i], self.y[i], self.z[i]))
            f.write("\n")
            f.write("Bonds\n\n")
            for i in range(1, self.nbonds+1):
                f.write("%d %d %d %d %d\n" % (i, self.bond_type[i], self.bond_particle1[i], self.bond_particle2[i], self.bond_length[i]))
            f.write("\n")
            f.write("Angles\n\n")
            for i in range(1, self.angles+1):
                f.write("%d %d %d %d %d\n" % (i, self.bond_type[i], self.bond_particle1[i], self.bond_particle2[i], self.bond_length[i]))
            f.write("\n")
            f.write("Dihedrals\n\n")
            for i in range(1, self.ndihedrals+1):
                f.write("%d %d %d %d %d\n" % (i, self.bond_type[i], self.bond_particle1[i], self.bond_particle2[i], self.bond_length[i]))
            f.write("\n")

# 创建一个 System 对象
system = System()

# 创建一个 IO 对象
io = IO()

# 读取文件，并使用读取到的数据更新 system 对象
io.read_file('capsule.data', system)

# 写入文件，并把system对象输出到文件中
io.output('init.data',system)