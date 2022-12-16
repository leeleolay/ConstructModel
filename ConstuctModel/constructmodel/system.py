import sys
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

        self.atomtypes = 0
        self.bondtypes = 0
        self.anglestypes = 0
        self.dihedraltypes = 0

        self.masstypes = []

        # Initialize atom, bond, angle and dihedral properties
        self.atoms = [Atom]
        self.bonds = [Bond]
        self.angles = [Angle]
        self.dihedrals = [Dihedral]      

class Atom:
    def __init__(self, idx = 0, type = 0, x = 0.0, y = 0.0, z = 0.0, mass = 0.0):
        self.idx = idx
        self.type = type
        self.x = x
        self.y = y
        self.z = z
        self.mass = mass  

class Bond:
    def __init__(self, idx, type, atom1, atom2, length):
        self.idx = idx
        self.type = type
        self.atom1 = atom1
        self.atom2 = atom2
        self.length = length

class Angle:
    def __init__(self, idx, type, atom1, atom2, atom3, area):
        self.idx = idx
        self.type = type
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.area = area

class Dihedral:
    def __init__(self, idx, type, atom1, atom2, atom3, atom4):
        self.idx = idx
        self.type = type
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4

class Molecule(Atom, Bond, Angle, Dihedral):
    def __init__(self, system:System, molidx = 0):
        self.atoms = system.atoms
        self.bonds = system.bonds
        self.angles = system.angles
        self.dihedrals = system.dihedrals
        self.molidx = 1000 + molidx

    @property
    def centroid(self):
        total_mass = 0.0
        centroid_x = 0.0
        centroid_y = 0.0
        centroid_z = 0.0

        for atom in self.atoms:
            total_mass += atom.mass
            centroid_x += atom.x * atom.mass
            centroid_y += atom.y * atom.mass
            centroid_z += atom.z * atom.mass

        centroid_x /= total_mass
        centroid_y /= total_mass
        centroid_z /= total_mass

        return (centroid_x, centroid_y, centroid_z)

    @property
    def geometric_center(self):
        min_x = float('inf')
        min_y = float('inf')
        min_z = float('inf')
        max_x = float('-inf')
        max_y = float('-inf')
        max_z = float('-inf')

        for atom in self.atoms:
            min_x = min(min_x, atom.x)
            min_y = min(min_y, atom.y)
            min_z = min(min_z, atom.z)
            max_x = max(max_x, atom.x)
            max_y = max(max_y, atom.y)
            max_z = max(max_z, atom.z)

        geometric_center_x = (min_x + max_x) / 2
        geometric_center_y = (min_y + max_y) / 2
        geometric_center_z = (min_z + max_z) / 2

        return (geometric_center_x, geometric_center_y, geometric_center_z)

    @property
    def min_x(self):
        min_x = float('inf')
        for atom in self.atoms:
            min_x = min(min_x, atom.x)
        return min_x
    
    @property
    def min_y(self):
        min_y = float('inf')
        for atom in self.atoms:
            min_y = min(min_y, atom.y) 
        return min_y
    
    @property
    def min_z(self):
        min_z = float('inf')
        for atom in self.atoms:
            min_z = min(min_z, atom.z)
        return min_z
    
    @property
    def max_x(self):
        max_x = float('inf')
        for atom in self.atoms:
            max_x = max(max_x, atom.x)
        return max_x
    
    @property
    def max_y(self):
        max_y = float('inf')
        for atom in self.atoms:
            max_y = max(max_y, atom.y)
        return max_y

    @property
    def max_z(self):
        max_z = float('inf')
        for atom in self.atoms:
            max_z = max(max_z, atom.z)
        return max_z

class IO:
    def __init__(self, input_file_name:str, output_file_name:str):
        self.input_file_name = input_file_name
        self.output_file_name = output_file_name
        self._label_tag_atom = []
        self._label_tag_bond = []
        self._label_tag_angle = []
        self._label_tag_dihedral = []

    def read_file(self,system:System):
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
    
    def output(self,system:System):
        with open(self.output_file_name, 'w') as f:
            f.write("LAMMPS Description\n\n")
            f.write("%d atoms\n" % system.natoms)
            f.write("%d bonds\n" % system.nbonds)
            f.write("%d angles\n" % system.nangles)
            f.write("%d dihedrals\n" % system.ndihedrals)
            f.write("\n")
            f.write("%d atom types\n" % system.atomtypes)
            f.write("%d bond types\n" % system.bondtypes)
            f.write("%d angle types\n" % system.anglestypes)
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

# 创建一个 System 对象
system = System()

# 创建一个 IO 对象
io = IO()

# 读取文件，并使用读取到的数据更新 system 对象
io.read_file('capsule.data', system)

# 创建cell分子
cell = Molecule(system, 1)

# 设置系统atom，bond，angle，dihedral类型
system.atomtypes = 3
system.bondtypes = 2
system.anglestypes = 2
system.dihedraltypes = 2

# 设置系统质量类型和参数
system.masstypes = [1.0 for _ in range(system.atomtypes)]

# 写入文件，并把system对象输出到文件中
io.output('init.data',system)