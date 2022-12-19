from .system import System
from typing import List

class Atom:
    def __init__(self, idx:int = 0, type:int = 0, x:float = 0.0, y:float = 0.0, z:float = 0.0, mass:float = 0.0, molidx:int = None):
        self.idx:int = idx
        self.type = type
        self.x = x
        self.y = y
        self.z = z
        self.mass:float = mass
        self.molidx = 1000 + molidx

class Bond:
    def __init__(self, idx:int = 0, type:int = 0, atom1:int = 0, atom2:int = 0, length:float = 0):
        self.idx = idx
        self.type = type
        self.atom1 = atom1
        self.atom2 = atom2
        self.length = length

class Angle:
    def __init__(self, idx:int = 0, type:int = 0, atom1:int = 0, atom2:int = 0, atom3:int = 0, area:float = 0.0):
        self.idx = idx
        self.type = type
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.area = area

class Dihedral:
    def __init__(self, idx:int = 0, type:int = 0, atom1:int = 0, atom2:int = 0, atom3:int = 0, atom4:int = 0):
        self.idx = idx
        self.type = type
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4

class Model(System):
    def __init__(self, molidx = 0):
        self.atoms: List[Atom] = []

        self.bonds: List[Bond] = []
        self.angles: List[Angle] = []
        self.dihedrals: List[Dihedral] = []
        self.molidx = 1000 + molidx

    def creat_from_system(self, system:System):
        self.atoms: List[Atom] = system.atoms
        self.bonds: List[Bond] = system.bonds
        self.angles: List[Angle] = system.angles
        self.dihedrals: List[Dihedral] = system.dihedrals


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