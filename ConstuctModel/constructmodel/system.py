from .box import Box
from .data import *

class System:
    def __init__(self, box:Box = None):
        # Initialize system properties
        self.box = box
        self.natoms = 0
        self.nbonds = 0
        self.nangles = 0
        self.ndihedrals = 0

        self.atomtypes = 0
        self.bondtypes = 0
        self.angletypes = 0
        self.dihedraltypes = 0

        self.masstypes = 0

        # Initialize atom, bond, angle and dihedral properties
        self.atoms: List[Atom] = []
        self.bonds: List[Bond] = []
        self.angles: List[Angle] = []
        self.dihedrals: List[Dihedral] = []  
    
    def update_box(self, box:Box):
        self.box = box
    
    def update_system(self, model:Model):
        self.update_atom(model.atoms)
        self.update_bond(model.bonds)
        self.update_angle(model.angles)
        self.update_dihedral(model.dihedrals)

    def update_atom(self, model:List[Atom]):
        for i in range(model.__len__):
            item:Atom = model[i]
            item.idx += self.natoms
            item.type += self.atomtypes
            self.atoms.append(item)
            self.natoms += 1
        self.masstypes += 1
        self.atomtypes += 1
    
    def update_bond(self, model:List[Bond]):
        for i in range(model.__len__):
            item:Bond = model[i]
            item.idx += self.nbonds
            item.type += self.bondtypes
            self.bonds.append(item)
            self.nbonds += 1
        self.bondtypes += 1

    def update_angle(self, model:List[Angle]):
        for i in range(model.__len__):
            item:Angle = model[i]
            item.idx += self.nangles
            item.type += self.angletypes
            self.angles.append(item) 
            self.nangles += 1
        self.angletypes += 1

    def update_dihedral(self, model:List[Dihedral]):
        for i in range(model.__len__):
            item:Dihedral = model[i]
            item.idx += self.ndihedrals
            item.type += self.dihedraltypes
            self.dihedrals.append(item)
            self.ndihedrals += 1
        self.dihedraltypes += 1
