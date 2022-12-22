from box import Box
from data import *

class System:
    def __init__(self):
        self.box = Box()
        self.natoms = 0
        self.nbonds = 0
        self.nangles = 0
        self.ndihedrals = 0

        self.nmolecules = 0

        self.atomtypes = 0
        self.bondtypes = 0
        self.angletypes = 0
        self.dihedraltypes = 0

        self.masses = []

        self.atoms: List[Atom] = []
        self.bonds: List[Bond] = []
        self.angles: List[Angle] = []
        self.dihedrals: List[Dihedral] = []  
    
    def update_box(self, box:Box):
        self.box = box
    
    def update_molecule(self, molecule:Molecule):
        self.update_atom(molecule.atoms)
        self.update_bond(molecule.bonds)
        self.update_angle(molecule.angles)
        self.update_dihedral(molecule.dihedrals)

    def update_atom(self, model:List[Atom]):
        for i in range(len(model)):
            item:Atom = model[i]
            item.idx += self.natoms
            item.type += self.atomtypes
            item.molidx += self.nmolecules
            self.atoms.append(item)
            self.natoms += 1
        self.masses.append(model[-1].mass)
        self.atomtypes += 1
    
    def update_bond(self, model:List[Bond]):
        for i in range(len(model)):
            item:Bond = model[i]
            item.idx += self.nbonds
            item.type += self.bondtypes
            self.bonds.append(item)
            self.nbonds += 1
        self.bondtypes += 1

    def update_angle(self, model:List[Angle]):
        for i in range(len(model)):
            item:Angle = model[i]
            item.idx += self.nangles
            item.type += self.angletypes
            self.angles.append(item) 
            self.nangles += 1
        self.angletypes += 1

    def update_dihedral(self, model:List[Dihedral]):
        for i in range(len(model)):
            item:Dihedral = model[i]
            item.idx += self.ndihedrals
            item.type += self.dihedraltypes
            self.dihedrals.append(item)
            self.ndihedrals += 1
        self.dihedraltypes += 1

    def update_box(self, box:Box):
        self.box = box
