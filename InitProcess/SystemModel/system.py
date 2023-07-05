import copy

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
    
    def update_molecules(self, molecule:Molecule):
        self.update_atoms(molecule.atoms)
        self.update_bonds(molecule.bonds)
        self.update_angles(molecule.angles)
        self.update_dihedrals(molecule.dihedrals)

    def update_atoms(self, model: List[Atom]) -> None:
        atom_index = copy.deepcopy(self.natoms)
        atom_styles = copy.deepcopy(self.atomtypes)
        atom_molidx = copy.deepcopy(self.nmolecules)

        nummol = set()

        for atom in copy.deepcopy(model):
            atom.idx += atom_index
            atom.type += atom_styles
            atom.molidx += atom_molidx

            self.atoms.append(atom)
            self.natoms += 1
            nummol.add(atom.molidx)

        self.masses.append(model[-1].mass)
        self.atomtypes += 1
        self.nmolecules += len(nummol)
    
    def update_bonds(self, model:List[Bond]):

        bond_index = copy.deepcopy(self.nbonds)
        bond_styles = copy.deepcopy(self.bondtypes)

        for bond in copy.deepcopy(model):

            bond.idx += bond_index
            bond.type += bond_styles

            self.bonds.append(bond)
            self.nbonds += 1

        self.bondtypes += 1

    def update_angles(self, model:List[Angle]):

        angle_index = copy.deepcopy(self.nangles)
        angle_styles = copy.deepcopy(self.angletypes)

        for angle in copy.deepcopy(model):

            angle.idx += angle_index
            angle.type += angle_styles

            self.angles.append(angle) 
            self.nangles += 1
        
        self.angletypes += 1

    def update_dihedrals(self, model:List[Dihedral]):

        dihedral_index = copy.deepcopy(self.ndihedrals)
        dihedral_styles = copy.deepcopy(self.dihedraltypes)

        for dihedral in copy.deepcopy(model):

            dihedral.idx += dihedral_index
            dihedral.type += dihedral_styles

            self.dihedrals.append(dihedral)
            self.ndihedrals += 1

        self.dihedraltypes += 1
