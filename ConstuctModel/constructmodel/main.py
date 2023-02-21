import sys
import scipy.spatial as ss
import numpy as np
from typing import Dict
import copy

from system import System
from myio import MyIO
from data import *
from box import Box

def create_box(params:Dict):
    box = Box()
    box.xlo = params['xlo']
    box.xhi = params['xhi']
    box.ylo = params['ylo']
    box.yhi = params['yhi']
    box.zlo = params['zlo']
    box.zhi = params['zhi']
    return box

def create_mol(io:MyIO):
    system_blank = System()
    moleucle = Molecule()
    io.read_LAMMPS(system_blank)
    cell = moleucle.creat_from_system(system_blank)
    return cell

def create_mols(input_model: Molecule, params: Dict) -> Molecule:
    """
    Create a new molecule model by copying and modifying the original model multiple times.

    Parameters:
    - input_model: the original molecule model to be copied and modified
    - params: a dictionary containing the configuration parameters for the function

    Returns:
    - a new molecule model created by copying and modifying the original model
    """
    num = params['num_cell_x'] * params['num_cell_y']
    output_model = Molecule()

    for i in range(num):
        # Create a copy of the original model
        cell = copy.deepcopy(input_model)

        # Modify the xyz coordinates and indices of the atoms in the copied model
        dist_trans_x = params['x_center'][i]
        dist_trans_y = params['y_center'][i]
        dist_trans_z = params['z_center'][i]

        for atom in cell.atoms:
            atom.x += dist_trans_x
            atom.y += dist_trans_y
            atom.z += dist_trans_z
            atom.molidx += i
            atom.idx += len(cell.atoms)*i

        # Modify the indices of the bonds, angles and dihedrals in the copied model
        for bond in cell.bonds:
            bond.idx += len(cell.bonds)*i
            bond.atom1 += len(cell.atoms)*i
            bond.atom2 += len(cell.atoms)*i 

        for angle in cell.angles:
            angle.idx += len(cell.angles)*i
            angle.atom1 += len(cell.atoms)*i
            angle.atom2 += len(cell.atoms)*i
            angle.atom3 += len(cell.atoms)*i

        for dihedral in cell.dihedrals:
            dihedral.idx += len(cell.dihedrals)*i
            dihedral.atom1 += len(cell.atoms)*i
            dihedral.atom2 += len(cell.atoms)*i
            dihedral.atom3 += len(cell.atoms)*i
            dihedral.atom4 += len(cell.atoms)*i

        # Add the atoms, bonds, angles and dihedrals from the copied model to the final model
        output_model.atoms.extend(cell.atoms)
        output_model.bonds.extend(cell.bonds)
        output_model.angles.extend(cell.angles)
        output_model.dihedrals.extend(cell.dihedrals)

    return output_model

def creat_wall(params:Dict) -> List[Atom]:
    wall:List[Atom] = []
    atom = Atom()
    atom.idx = 1
    atom.mass = 2.0
    atom.molidx = 1001
    atom.type = 1
    for i in np.arange(params['wall_up_xlo'],params['wall_up_xhi'],params['wall_rho']):
        for j in np.arange(params['wall_up_ylo'],params['wall_up_yhi'],params['wall_rho']):
            for k in np.arange(params['wall_up_zlo'],params['wall_up_zhi'],params['wall_rho']):
                atom.x = i
                atom.y = j
                atom.z = k
                wall.append(copy.copy(atom))
                atom.idx += 1
    atom.molidx += 1
    for i in np.arange(params['wall_down_xlo'],params['wall_down_xhi'],params['wall_rho']):
        for j in np.arange(params['wall_down_ylo'],params['wall_down_yhi'],params['wall_rho']):
            for k in np.arange(params['wall_down_zlo'],params['wall_down_zhi'],params['wall_rho']):
                atom.x = i
                atom.y = j
                atom.z = k
                wall.append(copy.copy(atom))
                atom.idx += 1
    return wall

def creat_particle(params:Dict):
    particles:List[Atom] = []
    atom = Atom()
    atom.idx = 1
    atom.mass = 3.0
    atom.molidx = 1001
    atom.type = 1
    for i in np.arange(params['particle_xlo'],params['particle_xhi'],params['particle_rho']):
        for j in np.arange(params['particle_ylo'],params['particle_yhi'],params['particle_rho']):
            for k in np.arange(params['particle_zlo'],params['particle_zhi'],params['particle_rho']):
                atom.x = i
                atom.y = j
                atom.z = k
                particles.append(copy.copy(atom))
                atom.idx += 1
    return particles

def creat_particle_single(params:Dict):
    particles:List[Atom] = []
    atom = Atom()
    atom.idx = 1
    atom.mass = 3.0
    atom.molidx = 1001
    atom.type = 1
    for i in range(len(params['x'])):
        atom.x = params['x'][i]
        atom.y = params['y'][i]
        atom.z = params['z'][i]
        particles.append(copy.copy(atom))
        atom.idx += 1
    return particles

def init_io(io:MyIO):
    params = dict()
    args = MyIO.parse_args()
    params:Dict = MyIO.load_json(args.jsonfile)
    process_params(params)
    io.input_file_name = args.input
    io.output_file_name = args.output
    return params


def process_params(params):
    # 处理box参数
    box_length_x = params['cell']['num_cell_x'] * (params['cell']['gap_of_cells']+params['cell']['length_x'])
    box_length_y = params['cell']['num_cell_y'] * (params['cell']['gap_of_cells']+params['cell']['length_y'])
    box_xlo = -box_length_x/2
    box_xhi = box_length_x/2
    box_ylo = -box_length_y/2
    box_yhi = box_length_y/2  
    params['box'].update(
        xlo = box_xlo,
        xhi = box_xhi,
        ylo = box_ylo,
        yhi = box_yhi
    )

    # 处理wall参数
    wall_up_xlo, wall_down_xlo = [params['box']['xlo']] * 2
    wall_up_xhi, wall_down_xhi = [params['box']['xhi']] * 2
    wall_up_ylo, wall_down_ylo = [params['box']['ylo']] * 2
    wall_up_yhi, wall_down_yhi = [params['box']['yhi']] * 2
    wall_up_zlo = params['box']['zhi']-params['wall']['dist_from_box']-params['wall']['wall_rho']
    wall_up_zhi = params['box']['zhi']-params['wall']['dist_from_box']
    wall_down_zlo = params['box']['zlo']+params['wall']['dist_from_box']
    wall_down_zhi = params['box']['zlo']+params['wall']['dist_from_box']+params['wall']['wall_rho']
    params['wall'].update(
        wall_up_xlo = wall_up_xlo,
        wall_down_xlo = wall_down_xlo,
        wall_up_xhi = wall_down_xhi,
        wall_down_xhi = wall_down_xhi,
        wall_up_ylo = wall_up_ylo,
        wall_down_ylo = wall_down_ylo,
        wall_up_yhi = wall_up_yhi,
        wall_down_yhi = wall_down_yhi,
        wall_up_zlo = wall_up_zlo,
        wall_down_zlo = wall_down_zlo,
        wall_up_zhi = wall_up_zhi,
        wall_down_zhi = wall_down_zhi
    )

    # 处理cell参数
    x_center = []
    y_center = []
    z_center = []
    for i in range(params['cell']['num_cell_x']):
        for j in range(params['cell']['num_cell_y']): 
            for k in range(params['cell']['num_cell_z']):
                x_center.append(params['box']['xlo'] + (params['cell']['gap_of_cells']+params['cell']['length_x']) * (0.5 + i) )
                y_center.append(params['box']['ylo'] + (params['cell']['gap_of_cells']+params['cell']['length_y']) * (0.5 + j) )
                z_center.append(0.0)
    params['cell'].update(
        x_center = x_center,
        y_center = y_center,
        z_center = z_center
    )

    # 处理particle参数
    params['particle'].update(
        particle_xlo = params['box']['xlo'],
        particle_xhi = params['box']['xhi'],
        particle_ylo = params['box']['ylo'],
        particle_yhi = params['box']['yhi'],
    )

def main():
    # 初始化
    io = MyIO('capsule.data', 'init.data')
    params = init_io(io)

    # 构造box
    box = create_box(params['box'])
    # 构造cell模型
    cell = create_mol(io)
    cells = create_mols(cell, params['cell'])
    # 构造wall
    #wall = creat_wall(params['wall'])
    # 构造dpd粒子
    particle = creat_particle(params['particle'])
    #particle_single = creat_particle_single(params['particle_single'])

    # 构建系统
    system = System()
    system.update_box(box)
    system.update_molecule(cells)
    #system.update_atoms(wall)
    system.update_atoms(particle)
    #system.update_atoms(particle_single)

    # 写入文件，并把system对象输出到文件中
    io.write_LAMMPS(system)

if __name__ == "__main__":
    main()