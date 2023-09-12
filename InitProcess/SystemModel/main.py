import numpy as np
from typing import Dict, List
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

def creat_walls(params:Dict) -> List[Atom]:
    wall:List[Atom] = []
    atom = Atom()
    atom.idx = 1
    atom.mass = 2.0
    atom.molidx = 1001
    atom.type = 1
    for i in np.arange(params['wallup_xlo'],params['wallup_xhi'],params['wall_rho']):
        for j in np.arange(params['wallup_ylo'],params['wallup_yhi'],params['wall_rho']):
            for k in np.arange(params['wallup_zlo'],params['wallup_zhi'],params['wall_rho']):
                atom.x = i
                atom.y = j
                atom.z = k
                wall.append(copy.copy(atom))
                atom.idx += 1
    atom.molidx += 1
    for i in np.arange(params['walldown_xlo'],params['walldown_xhi'],params['wall_rho']):
        for j in np.arange(params['walldown_ylo'],params['walldown_yhi'],params['wall_rho']):
            for k in np.arange(params['walldown_zlo'],params['walldown_zhi'],params['wall_rho']):
                atom.x = i
                atom.y = j
                atom.z = k
                wall.append(copy.copy(atom))
                atom.idx += 1
    return wall

def creat_particles(params:Dict):
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

def generate_initial_configuration(num_particles, region, target_average_speed, target_average_direction, seed=None):
    if seed is not None:
        np.random.seed(seed)

    x_range = region[1] - region[0]
    y_range = region[3] - region[2]
    z_range = region[5] - region[4]

    positions = np.random.rand(num_particles, 3) * [x_range, y_range, z_range] + [region[0], region[2], region[4]]
    initial_velocities = np.random.normal(0, 1, (num_particles, 3))
    average_velocity = np.mean(initial_velocities, axis=0)
    
    initial_velocities -= average_velocity  # 让所有粒子的平均速度为0
    normalized_target_direction = target_average_direction / np.linalg.norm(target_average_direction)
    target_average_velocity = target_average_speed * normalized_target_direction
    initial_velocities += target_average_velocity  # 将所有粒子的平均速度设置为指定的值

    return positions, initial_velocities

def creat_particle_random(params: Dict):
    particles: List[Atom] = []
    
    num_particles = int(((params['particle_xhi'] - params['particle_xlo']) / params['particle_rho']) *
                        ((params['particle_yhi'] - params['particle_ylo']) / params['particle_rho']) *
                        ((params['particle_zhi'] - params['particle_zlo']) / params['particle_rho']))
    
    region = np.array([params['particle_xlo'], params['particle_xhi'],
                       params['particle_ylo'], params['particle_yhi'],
                       params['particle_zlo'], params['particle_zhi']])
    
    positions, velocities = generate_initial_configuration(num_particles, region, params['target_average_speed'], params['target_average_direction'], params['seed'])

    for i, (position, velocity) in enumerate(zip(positions, velocities)):
        atom = Atom()
        atom.idx = i + 1
        atom.mass = 3.0
        atom.molidx = 1001
        atom.type = 1
        atom.x, atom.y, atom.z = position
        atom.vx, atom.vy, atom.vz = velocity
        particles.append(atom)

    return particles

def create_bonds_between_cell_and_wall(params:Dict, system:System):
    cell_top_ids = params['cell_tom']
    cell_bottom_ids = params['cell_bottom']
    cell_ids = cell_top_ids + cell_bottom_ids
    cell_atoms = [atom for atom in system.atoms if atom.idx in cell_ids]
    wall_atoms = [atom for atom in system.atoms if atom.type == 2]

    bonds_to_add:List[Bond] = []

    bondidx = 1
    for cell_atom in cell_atoms:
        min_distance = float('inf')
        closest_wall_atom = None

        for wall_atom in wall_atoms:
            distance = ((cell_atom.x - wall_atom.x) ** 2 + 
                        (cell_atom.y - wall_atom.y) ** 2 + 
                        (cell_atom.z - wall_atom.z) ** 2) ** 0.5

            if distance < min_distance:
                min_distance = distance
                closest_wall_atom = wall_atom

        bond = Bond()
        bond.idx =  bondidx
        bond.type = 1
        bond.atom1 = cell_atom.idx
        bond.atom2 = closest_wall_atom.idx
        bond.length = min_distance
        bonds_to_add.append(bond)
        bondidx += 1
    
    system.update_bonds(bonds_to_add)

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
    wallup_xlo = params['box']['xlo']
    wallup_xhi = params['box']['xhi']
    wallup_ylo = params['box']['ylo']
    wallup_yhi = params['box']['yhi']
    wallup_zlo = params['box']['zhi']-params['wall']['dist_from_box_edge']
    wallup_zhi = params['box']['zhi']-params['wall']['dist_from_box_edge']+0.1

    walldown_xlo = params['box']['xlo']
    walldown_xhi = params['box']['xhi']
    walldown_ylo = params['box']['ylo']
    walldown_yhi = params['box']['yhi']
    walldown_zlo = params['box']['zlo']+params['wall']['dist_from_box_edge']-0.1
    walldown_zhi = params['box']['zlo']+params['wall']['dist_from_box_edge']

    params['wall'].update(
        wallup_xlo = wallup_xlo,
        wallup_xhi = wallup_xhi,
        wallup_ylo = wallup_ylo,
        wallup_yhi = wallup_yhi,
        wallup_zlo = wallup_zlo,
        wallup_zhi = wallup_zhi,

        walldown_xlo = walldown_xlo,        
        walldown_xhi = walldown_xhi,        
        walldown_ylo = walldown_ylo,        
        walldown_yhi = walldown_yhi,        
        walldown_zlo = walldown_zlo,       
        walldown_zhi = walldown_zhi
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
    params['particle']['target_average_direction'] = np.array(params['particle']['target_average_direction'])

if __name__ == "__main__":
    # 初始化
    io = MyIO()
    params = init_io(io)

    # 构建系统模型
    system = System()
    system.update_box(create_box(params['box']))
    system.update_molecules(create_mols(create_mol(io), params['cell']))
    system.update_atoms(creat_walls(params['wall']))
    if params['particle']['type'] == 'single':
        system.update_atoms(creat_particle_single(params['particle']))
    elif params['particle']['type'] == 'random':
        system.update_atoms(creat_particle_random(params['particle']))
    else:
        system.update_atoms(creat_particles(params['particle']))
    create_bonds_between_cell_and_wall(params['fixed_particle_space'], system)

    # 写入文件，并把system对象输出到文件中
    io.write_LAMMPS(system)