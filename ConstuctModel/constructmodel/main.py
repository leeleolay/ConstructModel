import sys
import scipy.spatial as ss
import numpy as np
from typing import Dict

from .system import System
from .myio import MyIO
from .data import *
from .box import Box

def create_mol():
    system = System()
    MyIO.read_LAMMPS('capsule.data', system)
    cell = Molecule.creat_from_system(system)
    return cell

def create_box(params:Dict):
    box = Box()
    box.xlo = params['xlo']
    box.xhi = params['xhi']
    box.ylo = params['ylo']
    box.yhi = params['yhi']
    box.zlo = params['zlo']
    box.zhi = params['zhi']
    return box

def init_io(io:MyIO):
    global parameters
    args = MyIO.parse_args()
    parameters:Dict = MyIO.load_json(args.jsonfile)
    io.input_file_name = args.input
    io.output_file_name = args.output
    return parameters


def process_params(params):
    # 处理wall参数
    xlo_wall_up = xlo_wall_down = params['box']['xlo']
    xhi_wall_up = xhi_wall_down = params['box']['xhi']
    ylo_wall_up = ylo_wall_down = params['box']['ylo']
    yhi_wall_up = yhi_wall_down = params['box']['yhi']
    zlo_wall_up = params['box']['zlo']-params['wall']['dist_from_box']-params['wall']['rho_wall']
    zhi_wall_up = params['box']['zhi']-params['wall']['dist_from_box']
    zlo_wall_down = params['box']['zlo']+params['wall']['dist_from_box']
    zhi_wall_down = params['box']['zhi']+params['wall']['dist_from_box']+params['wall']['rho_wall']
    params['wall']:Dict.update(
        xlo_wall_up = xlo_wall_up,
        xlo_wall_down = xlo_wall_down,
        xhi_wall_up = xhi_wall_down,
        xhi_wall_down = xhi_wall_down,
        ylo_wall_up = ylo_wall_up,
        ylo_wall_down = ylo_wall_down,
        yhi_wall_up = yhi_wall_up,
        yhi_wall_down = yhi_wall_down,
        zlo_wall_up = zlo_wall_up,
        zlo_wall_down = zlo_wall_down,
        zhi_wall_up = zhi_wall_up,
        zhi_wall_down = zhi_wall_down
    )
    return params

def creat_wall(params:Dict):
    wall:List[Atom] = []
    atom = Atom()
    atom.idx = 1
    atom.mass = 1.0
    atom.molidx = 1001
    atom.type = 1
    for i in range(params['xlo_wall_up'],params['xhi_wall_up'],params['rho_wall']):
        for j in range(params['ylo_wall_up'],params['yhi_wall_up'],params['rho_wall']):
            for k in range(params['zlo_wall_up'],params['zhi_wall_up'],params['rho_wall']):
                atom.x = i
                atom.y = j
                atom.z = k
                wall.append(atom)
                atom.idx += 1
    atom.molidx += 1
    for i in range(params['xlo_wall_down'],params['xhi_wall_down'],params['rho_wall']):
        for j in range(params['ylo_wall_down'],params['yhi_wall_down'],params['rho_wall']):
            for k in range(params['zlo_wall_down'],params['zhi_wall_down'],params['rho_wall']):
                atom.x = i
                atom.y = j
                atom.z = k
                wall.append(atom)
                atom.idx += 1
    return wall

def creat_particle(params:Dict):
    wall:List[Atom] = []
    atom = Atom()
    atom.idx = 1
    atom.mass = 1.0
    atom.molidx = 1001
    atom.type = 1
    for i in range(params['xlo_particle'],params['xhi_particle'],params['rho_particle']):
        for j in range(params['ylo_particle'],params['yhi_particle'],params['rho_particle']):
            for k in range(params['zlo_particle'],params['zhi_particle'],params['rho_particle']):
                atom.x = i
                atom.y = j
                atom.z = k
                wall.append(atom)
                atom.idx += 1

def main():
    # 初始化
    io = MyIO('capsule.data', 'init.data')
    params = init_io(io)
    params = process_params(params)

    # 构造box
    box = create_box(params['box'])
    # 构造cell模型
    cell = create_mol(io)
    # 构造wall
    wall = creat_wall(params['wall'])
    # 构造dpd粒子
    particle = creat_particle(params['particle'])

    # 构建系统
    system = System()
    system.update_box(box)
    system.update_molecule(cell)
    system.update_atom(wall)
    system.update_atom(particle)

    # 写入文件，并把system对象输出到文件中
    io.write_LAMMPS(io.output_file_name,system)

if __name__ == "__main__":
    main()