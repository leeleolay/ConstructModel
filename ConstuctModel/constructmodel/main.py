import sys
import scipy.spatial as ss
import numpy as np

from .system import System
from .myio import MyIO
from .data import *

def creat_mol(io:MyIO):
    system = System()
    io.read_file_LAMMPS('capsule.data', system)
    cell = Model(system, molidx = 1)
    return cell

def init():
    # 创建一个 IO 对象
    global parameters
    args = MyIO.parse_args()
    parameters = MyIO.load_json(args.json_file)
    return parameters

def init_system(system:System, parameters):
    json = parameters.json_file
    # 设置系统atom，bond，angle，dihedral类型
    system.atomtypes = 3
    system.bondtypes = 2
    system.angletypes = 2
    system.dihedraltypes = 2
    # 设置系统质量类型和参数
    system.masstypes = [1.0 for _ in range(system.atomtypes)]

def main():
    # 初始化
    io = MyIO('capsule.data', 'init.data')
    system = System()
    params = init()

    # 构造cell模型
    cell = creat_mol(io)
    # 构造wall
    # 构造dpd粒子

    # 构建系统
    init_system(system, params)
    system.update_system(cell)

    # 写入文件，并把system对象输出到文件中
    io.output_LAMMPS('init.data',system)

if __name__ == "__main__":
    main()