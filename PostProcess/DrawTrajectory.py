import argparse
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt

def read_lammps_traj(file_name, atom_type, min_timestep, max_timestep):
    atom_positions = {}
    reading_atoms = False
    box_bounds = None
    timestep = None

    with open(file_name, 'r') as file:
        for line in tqdm(file, desc="Reading trajectory file", unit="line"):
            if 'ITEM: TIMESTEP' in line:
                reading_atoms = False
                timestep = int(next(file).strip())
                if timestep > max_timestep:
                    break

            if 'ITEM: BOX BOUNDS' in line:
                box_bounds = []
                for _ in range(3):
                    box_bounds.append(list(map(float, next(file).strip().split())))
                box_bounds = np.array(box_bounds)

            if reading_atoms:
                data = line.strip().split()
                atom_id, atom_type, xs, ys, zs = int(data[0]), int(data[1]), float(data[2]), float(data[3]), float(data[4])
                if atom_type == atom_type:
                    x, y, z = box_bounds[:, 0] + np.array([xs, ys, zs]) * (box_bounds[:, 1] - box_bounds[:, 0])
                    if atom_id not in atom_positions:
                        atom_positions[atom_id] = {}
                    atom_positions[atom_id][timestep] = np.array([x, y, z])

            if 'ITEM: ATOMS' in line:
                if timestep >= min_timestep:
                    reading_atoms = True
                else:
                    reading_atoms = False

    return atom_positions,box_bounds

def plot_particle_trajectory(particle_positions, boxbounds, save_file=None):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for particle_id, positions in particle_positions.items():
        x = [pos[0] for pos in positions.values()]
        y = [pos[1] for pos in positions.values()]
        z = [pos[2] for pos in positions.values()]

        ax.plot(x, y, z, label=f'Particle {particle_id}')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.legend()

    ax.set_xlim(boxbounds[0,0],boxbounds[0,1])
    ax.set_ylim(boxbounds[1,0],boxbounds[1,1])
    ax.set_zlim(boxbounds[2,0],boxbounds[2,1])

    if save_file:
        plt.savefig(save_file)
    else:
        plt.show()

def parse_arguments():
    parser = argparse.ArgumentParser(description='Process trajectory data.')
    parser.add_argument('--file_name', type=str, required=True, help='The trajectory file.')
    parser.add_argument('--min_timestep', type=int, default=0, help='The minimum timestep to consider.')
    parser.add_argument('--max_timestep', type=int, default=1000, help='The maximum timestep to consider.')
    parser.add_argument('--atom_type', type=int, required=True, help='The type of atom to analyze.')
    parser.add_argument('--particle_ids', nargs='+', type=int, required=True, help='The IDs of the particles to analyze.')
    
    args = parser.parse_args()

    return args

def main():
    args = parse_arguments()

    # Read atom positions from trajectory file
    atom_positions, boxbounds = read_lammps_traj(args.file_name, args.atom_type, args.min_timestep, args.max_timestep)

    # Select positions for the desired particles
    particle_positions = {}
    for particle_id in args.particle_ids:
        if particle_id in atom_positions:
            particle_positions[particle_id] = atom_positions[particle_id]

    # Plot particle trajectories
    save_file = "trajectory.png"
    plot_particle_trajectory(particle_positions, boxbounds, save_file)

if __name__ == "__main__":
    main()
