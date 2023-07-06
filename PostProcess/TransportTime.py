import argparse
import numpy as np
from tqdm import tqdm

def is_inside(position, region):
    """Check if a position is inside a region."""
    min_point, max_point = region
    return np.all(min_point <= position) and np.all(position <= max_point)

def calculate_time(atom_positions, particle_ids, start_region, end_region):
    particle_times = {}
    for particle_id in particle_ids:  # Loop over particle IDs
        particle_positions = atom_positions.get(particle_id, {})
        start_times = []
        end_times = []

        for timestep, position in particle_positions.items():
            if is_inside(position, start_region):
                start_times.append(timestep)
            elif is_inside(position, end_region):
                end_times.append(timestep)

        if start_times and end_times:
            particle_times[particle_id] = min(end_times) - min(start_times)
        else:
            particle_times[particle_id] = None
    return particle_times

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

    return atom_positions

def parse_arguments():
    parser = argparse.ArgumentParser(description='Process trajectory data.')
    parser.add_argument('--file_name', type=str, required=True, help='The trajectory file.')
    parser.add_argument('--atom_type', type=int, required=True, help='The type of atom to analyze.')
    parser.add_argument('--min_timestep', type=int, default=0, help='The minimum timestep to consider.')
    parser.add_argument('--max_timestep', type=int, default=1000, help='The maximum timestep to consider.')
    parser.add_argument('--particle_ids', nargs='+', type=int, required=True, help='The IDs of the particles to analyze.')
    parser.add_argument('--start_region', nargs='+', type=float, required=True, help='The start region for the particles.')
    parser.add_argument('--end_region', nargs='+', type=float, required=True, help='The end region for the particles.')

    args = parser.parse_args()
    args.start_region = np.array(args.start_region).reshape(2, -1)
    args.end_region = np.array(args.end_region).reshape(2, -1)

    return args

def main():
    args = parse_arguments()

    # Read atom positions from trajectory file
    atom_positions = read_lammps_traj(args.file_name, args.atom_type, args.min_timestep, args.max_timestep)

    # Calculate particle times with hardcoded parameters
    particle_times = calculate_time(atom_positions, args.particle_ids, args.start_region, args.end_region)

    # Process the calculated times
    times = []
    with open("transport_results.txt", "w") as f:
        for particle_id, time in particle_times.items():
            if time is None:
                print(f"Particle {particle_id} did not move from the start region to the end region.")
                f.write(f"Particle {particle_id} did not move from the start region to the end region.\n")
            else:
                print(f"Particle {particle_id} moved from the start region to the end region in {time} timesteps.")
                f.write(f"Particle {particle_id} moved from the start region to the end region in {time} timesteps.\n")
                times.append(time)

        if times:
            max_time = np.max(times)
            min_time = np.min(times)
            avg_time = np.mean(times)

            f.write(f"\nMax time: {max_time} timesteps.\n")
            f.write(f"Min time: {min_time} timesteps.\n")
            f.write(f"Avg time: {avg_time} timesteps.\n")

if __name__ == "__main__":
    main()

