import numpy as np
import matplotlib.pyplot as plt
import multiprocessing
from tqdm import tqdm
import argparse
import os
from datetime import datetime

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
                if atom_type == args.target_atom_type:
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

def calculate_msd_single(atom_pos):
    atom_id, positions = atom_pos
    msd_values = []

    for time1, pos1 in positions.items():
        for time2, pos2 in positions.items():
            if time1 < time2:
                msd_values.append((((pos1 - pos2) ** 2).sum(), time2 - time1))

    return msd_values

def calculate_msd_serial(atom_positions):
    msd_values = []
    ensemble_msd_values = {}

    msd_results = [calculate_msd_single(atom_pos) for atom_pos in atom_positions.items()]

    for msd_result in msd_results:
        for msd_value, time_delta in msd_result:
            msd_values.append((time_delta, msd_value))
            if time_delta not in ensemble_msd_values:
                ensemble_msd_values[time_delta] = {'sum': 0, 'count': 0, 'squared_sum': 0}
            ensemble_msd_values[time_delta]['sum'] += msd_value
            ensemble_msd_values[time_delta]['squared_sum'] += msd_value ** 2
            ensemble_msd_values[time_delta]['count'] += 1

    ensemble_msd_values = {k: {'mean': ensemble_msd_values[k]['sum'] / ensemble_msd_values[k]['count'],
                           'std': np.sqrt(ensemble_msd_values[k]['squared_sum'] / ensemble_msd_values[k]['count'] -
                                          (ensemble_msd_values[k]['sum'] / ensemble_msd_values[k]['count']) ** 2)}
                      for k in ensemble_msd_values}

    return msd_values, ensemble_msd_values

def calculate_msd(atom_positions):
    with multiprocessing.Pool(processes=args.num_processes) as pool:
        msd_results = list(tqdm(pool.imap_unordered(calculate_msd_single, atom_positions.items()), total=len(atom_positions), desc="Calculating MSD"))

    msd_values = []
    ensemble_msd_values = {}
    for msd_result in msd_results:
        for msd_value, time_delta in msd_result:
            msd_values.append((time_delta, msd_value))
            if time_delta not in ensemble_msd_values:
                ensemble_msd_values[time_delta] = {'sum': 0, 'count': 0, 'squared_sum': 0}
            ensemble_msd_values[time_delta]['sum'] += msd_value
            ensemble_msd_values[time_delta]['squared_sum'] += msd_value ** 2
            ensemble_msd_values[time_delta]['count'] += 1


    ensemble_msd_values = {k: {'mean': ensemble_msd_values[k]['sum'] / ensemble_msd_values[k]['count'],
                           'std': np.sqrt(ensemble_msd_values[k]['squared_sum'] / ensemble_msd_values[k]['count'] -
                                          (ensemble_msd_values[k]['sum'] / ensemble_msd_values[k]['count']) ** 2)}
                      for k in ensemble_msd_values}

    return msd_values, ensemble_msd_values

def calculate_diffusion_coefficient(msd_values, time_range=None):
    time_deltas, msd_averages = zip(*msd_values)

    time_deltas = np.array(time_deltas)
    msd_averages = np.array(msd_averages)

    if time_range is not None:
        mask = (time_deltas >= time_range[0]) & (time_deltas <= time_range[1])
        time_deltas = time_deltas[mask]
        msd_averages = msd_averages[mask]

    slope, _ = np.polyfit(time_deltas, msd_averages, 1)
    diffusion_coefficient = slope / 6

    return diffusion_coefficient

def plot_msd(msd_values, ensemble_msd_values, results_dir ,name):
    # Plot the original MSD values
    plt.figure(figsize=(8, 6))
    time_deltas = []
    msd_averages = []

    for time_delta, msd_value in msd_values:
        time_deltas.append(time_delta)
        msd_averages.append(msd_value)

    plt.scatter(time_deltas, msd_averages, label='MSD', alpha=0.5)
    plt.xlabel('Time Delta')
    plt.ylabel('MSD')
    plt.title('Original MSD')
    plt.legend()
    plt.savefig(os.path.join(results_dir, str('msd_plot_'+ name + '.png')))
    plt.show()

    # Plot the ensemble-averaged MSD values
    plt.figure(figsize=(8, 6))
    ensemble_time_deltas = np.array(list(ensemble_msd_values.keys()))
    ensemble_msd_averages = np.array([ensemble_msd_values[k]['mean'] for k in ensemble_msd_values])
    ensemble_msd_std = np.array([ensemble_msd_values[k]['std'] for k in ensemble_msd_values])
    plt.errorbar(ensemble_time_deltas, ensemble_msd_averages, yerr=ensemble_msd_std, label='Ensemble Average MSD', color='red', fmt='o', capsize=5)

    plt.xlabel('Time Delta')
    plt.ylabel('MSD')
    plt.title('Ensemble Averaged MSD')
    plt.legend()
    plt.savefig(os.path.join(results_dir, str('ensemble_averaged_msd_plot_' + name + '.png')))
    plt.show()

def plot_diffusion_coefficient(diffusion_coefficient, results_dir, name):
    plt.bar(["D"], [diffusion_coefficient])
    plt.ylabel('Diffusion Coefficient')
    plt.savefig(os.path.join(results_dir, str('diffusion_coefficient_plot_' + name + '.png')))
    plt.show()

def get_unique_filename(filename, ext, directory='.'):
    base_filename, _ = os.path.splitext(filename)
    index = 0
    candidate = f"{base_filename}{ext}"

    while os.path.exists(os.path.join(directory, candidate)):
        index += 1
        candidate = f"{base_filename}_{index}{ext}"

    return candidate

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calculate MSD and diffusion coefficient from LAMMPS trajectory file")
    parser.add_argument("traj_file", help="Path to the LAMMPS trajectory file")
    parser.add_argument("target_atom_type", type=int, help="Target atom type for MSD calculation")
    parser.add_argument("min_timestep", type=int, help="Minimum timestep to consider")
    parser.add_argument("max_timestep", type=int, help="Maximum timestep to consider")
    parser.add_argument("output", type=str, help="Output file location")
    parser.add_argument("--num_processes", type=int, default=1, help="Number of processes to use (optional)")

    args = parser.parse_args()
    traj_name = os.path.basename(args.traj_file)

    results_dir = f'{args.output}results_msd'
    os.makedirs(results_dir, exist_ok=True)

    atom_positions = read_lammps_traj(args.traj_file, args.target_atom_type, args.min_timestep, args.max_timestep)

    if args.num_processes > 1:
        msd_values, ensemble_msd_values = calculate_msd(atom_positions)
    else:
        msd_values, ensemble_msd_values = calculate_msd_serial(atom_positions)

    diffusion_coefficient = calculate_diffusion_coefficient(msd_values)

    print(f'Diffusion Coefficient: {diffusion_coefficient}')

    plot_diffusion_coefficient(diffusion_coefficient, results_dir, traj_name)
    plot_msd(msd_values, ensemble_msd_values, results_dir, traj_name)
    
    # Save the data to a txt file
    msd_data_filename = get_unique_filename('msd_data', '.txt', results_dir)
    with open(os.path.join(results_dir, str(msd_data_filename+'_'+traj_name)), 'w') as f:
        f.write("MSD Values:\n")
        for time_delta, msd_value in msd_values:
            f.write(f"{time_delta} {msd_value}\n")

        f.write("\nEnsemble MSD Values:\n")
        for time_delta, ensemble_msd_value in ensemble_msd_values.items():
            f.write(f"{time_delta} {ensemble_msd_value}\n")

        f.write(f"\nDiffusion Coefficient: {diffusion_coefficient}\n")
