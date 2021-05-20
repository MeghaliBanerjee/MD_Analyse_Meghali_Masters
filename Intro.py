# This little script is just meant to be an introduction to the type of scripting we do a lot
# on projects like these.
# If it seems unfamiliar, the first step is to just try and work out line by line what is happening!
# It's made to be slightly confusing in parts, but never unnecessarily so, just to introduce you to this
# way of writing code and some useful concepts!
# The focus here is on producing a 'data pipeline' that takes raw unprocessed simulation data
# and converts it into a more user friendly and useful format.


# Step 1 is to get a working python installation that can run this file!
import mdtraj as md
import time
import os
import pandas as pd


def find_files_by_extension(path, extension, verbose=True):
    list_of_paths = []
    for root, folder, files in os.walk(path):
        for file in files:
            if file.rsplit(".", 1)[-1] == extension:
                list_of_paths.append(os.path.join(root, file))
                if verbose:
                    print(f".{extension} file number {len(list_of_paths)} found at: {list_of_paths[-1]}")
    if verbose:
        print(f"{len(list_of_paths)} .{extension} files found")
    return list_of_paths


def find_number_of_atoms(trajectory, topology):
    t1 = time.time()
    traj = md.load(trajectory, top=topology)
    t2 = time.time()
    n_atoms = traj.n_atoms
    print(f"Traj: {trajectory} \n Topo: {topology} \n Number of atoms: {n_atoms}")
    print(f"Trajectory took {t2 - t1} to load\n")
    return n_atoms


def find_residue_sequence(traj):
    fasta_sequence = ""
    for residue in traj.top.residues:
        if residue.code is not None:
            fasta_sequence += residue.code
        else:
            fasta_sequence += "_"
    return fasta_sequence


def collect_traj_data(coords, top):
    traj = md.load(coords, top=top)
    # print(f"Using sim {name}")
    fasta_seq = find_residue_sequence(traj)  # find fasta
    sim_data = [traj.n_atoms, traj.n_frames, traj.n_chains, traj.n_residues, fasta_seq]

    # add new simulation to dictionary
    name = coords[2:].rsplit(".", 1, top)[0]
    sim_dictionary[name] = sim_data


if __name__ == "__main__":

    # First Create a list of topology files
    topologies = find_files_by_extension(".", "pdb")  # NAMD simulations
    # topologies += find_files_by_extension(".", "parm7")  # Amber simulations
    # topologies += find_files_by_extension(".", "gro")  # GROMACS simulations

    # Next lets get a list of coordinate data files
    #trajectories = find_files_by_extension(".", "dcd")
    # trajectories += find_files_by_extension(".", "xtc")

    # Now let's sort the lists so the directories match
    topologies = sorted(topologies)
    trajectories = sorted(trajectories[2:4])

    sim_dictionary = {}

    # Finds number of atoms, frames, chains, res, fasta
    # for coords, top in zip(trajectories, topologies):
    #    collect_traj_data(coords, top)

    # isolate amino acids in pdbs and dcds files
    for coords, top in zip(trajectories, topologies):
        name = coords[2:].rsplit(".", 1)[0]
        print(f"Using sim {name}")
        traj = md.load(coords, top=top)
        residue_indices = traj.topology.select("protein")
        isolate_traj = traj.atom_slice(residue_indices)
        print(f"Saving files for {name}")
        md.Trajectory.save_pdb(isolate_traj, f"{name}_isolated.pdb")
        md.Trajectory.save_dcd(isolate_traj, f"{name}_isolated.dcd")

    # write to dataframe
    # df = pd.DataFrame.from_dict(sim_dictionary, orient='index', columns=['Number of atoms', 'Number of frames', 'Number of chains', 'Number of residues', 'Fasta Sequence'])
    # write to csv
    # df.to_csv('Simulation_Data.csv', encoding='utf-8')
