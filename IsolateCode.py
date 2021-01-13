# script to isolate amino acid proteins and save pdb and dcd files

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


if __name__ == "__main__":

    # load in files
    topologies = find_files_by_extension(".", "pdb")  # NAMD simulations
    topologies += find_files_by_extension(".", "parm7")  # Amber simulations
    trajectories = find_files_by_extension(".", "dcd")
    topologies = sorted(topologies)
    trajectories = sorted(trajectories[2:4])

    # isolate amino acids in pdbs and dcds files
    for coords, top in zip(trajectories, topologies):
        name = coords[2:].rsplit(".", 1)[0]
        print(f"Using sim {name}")
        traj = md.load(coords, top=top)  # load in as trajectory

        residue_indices = traj.topology.select("protein")  # find indices for residue
        isolate_traj = traj.atom_slice(residue_indices)  # isolate these indices in trajectory

        print(f"Saving files for {name}")
        isolate_traj[0].save(f"{name}_isolated.pdb")
        isolate_traj.save(f"{name}_isolated.dcd")
