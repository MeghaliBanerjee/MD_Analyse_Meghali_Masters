# script to isolate amino acid proteins and save pdb and dcd files

import mdtraj as md
import numpy as np
import time
import os
import re


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


def fix_topology(traj):
    for i, res in enumerate(traj.top._residues):
        res.resSeq = i
        if res.name[:2] == "CY":
            res.name = "CYS"
        elif res.name[:2] == "HD" or res.name[:2] == "HE":
            res.name = "HIS"
    return traj


if __name__ == "__main__":

    # load in files
    topologies = find_files_by_extension(".", "gro")  # GROMACS simulations
    topologies += find_files_by_extension(".", "pdb")  # NAMD simulations
    topologies += find_files_by_extension(".", "parm7")  # Amber simulations

    trajectories = find_files_by_extension(".", "xtc")
    trajectories += find_files_by_extension(".", "dcd")

    topologies = sorted(topologies)
    trajectories = sorted(trajectories)

    for coords, top in zip(trajectories, topologies):
        name = coords.split(".")[1].split("\\")[2]

        print(f"Using top {top} with traj {coords}")

        traj = md.load(coords, top=top)  # load in trajectory

        # rename residues
        traj_fixed = fix_topology(traj)

        # isolate proteins
        protein_index = traj_fixed.top.select("protein")  # find indices for atoms part of residue
        isolated_traj = md.Trajectory.atom_slice(traj_fixed[0], protein_index)  # trajectory isolated

        # shorten topology to include only residues SER to ARG
        ser_index = [atom.index for atom in isolated_traj.top.atoms if (atom.residue.name == "SER")]
        arg_index = [atom.index for atom in isolated_traj.top.atoms if (atom.residue.name == "ARG")]

        shorter_traj = md.Trajectory.atom_slice(isolated_traj[0],
                                                range(ser_index[0], arg_index[-1]))  # trajectory isolated

        shorter_traj[0].save(f'IsolatedSpliced\\{name}_isolatedSpliced.pdb')
