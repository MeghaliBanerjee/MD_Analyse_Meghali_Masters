# Calculate number of contacts

import mdtraj as md
import numpy as np
import time
import os
import re
from contact_map import ContactFrequency, ContactDifference

import matplotlib.pyplot as plt


def fix_topology(traj):
    for i, res in enumerate(traj.top._residues):
        res.resSeq = i
        if res.name[:2] == "CY":
            res.name = "CYS"
        elif res.name[:2] == "HD" or res.name[:2] == "HE":
            res.name = "HIS"
    return traj


def rename_residues(traj):
    print("Renaming residues...")
    for i, res in enumerate(traj.top._residues):
        res.resSeq = i
        if res.name[:2] == "CY":
            res.name = "CYS"
        elif res.name[:2] == "HD" or res.name[:2] == "HE":
            res.name = "HIS"
    print("Process Complete")
    return traj


def reduce_system(traj):
    print("Reducing size of system to include atoms in residue SER to ARG...")
    protein_index = traj.top.select("protein")  # find indices for atoms part of residue
    isolated_traj = md.Trajectory.atom_slice(traj, protein_index)  # trajectory isolated

    # shorten topology to include only residues SER to ARG
    ser_index = [atom.index for atom in isolated_traj.top.atoms if (atom.residue.name == "SER")]
    arg_index = [atom.index for atom in isolated_traj.top.atoms if (atom.residue.name == "ARG")]
    print("Process Complete")
    return md.Trajectory.atom_slice(isolated_traj,
                                    range(ser_index[0], arg_index[-1]))


if __name__ == "__main__":
    # Load trajectory
    traj = "140mM\helicase_140mM.dcd"
    top = "140mM\helicase_140mM.parm7"

    print(f"Loading trajectory for top {top} with trajectory {traj}")
    loaded = md.load(traj, top=top)

    # every 1000th frame for shorter trajectory
    test_traj = loaded[::1000]

    # rename proteins, reduce system
    traj_renamed = rename_residues(test_traj)  # rename proteins
    traj_short = reduce_system(traj_renamed)

    #frame_contacts = ContactFrequency(traj_short[0])

    large_cutoff = ContactFrequency(trajectory=traj_short[0], cutoff=1.5)
    print(large_cutoff)

    """
    (fig, ax) = frame_contacts.residue_contacts.plot(cmap='seismic', vmin=-1, vmax=1)
    plt.title(f"Contact map 140mM")
    plt.xlabel("Residue")
    plt.ylabel("Residue")
    plt.show()
    """

