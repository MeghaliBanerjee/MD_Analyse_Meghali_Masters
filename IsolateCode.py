# script to isolate amino acid proteins and save pdb and dcd files

import mdtraj as md
import numpy as np
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
    trajectories = sorted(trajectories)

    for coords, top in zip(trajectories, topologies):
        name = coords[2:].rsplit(".", 1)[0]
        print(f"Using sim {name}")

        # load in as trajectory
        traj = md.load(coords, top=top)
        protein_index = traj.top.select("protein")  # find indices for residue

        # top to dataframe
        table, bonds = traj.top.to_dataframe()

        #table = table_all.loc[protein_index]

        # rename some proteins
        mask = table['resName'].str.contains("CY")
        table.loc[mask, 'resName'] = 'CYS'

        mask2 = table['resName'].str.contains("HD|HE")
        table.loc[mask2, 'resName'] = 'HIS'

        #shorten table to include only residues SER to ARG
        start_index = table[table['resName'] == 'SER'].index[0]  # find start index
        end_index = table[table['resName'] == 'ARG'].index[-1]  # find end index
        short_table = table.loc[start_index:end_index]  # splice data frame

        #reset index
        short_table.reset_index(drop=True, inplace=True)

        # topology back from frame
        top2 = md.Topology.from_dataframe(short_table)

        # shortened table:
        top2.save(f"{name}_isolatedSpliced.pdb")
