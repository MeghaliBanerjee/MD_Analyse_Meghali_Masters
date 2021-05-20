import mdtraj as md
import numpy as np
import pandas as pd
from sklearn.decomposition import IncrementalPCA
from sklearn.preprocessing import StandardScaler
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns;

sns.set()


def draw_vector(v0, v1, ax=None):
    ax = ax or plt.gca()
    arrowprops = dict(arrowstyle='->',
                      linewidth=2,
                      shrinkA=0, shrinkB=0)
    ax.annotate('', v1, v0, arrowprops=arrowprops)


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

    # top to dataframe
    flattened = [list(np.array(row).flat) for row in traj_short.xyz]  # array of lists of coordinate data
    traj_df = pd.DataFrame(flattened)  # to dataframe
    print(f"The data frame is in shape {np.shape(traj_df)}")

    # Incremental PCA on data
    transformer = IncrementalPCA()
    transformer.partial_fit(traj_df)
    traj_transformed = transformer.fit_transform(traj_df)  # transform orig data
    orig_var_ratio = transformer.explained_variance_ratio_
    print("Incremental PCA process complete.")

    """
    # scale data?
    scaler = StandardScaler().fit(traj_df)
    traj_df_scaled = scaler.transform(traj_df)
    print("Original Data Scaled.")


    # PCA on scaled data
    transformer_scaled = IncrementalPCA()
    transformer_scaled.partial_fit(traj_df)
    traj_transformed_scaled = transformer.fit_transform(traj_df)  # transform orig data
    orig_var_ratio_scaled = transformer.explained_variance_ratio_
    print("Incremental PCA process on scaled data complete.")
    """

    """
    # plot cumulative variance explained
    plt.hist(orig_var_ratio, bins=range(len(orig_var_ratio)), label="Original")
    # plt.plot(np.cumsum(orig_var_ratio_scaled), label="Scaled")
    plt.title("Cumulative Variance by Principal Components")
    plt.legend()
    plt.xlabel('number of components')
    plt.ylabel('cumulative explained variance');
    plt.show()
    """

    plt.figure()
    traj_df = traj_df.to_numpy()
    X_new = transformer.inverse_transform(traj_transformed)
    plt.scatter(traj_df[:, 0], traj_df[:, 1], label="orig")
    #plt.scatter(X_new[:, 0], X_new[:, 1], label="PCA")
    plt.legend()
    plt.show()

"""
    plt.scatter(traj_df, alpha=0.2)
    for length, vector in zip(transformer.explained_variance_, transformer.components_):
        v = vector * 3 * np.sqrt(length)
        draw_vector(transformer.mean_, transformer.mean_ + v)
    plt.axis('equal'"""
