import mdtraj as md
import os
import glob


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

    # Collect new topologies and trajectories

    pathlist = glob.glob("C:/Users/megha/PycharmProjects/pythonProject/IsolatedSpliced/*")
    topologies = []

    for path in pathlist:
        topologies.append(path.split("/")[-1])

    trajectories = find_files_by_extension(".", "xtc")
    trajectories += find_files_by_extension(".", "dcd")

    # verify all topologies work with first trajectory
    traj = trajectories[0]
    print(f"Test all topologies with trajectory: {traj}.")

    for top in topologies:
        print(f"Testing topology {top}.")
        loaded = md.load(traj, top=top)  # load in trajectory
        print("Success.")

    # verify all trajectories run with first topology
    topology = topologies[0]

    for trajectory in trajectories:
        print(f"Testing trajectory {trajectory}.")
        loaded = md.load(trajectory, top=topology)  # load in trajectory
        print("Success.")
