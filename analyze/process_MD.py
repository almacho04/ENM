
from enm import *
from MDAnalysis import *
from MDAnalysis.analysis.encore.covariance import covariance_matrix
from MDAnalysis.analysis import align
from MDAnalysis import transformations
from multiprocessing import cpu_count
from scipy.stats import pearsonr

from multiprocessing import Pool
import numpy as np
import matplotlib.pyplot as plt
import os
import argparse
import time

n_cores = cpu_count()



def check_pdb_dcd_files(pdb_files, dcd_files, current_dir):
    if len(pdb_files) == 1:
        print(format_colored("Found pdb file: {}".format(pdb_files[0])))
    elif len(pdb_files) > 1:
        print(format_colored("ERROR: Found multiple pdb files: {}".format(pdb_files)))
        # throw error
        ValueError("Found multiple pdb files")
        return
    elif len(pdb_files) == 0:
        print(format_colored("ERROR: No pdb file found in directory {}".format(current_dir)))
        # throw error and stop execution
        ValueError("No pdb file found")
        return
    if len(dcd_files) > 0:
        print(format_colored("Found dcd files: {}".format(dcd_files)))
    else:
        print(format_colored("ERROR: No dcd files found in directory {}".format(current_dir)))
        # throw error and stop execution
        ValueError("No dcd files found")
        # quit the program without exit()
        return

def check_processed_files(output_directory, protein_name):
    result = {}

    # Define the expected file paths
    files_to_check = [
        os.path.join(output_directory, "trajectories", f"{protein_name}_coords.npy"),
        os.path.join(output_directory, "covariances", f"{protein_name}_cov.npy"),
        os.path.join(output_directory, "fluctuations", f"{protein_name}_flucts.npy"),
        os.path.join(output_directory, "universes", f"{protein_name}.pkl"),
    ]

    for file_path in files_to_check:
        if not os.path.exists(file_path):
            result[file_path] = "missing"

    return result

def load_md_data(pdb_file, dcd_files, current_dir):
    u = Universe(pdb_file, dcd_files, in_memory=True)
    print(format_colored("Made universe with trajectories:"), dcd_files)
    transform = transformations.NoJump(max_threads=n_cores, parallelizable=True, check_continuity=True)
    u.trajectory.add_transformations(transform)
    print(format_colored("Added transformations to trajectory."))
    aligner = align.AlignTraj(u, reference=u, select="protein", in_memory=True, max_cores=n_cores)
    aligner.run()
    print(format_colored("Aligned trajectory."))
    return u

def main(args):
    current_dir = args.dir
    output_directory = os.path.expanduser("~/Desktop/soldier_girl/everything/processed_MD/")
    
    print(current_dir)
    dcd_files = [os.path.join(current_dir, f) for f in os.listdir(current_dir) if os.path.isfile(os.path.join(current_dir, f)) and f.endswith(".dcd")]
    # sort dcd files by their name
    dcd_files.sort()
    pdb_files = [os.path.join(current_dir, f) for f in os.listdir(current_dir) if os.path.isfile(os.path.join(current_dir, f)) and f.endswith(".pdb")]
    mae_files = [os.path.join(current_dir, f) for f in os.listdir(current_dir) if os.path.isfile(os.path.join(current_dir, f)) and f.endswith(".mae")]
    psf_files = [os.path.join(current_dir, f) for f in os.listdir(current_dir) if os.path.isfile(os.path.join(current_dir, f)) and f.endswith(".psf")]
    print(pdb_files)
    check_pdb_dcd_files(pdb_files, dcd_files, current_dir)
    pdb_file = None
    if(len(pdb_files) == 0):
        if(len(psf_files) == 0):
            pdb_file = mae_files[0]
        elif(len(psf_files) == 1):
            pdb_file = psf_files[0]
    else:
        pdb_file = pdb_files[0]
    protein_name = os.path.splitext(os.path.basename(pdb_file))[0]         
    # check if the protein is already processed
    missing_files = check_processed_files(output_directory, protein_name)

    if missing_files:
        print("ERROR: Protein not fully processed. The following files are missing:")
        for file_path, status in missing_files.items():
            print("\n")
            print(format_colored(f"{file_path} - {status}", color_code='36'))
    else:
        print("\n\n\n\n")
        print(format_colored("Protein already processed.", color_code="36"))
        return
    
    u = load_md_data(pdb_file, dcd_files, current_dir)
    save_universe_as_pickle(u, protein_name+".pkl", output_directory + "universes/")
    C = get_coords(u) # returns a numpy array of shape (n_frames, n_atoms, 3)
    

    # ! Saving the coordinates, covariance matrix, and fluctuations
    # * -----------------------------------------------------------
    np.save(os.path.join(output_directory+"trajectories/", protein_name+"_coords"), C)
    print(format_colored("Saved coordinates as"), protein_name + ".npy in directory", output_directory+"trajectories/")
    mdcov = covariance_matrix(u, select="name CA")
    np.save(os.path.join(output_directory+"covariances/", protein_name+"_cov"), mdcov)
    print(format_colored("Saved covariance matrix as"), protein_name + ".npy in directory", output_directory+"covariances/")
    flucts = get_rmsf(C)
    np.save(os.path.join(output_directory+"fluctuations/", protein_name+"_flucts"), flucts)
    print(format_colored("Saved fluctuations as"), protein_name + ".npy in directory", output_directory+"fluctuations/")
    # * -----------------------------------------------------------


    print("\n\n\n\n")
    print(format_colored("Done!", color_code="36"))
    print(format_colored(u.trajectory, color_code="36"))
    return

def process_directory(directory):
    args = argparse.Namespace(dir=directory)
    main(args)

def process_directories(args):
    start_time = time.time()
    # Get a list of directories to process (e.g., from command line arguments)
    directories = args.dirs  # 'dirs' should be a list of directory paths

    # Create a multiprocessing pool to run 'main' for each directory
    N_CORES = min(n_cores, len(directories))
    print(f"Running on {N_CORES} cores")
    with Pool(processes=N_CORES) as pool:
        pool.map(process_directory, directories)
    print(format_colored("Done with " + str(len(directories)) + " directories!", color_code="36"))
    print("\n", np.array(directories))
    elapsed_time = str(time.time() - start_time) + " seconds"
    print(format_colored("Elapsed time: " + elapsed_time, color_code="36"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process MD data for multiple directories")
    parser.add_argument("--dirs", nargs="+", type=str, help="Directories where the pdb and dcd files are located")
    args = parser.parse_args()

    process_directories(args)