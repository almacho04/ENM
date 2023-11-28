from enm import *
from multiprocessing import Pool
import numpy as np
import os
import argparse
import time
import csv
from multiprocessing import cpu_count

output_directory = os.path.expanduser("~/Desktop/soldier_girl/everything/ENM_models/")

csv_path = os.path.expanduser("~/Desktop/soldier_girl/everything/protein_names.csv")
file_names_to_abstracted_names = {}
copy_filenames = []
with open(csv_path, "r") as csv_file:
    reader = csv.reader(csv_file)
    next(reader)
    for row in reader:
        file_names_to_abstracted_names[row[0]] = row[1]
        copy_filenames.append(row[0])

k_params = [{'kconst': k} for k in ['homo', 'exp', 'dENM', 'sENM10', 'sdENM', 'sENM13']]

edge_params = [{'adj': 'cutOff', 'cutOff': c} for c in range(8, 22, 2)] + [{'adj': 'delaunay'}]
n_cores = cpu_count()
abstract_names_to_file_names = {v: k for k, v in file_names_to_abstracted_names.items()}

def check_if_file_exists(file_path):
    print(format_colored(file_path))
    print(os.path.exists(file_path))
    # Dont calculate if already done.
    return os.path.exists(file_path)

def load_pdb(file_path):
    enm = ENM(file_path)
    enm.filter("CA")
    return enm


def get_all_pdb_files(dir_path):
    pdb_files = []
    for file in os.listdir(dir_path):
        if file.endswith(".pdb"):
            pdb_files.append(os.path.join(dir_path, file))
    return pdb_files


def process_enm_model(pdb_file, abstract_name):
    model = load_pdb(pdb_file)

    try:
        for ep in edge_params:
            for kp in k_params:
                ep_part = None
                if ep['adj'] == 'cutOff':
                    ep_part = str(ep['cutOff'])
                elif ep['adj'] == 'delaunay':
                    ep_part = 'delaunay'

                covariance_file_path = os.path.expanduser(f"~/Desktop/soldier_girl/everything/ENM_models/enm_covariances/{abstract_name}/" + f"{abstract_name}_{ep_part}_{kp['kconst']}_Cov.npy")
                if not check_if_file_exists(covariance_file_path):
                    H = model.getHessian(**ep, **kp)

                    print(f"{abstract_name}_{ep_part}_{kp['kconst']}")
                    np.save(os.path.join(output_directory + f"/enm_hessians/{abstract_name}", f"{abstract_name}_{ep_part}_{kp['kconst']}_H.npy"), H)
                    D, eigvals, eigvecs = get_enm_displacements_from_H(H)
                    print(format_colored(f"Finished calculating displacements for {abstract_name}_{ep_part}_{kp['kconst']}", color_code="32"))
                    np.save(os.path.join(output_directory + f"/enm_displacements/{abstract_name}", f"{abstract_name}_{ep_part}_{kp['kconst']}_Disp.npy"), D)
                    print(format_colored(f"Finished saving displacements for {abstract_name}_{ep_part}_{kp['kconst']}", color_code="32"))
                    np.save(os.path.join(output_directory + f"/enm_eigvalues/{abstract_name}", f"{abstract_name}_{ep_part}_{kp['kconst']}_Eigvals.npy"), eigvals)
                    np.save(os.path.join(output_directory + f"/enm_eigvectors/{abstract_name}", f"{abstract_name}_{ep_part}_{kp['kconst']}_Eigvecs.npy"), eigvecs)
                    print(format_colored(f"Finished saving eigenvalues and eigenvectors for {abstract_name}_{ep_part}_{kp['kconst']}", color_code="32"))
                    Cov = get_cov(H, eigvals, eigvecs)
                    np.save(os.path.join(output_directory + f"/enm_covariances/{abstract_name}", f"{abstract_name}_{ep_part}_{kp['kconst']}_Cov.npy"), Cov)
                    print("\n")
                    print(format_colored(f"Done with: {abstract_name}_{ep_part}_{kp['kconst']} \n {abstract_names_to_file_names[abstract_name]}", color_code="32"))
                else:
                    print(format_colored(f"Already done with: {abstract_name}_{ep_part}_{kp['kconst']} \n {abstract_names_to_file_names[abstract_name]}", color_code="32"))
    except Exception as e:
        print(f"Error in {abstract_name}: {str(e)}")
    print(format_colored("FINISHED PROCESSING: " + abstract_name))


def process_enm_model_wrapper(args):
    pdb_file, abstract_name = args
    process_enm_model(pdb_file, abstract_name=abstract_name)


if __name__ == "__main__":
    start = time.time()
    parser = argparse.ArgumentParser(description="process ENM data")
    parser.add_argument("--dir", type=str, help="directory containing pdb")
    args = parser.parse_args()
    pdb_files = get_all_pdb_files(args.dir)
    N_CORES = min(n_cores, len(pdb_files))

    # Create a pool of worker processes
    with Pool(processes=N_CORES) as pool:
        pool.map(process_enm_model_wrapper, [(pdb_file, file_names_to_abstracted_names[extract_filename(pdb_file)]) for pdb_file in pdb_files])

    end = time.time()
    stringTime = "Time taken: " + str(end - start)
    print(format_colored(stringTime, color_code="36"))
