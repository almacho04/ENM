from enm import *
import numpy as np 
import csv
from scipy.stats import pearsonr
import glob
from multiprocessing import cpu_count
from multiprocessing import Pool
import time
import os



output_directory = os.path.expanduser("~/Desktop/soldier_girl/everything/")

csv_path = os.path.expanduser("~/Desktop/soldier_girl/everything/protein_names.csv")


total_cores = cpu_count()

def check_if_entry_exists(protein_name):
    # check if the csv file exists
    if not os.path.exists(output_csv_path):
        return False
    splitted_name = protein_name.split("_")
    with open(output_csv_path, "r") as csvfile:
        reader = csv.reader(csvfile)
        next(reader)
        for row in reader:
            if row[0] == splitted_name[0] + "_" + splitted_name[1] and row[1] == splitted_name[2] and row[2] == splitted_name[3]:
                return True
    return False
    

def load_enm_data(abstract_name):
    try: 
        enm_covs = {}
        enm_displacements = {}

        for file in os.listdir(output_directory + f"/ENM_models/enm_covariances/{abstract_name}"):
            if file.endswith(".npy"):
                name_of_file = extract_filename(file)
                dir_of_enm_covs[name_of_file] = np.load(output_directory + f"ENM_models/enm_covariances/{abstract_name}/{file}")
                dir_of_enm_displacements[name_of_file[:-3] + "Disp"] = np.load(output_directory + f"ENM_models/enm_displacements/{abstract_name}/{file[:-7] + 'Disp.npy'}")
                print(format_colored(f"Loaded ENM covariance matrices and displacements for {name_of_file}", color_code = "36"))

        

        print(format_colored(f"Loaded ENM covariance matrices and displacements for {abstract_name}"))
        
        return enm_covs, enm_displacements
    except Exception as e:
        print(format_colored(f"Error loading ENM data for {abstract_name}: {str(e)}", color_code="31"))
        raise

def write_to_csv(protein_name, enm_cov, enm_displacement):
    if(protein_name.startswith("protein_019")):
        print(format_colored(f"SKIPPING {protein_name}"))
        print("\n")
        return
    # before writing to csv, check if the entry already exists in the csv file
    if check_if_entry_exists(protein_name):
        print(format_colored(f"Entry {protein_name} already exists in the csv file"))
        return
    splitted_name = protein_name.split("_")
    results = []
    md_cov = dir_of_md_covs[splitted_name[0] + "_" + splitted_name[1]]
    md_fluct = dir_of_md_rmsfs[splitted_name[0] + "_" + splitted_name[1]]
    std_enm_cov = standardize(enm_cov)
    std_md_cov = standardize(md_cov)
    idx = np.triu_indices(std_md_cov.shape[0], k=1)
    if(splitted_name[1] == "019"):
        print(format_colored(f"ENM SHAPE: {std_enm_cov.shape}"))
        print(format_colored(f"MD SHAPE: {std_md_cov.shape}"))
    
    try:
        pCov = pearsonr(std_md_cov[idx], std_enm_cov[idx])
        pDisp = pearsonr(md_fluct, np.linalg.norm(enm_displacement.real.reshape(-1,3), axis=1))
        bhtCov = bht_dist(std_md_cov, std_enm_cov.real)
        result_entry = {
            "protein_name" : splitted_name[0] + "_" + splitted_name[1],
            "adj" : splitted_name[2],
            "Kconst" : splitted_name[3],
            "pearsonCovariance" : pCov[0],
            "pearsonDisplacement" : pDisp[0],
            "bhtDistance" : bhtCov,
            "bhtCoef" : np.exp(-bhtCov),
            "isComplex" : "False"
        }
        results.append(result_entry)

        if pCov[0] > 0.7:
            print(format_colored("GOOD RESULT", color_code="32"))
            print(format_colored(f"{protein_name:40} {pCov}", color_code="32"))
        elif pCov[0] > 0.4:
            print(format_colored(f"{protein_name:40} {pCov}", color_code="36"))
        else:
            print(format_colored(f"{protein_name:40} {pCov}", color_code="19"))
    except Exception as e:
        # print(format_colored("{" + f" Doesn't work with {protein_name}"))
        # print(e ,format_colored("}"))
        try:
            pCov = pearsonr(std_md_cov[idx], std_enm_cov[idx].real)
            pDisp = pearsonr(md_fluct, np.linalg.norm(enm_displacement.real.reshape(-1,3), axis=1)) # refer to test_flucts.py
            bhtCov = bht_dist(std_md_cov, std_enm_cov.real)
            result_entry = {
                "protein_name" : splitted_name[0] + "_" + splitted_name[1],
                "adj" : splitted_name[2],
                "Kconst" : splitted_name[3],
                "pearsonCovariance" : pCov[0],
                "pearsonDisplacement" : pDisp[0],
                "bhtDistance" : bhtCov,
                "bhtCoef" : np.exp(-bhtCov),
                "isComplex" : "True"
            }
            results.append(result_entry)
            if pCov[0] > 0.7:
                print(format_colored("GOOD RESULT", color_code="32"))
                print(format_colored(f"{protein_name:40} {pCov}", color_code="32"))
            elif pCov[0] > 0.4:
                print(format_colored(f"{protein_name:40} {pCov}", color_code="36"))
            else:
                print(format_colored(f"{protein_name:40} {pCov}", color_code="19"))
        except Exception as e:
            print(format_colored("\n{" + f" Doesn't work with {protein_name}" + "}"))
            print(format_colored(f"ENM SHAPE: {std_enm_cov[idx].shape}"))
            print(format_colored(f"MD SHAPE: {std_md_cov[idx].shape}"))

    if results:
        # fields = ["protein_name", "adj", "Kconst", "pearsonCovariance", "isComplex"]
        fields = ["protein_name", "adj", "Kconst", "pearsonCovariance", "pearsonDisplacement", "bhtDistance", "bhtCeof", "isComplex"]
        file_exists = os.path.exists(output_csv_path)
        with open(output_csv_path, "a") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fields)
            if not file_exists:
                writer.writeheader()
            writer.writerows(results)
    
    return

if __name__ == "__main__":
    start = time.time()
    # save all the pearson results in a csv file, but not only one protein, but all of them from the abstracte_names_to_file_names

    file_names_to_abstracted_names = {}
    copy_filenames = []
    with open(csv_path, "r") as csv_file:
        reader = csv.reader(csv_file)
        next(reader)
        for row in reader:
            file_names_to_abstracted_names[row[0]] = row[1]
            copy_filenames.append(row[0])
    abstract_names_to_file_names = {v: k for k, v in file_names_to_abstracted_names.items()}
    
    N_CORES = min(total_cores, len(abstract_names_to_file_names))
    
    abstract_names = list(abstract_names_to_file_names.keys())

    dir_of_enm_covs = {}
    dir_of_enm_displacements = {}
    for abstract_name in abstract_names:
        enm_covs, enm_displacements = load_enm_data(abstract_name)
        dir_of_enm_covs.update(enm_covs)
        dir_of_enm_displacements.update(enm_displacements)
    
    
    # print(format_colored(f"Loaded ENM covariance matrices and displacements for {len(dir_of_enm_covs)} proteins", color_code = "36"))
    # print(f"Example: {list(dir_of_enm_covs.keys())[0]}")
    # print(f"Example: {list(dir_of_enm_displacements.keys())[0]}")
    # print("\n")
    
    # abstract_name = "protein_001"
    # for file in os.listdir(output_directory + f"/ENM_models/enm_covariances/{abstract_name}"):
    #     if file.endswith(".npy"):
    #         name_of_file = extract_filename(file)
            # dir_of_enm_covs[name_of_file] = np.load(output_directory + f"ENM_models/enm_covariances/{abstract_name}/{file}")
            # dir_of_enm_displacements[name_of_file[:-3] + "Disp"] = np.load(output_directory + f"ENM_models/enm_displacements/{abstract_name}/{file[:-7] + 'Disp.npy'}")
            # print(format_colored(f"Loaded ENM covariance matrices and displacements for {name_of_file}", color_code = "36"))
    

    dir_of_md_covs = {}
    for file in os.listdir(output_directory + f"/processed_MD/covariances"):
        dir_of_md_covs[file_names_to_abstracted_names[extract_filename(file)[:-4]]] = np.load(output_directory + f"/processed_MD/covariances/{file}")
    print(format_colored(f"Loaded MD covariance matrices for {len(dir_of_md_covs)} proteins", color_code = "36"))

    dir_of_md_rmsfs = {}
    for file in os.listdir(output_directory + f"/processed_MD/fluctuations"):
        dir_of_md_rmsfs[file_names_to_abstracted_names[extract_filename(file)[:-7]]] = np.load(output_directory + f"/processed_MD/fluctuations/{file}")
    print(format_colored(f"Loaded MD rmsfs for {len(dir_of_md_rmsfs)} proteins", color_code = "36"))

    # before writing to csv, check if the entry already exists in the csv file


    output_csv_path = os.path.expanduser("~/Desktop/soldier_girl/everything/analsys/results.csv")
    protein_names = list(dir_of_enm_covs.keys())
    
    
    for protein_name in protein_names:
        write_to_csv(protein_name, dir_of_enm_covs[protein_name], dir_of_enm_displacements[protein_name[:-3] + "Disp"])
    end = time.time()
    stringTime = "Time taken: " + str(end - start)
    print(format_colored(stringTime, color_code="36"))





#! END of the script










































# def print_all_pearsons(md_cov, dir_of_cov_matrices):
#     std_md_cov = standardize(md_cov)

#     for protein_name, cov in dir_of_cov_matrices.items():
#         std_cov = standardize(cov)
#         idx = np.triu_indices(std_md_cov.shape[0], k=1)
#         if np.iscomplex(std_cov[idx]).any():
#             print("\n")
#             print(format_colored(f"{protein_name:40} Data is complex", color_code="31"))
#             # get the hessian of this abstracted protein from /ENM_models/enm_hessians/{abstract_name}
#             H = np.load(output_directory + f"/ENM_models/enm_hessians/{abstract_name}/{protein_name[:-3]}H.npy")

#             if np.iscomplex(H).any():
#                 print(format_colored(f"{protein_name:40} H is complex", color_code="31"))
#                 continue
#             else:
#                 print(format_colored(f"{protein_name:40} H is real", color_code="32"))

#         try:
#             p = pearsonr(std_md_cov[idx], std_cov[idx])
#             # if pearson result is higher than 0.5 print with color code 32
#             if p[0] > 0.7:
#                 print(format_colored("GOOD RESULT", color_code="32"))
#                 print(format_colored(f"{protein_name:40} {p}", color_code="32"))
#             elif p[0] > 0.4:
#                 print(format_colored(f"{protein_name:40} {p}", color_code="36"))
#             else:
#                 print(format_colored(f"{protein_name:40} {p}", color_code="19"))

        

#         except Exception as e:

#             # print(format_colored("{" + f" Doesn't work with {protein_name}"))
#             # print(e ,format_colored("}"))
#             p = pearsonr(std_md_cov[idx], std_cov[idx].real)
#             if p[0] > 0.7:
#                 print(format_colored("GOOD RESULT", color_code="32"))
#                 print(format_colored(f"{protein_name:40} {p}", color_code="32"))
#             elif p[0] > 0.4:
#                 print(format_colored(f"{protein_name:40} {p}", color_code="36"))
#             else:
#                 print(format_colored(f"{protein_name:40} {p}", color_code="19"))
#             print("\n")