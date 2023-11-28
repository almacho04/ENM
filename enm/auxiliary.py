import os
import pickle
import numpy as np
from scipy.signal import convolve2d
from sklearn.decomposition import PCA
from scipy.linalg import sqrtm, eigh
from scipy.stats import chi2

# ! All functions are in alphabetical order ! #

def bht_dist(covariance_matrix1, covariance_matrix2):
    avg_covariance_matrix = (covariance_matrix1 + covariance_matrix2) / 2.0
    eigenvalues, eigenvectors = eigh(avg_covariance_matrix)
    
    idx = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]
    
    cumulative_variance = np.cumsum(eigenvalues) / np.sum(eigenvalues)
    q = np.argmax(cumulative_variance >= 0.95) + 1
    principal_components = eigenvectors[:, :q]
    
    
    proj_covariance_matrix1 = np.dot(np.dot(principal_components.T, covariance_matrix1), principal_components)
    proj_covariance_matrix2 = np.dot(np.dot(principal_components.T, covariance_matrix2), principal_components)
    
    
    avg_covariance_matrix_proj = (proj_covariance_matrix1 + proj_covariance_matrix2) / 2.0
    
    
    det_term1 = np.linalg.det(avg_covariance_matrix_proj)
    det_term2 = np.sqrt(np.abs(np.linalg.det(proj_covariance_matrix1) * np.linalg.det(proj_covariance_matrix2)))
    print(f"det_term1: {det_term1}")
    print(f"det_term2: {det_term2}")
    
    bhattacharyya_dist = 0.5 * np.log(det_term1 / det_term2)
    
    return bhattacharyya_dist

def calculate_covariance_matrix(hessian):
    # Returns covariance matrix from hessian matrix
    inverse_hessian = np.linalg.pinv(hessian)
    variances = np.diag(inverse_hessian)
    covariance_matrix = inverse_hessian.copy()
    print("\033[1;31mCalculated covariance matrix\033[0m")
    return covariance_matrix

def create_directory(directory_name):
    # Creates a directory with the name of the pdb file
    filename = extract_filename(directory_name)
    
    if not os.path.exists(filename):
        print("No directory found named", filename, "\nCreating directory...")
        os.mkdir(filename)



def extract_filename(file_path):
    # Returns the filename without the extension
    filename = os.path.basename(file_path)
    filename_without_extension, _ = os.path.splitext(filename)
    return filename_without_extension


def format_colored(*args, color_code="31"):
    # Returns a string with the arguments formatted with the color code
    formatted_args = [f"\033[{color_code}m{arg}\033[0m" for arg in args]
    return " ".join(formatted_args)



def get_coords(u, nmax=0):
    # Returns the coordinates of the universe object
    # in a form of a numpy array (n_frames, n_atoms, 3)
    print(format_colored("Getting coordinates..."))
    coords = []
    for ts in u.trajectory:
        # print(ts)
        if nmax > 0:
            if ts.frame >= nmax:
                break
        # coords.append(u.select_atoms('protein and name CA').positions)
        coords.append(u.select_atoms('name CA').positions)
    return np.array(coords)

def get_cov(H, val, vec):
    print(format_colored("Started calculating covariance matrix out of Hessian matrix, eigval, eigvec...", color_code="32"))
    val_inv = 1 / val
    print(format_colored("Finished calculating inverse of eigenvalues", color_code="32"))
    val_inv[np.abs(val_inv)>1000] = 0
    print(format_colored("Finished setting eigenvalues to 0 if they are too large", color_code="32"))
    diag_inv = np.zeros(H.shape)
    print(format_colored("Finished creating diagonal matrix", color_code="32"))
    np.fill_diagonal(diag_inv, val_inv)
    print(format_colored("Finished filling diagonal matrix", color_code="32"))
    return vec @ diag_inv @ np.linalg.inv(vec)

def get_enm_displacement(model, **kwargs):
    H = model.getHessian(**kwargs)
    val, vec = np.linalg.eig(H) # eigenvalues and eigenvectors
    idx = np.argsort(val)[6:] # exclude the first 6 modes
    disp = np.sum(vec[:,idx]**2 / val[idx]**2, axis=1) # displacement
    return disp


def get_enm_displacements_from_H(H):
    # Takes in a Hessian matrix and returns the displacements, eigenvalues, and eigenvectors
    val, vec = np.linalg.eig(H) # eigenvalues and eigenvectors
    idx = np.argsort(val)[6:] # exclude the first 6 modes
    disp = np.sum(vec[:,idx]**2 / val[idx]**2, axis=1) # displacement
    return disp, val, vec

def get_ranges(cov):
    # Returns the ranges for the colorbar
    # in a form of [x,y], where x is the minimum and y is the maximum
    ranges = np.quantile(cov.ravel(), [0.1, 0.9])
    # now we want to make sure to round the ranges, so that if ranges[0] is 0.07 make it 0.1 or if it's 6.7 make it 10, etc.
    print("BEFORE ROUNDING: ", ranges)
    ranges = np.around(ranges, decimals=1)

    print("RANGES: ", ranges)
    
    return ranges

def get_rmsf(coord_traj):
    print(format_colored("Getting RMSF..."))
    mean_coord = np.mean(coord_traj, axis=0)
    disp_xyz = coord_traj - mean_coord.reshape((1,) + mean_coord.shape)
    rmsf = np.sum(disp_xyz**2, axis=(0,2)) / (len(disp_xyz) - 1)
    return rmsf

def norm_cov(cov):
    trace_cov = np.trace(cov)
    return cov / trace_cov

def regularize_covariance(cov):
    cov = convolve2d(cov, np.ones(9).reshape(3,3), mode='valid')[::3,::3]
    diag = cov.diagonal()
    return cov / np.outer(diag, diag)**0.5


def save_universe_as_pickle(u, filename, directory):
    print(format_colored("Saving universe object as", filename), "in directory", directory)
    # Saves the universe object as a pickle file
    if not os.path.exists(directory):
        print(format_colored("ERROR: Directory does not exist"), directory)
        exit(1)
    with open(os.path.join(directory, filename), "wb") as f:
        pickle.dump(u, f)
    print(format_colored("Successfully saved universe object as", filename))

def standardize(matrix):
    # Returns a standardized matrix
    matrix = np.triu(matrix, k=1)
    mean = np.mean(matrix)
    std_dev = np.std(matrix)
    standardized_matrix = (matrix - mean) / std_dev
    return standardized_matrix
