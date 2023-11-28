
# ! Just an example code for one protein 
from enm import *
from prody import *
from pylab import *
from MDAnalysis import *
from MDAnalysis.analysis.encore.covariance import covariance_matrix
from MDAnalysis.analysis import align
from MDAnalysis import transformations
from multiprocessing import cpu_count

import numpy as np
import matplotlib.pyplot as plt
import argparse
import os


def load_md_data(topology, trajectory, pdb_file):
    print("\033[1;31mLoading MD data...\033[0m")
    n_cores = cpu_count()

    reference = Universe(pdb_file)
    u = Universe(topology, trajectory)
    prot = u.select_atoms("protein")

    transform = transformations.NoJump(max_threads=n_cores, parallelizable=True, check_continuity=True)
    u.trajectory.add_transformations(transform)
    reference_frame = u.trajectory[0]
    aligner = align.AlignTraj(u, reference=u, select="protein", in_memory=True, max_cores=n_cores)
    aligner.run()
    calphasMD = u.select_atoms("protein")

    

    mdcov = covariance_matrix(u, select="name CA")

    znw = parsePDB(pdb_file)
    calphasANM = znw.select("calpha")
    anm = ANM('ANM analysis')
    anm.buildHessian(calphasANM, cutoff=11.)

    anm.calcModes()
    print(format_colored("Successfully loaded MD data", "!!"))
    return u, calphasMD, mdcov, anm






def save_plot(plot, directory, filename):
    plot.savefig(os.path.join(directory, f"{filename}.png"))
    print(format_colored("Successfully saved plot as"), filename + ".png in directory", directory)

def main(args):
    topology = args.topology
    trajectory = args.trajectory
    pdb_file = args.pdb

    u, calphasMD, mdcov, anm = load_md_data(topology, trajectory, pdb_file)
    create_directory(u.filename)
    # save_universe_as_pickle(u, extract_filename(u.filename) + ".pkl", extract_filename(u.filename))
    
    enm = ENM(pdb_file)
    
    
    enm.filter('CA')
    b = enm.getCoords()
    H = enm.getHessian(adj = 'cutOff')
    enmcov = calculate_covariance_matrix(H)

    std_anm_cov, std_md_cov = standardize(np.triu(calcCovariance(anm[:3]), k=1)), standardize(np.triu(mdcov, k=1))
    std_enm_cov = standardize(np.triu(enmcov, k=1))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    cax = ax.matshow(std_md_cov, vmin=-0.1, vmax=0.1)
    fig.legend(loc="upper right", title="MD")
    fig.colorbar(cax)
    save_plot(fig, extract_filename(u.filename), "MD")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    cax = ax.matshow(std_enm_cov, vmin=-0.4, vmax=0.4)
    fig.colorbar(cax)
    fig.legend(loc="upper right", title="ENM")
    save_plot(fig, extract_filename(u.filename), "ENM")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    cax = ax.matshow(std_anm_cov, vmin=-0.4, vmax=0.4)
    fig.colorbar(cax)
    fig.legend(loc="upper right", title="ANM")
    save_plot(fig, extract_filename(u.filename), "ANM")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="MD Trajectory, ANM and ENM Analysis Script")
    parser.add_argument("--topology", required=True, help="Path to the topology file (e.g., tpr file)")
    parser.add_argument("--trajectory", required=True, help="Path to the trajectory file (e.g., xtc file)")
    parser.add_argument("--pdb", required=True, help="Path to the PDB file for ANM analysis")

    args = parser.parse_args()
    main(args)
