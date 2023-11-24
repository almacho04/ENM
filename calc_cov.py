from pathlib import Path
import time

from MDAnalysis import Universe
from MDAnalysis.analysis import encore
from MDAnalysis import transformations
import MDAnalysis as mda
import numpy as np
from scipy.signal import convolve2d


def unwrap_and_align(u, trans=['nojump', 'fit_rot_trans']):
    ag = u.select_atoms("protein and name CA")
    ag_ref = u.copy().select_atoms("protein and name CA")

    possible_trans = {'nojump':transformations.nojump.NoJump(ag),
                      'fit_rot_trans': transformations.fit_rot_trans(ag, ag_ref)}
    workflow = [possible_trans[t] for t in trans]
    u.trajectory.add_transformations(*workflow)
    return u


def get_coords(u, nmax=0):
    coords = []
    for ts in u.trajectory:
        print(ts)
        if nmax > 0:
            if ts.frame >= nmax:
                break
        coords.append(u.select_atoms('protein and name CA').positions)
        print(u.select_atoms('protein and name CA').positions)
    return np.array(coords)


def run_one(path):
#   if Path(f"../Data/Covariance/{path.stem}_noalign.npy").exists():
#       return
    U = unwrap_and_align(Universe(path.joinpath("md_0_1.tpr"), path.joinpath("md_0_1.xtc")))
    C = get_coords(U)
    path_out = f"../Data/Trajectories/{path.stem}_align.npy"
    np.save(path_out, C)

#   cov = mda.analysis.encore.covariance.covariance_matrix(U, select='name CA')
#   path_out = f"../Data/Covariance/{path.stem}_align.npy"
#   np.save(path_out, cov)

    U = Universe(path.joinpath("md_0_1.tpr"), path.joinpath("md_0_1.xtc"), trans=['nojump'])
    C = get_coords(U)
    path_out = f"../Data/Trajectories/{path.stem}_noalign.npy"
    np.save(path_out, C)

#   cov = mda.analysis.encore.covariance.covariance_matrix(U, select='name CA')
#   path_out = f"../Data/Covariance/{path.stem}_noalign.npy"
#   np.save(path_out, cov)



def run_all():
    path_list = sorted([p for p in Path("../MD/1ZNW").glob("*") if p.is_dir()])
    for path in path_list:
        print(path)
        try:
            run_one(path)
        except:
            pass


def regularize_covariance(cov):
    cov = convolve2d(cov, np.ones(9).reshape(3,3), mode='valid')[::3,::3]
    diag = cov.diagonal()
    return cov / np.outer(diag, diag)**0.5


def get_rmsf(coord_traj):
    mean_coord = np.mean(coord_traj, axis=0)
    disp_xyz = coord_traj - mean_coord.reshape((1,) + mean_coord.shape)
    rmsf = np.sum(disp_xyz**2, axis=(0,2)) / (len(disp_xyz) - 1)
    return rmsf


if __name__ == "__main__":
    run_all()
    

    
