from enm import *
from prody import *
from pylab import *

from MDAnalysis import *
from MDAnalysis.analysis.gnm import GNMAnalysis as GNM
from MDAnalysis.analysis.encore.covariance import covariance_matrix
from MDAnalysis import transformations
from multiprocessing import cpu_count

n = cpu_count()
transformations.NoJump
transformations.unwrap
import numpy as np
from scipy.spatial.distance import mahalanobis
from scipy.linalg import eig
import matplotlib.pyplot as plt




topology = "./../testvmd/md_0_1.gro"
trajectory = "./../testvmd/md_0_1.xtc"
u = Universe(topology, trajectory)

calphasMD = u.select_atoms("protein")
print("NUMBER OF ATOMS",calphasMD.atoms)
print("ENDED")
mdcov = covariance_matrix(u, select="name CA")
print("SHAPE OF MDCOV", mdcov.shape)
print("RESMD", mdcov)

ion()
a = ENM('./../testvmd/AF-O53675-F1-model_v4.pdb')
a = a.filter('CA')
b = ENM(a.getCoords())
H = a.getHessian(adj = 'delaunay')
# anm = ANM("Using external Hessian")

inverse_hessian = np.linalg.inv(H)
variances = np.diag(inverse_hessian)
covariances = np.triu(inverse_hessian, k=1) + np.tril(inverse_hessian, k=-1).T
print(covariances)


# znw = parseGromacsModes("./../testvmd/md_0_1.gro")
znw = parsePDB("./../testvmd/AF-O53675-F1-model_v4.pdb")
calphas = znw.select("calpha")
anm = ANM('znw ANM analysis')
anm.buildHessian(calphas, cutoff = 15.)

# anm.setHessian(H)
anm.calcModes()

print(anm.getHessian().shape)

# slowest_mode = anm[0]
print("NUMBER OF ATOMS",anm.numAtoms()) # returns 118
anmcov = calcCovariance(anm[:3]) # returns (354, 354)
print("RESCOV",anmcov)
print(anmcov.shape)
# writeNMD("testing_water.nmd", anm[:3], calphas)
cov_matrix = mdcov

# now plot a heatmap of the covariance matrix
fig = plt.figure()
ax = fig.add_subplot(111)
cax = ax.matshow(cov_matrix, interpolation='nearest')
fig.colorbar(cax)
plt.show()
