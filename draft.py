import numpy as np
from scipy.spatial.distance import cdist



#####################################################
### Basic Elastic Network Model

### Adjacancy matrix
def get_A(ca, cut=11):
    return cdist(ca, ca) <= cut


### Laplace operator
def get_grad(A):
    Nb = int(np.sum(A==1)/2)
    grad = np.zeros((Nb, len(A)), int)
    bonds = np.array([(i, j) for i, j in zip(*np.where(A)) if j > i])
    for k, (i, j) in enumerate(zip(*bonds.T)):
        grad[k,i] = 1
        grad[k,j] = -1
    return grad


### Kirchoff matrix with spring constants set by sequence
def get_K(A, kmat, seq):
    seq = np.array(list(seq), int) - 1
    bonds = np.array([(i, j) for i, j in zip(*np.where(A)) if i > j])
    Nb = len(bonds)
    K = np.zeros((Nb, Nb), float)
    np.fill_diagonal(K, [kmat[seq[i],seq[j]] for i, j in zip(*bonds.T)])
    return K


### Kirchoff matrix with homogeneous spring constants
def get_K_homo(A, k):
    bonds = np.array([(i, j) for i, j in zip(*np.where(A)) if i > j])
    return np.eye(len(bonds))


### Direction vector between connected nodes
def get_nij(coord, d=3):
    N = len(coord)
    nij = np.zeros((N, N, d), float)
    for i, j in zip(*np.tril_indices(N, k=-1)):
        rij = coord[i] - coord[j]
        rij = rij / np.linalg.norm(rij)
        nij[i,j] = rij
        nij[j,i] = -rij
    return nij


### D = grad X nij
### This code can proabably be simplied
def get_D(grad, nij, d=3):
    Nb = np.sum(grad==1)
    Na = nij.shape[1]
    D = np.zeros((Nb, int(Na*d)), float)
    for i in range(Nb):
        d0 = np.zeros((Na,d), float)
        j, k = np.where(grad[i] != 0)[0]
        d0[j] = grad[i,j] * nij[j,k]
        d0[k] = grad[i,k] * nij[j,k]
        D[i] = d0.reshape(Na*d)
    return D


### Hessian matrix
### Input = atomic coordinates (N x d=3)
def get_H(ca):
    adj = get_A(ca)
    K = get_K_homo(adj, 1)
    nij = get_nij(ca)
    grad = get_grad(adj)
    D = get_D(grad, nij)
    H = D.T.dot(K).dot(D)
    return H


#####################################################
### Extras




### Green's function
def get_G(ca):
    H = get_H(ca)
    val, vec = np.linalg.eig(H)

    val_inv = 1 / val
    val_inv[np.abs(val_inv)>1000] = 0
    diag_inv = np.zeros(H.shape)
    np.fill_diagonal(diag_inv, val_inv)
    return vec @ diag_inv @ np.linalg.inv(vec)
