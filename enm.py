import numpy as np

from Bio.PDB import PDBParser, parse_pdb_header, Structure, Model, Chain, Residue, Atom
from Bio.PDB.PDBParser import PDBParser

from scipy.spatial.distance import cdist
from scipy.spatial import Delaunay

from .node import Node, MutableObject
from .selection import Selection

from .parser import ConvertToChainList

np.set_printoptions(suppress = True)

class ENM:
	def __init__(self, pdb_obj = None, **kwargs):
		self.atoms = []
		self.info = ''

		# Indicates data type of a given pdb_obj

		# Build from ENM instance
		if isinstance(pdb_obj, ENM):
			self.atoms = pdb_obj.atoms.copy()
			self.info = pdb_obj.info + 'New object\n'
			return None

		pdb_obj = ConvertToChainList(self, pdb_obj, **kwargs)

		try:
			for chain in pdb_obj:
				for residue in chain:
					for atom in residue:
						self.atoms.append(Node(atom))
		except:
			raise TypeError(f'Invalid type --> {type(pdb_obj)}')

	def __repr__(self):
		return self.info

	def __getitem__(self, key):
		return self.atoms[key]

	def __setitem__(self, key, nval):
		raise NotImplementedError(f'Please modify via ENM.atoms[{key}]')

	# Easy to create other ENM class
	def select(self, filt):
		self.filter(filt)
		return self

	def filter(self, filt, cache = None):
		self.info += f' -> Filtering with script "{filt}"\n'
		if not callable(filt):
			filt = Selection(filt, cache)
		self.atoms = [node for node in self.atoms if filt(node, cache)]

	def getCoords(self):
		return np.array([node.coord for node in self.atoms])

	def getHessian(self, **kwargs):
		return self._buildHessian(self.getCoords(), **kwargs)

	def _buildHessian(self, coordinateArray, **kwargs):
		adj_type = kwargs.get('adj', 'cutOff')
		if adj_type == 'delaunay':
			adj = self._adjacancyMatrixDelaunay(coordinateArray)
		elif adj_type == 'cutOff':
			adj = self._adjacancyMatrix(coordinateArray, kwargs.get('cutOff', 11))
		self.edges = adj
		K = self._get_K_homo(adj, 1)
		nij = self._get_nij(coordinateArray)
		grad = self._get_grad(adj)
		D = self._get_D(grad, nij)
		H = D.T.dot(K).dot(D)
		return H

	def _adjacancyMatrix(self, coordinateArray, cutOff):
		return cdist(coordinateArray, coordinateArray) <= cutOff

	def _adjacancyMatrixDelaunay(self, coordinateArray):
		tri = Delaunay(coordinateArray)
		vec = tri.simplices
		n = coordinateArray.shape[0]
		adjacency_matrix = np.zeros((n, n), dtype = coordinateArray.dtype)
		for row in vec:
			for x in row:
				for y in row:
					adjacency_matrix[x, y] = 1
		return adjacency_matrix

	def _get_grad(self, adj):
		Nb = int(np.sum(adj == 1) / 2)
		grad = np.zeros((Nb, len(adj)), int)
		bonds = np.array([(i, j) for i, j in zip(*np.where(adj)) if j > i])
		for k, (i, j) in enumerate(zip(*bonds.T)):
			grad[k,i] = 1
			grad[k,j] = -1
		return grad

	def _get_K_homo(self, adj, k):
		bonds = np.array([(i, j) for i, j in zip(*np.where(adj)) if i > j])
		return np.eye(len(bonds))

	def _get_nij(self, coord, d = 3):
		N = len(coord)
		nij = np.zeros((N, N, d), float)
		for i, j in zip(*np.tril_indices(N, k = -1)):
			rij = coord[i] - coord[j]
			rij = rij / np.linalg.norm(rij)
			nij[i,j] = rij
			nij[j,i] = -rij
		return nij

	def _get_D(self, grad, nij, d = 3):
		Nb = np.sum(grad == 1)
		Na = nij.shape[1]
		D = np.zeros((Nb, int(Na*d)), float)
		for i in range(Nb):
			d0 = np.zeros((Na,d), float)
			j, k = np.where(grad[i] != 0)[0]
			d0[j] = grad[i,j] * nij[j,k]
			d0[k] = grad[i,k] * nij[j,k]
			D[i] = d0.reshape(Na*d)
		return D
