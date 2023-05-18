import logging
import numpy as np

from Bio.PDB import PDBParser, Structure, Model, Chain, Residue, Atom
from Bio.PDB.PDBParser import PDBParser
from scipy.spatial.distance import cdist
from scipy.spatial import Delaunay

from .edge import Edge
from .node import Node
from .selection import Selection

np.set_printoptions(suppress = True)

logging.basicConfig(encoding = 'utf-8', level = logging.DEBUG)

class ENM:
	def __init__(self, pdb_obj, **kwargs):
		self.parse(pdb_obj, **kwargs)

	def _reset(self):
		self.selection = Selection('all')
		self.atoms = []
		self.nodes = None
		self.hessian = None
		self.cutOff = None
		self.gamma = None
		self.info = None

	def parse(self, pdb_obj, **kwargs):
		self._reset()
		self.roundCoord = kwargs.get('roundCoord', None)
		self.info = ''
		if isinstance(pdb_obj, str):
			# Default parsing method
			parser = PDBParser(PERMISSIVE = 1)
			structure_id = pdb_obj
			filename = kwargs.get('filename', None)
			if filename is None:
				filename = pdb_obj + '.pdb'
			pdb_obj = parser.get_structure(structure_id, filename)
			logging.info(f'Parsed file "{filename}"')
			logging.info(f'Selected structure "{pdb_obj.id}"')
			self.info += f' -> FILE:{filename}'
			self.info += f' -> STRUCTURE:{pdb_obj.id}'

		if isinstance(pdb_obj, Structure.Structure):
			pdb_obj = next(iter(pdb_obj))
			logging.info(f'Selected model "{pdb_obj.id}"')
			self.info += f' -> MODEL:{pdb_obj.id}'

		if isinstance(pdb_obj, Model.Model):
			pdb_obj = next(iter(pdb_obj))
			logging.info(f'Selected chain "{pdb_obj.id}"')
			self.info += f' -> CHAIN:{pdb_obj.id}'

		if isinstance(pdb_obj, Chain.Chain):
			for residue in pdb_obj:
				for atom in residue:
					#####################################################
					atom.coord = np.array(atom.coord, dtype = np.float64)
					if self.roundCoord is not None:
						atom.coord = atom.coord.round(self.roundCoord)
					self.atoms.append(Node(atom.coord, atom))
		elif isinstance(pdb_obj, np.ndarray):
			if isinstance(pdb_obj[0][0], Atom.Atom):
				for atom in pdb_obj:
					atom.coord = np.array(atom.coord, dtype = np.float64)
					if self.roundCoord is not None:
						atom.coord = atom.coord.round(self.roundCoord)
					self.atoms.append(Node(atom.coord, atom))
			else:
				for coord in pdb_obj:
					coord = np.array(coord, dtype = np.float64)
					if self.roundCoord is not None:
						coord = coord.round(self.roundCoord)
					self.atoms.append(Node(coord))
		else:
			raise TypeError(f'Invalid type --> {pdb_obj}')
		
	def __repr__(self):
		return self.info

	def selector(self, value, lambda_func = False):
		self.selection = Selection(value, lambda_func)
		self.nodes = [node for node in self.atoms if self.selection._separator(node)]

	def getCoords(self):
		if self.nodes is None:
			self.selector('all')
		return np.array([node.coord for node in self.nodes])

	"""
	def getEdges(self, **kwargs):
		assert False
		H = self.getHessian(**kwargs)
		n = H.shape[0]
		assert H.shape[1] == n
		edges = np.empty(shape = H.shape, dtype = Edge)
		print(len(self.nodes))
		print(H.shape)
		# return [[Edge(x, self.gamma, [self.nodes[], self.nodes[]]) for x in y] for y in H]
		for i in range(n):
			for j in range(n):
				edges[i][j] = Edge(H[i][j], self.gamma, [self.nodes[i], self.nodes[j]])
		return edges
	"""

	def getHessian(self, **kwargs):
		if self.hessian is None or len(kwargs) != 0:
			self.hessian = self._buildHessian(self.getCoords(), **kwargs)
		return self.hessian
	
	######################PRIVATE FUNCTIONS

	def _buildHessian(self, coordinateArray, **kwargs):
		self.cutOff = kwargs.get('cutOff', 11)
		self.gamma = kwargs.get('gamma', 1)
		adj = self._adjacancyMatrix(coordinateArray)
		K = self._get_K_homo(adj, 1)
		nij = self._get_nij(coordinateArray)
		grad = self._get_grad(adj)
		D = self._get_D(grad, nij)
		return D.T.dot(K).dot(D)

	def _adjacancyMatrix(self, coordinateArray):
		return cdist(coordinateArray, coordinateArray) <= self.cutOff

	def _adjacancyMatrixWithDelaunay(self, coordinateArray):
		tri = Delaunay(coordinateArray)
		matrix = tri.simplices

		num_nodes = len(coordinateArray)

		# Create an empty adjacency matrix
		adjacency_matrix = np.zeros((num_nodes, num_nodes), dtype=int)

		for row in matrix:
			# Iterate over the four elements (nodes) in the row
			for i in range(4):
				for j in range(4):
					# Mark the corresponding entries in the adjacency matrix as 1
					adjacency_matrix[row[i], row[j]] = 1

		print(matrix)

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
