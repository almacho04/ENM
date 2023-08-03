from turtle import distance

import numpy as np

from Bio.PDB import PDBIO, parse_pdb_header, Structure, Model, Chain, Residue, Atom 

from scipy.spatial.distance import cdist
from scipy.spatial import Delaunay

from .node import Node, MutableObject
from .selection import Selection

from .parser import ConvertToChainList

import re

class ENM:
	def __init__(self, pdb_obj = None, **kwargs):
		self.dist_mat = None
		self.atoms = []
		self.info = ' -> ENM Object\n'

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

		self.dist_mat = self._calculate_dist_mat()

	def _calculate_dist_mat(self):
		coordinates = self.getCoords()
		num_atoms = len(coordinates)
		dist_mat = np.zeros((num_atoms, num_atoms))
		for i in range(num_atoms):
			for j in range(num_atoms):
				dist_mat[i, j] = self.distance_function(coordinates[i], coordinates[j])
		return dist_mat

	def __repr__(self):
		return self.info

	def __getitem__(self, key):
		return self.atoms[key]

	def __setitem__(self, key, nval):
		raise NotImplementedError(f'Please modify via ENM.atoms[{key}]')

	# Easy to create other ENM class
	def filter(self, filt, cache = {}):
		self.info += f' -> Filtering with script "{filt}"\n'
		if not callable(filt):
			filt = Selection(filt)
		ans = list()
		for node in self.atoms:
			if filt(node, cache):
				ans.append(node)
			else:
				if node.parent is not None:
					node.parent.detach_child(node.id)
					node = node.parent
				while node is not None and len(node.child_dict) == 0 and node.parent is not None:
					node.parent.detach_child(node.id)
					node = node.parent
		self.atoms = ans
		return self

	def getCoords(self):
		return np.array([node.coord for node in self.atoms])

	def save(self, filename, save_infos = True, rebuild = False):
		io = PDBIO()

		items = self.atoms

		if rebuild:
			# restore Residues
			last = 0
			for item in items:
				if not hasattr(item, 'parent') or item.parent is None:
					res = Residue.Residue((' ', last + 1, ' '), 'ALIEN', ' ')
					res.add(item)
				last = item.parent.id[1]

			# restore Chains
			prev = Chain.Chain('Z')
			for item in items:
				if not hasattr(item.parent, 'parent') or item.parent.parent is None:
					prev.add(item.parent)
				prev = item.parent.parent

		# Atoms -> Residues -> Chains
		items = dict((item.parent.parent.id ,item.parent.parent) for item in items).values()

		model = Model.Model(0)
		for item in items:
			model.add(item)
		structure = Structure.Structure('0')
		structure.add(model)
		io.set_structure(structure)
		io.save(filename)
		if save_infos:
			with open(filename, 'r') as file:
				content = file.read()
			with open(filename, 'w') as file:
				for remark in self.info.split('\n')[:-1]:
					file.write('REMARK ' + remark[4:] + '\n')
				file.write(content)
				file.close()

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

		#Example to create Kirhoff matrix with exp dis dependence
		K_exp_deis = self._get_K_exp_dist(adj,1)

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

	def distance_function(self,point1, point2):
		# Calculate Euclidean distance between two points
		return np.sqrt((point1[0] - point2[0]) ** 2 + (point1[1] - point2[1]) ** 2 + (point1[2] - point2[2]) ** 2)

	def _get_K_exp_dist(self, adj, k0):
		bonds = np.array([k0 * np.exp(-self.dist_mat[i][j]) for i, j in zip(*np.where(adj)) if i > j])
		K_distance = np.zeros((len(bonds), len(bonds)))
		np.fill_diagonal(K_distance, bonds)
		return K_distance

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
