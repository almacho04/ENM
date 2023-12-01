import numpy as np
import pandas as pd

from Bio.PDB import PDBIO, parse_pdb_header, Structure, Model, Chain, Residue, Atom

from scipy.spatial.distance import cdist
from scipy.spatial import Delaunay

from .node import Node, MutableObject
from .selection import Selection

from .parser import ConvertToChainList

import Bio.Data

import Bio.Data.IUPACData as bdi

import re
import json
import os

class ENM:
    def __init__(self, pdb_obj = None, **kwargs):
        self.dist_mat = None
        self.res_3_letter = []
        self.res_1_letter = []
        self.atoms = []
        self.info = ' -> ENM Object\n'

        current_directory = os.path.dirname(__file__)

        # Construct the path to the DB_1 folder and dENM.json file
        db_1_folder = os.path.join(current_directory, '..', 'DB_1')
        json_file_path_1 = os.path.join(db_1_folder, 'dENM.json')
        json_file_path_2 = os.path.join(db_1_folder, 'sENM10.json')
        json_file_path_3 = os.path.join(db_1_folder, 'sENM13.json')
        json_file_path_4 = os.path.join(db_1_folder, 'sdENM.json')

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
                    if not isinstance(residue, list):  # Check if the element is not a list
                        residue_3_letter_code = residue.resname  # Get the 3-letter code of the residue
                    for atom in residue:
                        if not isinstance(residue, list):
                            sr = residue_3_letter_code
                            self.res_3_letter.append(sr)

                            if atom.name == 'CA':
                                # ! This is not a good way to do this
                                # ! Find a better way to do this!
                                if(sr.title() == "Cym" or sr.title() == "Cyx" or sr.title() =="Ocs"):
                                    self.res_1_letter.append(bdi.protein_letters_3to1["Cys"])
                                elif(sr.title() == "Hid"):
                                    self.res_1_letter.append(bdi.protein_letters_3to1["His"])
                                elif(sr.title() == "Nma" or sr.title() == "Nlm" or sr.title() == "Nag" or sr.title() == "Nbg"): # ! Especially this one
                                    self.res_1_letter.append("X")
                                else:
                                    self.res_1_letter.append(bdi.protein_letters_3to1[sr.title()]) # * This is the correct way to do this, but it doesn't have all the residues
                                
                        self.atoms.append(Node(atom))
        except:
            raise TypeError(f'Invalid type --> {type(pdb_obj)}')

        with open(json_file_path_1, 'r') as json_file:
            self.denm_data = json.load(json_file)

        with open(json_file_path_2, 'r') as json_file:
            self.senm10_data = json.load(json_file)

        with open(json_file_path_3, 'r') as json_file:
            self.senm13_data = json.load(json_file)

        self.sdenm = pd.read_json(json_file_path_4)
        self.amino_acid_kappa_10 = {(entry['Amino Acid 1'], entry['Amino Acid 2']): entry['kappa'] for entry in self.senm10_data}
        self.amino_acid_kappa_13 = {(entry['Amino Acid 1'], entry['Amino Acid 2']): entry['kappa'] for entry in self.senm13_data}

    def _calculate_dist_mat(self):
        return cdist(self.getCoords(), self.getCoords(), 'euclidean')

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
    def numAtoms(self):
        return len(self.atoms)

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
        # Assign edges
        adj_type = kwargs.get('adj', 'cutOff')
        if adj_type == 'delaunay':
            self._adjacencyMatrixDelaunay(coordinateArray)
        elif adj_type == 'cutOff':
            self._adjacencyMatrix(coordinateArray, kwargs.get('cutOff', 11))
        elif adj_type == 'all':
            self.adjacencyMatrix = np.ones((len(coordinateArray), len(coordinateArray)))
            np.fill_diagonal(self.adjacencyMatrix, 0)

        # Assign spring constants
        k_type = kwargs.get("kconst", 'homo')
        if k_type == 'homo':
            # Homogeneous spring constant
            K = self._get_K_homo(self.adjacencyMatrix, 1)

        elif k_type == 'exp':
            #Example to create Kirhoff matrix with exp dis dependence
            K = self._get_K_exp_dist(self.adjacencyMatrix,1)

        elif k_type == 'dENM':
            # Example to create Kirhoff matrix with distance dependence
            K = self._get_K_dENM(self.adjacencyMatrix,1)

        elif k_type == 'sENM10':
            #Example to create Kirhoff matrix with amino acid table
            K = self._get_K_sENM10(self.adjacencyMatrix,1)

        elif k_type == 'sENM13':
            # Example to create Kirhoff matrix with amino acid table
            K = self._get_K_sENM13(self.adjacencyMatrix, 1)

        elif k_type == 'sdENM':
            # Example to create Kirhoff matrix with amino acid and distance table
            K = self._get_K_sdENM(self.adjacencyMatrix, 1)

        nij = self._get_nij(coordinateArray)
        grad = self._get_grad(self.adjacencyMatrix)
        D = self._get_D(grad, nij)
        H = D.T.dot(K).dot(D)
        return H


    ### Adjacency matrix algorithms

    def _adjacencyMatrix(self, coordinateArray, cutOff):
        if isinstance(self.dist_mat, type(None)):
            self.dist_mat = cdist(coordinateArray, coordinateArray)
        elif len(self.dist_mat) != len(coordinateArray):
            self.dist_mat = cdist(coordinateArray, coordinateArray)
        self.adjacencyMatrix = self.dist_mat <= cutOff

    def _adjacencyMatrixDelaunay(self, coordinateArray):
        tri = Delaunay(coordinateArray)
        vec = tri.simplices
        n = coordinateArray.shape[0]
        self.adjacencyMatrix = np.zeros((n, n), dtype = coordinateArray.dtype)
        for row in vec:
            for x in row:
                for y in row:
                    self.adjacencyMatrix[x, y] = 1


    ### Spring constant assignment algorithms

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

    def find_kappa_by_distance(self, distance):
        # refer to ENM/get_median_kappa.py
        default = 0.033
        for entry in self.denm_data:
            interval_start, interval_end = entry['distance_interval']
            if interval_start <= distance < interval_end:
                return entry['kappa']
        return default

    def _get_K_dENM(self, adj, k0):
        bonds = np.array([self.find_kappa_by_distance(self.dist_mat[i][j]) for i, j in zip(*np.where(adj)) if i > j])
        K_distance = np.zeros((len(bonds), len(bonds)))
        np.fill_diagonal(K_distance, bonds)
        return K_distance

    def _get_K_sENM10(self, adj, k0):
        default = 0.911 # refer to ENM/get_median_kappa.py
        bonds = np.array([self.amino_acid_kappa_10.get((min(self.res_1_letter[i], self.res_1_letter[j]), max(self.res_1_letter[i], self.res_1_letter[j])), default) for i, j in zip(*np.where(adj)) if i > j])
        K_distance = np.zeros((len(bonds), len(bonds)))
        np.fill_diagonal(K_distance, bonds)
        return K_distance

    def _get_K_sENM13(self, adj, k0):
        default  = 0.913 # refer to ENM/get_median_kappa.py
        bonds = np.array([self.amino_acid_kappa_13.get((min(self.res_1_letter[i], self.res_1_letter[j]), max(self.res_1_letter[i], self.res_1_letter[j])), default) for i, j in zip(*np.where(adj)) if i > j])
        K_distance = np.zeros((len(bonds), len(bonds)))
        np.fill_diagonal(K_distance, bonds)
        return K_distance

    def find_kappa_by_distance_and_amino_acid(self, distance, res_i, res_j):
        # d1 = self.sdenm.loc[(self.sdenm['Amino Acid 1'] == res_i) & (self.sdenm['Amino Acid 2'] == res_j)]

        # for index, entry in d1.iterrows():
        #     interval_start, interval_end = entry['distance_interval']
        #     if interval_start <= distance < interval_end:
        #         return entry['kappa']
        # return None
        d1 = self.sdenm.loc[(self.sdenm['Amino Acid 1'] == res_i) & (self.sdenm['Amino Acid 2'] == res_j)]
        if not d1.empty:
            distance_mask = (d1['distance_interval'].apply(lambda x: x[0] <= distance < x[1]))
            if distance_mask.any():
                return d1.loc[distance_mask, 'kappa'].values[0]
        return 0.014 # refer to ENM/get_median_kappa.py
        

    def _get_K_sdENM(self, adj, k0):
        
        bonds = np.array([self.find_kappa_by_distance_and_amino_acid(self.dist_mat[i][j], min(self.res_1_letter[i], self.res_1_letter[j]), max(self.res_1_letter[i], self.res_1_letter[j])) for i, j in zip(*np.where(adj)) if i > j])
        
        K_distance = np.zeros((len(bonds), len(bonds)))
        np.fill_diagonal(K_distance, bonds)
        return K_distance
        
    
        
        


    ### Functions for building the Hessian

    def _get_grad(self, adj):
        Nb = int(np.sum(adj == 1) / 2)
        grad = np.zeros((Nb, len(adj)), int)
        bonds = np.array([(i, j) for i, j in zip(*np.where(adj)) if j > i])
        for k, (i, j) in enumerate(zip(*bonds.T)):
            grad[k,i] = 1
            grad[k,j] = -1
        return grad


    def _get_nij(self, coord, d = 3):
        N = len(coord)
        I, J = np.where(self.adjacencyMatrix)
        nij = np.zeros((N, N, d), float)
        rij = (coord[J] - coord[I])
        rij = rij / np.linalg.norm(rij, axis=1).reshape(-1, 1)
        nij[I,J] = rij
        nij[J,I] = - rij
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


