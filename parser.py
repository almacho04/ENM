import numpy as np

from Bio.PDB import PDBParser, parse_pdb_header, Structure, Model, Chain, Residue, Atom
from Bio.PDB.PDBParser import PDBParser

from .node import Node, MutableObject
def ConvertToChainList(self, pdb_obj, **kwargs):
	# pdb_obj can be:
	#				str (filepath)
	#				Bio.PDB.Structure
	#				Bio.PDB.Model
	#				Bio.PDB.Chain
	#				List of Bio.PDB.Chain
	#				Bio.PDB.Residue
	#				List of Bio.PDB.Residue
	#				Bio.PDB.Atom
	#				List of Bio.PDB.Atom
	#				List of coord tuple -> ((,3) ndarray)



	# file(pdb_obj)_to_structure
	if isinstance(pdb_obj, str) or 'filename' in kwargs:
		filename = pdb_obj or kwargs.get('filename')
		structure_id = kwargs.get('structure') or parse_pdb_header(filename)['idcode'] or 0

		# Needs to add **kwargs
		parser = PDBParser(kwargs.get('QUIET', True))
		pdb_obj = parser.get_structure(structure_id, filename)
		self.info += f' -> File:{filename}\n'
		self.info += f' -> Structure:{pdb_obj.id}\n'

	# structure_to_model
	if isinstance(pdb_obj, Structure.Structure):
		structure = pdb_obj

		pdb_obj = structure[kwargs.get('model', 0)]
		self.info += f' -> Model:{pdb_obj.id}\n'

	# model_to_chains
	if isinstance(pdb_obj, Model.Model):
		model = pdb_obj

		if 'chains' not in kwargs:
			pdb_obj = model.get_chains()
		else:
			pdb_obj = filter(lambda chain:True if chain.get_id() in kwargs.get('chains') else False, model.get_chains())
		pdb_obj = list(pdb_obj)
		self.info += f' -> Chains:{[chain.id for chain in pdb_obj]}\n'

	if isinstance(pdb_obj, np.ndarray) and pdb_obj.shape[1] == 3:
		pdb_obj = [[[MutableObject(coord = node, id = 'CoordAtom') for node in pdb_obj]]]
	elif isinstance(pdb_obj, (Atom.Atom, Node, MutableObject)):
		pdb_obj = [[[pdb_obj]]]
	elif isinstance(pdb_obj, list) and isinstance(pdb_obj[0], (Atom.Atom, Node, MutableObject)):
		pdb_obj = [[pdb_obj]]
	elif isinstance(pdb_obj, Residue.Residue):
		pdb_obj = [[pdb_obj]]
	elif isinstance(pdb_obj, list) and isinstance(pdb_obj[0], Residue.Residue):
		pdb_obj = [pdb_obj]
	elif isinstance(pdb_obj, Chain.Chain):
		pdb_obj = [pdb_obj]
	return pdb_obj