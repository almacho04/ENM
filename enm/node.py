from Bio.PDB.Atom import Atom

class MutableObject:
	def __init__(self, **kwargs):
		for x, y in kwargs.items():
			setattr(self, x, y)

class Node(Atom):
	def __init__(self, atom):
		for x, y in atom.__dict__.items():
			setattr(self, x, y)
	@property
	def X(self):
		return self.coord[0]
	@property
	def Y(self):
		return self.coord[1]
	@property
	def Z(self):
		return self.coord[2]
	@property
	def residue(self):
		return self.parent
	@property
	def resname(self):
		return self.parent.resname
	@property
	def resid(self):
		return self.parent.id[1]
	@property
	def chain(self):
		return self.full_id[2]