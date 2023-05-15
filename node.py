from . import edge

class Node:
	def __init__(self, coord, atom = None):
		self.atom = atom
		self.coord = coord
		self.mass = 1.0

	@property
	def aminoAcid(self):
		if self.atom is not None:
			return self.atom.get_parent()
		else:
			raise ValueError('Atom is not defined')

	def __repr__(self):
		return self.atom.__repr__()