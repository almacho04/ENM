#Spring
class Edge():
	def __init__(self, D, K, pairs = [None, None]):
		self.K = K
		self.D = D
		self._pairs = pairs

	#ENM class
	#@property
	#def model(self):
	#	return self._model

	#Ends of bonding (returns pair of nodes)
	@property
	def pairs(self):
		return self._pairs

	def __repr__(self):
		# return f'K={self.K}, D={self.D}'
		return f'{self.D / self.K}'