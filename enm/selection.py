import re
import ast
from .EXTENDER import EXTENDER
DEFAULT_VALUE = False

class Selection:
	def __init__(self, script, SHOW_SCRIPT = False):
		# Precompiling
		# Convert all variables that in EXTENDER to its definition
		script = self.translator(script)

		# Define none existing variables to False
		# Go throw variables
		tree = ast.parse(script)
		self.vars = list()
		def traverse(node):
			if isinstance(node, ast.Name):
				self.vars.append(node.id)
			for child_node in ast.iter_child_nodes(node):
				traverse(child_node)
		traverse(tree)
		self.func = eval('lambda ' + ','.join(self.vars) + ': ' + script)
		if SHOW_SCRIPT:
			print('Filter --> lambda ' + ','.join(self.vars) + ': ' + script)

	def translator(self, script):
		tree = ast.parse(script)
		def traverse(node):
			if isinstance(node, ast.Name) and node.id in EXTENDER:
				node.id = '(' + EXTENDER[node.id] + ')'
			for child_node in ast.iter_child_nodes(node):
				traverse(child_node)
		traverse(tree)
		return ast.unparse(tree)

	def __call__(self, node, cache):
		flags = dict()

		# Make Node, Atom attributes accessible from global
		flags['node'] = node
		flags.update(node.__dict__)

		# Assigns to True atom elements, Ex: 'C or water' it will look 'True or (HOH)'
		# Or you can choose whole Carbon atoms using _C
		# Assign _"Element name" variable to True
		# Assign "Atom name with extension" variable to True
		if hasattr(node, 'element'):
			flags['_' + node.element] = True
		if hasattr(node, 'name'):
			flags[node.name] = True

		if node.parent:
			# Make Residue attributes accessible from global
			flags['residue'] = node.parent
			flags.update(flags['residue'].__dict__)

			# Assign "Residue name" variable to True
			if flags['resname']:
				flags[flags['resname']] = True
		# ans = bool(eval(self.script, flags))
		args = list(map(lambda id: flags.get(id, DEFAULT_VALUE), self.vars))
		ans = self.func(*args)
		cache['idx_elem' + node.element] = cache.get('idx_elem' + node.element, 0) + 1
		cache['idx_' + node.name] = cache.get('idx_' + node.name, 0) + 1
		cache['idx_' + flags['resname']] = cache.get('idx_' + flags['resname'], 0) + 1
		if ans:
			cache['num_elem' + node.element] = cache.get('num_elem' + node.element, 0) + 1
			cache['num_' + node.name] = cache.get('num_' + node.name, 0) + 1
			cache['num_' + flags['resname']] = cache.get('num_' + flags['resname'], 0) + 1
		return ans