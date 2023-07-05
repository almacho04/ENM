import re
import ast
from .EXTENDER import EXTENDER

class Selection:
	def __init__(self, script, cache = None):
		# def_flags used to assign variable default value to False
		# avoids ValueError for variables which are used in a script
		self.def_flags = dict()

		# Precompiling
		# Convert all variables that in EXTENDER to its definition
		self.script = self.translator(script)

		# Define none existing variables to False
		# Go throw variables
		tree = ast.parse(self.script)
		def traverse(node):
			if isinstance(node, ast.Name):
				self.def_flags[node.id] = False
			for child_node in ast.iter_child_nodes(node):
				traverse(child_node)
		traverse(tree)

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
		flags = {**self.def_flags}

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

		if node.get_parent():
			# Make Residue attributes accessible from global
			flags['residue'] = node.get_parent()
			flags.update(flags['residue'].__dict__)

			# Assign "Residue name" variable to True
			if hasattr(flags['residue'], 'resname'):
				flags[flags['residue'].resname] = True
		return bool(eval(self.script, flags))