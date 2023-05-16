import re

# var_name = r'([^a-zA-Z_\"\'])([a-zA-Z_][a-zA-Z0-9_]*)([^a-zA-Z0-9_\"\'])'
var_name = r'[a-zA-Z_][a-zA-Z0-9_]*'
keywords = ['and', 'as', 'assert', 'break', 'class', 'continue', 'def', 'del', 'elif', 'else', 'except', 'False', 'finally', 'for', 'from', 'global', 'if', 'import', 'in', 'is', 'lambda', 'None', 'nonlocal', 'not', 'or', 'pass', 'raise', 'return', 'True', 'try', 'while', 'with', 'yield']

EXTENDER = {
	'all': '(True)',
	'water': '(HOH)',
	'heavy': '(not H)'
}

class Selection:
	def __init__(self, s, lambda_func = False):
		if lambda_func:
			self._separator = s
			return
		self.def_flags = {}
		s = self.translate(s)
		self.script = s
		for i in re.finditer(var_name, s):
			a = i.group()
			x, y = i.span()
			if s[x - 1] == '\"' or s[x - 1] == '\'' or s[y] == '\"' or s[y] == '\'':
				continue
			if a in keywords:
				continue
			self.def_flags[a] = False

	def __repr__(self):
		return f'Node subset of $({s}) atoms'

	def translate(self, s):
		s = '(' + s + ')'
		script = list(s)
		for i in re.finditer(var_name, s):
			a = i.group()
			x, y = i.span()
			if s[x - 1] == '\"' or s[x - 1] == '\'' or s[y] == '\"' or s[y] == '\'':
				continue
			if a in keywords:
				continue
			if a in EXTENDER:
				for j in range(x, y):
					script[j] = ''
				script[x] = EXTENDER[a]
		return ''.join(script)

	def _separator(self, node):
		flags = {**self.def_flags}
		flags['X'], flags['Y'], flags['Z'] = node.coord
		flags['node_mass'] = node.mass

		atom = node.atom
		if atom is not None:
			flags[atom.name] = True
			flags['_' + atom.element] = True
			flags['atom'] = atom

			flags['name'] = atom.name
			flags['element'] = atom.element
			flags['bfactor'] = atom.bfactor
			flags['occupancy'] = atom.occupancy
			flags['mass'] = atom.mass
			residue = atom.get_parent()
			if residue is not None:
				flags[residue.resname] = True
				flags['aminoAcid'] = residue
				flags['resname'] = residue.resname
				flags['hetflag'], flags['resseq'], flags['icode'] = residue.full_id[3]
		return bool(eval(self.script, flags))
