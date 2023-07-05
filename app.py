from enm import *
a = ENM('1znw.pdb')
a = a.filter('CA')

print(a.atoms)