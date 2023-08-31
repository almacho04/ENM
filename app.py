from enm import *
a = ENM('1znw.pdb')
a = a.filter('CA')
b = ENM(a.getCoords())
H = a.getHessian(adj = 'delaunay')


print(a.atoms)