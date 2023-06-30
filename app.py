from enm import *
a = ENM('1znw.pdb')
a = a.filter('CA')

capture = ['coord', 'id', 'name', 'altloc', 'disordered_flag', 'element', 'fullname', 'bfactor', 'occupancy']
print(*capture)

b = list()
for atom in a.atoms:
    kwargs = dict()
    for name in capture:
        kwargs[name] = getattr(atom, name)
    b.append(MutableObject(**kwargs))

b = ENM(b)
b.save('output.pdb', 1, 1)