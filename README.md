<br>

## About The Project

The Elastic Network Model is a simplified computational model used in the field of structural biology to study the dynamics and flexibility of proteins and other biomolecules. It approximates the system as a network of interconnected beads representing atoms or groups of atoms, connected by springs. The model assumes that the dominant motions of the biomolecule can be described by the harmonic oscillation of these springs.

ENM provides insights into the collective motions and fluctuations of biomolecules, such as protein domains, without explicitly simulating the motion of every atom. By simplifying the system, ENM allows for computationally efficient analysis of large biomolecular systems and can provide valuable information about the global dynamics, conformational changes, and functional motions of biomolecules.

ENM has been widely used in structural biology and biophysics research to understand the relationship between structure, dynamics, and function of biomolecules. It has proven to be a valuable tool for studying protein folding, protein-protein interactions, allosteric regulation, and other aspects of biomolecular behavior.

## Prerequisites

```sh
Python    >= 3.10.6
biopython >= 1.79
numpy     >= 1.23.5
scipy     >= 1.10.1
```

## Examples

### Parsing
```python
from enm import *
a = ENM('1znw.pdb')
```
```python
from enm import *
a = ENM('1znw.pdb')
a = a.filter('CA')
b = ENM(a.getCoords())
H = b.getHessian(adj = 'delaunay')
```
### Filtering
```python
from enm import *
a = ENM('1znw.pdb')
b = ENM(a.filter('CA')) # Selection of Carbon Alpha atoms
b = ENM(a.filter('_C')) # Selection of all Carbon atoms
b = ENM(a.filter('_C and _N')) # Selection of all Carbon, Nitrogen atoms
b = ENM(a.filter('CA and water')) # Selection of all Carbon atoms, water (HOH) molecules
b = ENM(a.filter('VAL')) # Selection of all Valine residue atoms
b = ENM(a.filter('_C and num < 10')) # Selection of all Carbon atoms, limited number of atoms to num
b = ENM(a.filter('resid < 15')) # Selection of first 15 residue atoms
```
### Hessian
```python
from enm import *
a = ENM('1znw.pdb')
a = a.filter('CA')
H1 = a.getHessian(cutOff = 15) #(adj = 'cutOff', cutOff = 15)
H2 = a.getHessian(adj = 'delaunay')
```
### Exporting
```python
from enm import *
a = ENM('1znw.pdb')
a = a.filter('heavy')
a.save('output.pdb')
```

```python
# Parsing via MutableObject
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
```

## API

```sh
ENM Class
initializing:
    ENM(pbj_obj, **kwargs):
        pdb_obj can be:
                   str (filepath)
                   Bio.PDB.Structure
                   Bio.PDB.Model
                   Bio.PDB.Chain
                   List of Bio.PDB.Chain
                   Bio.PDB.Residue
                   List of Bio.PDB.Residue
                   Bio.PDB.Atom
                   List of Bio.PDB.Atom
                   List of coord tuple -> ((N, 3) ndarray)
        **kwargs(Optional):
            filename - str, If pbj_obj is not Defined (default from str(pdb_odb))
            structure - str , id of Structure (default from pdb_file or 0)
            model - str , id of model (default 0)
            chains - str , Ex: 'AB', 'A', 'C', 'ABC'

            For Bio.PDB.PDBParser
            PERMISSIVE = True,
            QUIET = True,
            Arguments:
                - PERMISSIVE - bool, if false, exceptions in
                constructing the SMCRA data structure are fatal. If true (DEFAULT),
                the exceptions are caught, but some residues or atoms will be missing.
                - QUIET - bool, if true (DEFAULT), warnings issued in constructing
                the SMCRA data will be suppressed. If false, they will be shown.

attributes:
    atoms - List of Bio.PDB.Atom
    info  - Information about Protein object
    
    [idx] = atoms[idx]
    
    Optional:
    edges - Adjacency matrix cache
```
```sh 
Node Class
Flexible Bio.PDB.Atom class
* All methods, attributes of Bio.PDB.Atom available
Extra property attributes:
    X - X axis of coordinate
    Y - Y axis of coordinate
    Z - Z axis of coordinate
    residue - Residue object
    resname - Residue name
    residx - Residue id
    chain - Str chain
```
```sh
ENM.Filter method:
    All attributes of Node are available
    All attributes of Residue are available
    Attribute _%ELEMENT_NAME% is setted to true
        Ex: _C --> For Carbon Node (CA, CB, CC, ...) it will be True
    Attribute %ATOM_NAME% is setted to true
        Ex: CA

    num_elem_%ELEMENT_NAME%
    num_%ATOM_NAME%
    num_%RESIDUE_NAME%
        :Number of already chosen

    idx_elem_%ELEMENT_NAME%
    idx_%ATOM_NAME%
    idx_%RESIDUE_NAME%
        :Number of objects viewed before,
            same as num, but counts even not included
```
```sh
Selection Class
    Used for creating lambda function for primitive script
    To view lambda script use .filter('CA or water', SHOW_SCRIPT = True)
```
```sh 
EXTENDER.py
    Add your own expressions to shorten script
    Ex:
        'A_letter_residue' : 'ALA or ARG or ASN or ASP'

```
```sh
ENM.save:
    filepath (Required) - filename for pdb file
    save_infos (Optional) - Default True, bool to log pdb information
            (Slow operation due to rewrite of file for writing log to the beggining of the file)
    rebuild (Optional) - Default False, bool to necessity to rebuild structure
                                (necessary if ENM builded from NOT pdb file)

    to build/save from MutableObject it should have attributes --> :
            'coord', 'id', 'name', 'altloc', 'disordered_flag', 'element', 'fullname', 'bfactor', 'occupancy'
```