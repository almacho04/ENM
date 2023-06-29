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
### Filtering
### Hessian
### Exporting

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
            is_pqr = False,
            Arguments:
                - PERMISSIVE - bool, if false, exceptions in
                constructing the SMCRA data structure are fatal. If true (DEFAULT),
                the exceptions are caught, but some residues or atoms will be missing.
                - QUIET - bool, if true (DEFAULT), warnings issued in constructing
                the SMCRA data will be suppressed. If false, they will be shown.
                - is_pqr - bool, specifies the type of file to be parsed.
                If false (DEFAULT) a .pdb file format is assumed. Set it to true if you
                want to parse a .pqr file instead.
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
    To view lambda script use .filter('CA and water', SHOW_SCRIPT = True)
```
```sh 
EXTENDER.py
    Add your own expressions to shorten script
    Ex:
        'A_letter_residue' : 'ALA or ARG or ASN or ASP'

```
