from pdbfixer import PDBFixer
from openmm.app import PDBFile
from openmm.app import Modeller
from openmm.unit import *
from collections import Counter
from openmm.vec3 import *
from openmm import unit as openmm_unit
import os.path
import string
import time
import sys

def pre_assemble(pdbs):
    tic = time.time()
    for i, pdb in enumerate(pdbs):
        mol = PDBFile(os.path.join(os.path.dirname(__file__),pdb))
        if i == 0: 
            modeller = Modeller(mol.topology, mol.positions)
        else:
            modeller.add(mol.topology, mol.positions)
    '''        
    pdb0 = PDBFile(os.path.join(os.path.dirname(__file__),pdbid,pdbs[0]))
    modeller = Modeller(pdb0.topology, pdb0.positions)
    ## add DNA ##
    pdb1 = PDBFile(os.path.join(os.path.dirname(__file__),pdbid,pdbs[1]))
    modeller.add(pdb1.topology, pdb1.positions)
    pdb2 = PDBFile(os.path.join(os.path.dirname(__file__),pdbid,pdbs[2]))
    modeller.add(pdb2.topology, pdb2.positions)
    '''
    PDBFile.writeFile(modeller.topology,modeller.positions, open(os.path.join(os.path.dirname(__file__),pdbid,'pre_assemble.pdb'), 'w'))
    toc = time.time()
    dur = toc-tic
    print("assemble structures {:4.1f}ms".format(dur*1000))


def solv(conc):
    tic = time.time()
    fixer = PDBFixer(os.path.join(os.path.dirname(__file__),pdbid,'pre_assemble.pdb'))
    c=Quantity(conc, molar)
    #CUBIC box
    # if water box is too big :    chains[-1].id = chr(ord(chains[-2].id)+1)  TypeError: ord() expected a character, but string of length 2 found
    maxSize = max(max((pos[i] for pos in fixer.positions))-min((pos[i] for pos in fixer.positions)) for i in range(3))
    maxSize = maxSize + Quantity(1, nanometer)
    boxSize = maxSize*Vec3(1, 1, 1)
    print(boxSize)
    fixer.addSolvent(boxSize,positiveIon='K+', negativeIon='Cl-', ionicStrength=c)
    '''
    fixer.addSolvent(padding=1.0 * openmm_unit.nanometer,positiveIon='K+', negativeIon='Cl-', ionicStrength=c)
    '''

    PDBFile.writeFile(fixer.topology, fixer.positions, open(os.path.join(os.path.dirname(__file__),pdbid,'pre_input.pdb'), 'w'))
    toc = time.time()
    dur = toc-tic
    print("solvation in {:4.3f}mM KCl {:4.1f}ms".format(conc,dur*1000))


## pre assemble only for protein and DNA molecules

pdbid=sys.argv[1]
n=len(sys.argv)
print(pdbid)
print(n)
print(sys.argv[2:])
pre_assemble(sys.argv[2:])
solv(0.15)

