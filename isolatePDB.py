from pdbfixer import PDBFixer
from openmm.app import PDBFile, PDBxFile, Modeller
from openmm.unit import *
from collections import Counter
from openmm.vec3 import *
import string
import time
import sys
import os

def isolate_chain(pdbid,fill):
    tic = time.time()
    #fixer = PDBFixer(pdbid=pdbid)
    fixer = PDBFixer(pdbid+'.pdb')
    #fixer.removeHeterogens(False) # keep waters

    ## isolate chains
    chains={}
    for chain in fixer.topology.chains():
        #chains[chain.id] = chain
        chains[chain] = chain.id
    print(chains)
    solvent=['HOH','EPE']
    for i in range(len(chains)):
        #fixer = PDBFixer(pdbid=pdbid)
        fixer = PDBFixer(pdbid+'.pdb')
        #fixer.removeHeterogens(False) # keep waters
        c = [key for key in range(len(chains)) if key != i]
        #print(c)
        fixer.removeChains(c)
        if fill:
            fixer.findMissingResidues()
            if len(list(fixer.missingResidues)) > 0: 
                print("Chain {}: missing Residues: {}".format(i, fixer.missingResidues))
            fixer.findNonstandardResidues()
            fixer.replaceNonstandardResidues()
            fixer.findMissingAtoms()
            fixer.addMissingAtoms()
        else:
            fixer.findNonstandardResidues()
            fixer.replaceNonstandardResidues()

        fixer.addMissingHydrogens(7.0)
        # list molecules 
        chainId = list(chains.values())[i]
        for chain in fixer.topology.chains():
            residues = list(chain.residues())
            ids = [int(r.id) for r in residues]
            resnames = [r.name for r in residues]
            nW = 0
            nS = 0
            for res in resnames:
                if res == 'HOH': 
                    nW+=1
                if res in solvent: 
                    nS+=1
            print("Chain {}: chain ID: {}, numRes: {}, numWater: {}, numSolvent: {}, savePDB: {}".format(i,chainId,len(ids),nW,nS,('chain'+str(i)+list(chains.values())[i]+'.pdb')))
            print("Chain {}: {}".format(i,Counter(map(str,resnames))))
            PDBFile.writeFile(fixer.topology,fixer.positions, open(os.path.join(os.path.dirname(__file__),pdbid,'chain'+str(i)+list(chains.values())[i]+'.pdb'), 'w'), keepIds=True)
    toc = time.time()
    dur = toc-tic
    print("isolate chains {:4.1f}ms".format(dur*1000))



pdbid=sys.argv[1]
fill = sys.argv[2] # True or False
if fill: 
    directory=pdbid+'_fill-missing'
else:
    directory=pdbid

if not os.path.exists(directory):
    os.makedirs(directory)
# fill missing residues
#fill = False  

isolate_chain(pdbid,fill)

