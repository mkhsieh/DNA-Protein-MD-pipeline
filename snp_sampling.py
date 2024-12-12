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
from test_run import mdrun

def fixing_pdb(pdbid,mut):
    tic = time.time()
    #fixer = PDBFixer(pdbid=pdbid)
    fixer = PDBFixer(pdbid+'.pdb')
    ## detect missing residues: should be filled in isolation step.
    fixer.findMissingResidues()
    if len(fixer.missingResidues) > 0:
        print(fixer.missingResidues)
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(True) # True: keep waters, protien, DNA
    #fixer.removeHeterogens(False) # False: keep protien, DNA 

    for i in range(len(mut)):
        fixer.applyMutations([mut[i]], 'A')
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    fixer.addMissingHydrogens(7.0)

    PDBFile.writeFile(fixer.topology, fixer.positions, open(os.path.join(os.path.dirname(__file__),project_name,'assemble_tmp.pdb'), 'w'))

    toc = time.time()
    dur = toc-tic
    print("mutated strucutre is saved to assemble_tmp.pdb file {:4.1f}ms".format(dur*1000))

def addResfromPdb(hetatm_pdb,keepRes):
    tic = time.time()
    ###
    ### pretreatment for the exsiting small molecule pdb 
    ###

    # create a residue list that compomnets will added directly withot add hydrogen process
    #skipaddH=['MG','NA','MN']
    skipaddH=['MG','NA','MN','DTP']

    # remove bonding list i.e., ions
    rmBond=['MG','NA','MN']

    deleteR = []
    deleteB = []
    for i in range(len(list(keepRes.keys()))):
        print("iteration for keys:", i,  list(keepRes.keys())[i])
        fixer = PDBFixer(hetatm_pdb)
        modeller=Modeller(fixer.topology,fixer.positions)

        # delete undesired residues
        res = [ r for r in modeller.topology.residues() if r.name != list(keepRes.keys())[i] ]
        #print("residue to delete:",res)
        for residue in res:
            deleteR.append(residue)
        modeller.delete(deleteR)

        # remove bonds if keepRes is ion
        if list(keepRes.keys())[i] in rmBond:
            #print("residue for remove bond:",list(keepRes.keys())[i])
            for bond in modeller.topology.bonds():
                deleteB.append(bond)
            #print(deleteB) 
            modeller.delete(deleteB)

        print(f"save residue {[r.name for r in modeller.topology.residues()]} to {list(keepRes.keys())[i]}.pdb")
        PDBFile.writeFile(modeller.topology,modeller.positions, open(os.path.join(os.path.dirname(__file__),project_name,list(keepRes.keys())[i]+'.pdb'), 'w'))

    pdblist = []
    deleteH = []
    addfromXML = True 
    for i in range(len(list(keepRes.keys()))):
        # add hydrogens
        #print(list(keepRes.keys())[i], skipaddH)
        if list(keepRes.keys())[i] not in skipaddH:

            if addfromXML==True:
                print("add H to the pdb",list(keepRes.keys())[i]+'.pdb')
                pdb = PDBFile(os.path.join(os.path.dirname(__file__),project_name,list(keepRes.keys())[i]+'.pdb'))
                modeller=Modeller(pdb.topology,pdb.positions)
                for r in modeller.topology.residues():
                    for atom in r.atoms():
                        #print(atom.name,atom.element, atom.element.symbol, atom.index, atom.residue, atom.id)
                        if atom.element.symbol  == 'H':
                            deleteH.append(atom)
                #print(deleteH)
                modeller.delete(deleteH)
                modeller.loadHydrogenDefinitions(os.path.join(os.path.dirname(__file__),'ligandhydrogens.xml'))
                modeller.addHydrogens() 
                PDBFile.writeFile(modeller.topology,modeller.positions, open(os.path.join(os.path.dirname(__file__),project_name,list(keepRes.keys())[i]+'_mod.pdb'), 'w'))
                pdb_mod = os.path.join(os.path.dirname(__file__),project_name,list(keepRes.keys())[i]+'_mod.pdb')
                pdblist.append(pdb_mod)

        else:
            print("residue does not add hydrogens:",list(keepRes.keys())[i])
            pdb = os.path.join(os.path.dirname(__file__),project_name,list(keepRes.keys())[i]+'.pdb')
            pdblist.append(pdb)
    for i, pdb in enumerate(pdblist):
        addpdb = PDBFile(pdb)
        if i == 0: 
            modeller = Modeller(addpdb.topology,addpdb.positions)  
        else: 
            modeller.add(addpdb.topology,addpdb.positions)
    toc = time.time()
    dur = toc-tic
    print("regenerate ligand and co-factor {:4.1f}ms".format(dur*1000))
    return modeller.topology,modeller.positions

def assemble_sys(hetatm_pdbs,keepRes):
    tic = time.time()
    pdb0 = PDBFile(os.path.join(os.path.dirname(__file__),project_name,'assemble_tmp.pdb'))
    modeller = Modeller(pdb0.topology, pdb0.positions)
    if hetatm_pdbs != None and keepRes != None:
        #print(list(keepRes.keys())[0])
        #print(list(keepRes.values())[0])
        #print(hetatm_pdbs)
        for hetpdb in hetatm_pdbs:
            #print(hetpdb)
            top, pos = addResfromPdb(hetpdb,keepRes)
            modeller.add(top, pos)

    PDBFile.writeFile(modeller.topology,modeller.positions, open(os.path.join(os.path.dirname(__file__),project_name,'assemble.pdb'), 'w'))
    toc = time.time()
    dur = toc-tic
    print("assemble structures {:4.1f}ms".format(dur*1000))

def solv(conc):
    tic = time.time()
    fixer = PDBFixer(os.path.join(os.path.dirname(__file__),project_name,'assemble.pdb'))
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

    PDBFile.writeFile(fixer.topology, fixer.positions, open(os.path.join(os.path.dirname(__file__),project_name,'input.pdb'), 'w'))
    toc = time.time()
    dur = toc-tic
    print("solvation in {:4.3f}mM KCl {:4.1f}ms".format(conc,dur*1000))



inputf = sys.argv[1]
gpuid = sys.argv[2]

#f=open(os.path.join(os.path.dirname(__file__), 'input.txt'))
f=open(os.path.join(os.path.dirname(__file__), inputf))
'''
input.txt
  pdbid, #_numtation, [mut1, mut2, ...], conc, method, [pdb1, residue1, smile1,]  
'''
lines=f.readlines()
lineacc=0
for line in lines:
    lineacc+=1
    word = line.split()
    print(word)
    n=len(word)

    pdbid=word[0]
    nmut=int(word[1])
    mut=[]
    #project_name=pdbid
    if nmut > 0:
        for i in range(nmut):
            mut.append(word[2+i])
            project_name=word[2+i]
    conc=float(word[2+nmut])
    method=int(word[3+nmut])

    directory=project_name
    if not os.path.exists(directory):
        os.makedirs(directory)

    directory=project_name+"/"+str(method)
    if not os.path.exists(directory):
        os.makedirs(directory)

    fixing_pdb(pdbid,mut) # read pdb file
    
    d = (n-4-nmut) // 3
    keep={}
    hetatm_pdbs=[]
    # no ligand
    if d == 0:
        #print("in:", lineacc, pdbid,conc,mut)
        assemble_sys(hetatm_pdbs=None,keepRes=None)
        #solv(conc)

    if d > 0:
        for i in range(4+nmut, n, 3):
            #print(i,word[i],word[i+1],word[i+2])
            hetatm_pdbs.append(os.path.join(os.path.dirname(__file__),word[i]+'.pdb'))
            keep[word[i+1]]=word[i+2]
        hetpdbset=set(hetatm_pdbs)
        hetatm_pdbs=list(hetpdbset)
        #print("in:", lineacc, pdbid,conc,mut,hetatm_pdbs,keep)

        assemble_sys(hetatm_pdbs,keep)
        #solv(conc)
    
    #inPBD=os.path.join(os.path.dirname(__file__),project_name,'input.pdb') # if regenerate hydrateion step
    inPBD=os.path.join(os.path.dirname(__file__),project_name,'assemble.pdb')
    mdrun(inPBD,list(keep.values())[0],project_name,method,gpuid)

