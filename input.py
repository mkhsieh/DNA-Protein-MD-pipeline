#from PDAnalysis import pdb_parser, utils, protein, deformation 
from Bio.PDB import * 
from pdbfixer import PDBFixer
from openmm.app import PDBFile
import os.path
import sys


def read_struct2seq(pdb,chain=''):
    # You can use a dict to convert three letter code to one letter code
    d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
    'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
    'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
    'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}


    # run parser
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('struct', pdb)   
    # iterate each model, chain, and residue
    # printing out the sequence for each chain
    sequence =''
    for model in structure:
        for c in model:
            if c.id==chain:
                seq = []
                for residue in c:
                    seq.append(d3to1[residue.resname])
                #print('>some_header\n',''.join(seq))
                sequence=''.join(seq)
    return sequence 


def struct2snp(pdb,chain=''):
    # You can use a dict to convert three letter code to one letter code
    d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
    'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
    'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
    'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}


    # run parser
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('struct', pdb)   
    # iterate each model, chain, and residue
    # printing out the sequence for each chain
    sequence =read_struct2seq(pdb,chain)
    #fasta = open('./test.fasta','w')
    print(sequence)
    snp={}
    a=list(d3to1.values())
    aaa=list(d3to1.keys())
    for i, L in enumerate(sequence):
        #print(L)
        Lto3 = aaa[a.index(L)]
        for S in a:
            if L != S:
                Sto3 = aaa[a.index(S)]
                snp[f'{Lto3}-{i+1}-{Sto3}']=sequence[:i]+S+sequence[i+1:]
                #snp[f'{L}-{i+1}-{S}']=sequence[:i]+S+sequence[i+1:]
                #print(f'{L}-{i+1}-{S}')
                #print(f'{Lto3}-{i+1}-{Sto3}')
                #fasta.write(f'>{Lto3}-{i+1}-{Sto3}\n')
                #fasta.write(f'{sequence[:i]+S+sequence[i+1:]}\n')
    return snp

pdb = './5OMV/pre_assemble.pdb'
#seq1=read_struct2seq(pdb, chain='A')
#print(seq1)
snp=struct2snp(pdb, chain='A')

n = 1500
totn=len(list(snp.keys()))
print(f'totn: {totn}')
nf = int(totn/n)
print(f'nf: {nf}')
endloop = totn % n
print(f'endloop: {endloop}')
ref_struc = '5OMV/pre_input'
nmut = 1
mut =""
conc=0.15
method=1
## change DTP to ATP in 5OMV/chain3A and save as 5OMV/chain3Am
lig_pdb = [0]*3 
lig_res = [0]*3 
lig_smiles = [0]*3 

lig_pdb[0]='5OMV/chain3A'
lig_res[0]='DTP'
lig_smiles[0]='P(=O)([O-])(OP(=O)([O-])OP(=O)([O-])([O-]))OCC(O1)C(O)CC1(N2C=NC3=C(N)N=CN=C32)'

lig_pdb[1]='5OMV/chain3A'
lig_res[1]='MG'
lig_smiles[1]='Mg+2'

lig_pdb[2]='5OMV/chain3A'
lig_res[2]='MN'
lig_smiles[2]='MN+2'

for i in range(nf):
    fafile=f'./snp_{i}.fa'
    inputfile=f'./snp_input_{i}.txt'
    fasta = open(fafile,'w')
    inputf = open(inputfile,'w')
    for j in range(n):
        fasta.write(f'>{i*n+j} {list(snp.keys())[i*n+j]}\n')
        fasta.write(f'>{list(snp.values())[i*n+j]}\n')
        inputf.write(f'{ref_struc} {nmut} {list(snp.keys())[i*n+j]} {conc} {method}') 
        for i in range(len(lig_pdb)):
            inputf.write(f' {lig_pdb[i]} {lig_res[i]} {lig_smiles[i]}')
        inputf.write('\n')
    fasta.close()
    inputf.close()

if endloop > 0:
    fafile=f'./snp_{nf}.fa'
    inputfile=f'./snp_input_{nf}.txt'
    fasta = open(fafile,'w')
    inputf = open(inputfile,'w')
    for j in range(endloop):
        fasta.write(f'>{nf*n+j} {list(snp.keys())[nf*n+j]}\n')
        fasta.write(f'>{list(snp.values())[nf*n+j]}\n')
        inputf.write(f'{ref_struc} {nmut} {list(snp.keys())[nf*n+j]} {conc} {method}')
        for i in range(len(lig_pdb)):
            inputf.write(f' {lig_pdb[i]} {lig_res[i]} {lig_smiles[i]}')
        inputf.write('\n')
    fasta.close()
    inputf.close()
