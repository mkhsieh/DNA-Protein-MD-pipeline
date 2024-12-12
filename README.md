# DNA-Protein-MD-pipeline

Min-Kang Hsieh

This is a molcular dynamics simulation pipeline for a DNA-protein-ligands complex to probe single-amino acid mutation effects on structure of complex that associated with its enzymatic functions.

# Usage:
A showcase for archaeal 9N DNA polymerase (B-family; PDB ID: 5OMV):

**First**

Check and isolate molecules from original pdb file:

python isolatePDB.py False

python isolatePDB.py True

Above, True means to fill the missing residues; False means not to fill missing residues, and the isolated Protein (chain0A.pdb), DNA (chain1T.pdb and chain2P.pdb) are generated in 5OMV and 5OMV_full-missing folders, respectively. A few missing residues at the end C terminus are decided not to be filled but a few DNA residues either missing or incorrect.

**Second**

Pre-assemble protein and DNA molecules, followed by hydration in 150 mM KCl to form pre_assemble.pdb (without solvents) and pre_input.pdb (with solvents) in 5OMV folder.

python pre_assemble.py 5OMV 5OMV/chain0A.pdb 5OMV_full-missing/chain1T.pdb 5OMV_full-missing/chain2P.pdb

**Third**

Prepare a input configuation file for setting up modeling as snp_input_test.txt. Or it is able to generate a list of mutation samples by using input.py.

(ref_structure) (number of mutation) (mutation) (KCl concentration) [ligand1:(pdb source) (residue name) (smiles code)] [ligand2:(pdb source) (residue name) (smiles code)] [ligand3:(pdb source) (residue name) (smiles code)] 

Example:

5OMV/pre_input 1 MET-1-CYS 0.15 1 5OMV/chain3A DTP P(=O)([O-])(OP(=O)([O-])OP(=O)([O-])([O-]))OCC(O1)C(O)CC1(N2C=NC3=C(N)N=CN=C32) 5OMV/chain3A MG Mg+2 5OMV/chain3A MN MN+2

**Fourth**

Sampling for single-amino acid mutations by executing the following.

python snp_sampling.py snp_input_test.txt 

**The project is currently at the data recruiting stage...**




