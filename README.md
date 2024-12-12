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
Pre-assemble protein and DNA molecules, followed by hydration in 150 mM KCl to form pre_assemble.pdb in 5OMV folder.
python pre_assemble.py 5OMV 5OMV/chain0A.pdb 5OMV_full-missing/chain1T.pdb 5OMV_full-missing/chain2P.pdb

**Third**
Prepare a input configuation file for setting up modeling.





