# DNA-Protein-MD-pipeline

Min-Kang Hsieh

This is a molcular dynamics simulation pipeline for a DNA-protein-ligands complex to probe single-amino acid mutation effects on structure of complex that associated with its enzymatic functions.
# Usage:
Showcase for 9N DNA polymerase (archaeal B-family; PDB ID: 5OMV):
**First**
check and isolate molecules from original pdb file:

python isolatePDB.py False
python isolatePDB.py True

Above, True will fill the missing residues; False will not, and the isolated pdb files will save in 5OMV and 5OMV_full-missing folders, respectively.




