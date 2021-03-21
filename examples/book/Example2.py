"""Example 2."""

from Definitions import *

# . Header.
logFile.Header ( )

# . Define the molecule.
moleculeName = "bala_c7eq"
smiles       = "CC(=O)NC(C)C(=O)NC"

# . Read the molecule from three different files.
molecules = [ ImportSystem ( os.path.join ( dataPath, tag, moleculeName + "." + tag ) ) for tag in ( "mol", "pdb", "xyz" ) ]

# . Generate the molecule from a SMILES string.
molecules.append ( SMILESReader.StringToSystem ( smiles ) )

# . Print summaries of the molecules.
for molecule in molecules:
    molecule.Summary ( )

# . Footer.
logFile.Footer ( )
