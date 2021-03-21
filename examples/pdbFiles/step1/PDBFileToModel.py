"""Read a PDB file and convert it to a PDB model."""

import os, os.path, sys
sys.path.append ( os.path.join ( os.path.dirname ( os.path.realpath ( __file__ ) ), "..", "data" ) )

from Definitions import dataPath      , \
                        outPath
from pBabel      import PDBFileReader , \
                        PDBModel
from pCore       import logFile

# . Header.
logFile.Header ( )

# . File names.
fileName = os.path.join ( dataPath, "1CDK.pdb" )

# . Read the PDB file and create a PDB model.
model = PDBFileReader.PathToPDBModel ( fileName )
model.Summary ( )

# . Write out the PDB model.
( head, tail ) = os.path.split ( fileName )
outName        = os.path.join ( outPath, tail[0:-4] + ".model" )
model.ToModelFile ( outName )

# . Footer.
logFile.Footer ( )
