"""Read a PDB model file and convert it to a system using atom data from the original PDB file."""

import os, os.path, sys
sys.path.append ( os.path.join ( os.path.dirname ( os.path.realpath ( __file__ ) ), "..", "data" ) )

from Definitions import dataPath      , \
                        outPath       , \
                        pdbDataPath   , \
                        step2Path
from pBabel      import ExportSystem  , \
                        PDBFileReader , \
                        PDBModel
from pCore       import logFile       , \
                        Pickle

# . Header.
logFile.Header ( )

# . Set the aliases - new/old.
_Aliases = [ ( "A:ATP.400:O3B", "A:ANP.400:N3B" ) ,
             ( "A:ATP.400:*"  , "A:ANP.400:*"   ) ,
             ( "A:MG.*:MG"    , "A:MN.*:MN"     ) ,
             ( "I:SER.17:*"   , "I:ALA.17:*"    ) ]

# . Get the model.
model = PDBModel.FromModelFile ( os.path.join ( step2Path, "step2.model" ) )
model.Summary ( )

# . Get its raw counterpart.
rawModel = PDBFileReader.PathToPDBModel ( os.path.join ( dataPath, "1CDK.pdb" ) )

# . Make the atomic model.
model.MakeAtomicModelFromComponentLibrary ( libraryPaths = [ pdbDataPath ] )
model.ExtractAtomData ( rawModel, aliases = _Aliases )
model.Summary ( )

# . Make a system.
system = model.MakeSystem ( )
system.Summary ( )

# . Save the system.
Pickle ( os.path.join ( outPath, "step3.pkl" ), system )

# . Write out a PDB file as a check as it can be useful to inspect this
# . file to ensure that everything has worked as planned. Something like
# . 9999.0 will appear in records for atoms with undefined coordinates.
ExportSystem ( os.path.join ( outPath, "step3.pdb" ), system )

# . Print out the total charge of the system as a check.
charge = 0
for atom in system.atoms:
    charge += getattr ( atom, "formalCharge", 0 )
logFile.Paragraph ( "Total System Charge = {:d}.".format ( charge ) )

# . Footer.
logFile.Footer ( )
