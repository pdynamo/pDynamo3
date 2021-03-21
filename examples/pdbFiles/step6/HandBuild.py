"""Hand-build coordinates for PKA."""

import os, os.path, sys
sys.path.append ( os.path.join ( os.path.dirname ( os.path.realpath ( __file__ ) ), "..", "data" ) )

from Definitions       import outPath
from pCore             import logFile       , \
                              Pickle        , \
                              Selection     , \
                              Unpickle
from pMolecule.NBModel import NBModelCutOff
from pSimulation       import BuildHydrogenCoordinates3FromConnectivity

# . Header.
logFile.Header ( )

# . Define the atoms IDs and the parameters.
_ToBuild = ( ( "A:LYS.8:CB" , "A:LYS.8:CA" , "A:LYS.8:N"  , "A:LYS.8:CG"  , 1.54, 109.5, 0.0 ) ,
             ( "I:SER.17:OG", "I:SER.17:CB", "I:SER.17:CA", "A:ATP.400:PG", 1.40, 109.5, 0.0 ) )

# . Get the system.
system = Unpickle ( os.path.join ( outPath, "step4.pkl" ) )
system.Summary ( )

# . Build the coordinates of the missing heavy atoms.
for ( id1, id2, id3, id4, r, theta, phi ) in _ToBuild:
    i = system.sequence.AtomIndex ( id1 )
    j = system.sequence.AtomIndex ( id2 )
    k = system.sequence.AtomIndex ( id3 )
    l = system.sequence.AtomIndex ( id4 )
    system.coordinates3.BuildPointFromDistanceAngleDihedral ( i, j, k, l, r, theta, phi )

# . Build the hydrogen atom coordinates.
BuildHydrogenCoordinates3FromConnectivity ( system )

# . Save the system.
Pickle ( os.path.join ( outPath, "step6.pkl" ), system )

# . Calculate an energy if all atoms now have coordinates.
if system.coordinates3.numberUndefined <= 0:
    system.DefineNBModel ( NBModelCutOff.WithDefaults ( ) )
    system.Energy ( doGradients = True )

# . Footer.
logFile.Footer ( )
