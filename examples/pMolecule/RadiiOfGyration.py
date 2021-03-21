"""Radii of gyration of various bALA structures."""

import glob, os

from Definitions        import dataPath
from pBabel             import ImportCoordinates3  , \
                               ImportSystem
from pCore              import Align               , \
                               logFile             , \
                               TestScriptExit_Fail
from pScientific.Arrays import Array

#===================================================================================================================================
# . Script.
#===================================================================================================================================
# . Header.
logFile.Header ( )

# . Paths.
molPath  = os.path.join ( dataPath, "mol" )
xyzPath  = os.path.join ( dataPath, "bAlaConformations" )

# . Conformations.
xyzFiles = glob.glob ( os.path.join ( xyzPath, "*.xyz" ) )
xyzFiles.sort ( )

# . Generate the molecule.
molecule = ImportSystem ( os.path.join ( molPath, "bAla_c7eq.mol" ) )
molecule.Summary ( )

# . Translate to principal axes using masses as weights.
masses = Array.FromIterable ( [ atom.mass for atom in molecule.atoms ] )
molecule.coordinates3.ToPrincipalAxes ( weights = masses )

# . Loop over structures and output at the same time.
table = logFile.GetTable ( columns = [ 20, 10 ] )
table.Start ( )
table.Title ( "Radii Of Gyration" )
table.Heading ( "Conformation" )
table.Heading ( "Value"        )
for xyzFile in xyzFiles:
    conformation     = os.path.split ( xyzFile )[-1].split ( "_", 1 )[-1][0:-4]
    coordinates3     = ImportCoordinates3 ( xyzFile )
    radiusOfGyration = coordinates3.RadiusOfGyration ( weights = masses )
    table.Entry ( conformation, align = Align.Left )
    table.Entry ( "{:.2f}".format ( radiusOfGyration ) )
table.Stop ( )

# . Footer.
logFile.Footer ( )
if False: TestScriptExit_Fail ( )
