"""Find CIP labels for several bALA configurations."""

import glob, os

from Definitions import dataPathM
from pBabel      import ImportCoordinates3  , \
                        ImportSystem
from pCore       import logFile             , \
                        TestScriptExit_Fail
from pSimulation import CIPLabelFinder

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Paths.
molPath  = os.path.join ( dataPathM, "mol" )
xyzPath  = os.path.join ( dataPathM, "bAlaConformations" )

#===================================================================================================================================
# . Script.
#===================================================================================================================================
# . Header.
logFile.Header ( )

# . Conformations.
xyzFiles = glob.glob ( os.path.join ( xyzPath, "*.xyz" ) )
xyzFiles.sort ( )

# . Generate the system.
system = ImportSystem ( os.path.join ( molPath, "bAla_c7eq.mol" ) )
system.Summary ( )

# . Initialization.
isOK = True

# . Loop over the structures in the xyz files.
for xyzFile in xyzFiles:
    system.coordinates3 = ImportCoordinates3 ( xyzFile )
    conformation = os.path.split ( xyzFile )[-1].split ( "_", 1 )[-1][0:-4]
    logFile.Heading ( "bALA Configuration " + conformation, includeBlankLine = True )
    results = CIPLabelFinder ( system )
    if results is None:
        localIsOK = False
    else:
        ( ( tCenters, rtCenters, stCenters, utCenters ), ( dCenters, edCenters, zdCenters, udCenters ) ) = results
        localIsOK = ( len ( tCenters  ) == 4 ) and ( len ( dCenters  ) == 0 ) and \
                    ( len ( stCenters ) == 1 ) and ( len ( utCenters ) == 3 )
    isOK = ( isOK and localIsOK )

# . Footer.
logFile.Footer ( )
if not isOK: TestScriptExit_Fail ( )
