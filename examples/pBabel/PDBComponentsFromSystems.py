"""Interconvert PDB components and systems."""

import glob, os, os.path

from Definitions import dataPathM
from pBabel      import ImportSystem, PDBComponent
from pCore       import logFile, LogFileActive, TestScriptExit_Fail

#===================================================================================================================================
# . Script.
#===================================================================================================================================
# . Header.
logFile.Header ( )

# . Initialization.
isOK         = True
numberErrors = 0

# . Output setup.
dataPath = os.path.join ( dataPathM, "mol" )

# . Get the files.
molFiles = glob.glob ( os.path.join ( dataPath, "*.mol" ) )
molFiles.sort ( )

# . Loop over the files.
for ( i, molFile ) in enumerate ( molFiles ):

    try:

        # . Get the system.
        system = ImportSystem ( molFile )
        system.Summary ( )

        # . Convert to PDB component.
        component = PDBComponent.FromSystem ( system, label = ( "{:03d}".format ( i ) ) )
        component.Summary ( )

        # . Convert back to a system.
        system = component.ToSystem ( )
        system.Summary ( )

        # . Finish up.
        logFile.Separator ( )

    # . Error.
    except Exception as e:
        numberErrors += 1
        logFile.Paragraph ( "Error occurred> " +  e.args[0] )

# . Footer.
logFile.Footer ( )
if ( not isOK ) or ( numberErrors != 0 ): TestScriptExit_Fail ( )
