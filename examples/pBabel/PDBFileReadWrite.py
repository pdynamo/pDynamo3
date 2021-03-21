"""Read and write PDB files."""

import glob, os, os.path

from Definitions      import dataPath            , \
                             outPath
from pBabel           import ExportSystem        , \
                             ImportSystem
from pCore            import logFile             , \
                             LogFileActive       , \
                             TestScriptExit_Fail
from pSimulation      import PrintComponentData

#===================================================================================================================================
# . Script.
#===================================================================================================================================
# . Header.
logFile.Header ( )

# . Initialization.
isOK         = True
numberErrors = 0

# . Paths.
sourcePath = os.path.join ( dataPath, "pdb" )
outPath    = os.path.join ( outPath , "pdb" )
if not os.path.exists ( outPath ): os.mkdir ( outPath )
outFiles = glob.glob ( os.path.join ( outPath, "*.pdb" ) )
for outFile in outFiles: os.remove ( outFile )

# . Loop over the files.
for pdbFile in sorted ( glob.glob ( os.path.join ( sourcePath, "*.pdb" ) ) ):

    # . Process file name.
    ( head, fileName ) = os.path.split ( pdbFile )

    # . Header.
    logFile.Paragraph ( "Processing " + fileName + ":" )

    try:

        # . Reading.
        system       = ImportSystem ( pdbFile )
        system.label = os.path.split ( pdbFile )[-1]
        system.Summary ( )
        PrintComponentData ( system.sequence, doFrequencies = True )

        # . Writing.
        ExportSystem ( os.path.join ( outPath, fileName ), system )

    # . Error.
    except Exception as e:
        numberErrors += 1
        logFile.Paragraph ( "Error occurred> " +  e.args[0] )

# . Footer.
logFile.Footer ( )
if ( not isOK ) or ( numberErrors != 0 ): TestScriptExit_Fail ( )
