"""Read and write mmCIF files."""

import glob, os, os.path

from Definitions import dataPath            , \
                        outPath
from pBabel      import ExportSystem        , \
                        ImportSystem
from pCore       import logFile             , \
                        LogFileActive       , \
                        TestScriptExit_Fail

#===================================================================================================================================
# . Script.
#===================================================================================================================================
# . Header.
logFile.Header ( )

# . Initialization.
isOK         = True
numberErrors = 0

# . Output setup.
sourcePath = os.path.join ( dataPath, "mmcif" )
outPath    = os.path.join ( outPath , "mmcif" )
if outPath is not None:
    if not os.path.exists ( outPath ): os.mkdir ( outPath )
    outFiles = glob.glob ( os.path.join ( outPath, "*.cif" ) )
    for outFile in outFiles: os.remove ( outFile )

# . File names.
cifFiles = glob.glob ( os.path.join ( sourcePath, "*.cif" ) )
cifFiles.sort ( )

# . Loop over the files.
for cifFile in cifFiles:

    # . Process file name.
    ( head, fileName ) = os.path.split ( cifFile )

    # . Header.
    logFile.Paragraph ( "Processing " + fileName + ":" )

    try:

        # . Reading.
        system = ImportSystem ( cifFile, format = "mmcif" )
        system.Summary ( )

        # . Writing.
        if outPath is not None:
            ExportSystem ( os.path.join ( outPath, fileName ), system, format = "mmcif" )

    # . Error.
    except Exception as e:
        numberErrors += 1
        logFile.Paragraph ( "Error occurred> " +  e.args[0] )

# . Footer.
logFile.Footer ( )
if ( not isOK ) or ( numberErrors != 0 ): TestScriptExit_Fail ( )
