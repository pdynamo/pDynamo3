"""Import and export systems."""

import glob, os, os.path

from Definitions import dataPath               , \
                        dataPathM              , \
                        outPath
from pBabel      import ExportImportPathOutput , \
                        ExportOptions          , \
                        ExportSystem           , \
                        ImportOptions          , \
                        ImportSystem
from pCore       import logFile                , \
                        LogFileActive          , \
                        TestScriptExit_Fail

#===================================================================================================================================
# . Script.
#===================================================================================================================================
# . Header.
logFile.Header ( )

# . Export and import options.
ExportOptions ( )
ImportOptions ( )

# . Files to test.
toTest = []
for ( rootPath, path ) in ( ( dataPath , "pdb" ) ,
                            ( dataPathM, "mol" ) ,
                            ( dataPathM, "xyz" ) ):
    for tag in ( "bala_c7eq", "water" ):
        toTest.extend ( glob.glob ( os.path.join ( rootPath, path, tag + ".*" ) ) )
toTest.sort ( )

# . Output setup.
outPath = os.path.join ( outPath, "importExport" )
if outPath is not None:
    if not os.path.exists ( outPath ): os.mkdir ( outPath )
    outFiles = glob.glob ( os.path.join ( outPath, "*.*" ) )
    for outFile in outFiles: os.remove ( outFile )

# . Loop over files.
failures = 0
for inPath in toTest:
    logFile.Paragraph ( "Processing \"" + ExportImportPathOutput ( inPath ) + "\":" )
    ( head, tail ) = os.path.split ( inPath )
    try:
        system = ImportSystem ( inPath )
        system.Summary ( )
        ExportSystem ( os.path.join ( outPath, tail ), system )
    except Exception as e:
        failures += 1
        logFile.Paragraph ( "Error occurred> " +  e.args[0] )

# . Footer.
logFile.Footer ( )
if failures != 0: TestScriptExit_Fail ( )
