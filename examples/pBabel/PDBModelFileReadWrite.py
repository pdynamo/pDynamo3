"""Read and write PDB model files."""

import glob, os, os.path

from Definitions import dataPath            , \
                        outPath
from pBabel      import PDBFileReader       , \
                        PDBModel
from pCore       import logFile             , \
                        TestScriptExit_Fail

#===================================================================================================================================
# . Script.
#===================================================================================================================================
# . Header.
logFile.Header ( )

# . Initialization.
isOK         = True
numberErrors = 0

# . Paths.
sourcePath = os.path.join ( dataPath, "pdb"      )
outPath    = os.path.join ( outPath , "pdbModel" )
if not os.path.exists ( outPath ): os.mkdir ( outPath )

# . File names.
pdbFiles = glob.glob ( os.path.join ( sourcePath, "*.pdb" ) )
pdbFiles.sort ( )

# . Loop over the files.
for pdbFile in pdbFiles:

    # . Process file name.
    ( head, fileName ) = os.path.split ( pdbFile )

    # . Header.
    logFile.Paragraph ( "Processing " + fileName + ":" )

    try:

       # . Reading.
        model = PDBFileReader.PathToPDBModel ( pdbFile )
        model.Summary ( title = "Summary for PDB Model from \"{:s}\"".format ( os.path.split ( pdbFile )[-1] ) )

        # . Writing.
        if outPath is not None:
            modelPath = os.path.join ( outPath, fileName[0:-4] + ".model" )
            model.ToModelFile ( modelPath )

            # . Reading again.
            model = PDBModel.FromModelFile ( modelPath )
            model.Summary ( title = "Summary for PDB Model from \"{:s}\"".format ( os.path.split ( modelPath )[-1] ) )

    # . Error.
    except Exception as e:
        numberErrors += 1
        logFile.Paragraph ( "Error occurred> " +  e.args[0] )

# . Footer.
logFile.Footer ( )
if ( not isOK ) or ( numberErrors != 0 ): TestScriptExit_Fail ( )
