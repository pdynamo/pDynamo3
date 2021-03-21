"""Test OPLS protein parameters."""

import glob, os

from Definitions       import structuresPath
from pBabel            import ImportSystem
from pCore             import logFile             , \
                              TestScriptExit_Fail
from pMolecule.MMModel import MMModelOPLS
from pMolecule.NBModel import NBModelFull

#===================================================================================================================================
# . Script.
#===================================================================================================================================
# . Header.
logFile.Header ( )

# . Paths.
dataPath = os.path.join ( structuresPath, "aminoAcids", "mol" )

# . Get all files.
molFiles = glob.glob ( os.path.join ( dataPath, "*.mol" ) )
molFiles.sort ( )

# . Read all mol files.
numberFailed = 0
for molFile in molFiles:

    logFile.Paragraph ( "Processing " + molFile + ":" )
    molecule = ImportSystem ( molFile )
    try:
        molecule.DefineMMModel ( MMModelOPLS.WithParameterSet ( "protein" ) )
        molecule.DefineNBModel ( NBModelFull.WithDefaults ( ) )
        molecule.Summary ( )
        molecule.Energy  ( doGradients = True )
    except Exception as e:
        numberFailed += 1
        logFile.Paragraph ( "Error occurred> " +  e.args[0] )

# . Summary of results.
logFile.SummaryOfItems ( [ ( "Successes", "{:d}".format ( len ( molFiles ) - numberFailed ) )   ,
                             ( "Failures" , "{:d}".format (                    numberFailed ) ) ] ,
                             order = False, title = "OPLS Protein Parameter Tests" )

# . Footer.
logFile.Footer ( )
if numberFailed != 0: TestScriptExit_Fail ( )
