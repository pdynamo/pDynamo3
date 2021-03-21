"""Test for hydrogen building."""

import glob, os

from Definitions               import dataPathM
from pBabel                    import ImportSystem
from pCore                     import logFile               , \
                                      Selection             , \
                                      TestScriptExit_Fail
from pMolecule.MMModel         import MMModelOPLS
from pMolecule.NBModel         import NBModelFull
from pScientific.RandomNumbers import RandomNumberGenerator
from pSimulation               import BuildHydrogenCoordinates3FromConnectivity

# . Header.
logFile.Header ( )

# . Paths.
dataPath = os.path.join ( dataPathM, "mol" )

# . Energy models.
mmModel = MMModelOPLS.WithParameterSet ( "protein" )
nbModel = NBModelFull.WithDefaults ( )

# . Get the files.
molFiles = glob.glob ( os.path.join ( dataPath, "*.mol" ) )

# . Initialization.
numberErrors   = 0
numberFailures = 0

# . Loop over the files.
for molFile in molFiles:

    try:

        # . Get the system.
        try:
            system = ImportSystem ( molFile )
            system.DefineMMModel ( mmModel )
            system.DefineNBModel ( nbModel )
            system.Summary ( )
        except:
            continue

        # . Calculate an energy.
        eBefore = system.Energy ( doGradients = True )

        # . Define all hydrogen positions as undefined.
        for ( i, atom ) in enumerate ( system.atoms ):
            if atom.atomicNumber == 1:
                system.coordinates3.FlagCoordinateAsUndefined ( i )

        # . Build as many undefined coordinates as possible.
        randomNumberGenerator = RandomNumberGenerator.WithSeed ( 957197 )
        BuildHydrogenCoordinates3FromConnectivity ( system, randomNumberGenerator = randomNumberGenerator )

        # . Calculate an energy if all coordinates have been defined.
        if system.coordinates3.numberUndefined > 0:
            numberFailures += 1
            logFile.Paragraph ( "Not all hydrogens have been rebuilt." )
        else:
            eAfter = system.Energy ( doGradients = True )
            logFile.Paragraph ( "Energy difference after rebuilding = {:.1f}.".format ( eAfter - eBefore ) )

    except Exception as e:
        label = os.path.split ( molFile )[-1][:-4]
        logFile.Paragraph ( "Error occurred for {:s}> {:s}.".format ( label, e.args[0] ) )
        numberErrors += 1

# . Footer.
logFile.Footer ( )
if ( ( numberErrors != 0 ) or ( numberFailures != 0 ) ): TestScriptExit_Fail ( )
