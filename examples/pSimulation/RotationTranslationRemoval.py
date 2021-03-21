"""Test that rotation and translation is correctly removed during a dynamics simulation."""

import glob, os

from Definitions               import dataPathM
from pBabel                    import ImportSystem
from pCore                     import Clone                  , \
                                      logFile                , \
                                      TestDataSet            , \
                                      TestScriptExit_Fail    , \
                                      TestReal
from pMolecule.MMModel         import MMModelOPLS
from pMolecule.NBModel         import NBModelFull
from pScientific.Arrays        import Array
from pScientific.RandomNumbers import NormalDeviateGenerator , \
                                      RandomNumberGenerator
from pSimulation               import LangevinDynamics_SystemGeometry

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Dynamics options.
_NSteps = 20000

# . Tolerances.
_RMSAbsoluteErrorTolerance = 1.0e-03

#===================================================================================================================================
# . Script.
#===================================================================================================================================
# . Header.
logFile.Header ( )

# . Paths.
dataPath = os.path.join ( dataPathM, "mol" )

# . Get the system.
molecule = ImportSystem ( os.path.join ( dataPath, "tyrosineDipeptide.mol" ) )
molecule.DefineMMModel ( MMModelOPLS.WithParameterSet ( "protein" ) )
molecule.DefineNBModel ( NBModelFull.WithDefaults ( ) )
molecule.Summary ( )
molecule.Energy  ( doGradients = True )

# . Save initial coordinates.
reference3 = Clone ( molecule.coordinates3 )

# . Do some dynamics.
normalDeviateGenerator = NormalDeviateGenerator.WithRandomNumberGenerator ( RandomNumberGenerator.WithSeed ( 247171 ) )
LangevinDynamics_SystemGeometry ( molecule                         ,
                                  collisionFrequency     =    25.0 ,
                                  logFrequency           =    1000 ,
                                  normalDeviateGenerator = normalDeviateGenerator ,
                                  steps                  = _NSteps ,
                                  temperature            =   300.0 ,
                                  timeStep               =   0.001 )

# . Check RMSs which should be the same as rotation and translation are removed.
masses = Array.FromIterable ( [ atom.mass for atom in molecule.atoms ] )
rms0   = molecule.coordinates3.RootMeanSquareDeviation ( reference3, weights = masses )
molecule.coordinates3.Superimpose ( reference3, weights = masses )
rms1 = molecule.coordinates3.RootMeanSquareDeviation ( reference3, weights = masses )

# . Get the observed and reference data.
observed      = { "RMS Deviation" : rms1 }
referenceData = TestDataSet.WithOptions ( label = "Rotation/Translation Removal" )
referenceData.AddDatum ( TestReal.WithOptions ( label = "RMS Deviation", value = rms0, parent = referenceData, absoluteErrorTolerance = _RMSAbsoluteErrorTolerance, toleranceFormat = "{:.3f}", valueFormat = "{:.3f}" ) )

# . Check for success/failure.
if len ( observed ) > 0:
    results = referenceData.VerifyAgainst ( observed )
    results.Summary ( )
    isOK = results.WasSuccessful ( )
else:
    isOK = True

# . Footer.
logFile.Footer ( )
if not isOK: TestScriptExit_Fail ( )
