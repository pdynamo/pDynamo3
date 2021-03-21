"""Testing for ion mobilities."""

import glob, math, os

from Definitions               import dataPathM                , \
                                      referenceDataPath
from pBabel                    import ImportSystem
from pCore                     import logFile                  , \
                                      TestDataSet              , \
                                      TestScriptExit_Fail      , \
                                      YAMLMappingFile_ToObject
from pScientific.RandomNumbers import RandomNumberGenerator
from pSimulation               import HardSphereIonMobilities

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Reference data path.
referenceDataPath = os.path.join ( referenceDataPath, "IonMobilities.yaml" )

# . Tolerance for acceptability (10%).
_percentErrorTolerance = 10.0

#===================================================================================================================================
# . Script.
#===================================================================================================================================
# . Header.
logFile.Header ( )

# . Define the molecule.
molecule = ImportSystem ( os.path.join ( dataPathM, "xyz", "serineOctamerZD4L4.xyz" ) )
molecule.Summary ( )

# . Define the randomNumberGenerator so as to have reproducible results.
randomNumberGenerator = RandomNumberGenerator.WithSeed ( 314159 )

# . Do the test - 250000+ trajectories are generally necessary.
observed = HardSphereIonMobilities ( molecule, nreflections = 30, ntrajectories = 100000, randomNumberGenerator = randomNumberGenerator )

# . Remove non-real values.
keys = list ( observed.keys ( ) )
for key in keys:
    if not isinstance ( observed[key], float ): del observed[key]

# . Verify the observed data against the reference data.
referenceData = YAMLMappingFile_ToObject ( referenceDataPath, TestDataSet )
results       = referenceData.VerifyAgainst ( observed )
isOK          = results.WasSuccessful ( )
results.Summary ( )

# . Footer.
logFile.Footer ( )
if not isOK: TestScriptExit_Fail ( )
