"""Gaussian basis set tests."""

import glob, math, os

from  collections                    import defaultdict
from pCore                           import logFile                  , \
                                            TestScriptExit_Fail      , \
                                            YAMLMappingFile_ToObject
from pMolecule.QCModel.GaussianBases import GaussianBasis            , \
                                            GaussianBasisType

#===================================================================================================================================
# . Script.
#===================================================================================================================================
# . Header.
logFile.Header ( )

# . Paths.
dataPath = os.path.join ( os.getenv ( "PDYNAMO3_PARAMETERS" ), "gaussianBasisSets" )
sources  = [ path for path in glob.glob ( os.path.join ( dataPath, "*" ) ) if os.path.isdir ( path ) ]

# . Initialization.
basisTypes       = defaultdict ( int )
maximumDeviation = 0.0
nDeviations      = 0
nFails           = 0
nSets            = 0
nSpherical       = 0
_Tolerance       = 5.0e-03

# . Loop over basis sets.
for source in sorted ( sources ):
    for path in sorted ( glob.glob ( os.path.join ( source, "*.yaml" ) ) ):
        basis     = YAMLMappingFile_ToObject ( path, GaussianBasis )
        report    = basis.Orthonormalize ( )
        deviation = report["Maximum Deviation"]
        nSets +=1
        if basis.isSpherical     : nSpherical  += 1
        if report["Status"] != 0 : nFails      += 1
        if deviation > _Tolerance:
            nDeviations += 1
            logFile.Paragraph ( "Tolerance exceeeded for {:s} ({:.5g} instead of {:.5g}).".format ( path, deviation, _Tolerance ) )
        basisTypes[basis.basisType.name] += 1
        maximumDeviation = max ( maximumDeviation, deviation )

# . Summarize the results.
items = [ ( "Number of Sets"          , "{:d}".format ( nSets              ) ) ,
          ( "Number of Failures"      , "{:d}".format ( nFails             ) ) ,
          ( "Cartesian Sets"          , "{:d}".format ( nSets - nSpherical ) ) ,
          ( "Spherical Harmonic Sets" , "{:d}".format ( nSpherical         ) ) ,
          ( "Number of Deviations"    , "{:d}".format ( nDeviations        ) ) ,
          ( "Maximum Deviation"       , "{:g}".format ( maximumDeviation   ) ) ]
for key in sorted ( basisTypes.keys ( ) ):
    n = basisTypes[key]
    if n > 0: items.append ( ( "{:s} Bases".format ( key ), "{:d}".format ( n ) ) )
logFile.SummaryOfItems ( items, order = False, title = "Gaussian Basis Sets Test Results" )

# . Footer.
logFile.Footer ( )
isOK = ( nDeviations == 0 ) and ( nFails == 0 )
if not isOK: TestScriptExit_Fail ( )
