"""Gaussian basis set tests."""

import glob, math, os

from pCore                           import logFile                  , \
                                            TestScriptExit_Fail      , \
                                            YAMLMappingFile_ToObject    
from pMolecule.QCModel               import MNDOParameters
from pMolecule.QCModel.GaussianBases import GaussianBasis            , \
                                            GaussianBasisType
from pScientific.LinearAlgebra       import OrthogonalizationMethod

#===================================================================================================================================
# . Script.
#===================================================================================================================================
# . Header.
logFile.Header ( )

# . Paths.
dataPath    = os.path.join ( os.getenv ( "PDYNAMO3_PARAMETERS" ), "mndoParameters" )
basisSource = os.path.join ( dataPath, "mndostong" )
sources     = set ( [ path for path in glob.glob ( os.path.join ( dataPath, "*" ) ) if os.path.isdir ( path ) ] )
sources.remove ( basisSource )

# . Initialization.
maximumDeviation = 0.0
missing          = set ( )
nDeviations      = 0
nFails           = 0
nOrbital         = 0
nSets            = 0
_Tolerance       = 1.0e-09

# . Loop over parameter sets.
for source in sorted ( sources ):
    for path in sorted ( glob.glob ( os.path.join ( source, "*.yaml" ) ) ):
        mndo       = YAMLMappingFile_ToObject ( path, MNDOParameters )
        basisLabel = mndo.orbitalBasisLabel
        if basisLabel is not None:
            basisPath = os.path.join ( basisSource, basisLabel + ".yaml" )
            if os.path.exists ( basisPath ):
                basis     = YAMLMappingFile_ToObject ( basisPath, GaussianBasis )
                for ( i, zeta ) in enumerate ( mndo.shellExponents ): basis.ScaleShellExponents ( i, zeta )
                report    = basis.Orthonormalize ( )
                deviation = report["Maximum Deviation"]
                nSets    +=1
                if report["Status"] != 0 : nFails      += 1
                if deviation > _Tolerance: nDeviations += 1
                if basis.basisType == GaussianBasisType.Orbital: nOrbital  += 1
                maximumDeviation = max ( maximumDeviation, deviation )
            else: missing.add ( basisLabel )

# . Summarize the results.
items =  [ ( "Number of Sets"          , "{:d}".format ( nSets            ) ) ,
           ( "Number of Failures"      , "{:d}".format ( nFails           ) ) ,
           ( "Number of Deviations"    , "{:d}".format ( nDeviations      ) ) ,
           ( "Maximum Deviation"       , "{:g}".format ( maximumDeviation ) ) ,
           ( "Number of Missing Bases" , "{:d}".format ( len ( missing )  ) ) ]
if nOrbital != nSets: items.append ( ( "Non-Orbital Sets"              , "{:d}".format ( nSets - nOrbital  ) ) )
logFile.SummaryOfItems ( items, order = False, title = "MNDO Parameter Basis Set Test Results" )
if len ( missing ) > 0: logFile.Paragraph ( "Missing bases: {:s}".format ( ", ".join ( sorted ( missing ) ) ) )

# . Footer.
logFile.Footer ( )
isOK = ( nDeviations == 0 ) and ( nFails == 0 ) and ( nOrbital == nSets )
if not isOK: TestScriptExit_Fail ( )
