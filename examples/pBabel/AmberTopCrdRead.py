"""Test for reading Amber Top and Crd Files."""

import glob, math, os

from Definitions       import dataPath                 , \
                              _FullVerificationSummary , \
                              outPath                  , \
                              referenceDataPath
from pBabel            import ImportCoordinates3       , \
                              ImportSystem
from pCore             import logFile                  , \
                              LogFileActive            , \
                              TestDataSet              , \
                              TestScriptExit_Fail      , \
                              TextLogFileWriter        , \
                              YAMLMappingFile_ToObject
from pMolecule         import SystemGeometryObjectiveFunction
from pMolecule.NBModel import NBModelFull

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Molecule names.
_SystemLabels = ( "dna", "glucose" )

# . Options.
_AbsoluteErrorTolerance         = 0.5
_GradientAbsoluteErrorTolerance = 1.0e-03
_MaximumAtoms                   = 100

# . Reference file name.
_ReferenceFileName = "AmberTopCrdRead.yaml"

#===================================================================================================================================
# . Script.
#===================================================================================================================================
# . Header.
logFile.Header ( )

# . Initialization.
observed = {}

# . Output setup.
dataPath = os.path.join ( dataPath, "amber" )

# . Loop over the systems.
for label in _SystemLabels:

    # . Read the data.
    system              = ImportSystem       ( os.path.join ( dataPath, label + ".top" ) )
    system.coordinates3 = ImportCoordinates3 ( os.path.join ( dataPath, label + ".crd" ) )

    # . Calculation.
    system.DefineNBModel ( NBModelFull.WithDefaults ( ) )
    system.Summary ( )
    logFile.Paragraph ( "Formula = " + system.atoms.FormulaString ( ) + "." )
    energy = system.Energy ( doGradients = True )

    # . Get the dictionary of energies.
    localObserved = dict ( system.scratch.energyTerms )

    # . Test the gradients.
    if len ( system.atoms ) <= _MaximumAtoms:
        of = SystemGeometryObjectiveFunction.FromSystem ( system )
        localObserved["Gradient Error"] = of.TestGradients ( )

    # . Accumulate observed data.
    observed[label] = localObserved

# . Verify the observed data against the reference data.
referenceDataPath = os.path.join ( referenceDataPath, _ReferenceFileName )
referenceData     = YAMLMappingFile_ToObject ( referenceDataPath, TestDataSet )
results           = referenceData.VerifyAgainst ( observed )
isOK              = results.WasSuccessful ( )
results.Summary ( fullSummary = _FullVerificationSummary )

# . Footer.
logFile.Footer ( )
if not isOK: TestScriptExit_Fail ( )
