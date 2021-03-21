"""Test for reading CHARMM Param and PSF files."""

import glob, math, os

from Definitions       import dataPath                          , \
                              _FullVerificationSummary          , \
                              outPath                           , \
                              referenceDataPath
from pBabel            import CHARMMParameterFileReader         , \
                              ExportSystem                      , \
                              ImportCoordinates3                , \
                              ImportSystem
from pCore             import Align                             , \
                              logFile                           , \
                              LogFileActive                     , \
                              TestDataSet                       , \
                              TestScriptExit_Fail               , \
                              TextLogFileWriter                 , \
                              YAMLMappingFile_ToObject
from pMolecule         import AtomSelection                     , \
                              SystemGeometryObjectiveFunction
from pMolecule.NBModel import NBModelFull
from pScientific       import Units
from pSimulation       import PrintComponentData

# . The energy values should be the same as CHARMM values to within 1.0e-3 kcal/mole except for non-bonding where the deviations
# . for the larger systems may reach 0.1 or so.

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Molecule names.
_SystemLabels = ( "ava", "citrate_synthase" )

# . Options.
_AbsoluteErrorTolerance         = 0.5
_GradientAbsoluteErrorTolerance = 1.0e-03
_MaximumAtoms                   = 100

# . Parameter sets.
_ParameterPaths = ( "par_all27_prot_na", "par_coa", "par_oaa" )

# . Patterns.
_Patterns = { "ava"              : ( "*:*:C", "*:VAL.2:*", "AAAA:*:*" ), \
              "citrate_synthase" : ( "*:*:C", "*:PRO.*:*", "", "AABW:*:*" ) }

# . Reference file name.
_ReferenceFileName = "CharmmPSFRead.yaml"

#===================================================================================================================================
# . Script.
#===================================================================================================================================
# . Header.
logFile.Header ( )

# . Initialization.
observed = {}

# . Output setup.
dataPath = os.path.join ( dataPath, "charmm" )
outPath  = os.path.join ( outPath , "pdb"    )
if not os.path.exists ( outPath ): os.mkdir ( outPath )

# . Get the parameters.
parameterPaths = []
for parameterPath in _ParameterPaths:
    parameterPaths.append ( os.path.join ( dataPath, parameterPath + ".prm" ) )
parameters = CHARMMParameterFileReader.PathsToParameters ( parameterPaths )

# . Generate systems.
for label in _SystemLabels:

    logFile.Heading ( "System \"{:s}\"".format ( label ), includeBlankLine = True )
    system              = ImportSystem ( os.path.join ( dataPath, label + ".psfx" ), isXPLOR = True, parameters = parameters )
    system.coordinates3 = ImportCoordinates3 ( os.path.join ( dataPath, label + ".chm" ) )
    system.label        = label
    system.DefineNBModel ( NBModelFull.WithDefaults ( ) )
    system.Summary ( )
    energy = system.Energy  ( )
    logFile.Paragraph ( "Energy (kcal/mole) = {:.4f}".format ( energy / Units.Energy_Kilocalories_Per_Mole_To_Kilojoules_Per_Mole ) )

    # . Get the dictionary of energies.
    localObserved = dict ( system.scratch.energyTerms )

    # . Test gradients.
    if len ( system.atoms ) <= _MaximumAtoms:
        of = SystemGeometryObjectiveFunction.FromSystem ( system )
        localObserved["Gradient Error"] = of.TestGradients ( )

    # . Write PDB file and do various sequence tests.
    if outPath is not None:
        ExportSystem ( os.path.join ( outPath, label + ".pdb" ), system, useSegmentEntityLabels = True )
    PrintComponentData ( system.sequence, doFrequencies = False )

    # . Selections.
    table = logFile.GetTable ( columns = [ 30, 6 ] )
    table.Start  ( )
    table.Title  ( "Selections" )
    for pattern in _Patterns[label]:
        table.Entry ( pattern, align = Align.Left )
        table.Entry ( "{:d}".format ( len ( AtomSelection.FromAtomPattern ( system, pattern ) ) ) )
    table.Stop ( )

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
