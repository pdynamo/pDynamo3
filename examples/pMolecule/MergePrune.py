"""Testing merging and pruning."""

import glob, os

from Definitions           import dataPath                 , \
                                  _FullVerificationSummary
from pBabel                import ImportSystem
from pCore                 import Clone                    , \
                                  logFile                  , \
                                  Selection                , \
                                  TestDataSet              , \
                                  TestReal                 , \
                                  TestScriptExit_Fail
from pMolecule.MMModel     import MMModelOPLS
from pMolecule.NBModel     import NBModelCutOff
from pScientific.Geometry3 import Vector3

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . The molecules to test.
_Molecules = ( "bAla_c7eq", "chlorideAnion", "water" )

# . The tests.
# . The number of molecules of each type to merge followed by the indices of the molecules to prune.
_Tests = ( ( ( 0, 2, 0 ), ( 0, ) ),
           ( ( 0, 2, 0 ), ( 1, ) ),
           ( ( 0, 4, 0 ), ( 1, 2 ) ),
           ( ( 0, 0, 2 ), ( 0, ) ),
           ( ( 0, 0, 2 ), ( 1, ) ),
           ( ( 0, 0, 4 ), ( 1, 2 ) ),
           ( ( 2, 0, 0 ), ( 0, ) ),
           ( ( 2, 0, 0 ), ( 1, ) ),
           ( ( 4, 0, 0 ), ( 1, 2 ) ),
           ( ( 1, 1, 1 ), ( 0, ) ),
           ( ( 1, 1, 1 ), ( 1, ) ),
           ( ( 1, 1, 1 ), ( 2, ) ),
           ( ( 2, 2, 2 ), ( 1, 3, 5 ) ) )

# . Tolerances.
_EnergyAbsoluteErrorTolerance = 1.0e-04

# . Translation.
_Displacement = 25.0

#===================================================================================================================================
# . Script.
#===================================================================================================================================
# . Header.
logFile.Header ( )

# . Initialization.
dataPath    = os.path.join ( dataPath, "mol" )
translation = Vector3.Uninitialized ( )

# . Get the individual systems.
energies  = []
molecules = []
for ( i, label ) in enumerate ( _Molecules ):
    molecule       = ImportSystem ( os.path.join ( dataPath, label + ".mol" ) )
    molecule.label = label
    molecule.DefineMMModel ( MMModelOPLS.WithParameterSet ( "bookSmallExamples" ) )
    molecule.DefineNBModel ( NBModelCutOff.WithDefaults ( ) )
    molecule.Summary ( )
    molecule.coordinates3.TranslateToCenter ( )
    energies.append  ( molecule.Energy  ( doGradients = True ) )
    molecules.append ( molecule )

# . Data initialization.
observed      = {}
referenceData = TestDataSet.WithOptions ( label = "Merge/Prune Energies" )

# . Loop over the tests.
for ( testIndex, ( moleculeFrequencies, moleculePruneIndices ) ) in enumerate ( _Tests ):

    # . Heading.
    logFile.Heading ( "Merge/Prune Test {:d}".format ( testIndex ), includeBlankLine = True )

    # . Initialization.
    mergedEnergy = 0.0
    prunedEnergy = 0.0
    translation.Set ( 0.0 )

    # . Gather items.
    index        = 0
    numberAtoms  = 0
    pruneIndices = []
    toMerge      = []
    for ( i, frequency ) in enumerate ( moleculeFrequencies ):
        molecule     = molecules[i]
        mergedEnergy += energies[i] * frequency
        for f in range ( frequency ):
            cloned = Clone ( molecule )
            translation[0] += _Displacement
            cloned.coordinates3.Translate ( translation )
            toMerge.append ( cloned )
            if index in moleculePruneIndices:
                pruneIndices.extend ( range ( numberAtoms, numberAtoms + len ( cloned.atoms ) ) )
                prunedEnergy += energies[i]
            index       += 1
            numberAtoms += len ( cloned.atoms )

    # . Merging.
    merged = toMerge[0].__class__.Merge ( toMerge )
    merged.DefineNBModel ( NBModelCutOff.WithDefaults ( ) )
    merged.Summary ( )
    eMerged = merged.Energy ( )

    # . Pruning.
    pruned = merged.Prune ( Selection.FromIterable ( pruneIndices ) )
    pruned.DefineNBModel ( NBModelCutOff.WithDefaults ( ) )
    pruned.Summary ( )
    ePruned = pruned.Energy ( )

    # . Get the observed and reference data.
    for ( tag, eObserved, eReference ) in ( ( "Merged Energy {:d}".format ( testIndex ), eMerged, mergedEnergy ), ( "Pruned Energy {:d}".format ( testIndex ), ePruned, prunedEnergy ) ):
        observed[tag] = eObserved
        referenceData.AddDatum ( TestReal.WithOptions ( label = tag, value = eReference, parent = referenceData, absoluteErrorTolerance = _EnergyAbsoluteErrorTolerance, toleranceFormat = "{:.3f}", valueFormat = "{:.3f}" ) )

# . Finish up.
logFile.Separator ( )

# . Check for success/failure.
if len ( observed ) > 0:
    results = referenceData.VerifyAgainst ( observed )
    results.Summary ( fullSummary = _FullVerificationSummary )
    isOK = results.WasSuccessful ( )
else:
    isOK = True

# . Footer.
logFile.Footer ( )
if not isOK: TestScriptExit_Fail ( )
