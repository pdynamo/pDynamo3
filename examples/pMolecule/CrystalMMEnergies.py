"""Crystal tests - MM.

The systems are small molecular crystals.

The AMBER files were set up using fDynamo and so the parameters may not be that good.
However, they are fine for testing.

All pure QC optimizations fail except three. This is due to protons in the primary image migrating to the secondary images.

Care must be taken when choosing a QCregion that the total MM charge is zero (or integer).
Likewise some selections (in this case only for GLYGLY) lead to exploded molecules.

For MNDO - results OK with existing parameters.
For DFT  - gradients worse. Require 1.0e-4 step for DFT Cartesian part but this gives bad crystal derivatives.
"""

import math, os.path

from Definitions                            import dataPath                        , \
                                                   _FullVerificationSummary
from pBabel                                 import ImportCoordinates3              , \
                                                   ImportSystem
from pCore                                  import Align                           , \
                                                   Clone                           , \
                                                   logFile                         , \
                                                   LogFileActive                   , \
                                                   Selection                       , \
                                                   TestDataSet                     , \
                                                   TestReal                        , \
                                                   TestScriptExit_Fail
from pMolecule                              import EnergyModelPriority             , \
                                                   SystemGeometryObjectiveFunction
from pMolecule.MMModel                      import MMModelOPLS
from pMolecule.NBModel                      import NBModelCutOff
from pMolecule.QCModel                      import ElectronicState
from pScientific.Geometry3                  import Transformation3                 , \
                                                   Transformation3Container
from pScientific.Symmetry                   import PeriodicBoundaryConditions
from pScientific.ObjectiveFunctionIterators import LBFGSMinimizer
from pSimulation                            import CrystalAnalyzeTransformations   , \
                                                   CrystalExpandToP1

#===================================================================================================================================
# . Parameters for system definitions.
#===================================================================================================================================
# . Paths.
_dataPath = os.path.join ( dataPath, "molecularCrystals" )

# . Molecule definitions.
_keywordLabels = ( "label", "canQCOptimize", "qcCharge", "qcmmCharge", "vacuumEnergy", "crystalEnergy", "spaceGroup", "qcSelection", "crystalParameters" )
_moleculeData  = \
 ( ( "ALAALA"  , False, 0,  0, -178.7095, -677.7298, "I4"     , [  6, 7, 8, 9 ]                                             , { "a" : 17.9850, "c" : 5.1540 } ),
   ( "ALAMET01", False, 0,  0, -135.0084, -652.9681, "P121/c1", [ 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26 ]              , { "a" : 13.089 , "b" : 5.329, "c" : 15.921, "beta" :  108.57 } ),
   ( "AQARUF"  , False, 0,  0, -144.5558, -608.3551, "P61"    , [  6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18 ]         , { "a" : 14.3720, "c" : 9.8282 } ),
   ( "BEVXEF01", False, 0,  0, -240.0604, -878.7555, "P212121", [ 23, 24, 25, 26, 27, 28 ]                                  , { "a" : 9.6590 , "b" : 9.6720, "c" : 10.7390 } ),
   ( "GLYALB"  , False, 0,  1,  -44.8402, -562.9905, "P212121", [ 0, 1, 2, 3, 4, 5, 6 ]                                          , { "a" : 9.6930 , "b" : 9.5240, "c" : 7.5370 } ),
   ( "GLYGLY"  , False, 0, -1, -184.8863, -687.4449, "P121/a1", [ 0, 1, 2, 3, 4, 5, 6 ]                                    , { "a" : 7.812  , "b" : 9.566, "c" : 9.410, "beta" :  124.60 } ),
   ( "GUFQON"  , False, 0,  0, -191.1490, -511.8117, "P212121", [ 13, 14, 15, 16, 17 ]                                      , { "a" : 7.2750 , "b" : 9.0970, "c" : 10.5070 } ),
   ( "HXACAN19", True , 0,  0,  157.7170,   24.6588, "P121/a1", [ 16, 17, 18, 19 ]                                          , { "a" : 12.8720, "b" : 9.3700, "c" : 7.0850, "beta" : 115.6200 } ),
   ( "IWANID"  , False, 0,  0, 4380.9282, 4070.8472, "P121"   , [ 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50 ]  , { "a" : 23.091 , "b" : 5.494, "c" : 17.510, "beta" : 117.88 } ),
   ( "LCDMPP10", True , 0,  0,  157.3277,   40.3629, "P1"     , [  4, 5, 6, 7]                                              , { "a" : 8.067  , "b" : 6.082, "c" : 5.155, "alpha" : 131.7, "beta" : 82.4, "gamma" : 106.6 } ),
   ( "WIRYEB"  , False, 0,  0, -155.7355, -612.6186, "P61"    , [  6, 7, 8, 9, 10, 11, 12, 13, 14, 15 ]                     , { "a" : 14.4240, "c" : 9.9960 } ),
   ( "WABZOO"  , True , 0,  0,  209.9697,  152.9967, "R3:R"   , [  4, 5, 6, 7, 8, 9, 14, 15, 16, 17, 25, 26, 27, 28 ]       , { "a" : 12.5940, "alpha" : 118.0300 } ) )

# . Options for expanding to P1 symmetry.
_defaultRange = range ( 1 ) #range ( -1, 2 )
_aRange       = _defaultRange
_bRange       = _defaultRange
_cRange       = _defaultRange
_numberCells  = len ( _aRange ) * len ( _bRange ) * len ( _cRange )

#===================================================================================================================================
# . Parameters for test.
#===================================================================================================================================
# . Geometry optimization.
_LogFrequency      =  1000
_NumberSteps       = 20000
_GradientTolerance = 0.1

# . Tolerances.
_EnergyTolerance                = 0.1
_GradientAbsoluteErrorTolerance = 5.0e-02

#===================================================================================================================================
# . Class for a crystal test system.
#===================================================================================================================================
class CrystalTestSystem:
    """Crystal test system."""

    def __init__ ( self, **options ):
        """Constructor."""
        for ( attribute, value ) in options.items ( ): setattr ( self, attribute, value )
        self.p1Factor = 1.0

    def GetSystem ( self, doQCMM = False, doQCQC = False, expandToP1 = False, log = logFile, qcModel = None, qcmmModels = None, useSymmetry = True ):
        """Get the system with the energy model defined."""
        # . Read the molecule.
        molecule              = ImportSystem       ( os.path.join ( _dataPath, self.label + ".top" ), mmModel = MMModelOPLS.WithDefaults ( ), log = log )
        molecule.coordinates3 = ImportCoordinates3 ( os.path.join ( _dataPath, self.label + ".crd" ), log = log )
        molecule.label        = self.label
        # . Set up symmetry.
        if useSymmetry:
            molecule.symmetry           = PeriodicBoundaryConditions.FromSpaceGroup ( self.spaceGroup )
            molecule.symmetryParameters = molecule.symmetry.MakeSymmetryParameters ( **self.crystalParameters )
            if expandToP1:
                self.p1Factor = float ( len ( molecule.symmetry.transformations ) * _numberCells )
                molecule      = CrystalExpandToP1 ( molecule, aRange = _aRange, bRange = _bRange, cRange = _cRange )

        # . Set up the QC model.
        if qcModel is not None:
            if doQCMM: charge = self.qcmmCharge
            else:      charge = self.qcCharge
            molecule.electronicState = ElectronicState.WithOptions ( charge = charge )
            if doQCMM:
                molecule.DefineQCModel ( qcModel, qcSelection = Selection.FromIterable ( self.qcSelection ) )
            elif doQCQC:
                molecule.DefineQCModel ( qcModel, qcSelection = Selection.FromIterable ( range ( len ( molecule.atoms ) ) ) )
            else:
                molecule.DefineQCModel ( qcModel )
        # . Set up the NB model.
        hasQCMMModels = ( doQCMM or doQCQC ) and ( qcmmModels is not None ) and ( len ( qcmmModels ) > 0 )
        if ( qcModel is None ) or doQCMM or doQCQC:
            molecule.DefineNBModel ( NBModelCutOff.WithDefaults ( ), assignQCMMModels = ( not hasQCMMModels ) )
            molecule.nbModel.useCentering = True
        # . QC/MM models.
        if hasQCMMModels:
            for key in sorted ( qcmmModels.keys ( ) ):
                molecule.AddEnergyModel ( key, qcmmModels[key], priority = EnergyModelPriority.QCMMModel )
        # . Summary.
        if LogFileActive ( log ):
            molecule.Summary ( log = log )
            log.Paragraph ( "Formula = " + molecule.atoms.FormulaString ( ) + "." )
        # . Finish up.
        return molecule

#===================================================================================================================================
# . Functions.
#===================================================================================================================================
def GeometryOptimizationSummary ( optimizationResults, log = logFile, useSymmetry = False ):
    """Output for geometry optimization results."""
    if LogFileActive ( log ) and ( len ( optimizationResults ) > 0 ):
        # . Parameters.
        tableEnergyWidth    = 20
        tableIntegerWidth   = 20
        tableNameWidth      = 30
        tableParameterWidth = 20
        # . System labels.
        labels = list ( optimizationResults.keys ( ) )
        labels.sort ( )
        # . Energies.
        # . Header.
        table = log.GetTable ( columns = [ tableNameWidth, tableEnergyWidth, tableEnergyWidth, tableEnergyWidth, tableIntegerWidth,tableIntegerWidth ] )
        table.Start   ( )
        table.Title   ( "Optimized Energies" )
        table.Heading ( "Label"           )
        table.Heading ( "Initial Energy"  )
        table.Heading ( "Final Energy"    )
        table.Heading ( "Energy Lowering" )
        table.Heading ( "Function Calls"  )
        table.Heading ( "Convergence"     )
        # . Data.
        for label in labels:
            data = optimizationResults[label]
            table.Entry ( label, align = Align.Left )
            table.Entry ( "{:20.4f}".format ( data["Initial Energy"]                          ) )
            table.Entry ( "{:20.4f}".format ( data["Final Energy"  ]                          ) )
            table.Entry ( "{:20.4f}".format ( data["Final Energy"  ] - data["Initial Energy"] ) )
            table.Entry ( "{:d}".format ( data["Function Calls"] ) )
            if data["Converged"]: table.Entry ( "T" )
            else:                 table.Entry ( "F" )
        table.Stop ( )
        # . Symmetry parameters.
        if useSymmetry:
            table = log.GetTable ( columns = [ tableNameWidth, tableParameterWidth, tableEnergyWidth, tableEnergyWidth, tableEnergyWidth ] )
            table.Start  ( )
            table.Title  ( "Optimized Symmetry Parameters" )
            table.Heading ( "Label"      )
            table.Heading ( "Parameter"  )
            table.Heading ( "Initial"    )
            table.Heading ( "Final"      )
            table.Heading ( "Difference" )
            for label in labels:
                # . Get the data.
                data      = optimizationResults[label]
                spInitial = data["Initial Parameters"]
                spFinal   = data["Final Parameters"  ]
                spLabels  = list ( spInitial.keys ( ) )
                spLabels.sort ( )
                for ( i, spLabel ) in enumerate ( spLabels ):
                    if i == 0: table.Entry ( label, align = Align.Left )
                    else:      table.Entry ( "" )
                    table.Entry ( spLabel, align = Align.Left )
                    table.Entry ( "{:.4f}".format ( spInitial[spLabel]                    ) )
                    table.Entry ( "{:.4f}".format ( spFinal  [spLabel]                    ) )
                    table.Entry ( "{:.4f}".format ( spInitial[spLabel] - spFinal[spLabel] ) )
            table.Stop ( )

def RunCrystalTest ( checkEnergies    = True     ,
                     dataSetTag       = "MM"     ,
                     doQCMM           = False    ,
                     doQCQC           = False    ,
                     expandToP1       = False    ,
                     geometryOptimize = True     ,
                     qcModel          = None     ,
                     qcmmModels       = None     ,
                     testGradients    = True     ,
                     useSymmetry      = True     ):
    """Run a crystal test."""
    # . Header.
    logFile.Header ( )
    # . Set up the systems.
    crystalTestSystems = []
    for values in _moleculeData:
        options = { key : value for ( key, value ) in zip ( _keywordLabels, values ) }
        crystalTestSystems.append ( CrystalTestSystem ( **options ) )
    # . Set up the optimizer.
    if geometryOptimize:
        optimizer = LBFGSMinimizer.WithOptions ( logFrequency         = _LogFrequency      ,
                                                 maximumIterations    = _NumberSteps       ,
                                                 rmsGradientTolerance = _GradientTolerance )
    # . Initialization.
    if checkEnergies:
        energyDifferences   = {}
    if geometryOptimize:
        optimizationResults = {}
        optimizer.Summary ( )
    if testGradients:
        maximumGradientDeviation = 0.0
    # . Loop over the molecules.
    numberEnergyFailures = 0
    for testSystem in crystalTestSystems:
        # . Get the molecule.
        molecule = testSystem.GetSystem ( doQCMM = doQCMM, doQCQC = doQCQC, qcModel = qcModel, qcmmModels = qcmmModels, useSymmetry = useSymmetry )
        # . Analyze the transformations.
        CrystalAnalyzeTransformations ( molecule )
        # . Calculate the energy and check it.
        try:
            energy   = molecule.Energy ( doGradients = True )
            energyOK = True
            if molecule.nbModel is not None: molecule.nbModel.StatisticsSummary ( molecule )
            if checkEnergies:
                if useSymmetry: energyDifferences[molecule.label] = math.fabs ( energy - testSystem.p1Factor * testSystem.crystalEnergy )
                else:           energyDifferences[molecule.label] = math.fabs ( energy -                       testSystem.vacuumEnergy  )
        except Exception as error:
            energyOK = False
            numberEnergyFailures += 1
            logFile.Paragraph ( "Error calculating initial energy: {:s}.".format ( repr ( error ) ) )
        # . Skip if not OK.
        if not energyOK: continue
        # . Test the gradients.
        # . The gradients are definitely more sensitive to the finite-difference step when using fractional coordinates.
        if testGradients:
            of = SystemGeometryObjectiveFunction.FromSystem ( molecule )
            of.IncludeSymmetryParameters ( )
            deviation = of.TestGradients ( delta = 5.0e-6, tolerance = 3.0e-4 )
            maximumGradientDeviation = max ( maximumGradientDeviation, deviation )
        # . Geometry optimize.
        if geometryOptimize and ( ( qcModel is None ) or ( doQCMM ) or ( doQCQC ) or ( ( qcModel is not None ) and ( testSystem.canQCOptimize ) ) ):
            # . Save the symmetry parameters for later.
            if useSymmetry:
                spInitial = molecule.symmetry.crystalSystem.GetUniqueSymmetryParameters ( molecule.symmetryParameters )
            of = SystemGeometryObjectiveFunction.FromSystem ( molecule )
            if useSymmetry: of.IncludeSymmetryParameters ( )
            state               = optimizer.Iterate ( of )
            finalEnergy         = molecule.Energy ( doGradients = True )
            optimizationResults[molecule.label] = { "Initial Energy" : energy, "Final Energy" : finalEnergy, "Function Calls" : state["Function Calls"], "Converged" : state["Converged"] }
            if useSymmetry:
                spFinal = molecule.symmetry.crystalSystem.GetUniqueSymmetryParameters ( molecule.symmetryParameters )
                optimizationResults[molecule.label].update ( { "Initial Parameters" : spInitial, "Final Parameters" : spFinal } )
    # . Results.
    # . Energy differences.
    if checkEnergies:
        # . System labels.
        labels = list ( energyDifferences.keys ( ) )
        labels.sort ( )
        # . Header.
        table = logFile.GetTable ( columns = [ 30, 20 ] )
        table.Start  ( )
        if useSymmetry: table.Title  ( "Differences in Crystal Energies" )
        else:           table.Title  ( "Differences in Vacuum Energies"  )
        table.Heading ( "System"     )
        table.Heading ( "Difference" )
        # . Data
        for label in labels:
            table.Entry ( label, align = Align.Left )
            table.Entry ( "{:20.4f}".format ( energyDifferences[label] ) )
        table.Stop ( )
    # . Optimizations.
    if geometryOptimize:
        GeometryOptimizationSummary ( optimizationResults, useSymmetry = useSymmetry )
    # . Get the observed and reference data.
    observed      = {}
    referenceData = TestDataSet.WithOptions ( label = "Crystal " + dataSetTag + " Energies" )
    if checkEnergies:
        observed.update ( energyDifferences )
        for label in energyDifferences:
            referenceData.AddDatum ( TestReal.WithOptions ( label = label, value = 0.0, parent = referenceData, absoluteErrorTolerance = _EnergyTolerance, toleranceFormat = "{:.3f}", valueFormat = "{:.3f}" ) )
    if testGradients:
        observed["Gradient Error"] = maximumGradientDeviation
        referenceData.AddDatum ( TestReal.WithOptions ( label = "Gradient Error", value = 0.0, parent = referenceData, absoluteErrorTolerance = _GradientAbsoluteErrorTolerance, toleranceFormat = "{:.3f}", valueFormat = "{:.3f}" ) )
    # . Check for success/failure.
    if len ( observed ) > 0:
        results = referenceData.VerifyAgainst ( observed )
        results.Summary ( fullSummary = _FullVerificationSummary )
        isOK    = results.WasSuccessful ( )
    else:
        isOK    = True
    # . Footer.
    logFile.Footer ( )
    isOK = isOK and ( numberEnergyFailures == 0 )
    if not isOK: TestScriptExit_Fail ( )

#===================================================================================================================================
# . Script.
#===================================================================================================================================
if __name__ == "__main__":
    RunCrystalTest ( checkEnergies    = False ,
                     dataSetTag       = "MM"  ,
                     doQCMM           = True  ,
                     doQCQC           = False ,
                     geometryOptimize = True  )
