"""MNDO CI QC/MM tests."""

import math, os.path

from Definitions        import _FullVerificationSummary
from pCore              import logFile                             , \
                               LogFileActive                       , \
                               TestDataSet                         , \
                               TestReal                            , \
                               TestScriptExit_Fail
from pMolecule          import SystemGeometryObjectiveFunction
from pMolecule.NBModel  import QCMMElectrostaticModelMultipoleFull , \
                               QCMMLennardJonesModelFull
from pMolecule.QCModel  import CIMethod                            , \
                               DIISSCFConverger                    , \
                               OccupancyType                       , \
                               QCModelMNDOCI
from pScientific.Arrays import ArrayPrint
from QCMMTestSystems    import qcmmTestSystems

#===================================================================================================================================
# . Parameters for system definitions.
#===================================================================================================================================
# . QC models.
_ElectronicStateKeys = ( "occupancyType", "isSpinRestricted", "numberFractionalHOOs", "permitRestrictedNonSinglet" )
_QCModelKeys         = ( "ciMethod", "activeElectrons", "activeOrbitals", "minimalMultiplicity", "multiplicity", "requiredRoot" )
_QCModels            = ( ( "Alanine Dipeptide"    , ( OccupancyType.Cardinal        , True, 0, True ) , ( CIMethod.SinglesDoubles,  8, 6, 1, 1, 0 ), ) ,
                         ( "Cyclohexane 6"        , ( OccupancyType.Cardinal        , True, 0, True ) , ( CIMethod.SinglesDoubles,  4, 4, 1, 1, 0 ), ) ,
                         ( "Cyclohexane 9"        , ( OccupancyType.Cardinal        , True, 0, True ) , ( CIMethod.SinglesDoubles,  6, 6, 1, 1, 1 ), ) ,
                         ( "Tyrosine Dipeptide 1" , ( OccupancyType.Cardinal        , True, 0, True ) , ( CIMethod.SinglesDoubles,  6, 6, 1, 3, 0 ), ) ,
                         ( "Tyrosine Dipeptide 2" , ( OccupancyType.FractionalFixed , True, 1, True ) , ( CIMethod.SinglesDoubles,  5, 6, 2, 2, 0 ), ) ,
                         ( "Water Dimer 1"        , ( OccupancyType.Cardinal        , True, 0, True ) , ( CIMethod.SinglesDoubles,  8, 6, 1, 1, 1 ), ) ,
                         ( "Water Dimer 2"        , ( OccupancyType.Cardinal        , True, 0, True ) , ( CIMethod.SinglesDoubles,  8, 6, 1, 3, 0 ), ) )

#===================================================================================================================================
# . Parameters for test.
#===================================================================================================================================
# . Options.
_MaximumAtoms  = 100
_TestGradients = True
_UseMultipoles = False

# . Tolerances.
_EnergyTolerance                = 0.1
_GradientAbsoluteErrorTolerance = 1.0e-02

# . Set up the systems.
_QCMMTestSystems = qcmmTestSystems

# . QC/MM models.
if _UseMultipoles:
    _MultipoleModel0 = QCMMElectrostaticModelMultipoleFull.WithDefaults ( )
    _MultipoleModel2 = QCMMElectrostaticModelMultipoleFull.WithOptions  ( multipoleOrder = 2 )
    _QCMMLJModel     = QCMMLennardJonesModelFull.WithDefaults ( )
    _QCMMModels      = { "qcmmElectrostatic" : _MultipoleModel0 ,
                         "qcmmLennardJones"  : _QCMMLJModel     }
else:
    _QCMMModels = None

#===================================================================================================================================
# . Script.
#===================================================================================================================================
# . Header.
logFile.Header ( )

# . Initialization.
numberErrors = 0
if _TestGradients:
    maximumGradientDeviation = 0.0

# . Loop over systems.
for ( label, esData, qcData ) in _QCModels:
    testSystem = _QCMMTestSystems[label]

    # . Get the electronic state options.
    esOptions = { key : value for ( key, value ) in zip ( _ElectronicStateKeys, esData ) }

    # . Define the QC model.
    qcOptions                    = { key : value for ( key, value ) in zip ( _QCModelKeys, qcData ) }
    qcOptions["hamiltonian"    ] = "am1"
    qcOptions["converger"      ] = DIISSCFConverger.WithOptions ( densityTolerance = 1.0e-14, maximumIterations = 250 )
    qcModel = QCModelMNDOCI.WithOptions ( **qcOptions )

    # . Get the molecule.
    molecule = testSystem.GetSystem ( electronicStateOptions = esOptions, qcModel = qcModel, qcmmModels = _QCMMModels )

    # . Energy.
    try:
        energy  = molecule.Energy ( doGradients = True )

        # . Charges.
        charges = molecule.qcModel.AtomicCharges ( molecule )
        ArrayPrint ( charges, itemFormat = "{:.3f}", title = "Charges" )
        logFile.Paragraph ( "Total Charge = {:.3f}".format ( sum ( charges ) ) )
        spins = molecule.qcModel.AtomicSpins ( molecule )
        ArrayPrint ( spins, itemFormat = "{:.3f}", title = "Spin Densities" )
        logFile.Paragraph ( "Total Spin Density = {:.3f}".format ( sum ( spins ) ) )

        # . Gradient testing.
        if _TestGradients and ( len ( molecule.atoms ) < _MaximumAtoms ):
            of = SystemGeometryObjectiveFunction.FromSystem ( molecule )
            gradientDeviation = of.TestGradients ( delta = 1.0e-06, tolerance = 1.0e-03 )
            maximumGradientDeviation = max ( maximumGradientDeviation, gradientDeviation )

    # . Error.
    except Exception as e:
        numberErrors += 1
        logFile.Paragraph ( "Error occurred> " +  e.args[0] )

# . Get the observed and reference data.
observed      = {}
referenceData = TestDataSet.WithOptions ( label = "MNDO CI QC/MM Energies" )
if _TestGradients:
    observed["Gradient Error"] = maximumGradientDeviation
    referenceData.AddDatum ( TestReal.WithOptions ( label = "Gradient Error", value = 0.0, parent = referenceData, absoluteErrorTolerance = _GradientAbsoluteErrorTolerance, toleranceFormat = "{:.3f}", valueFormat = "{:.3f}" ) )

# . Footer.
if len ( observed ) > 0:
    results = referenceData.VerifyAgainst ( observed )
    results.Summary ( fullSummary = _FullVerificationSummary )
    isOK = results.WasSuccessful ( )
else:
    isOK = True
logFile.Footer ( )
if ( not isOK ) or ( numberErrors > 0 ): TestScriptExit_Fail ( )
