"""MNDO CI tests."""

import math, os.path

from Definitions        import dataPath            , \
                               _FullVerificationSummary
from pBabel             import ImportSystem
from pCore              import logFile             , \
                               LogFileActive       , \
                               TestDataSet         , \
                               TestReal            , \
                               TestScriptExit_Fail
from pMolecule          import SystemGeometryObjectiveFunction
from pMolecule.QCModel  import CIDiagonalization   , \
                               CIMethod            , \
                               DIISSCFConverger    , \
                               ElectronicState     , \
                               OccupancyType       , \
                               QCModelMNDOCI
from pScientific.Arrays import ArrayPrint

#===================================================================================================================================
# . Parameters for system definitions.
#===================================================================================================================================
# . User-specified sets of micro-states.
_AllylSet1 = ( "1111000011100000", "0111100011100000", "0111010011100000", "0111001011100000" ,
               "0111000111100000", "1011100011100000", "1011010011100000", "1011001011100000" ,
               "1011000111100000", "1101100011100000", "1101010011100000", "1101001011100000" ,
               "1101000111100000", "1110100011100000", "1110010011100000", "1110001011100000" ,
               "1110000111100000", "1111000001110000", "1111000001101000", "1111000001100100" ,
               "1111000001100010", "1111000001100001", "1111000010110000", "1111000010101000" ,
               "1111000010100100", "1111000010100010", "1111000010100001", "1111000011010000" ,
               "1111000011001000", "1111000011000100", "1111000011000010", "1111000011000001" ,
               "1110100001110000", "1110010001110000", "1110001001110000", "1110000101110000" ,
               "1110100010110000", "1110010010110000", "1110001010110000", "1110000110110000" ,
               "1110100011010000", "1110010011010000", "1110001011010000", "1110000111010000" )
_AllylSet2 = ( "110100", "101100", "011100", "110010", "101010", "011010", "110001", "101001", "011001" )

# . Molecule data (general, electronic state, QC model options).
_OptionKeys          = ( "label", "directory", "electronicState", "qcModel" )
_ElectronicStateKeys = ( "charge", "numberFractionalHOOs", "numberFractionalLUOs", "occupancyType" )
_QCModelKeys         = ( "ciMethod", "activeElectrons", "activeOrbitals", "multiplicity", "minimalMultiplicity", "microStates", "requiredRoot" )
_MoleculeData        = ( ( "tyrosineDipeptide", "xyz"      , ( 0, 0, 0, OccupancyType.Cardinal        ) , ( CIMethod.SinglesDoubles , 10, 9, 1, 1, None      , 0 ) ) ,
                         ( "allyl"            , "radicals" , ( 0, 1, 0, OccupancyType.FractionalFixed ) , ( CIMethod.Singles        ,  7, 8, 2, 2, None      , 0 ) ) ,
                         ( "allyl"            , "radicals" , ( 0, 1, 0, OccupancyType.FractionalFixed ) , ( CIMethod.Doubles        ,  7, 8, 2, 2, None      , 0 ) ) ,
                         ( "allyl"            , "radicals" , ( 0, 1, 0, OccupancyType.FractionalFixed ) , ( CIMethod.SinglesDoubles ,  7, 8, 2, 2, None      , 0 ) ) ,
                         ( "allyl"            , "radicals" , ( 0, 1, 0, OccupancyType.FractionalFixed ) , ( CIMethod.User           ,  7, 8, 2, 2, _AllylSet1, 0 ) ) ,
                         ( "allyl"            , "radicals" , ( 0, 1, 0, OccupancyType.FractionalFixed ) , ( CIMethod.User           ,  3, 3, 2, 2, _AllylSet2, 0 ) ) ,
                         ( "methylene"        , "radicals" , ( 0, 1, 1, OccupancyType.FractionalFixed ) , ( CIMethod.Singles        ,  6, 6, 3, 3, None      , 0 ) ) ,
                         ( "methylene"        , "radicals" , ( 0, 1, 1, OccupancyType.FractionalFixed ) , ( CIMethod.Doubles        ,  6, 6, 3, 3, None      , 0 ) ) ,
                         ( "methylene"        , "radicals" , ( 0, 1, 1, OccupancyType.FractionalFixed ) , ( CIMethod.Full           ,  6, 6, 3, 1, None      , 0 ) ) ,
                         ( "methylene"        , "radicals" , ( 0, 1, 1, OccupancyType.FractionalFixed ) , ( CIMethod.SinglesDoubles ,  6, 6, 3, 3, None      , 0 ) ) ,
                         ( "methylene"        , "radicals" , ( 0, 1, 1, OccupancyType.FractionalFixed ) , ( CIMethod.SinglesDoubles ,  6, 6, 5, 3, None      , 0 ) ) ,
                         ( "methylene"        , "radicals" , ( 0, 1, 1, OccupancyType.FractionalFixed ) , ( CIMethod.Full           ,  6, 6, 5, 1, None      , 0 ) ) ,
                         ( "methylene"        , "radicals" , ( 0, 1, 1, OccupancyType.FractionalFixed ) , ( CIMethod.Full           ,  6, 6, 5, 5, None      , 0 ) ) ,
                         ( "water"            , "xyz"      , ( 0, 0, 0, OccupancyType.Cardinal        ) , ( CIMethod.SinglesDoubles ,  8, 6, 1, 1, None      , 0 ) ) )

_AlgorithmOptions    = { "ciDiagonalization"               : CIDiagonalization.Dense ,
                         "checkSparseDiagonalization"      : False ,
                         "numberOfStates"                  : 100   }
#                         "eigenvalueSolverIterations"      : 100000  ,
#                         "eigenvalueSolverPreconditioning" :  False  ,
#                         "eigenvalueSolverTolerance"       : 1.0e-12 }

#===================================================================================================================================
# . Parameters for test.
#===================================================================================================================================
# . Options.
_DoPrinting    = True
_MaximumAtoms  = 100
_TestGradients = True

# . Tolerances.
_EnergyTolerance                = 0.1
_GradientAbsoluteErrorTolerance = 1.0e-02

#===================================================================================================================================
# . Class for a CI test system.
#===================================================================================================================================
class CITestSystem:
    """CI test system."""

    def __init__ ( self, **options ):
        """Constructor."""
        for ( attribute, value ) in options.items ( ):
            setattr ( self, attribute, value )

    def GetSystem ( self, log = logFile ):
        """Get the system with the energy model defined."""
        # . Read the molecule.
        molecule       = ImportSystem  ( os.path.join ( dataPath, self.directory, self.label + ".xyz" ) )
        molecule.label = self.label
        # . Define the electronic state.
        esOptions                = { key : value for ( key, value ) in zip ( _ElectronicStateKeys, self.electronicState ) }
        molecule.electronicState = ElectronicState.WithOptions ( **esOptions )
        # . Define the QC model.
        qcOptions                = { key : value for ( key, value ) in zip ( _QCModelKeys, self.qcModel ) }
        qcOptions["hamiltonian"] = "am1"
        qcOptions["converger"  ] = DIISSCFConverger.WithOptions ( densityTolerance = 1.0e-14, maximumIterations = 250 )
        qcOptions.update ( _AlgorithmOptions )
        qcModel = QCModelMNDOCI.WithOptions ( **qcOptions )
        # . Finish up.
        molecule.DefineQCModel ( qcModel   )
        molecule.Summary       ( log = log )
        # . Finish up.
        return molecule

#===================================================================================================================================
# . Script.
#===================================================================================================================================
# . Header.
logFile.Header ( )

# . Set up the systems.
_CITestSystems = []
for values in _MoleculeData:
    options = { key : value for ( key, value ) in zip ( _OptionKeys, values ) }
    _CITestSystems.append ( CITestSystem ( **options ) )

# . Initialization.
numberErrors = 0
if _TestGradients:
    maximumGradientDeviation = 0.0

# . Loop over systems.
for testSystem in _CITestSystems:

    # . Get the molecule.
    molecule = testSystem.GetSystem ( )

    # . Determine an energy.
    try:
        energy  = molecule.Energy ( doGradients = True )

        # . Charges.
        charges = molecule.qcModel.AtomicCharges ( molecule )
        ArrayPrint ( charges, itemFormat = "{:.3f}", title = "Charges" )
        logFile.Paragraph ( "Total Charge = {:.3f}".format ( sum ( charges ) ) )
        spins = molecule.qcModel.AtomicSpins ( molecule )
        ArrayPrint ( spins, itemFormat = "{:.3f}", title = "Spin Densities" )
        logFile.Paragraph ( "Total Spin Density = {:.3f}".format ( sum ( spins ) ) )

        # . Printing.
        if _DoPrinting:
            molecule.qcModel.CIWavefunctionSummary ( molecule )
            molecule.qcModel.CIVectorsTable ( molecule, numberOfVectors = 20 )

        # . Gradient testing.
        if _TestGradients and ( len ( molecule.atoms ) < _MaximumAtoms ):
            of = SystemGeometryObjectiveFunction.FromSystem ( molecule )
            gradientDeviation = of.TestGradients ( delta = 5.0e-05, tolerance = 3.0e-04 )
            maximumGradientDeviation = max ( maximumGradientDeviation, gradientDeviation )

    # . Error.
    except Exception as e:
        numberErrors += 1
        logFile.Paragraph ( "Error occurred> " +  e.args[0] )

# . Get the observed and reference data.
observed      = {}
referenceData = TestDataSet.WithOptions ( label = "MNDO CI Energies" )
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
