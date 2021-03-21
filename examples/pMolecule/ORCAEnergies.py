"""Basic ORCA tests."""

# . Only limited accuracy can be expected in the gradients due to ORCA, not QC/MM.

import math, os, os.path

from Definitions        import _FullVerificationSummary
from pCore              import logFile                     , \
                               LogFileActive               , \
                               TestScriptExit_Fail         , \
                               TestScriptExit_NotInstalled , \
                               TestDataSet                 , \
                               TestReal
from pMolecule          import SystemGeometryObjectiveFunction
from pMolecule.NBModel  import NBModelORCA
from pMolecule.QCModel  import _ORCACommand                , \
                               QCModelORCA
from pScientific.Arrays import ArrayPrint
from QCMMTestSystems    import qcmmTestSystems

#===================================================================================================================================
# . Check for installation.
#===================================================================================================================================
command = os.getenv ( _ORCACommand )
if  ( command is None ) or not ( os.path.isfile ( command ) and os.access ( command, os.X_OK ) ):
    TestScriptExit_NotInstalled ( )

#===================================================================================================================================
# . Parameters for system definitions.
#===================================================================================================================================
# . QC models.
_qcModels = { "Cyclohexane 9" : [ [ "HF"          , "3-21G" , "SCFCONV10", "EXTREMESCF" , True  ] ] ,
              "Water Dimer 1" : [ [ "B3LYP"       , "6-31G*", "SCFCONV10", "EXTREMESCF" , True  ] ,
                                  [ "HF"          , "3-21G" , "SCFCONV10", "EXTREMESCF" , False ] ,
                                  [ "DFT-ENERGY+" ,           "SCFCONV10", "EXTREMESCF" , False ] ] }

#===================================================================================================================================
# . Parameters for test.
#===================================================================================================================================
# . Options.
_MaximumAtoms  = 100
_TestGradients = False

# . Tolerances.
_EnergyTolerance                = 0.1
_GradientAbsoluteErrorTolerance = 0.1 # . Not very precise.

# . Set up the systems.
_QCMMTestSystems = qcmmTestSystems

# . Define the QC models for each system.
_QCModels = {}
for label in sorted ( _qcModels.keys ( ) ):
    qcModels = []
    for options in _qcModels[label]:
        doQCMM = options.pop ( -1 )
        qcModels.append ( ( doQCMM, QCModelORCA.WithOptions ( keywords = options, randomJob = True ) ) )
    _QCModels[label] = qcModels

#===================================================================================================================================
# . Script.
#===================================================================================================================================
# . Header.
logFile.Header ( )

# . Initialization.
numberErrors = 0
if _TestGradients: maximumGradientDeviation = 0.0

# . Loop over systems and QC models.
for label in sorted ( _QCMMTestSystems.keys ( ) ):
    testSystem = _QCMMTestSystems[label]
    for ( doQCMM, qcModel ) in _QCModels.get ( label, [] ):
        # . Get the molecule.
        molecule = testSystem.GetSystem ( doQCMM = doQCMM, nbModel = NBModelORCA.WithDefaults ( ), qcModel = qcModel )
        # . Energy.
        try:
            energy  = molecule.Energy ( doGradients = True )
            # . Charges.
            charges = molecule.AtomicCharges ( )
            ArrayPrint ( charges, itemFormat = "{:.3f}", title = "Charges" )
            logFile.Paragraph ( "Total Charge = {:.3f}".format ( sum ( charges ) ) )
            spins = molecule.qcModel.AtomicSpins ( molecule )
            if spins is not None:
                ArrayPrint ( spins, itemFormat = "{:.3f}", title = "Spin Densities" )
                logFile.Paragraph ( "Total Spin Density = {:.3f}".format ( sum ( spins ) ) )
            # . Gradient testing.
            if _TestGradients and ( len ( molecule.atoms ) < _MaximumAtoms ):
                of = SystemGeometryObjectiveFunction.FromSystem ( molecule )
                gradientDeviation = of.TestGradients ( delta = 5.0e-04, tolerance = 1.0e-02 )
                maximumGradientDeviation = max ( maximumGradientDeviation, gradientDeviation )
        # . Error.
        except Exception as e:
            numberErrors += 1
            logFile.Paragraph ( "Error occurred> " +  e.args[0] )

# . Get the observed and reference data.
observed      = {}
referenceData = TestDataSet.WithOptions ( label = "ORCA Energies" )
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
