"""Simple ONIOM tests with either a QC(MNDO) or a QC(MNDO)/MM lower layer and a QC(ab initio) upper layer."""

import math, os.path

from Definitions       import _FullVerificationSummary
from pCore             import logFile                     , \
                              LogFileActive               , \
                              Selection                   , \
                              TestDataSet                 , \
                              TestReal                    , \
                              TestScriptExit_Fail         , \
                              TestScriptExit_NotInstalled
from pMolecule         import MultiLayerSystemGeometryObjectiveFunction
from pMolecule.NBModel import NBModelCutOff
from pMolecule.QCModel import DIISSCFConverger            , \
                              ElectronicState             , \
                              QCModelDFT                  , \
                              QCModelMNDO                 , \
                              QCModelORCA
from QCMMTestSystems   import qcmmTestSystems

#===================================================================================================================================
# . Parameters for test.
#===================================================================================================================================
# . Options.
_MaximumAtoms  = 100
_SkipORCATests = True
_TestGradients = True

# . Tolerances.
_EnergyTolerance                = 0.1
_GradientAbsoluteErrorTolerance = 1.0e-02

#===================================================================================================================================
# . Parameters for system definitions.
#===================================================================================================================================
# . Model definitions.
_converger   = DIISSCFConverger.WithOptions ( densityTolerance = 1.0e-12, maximumIterations = 250 )
_qcModelHF   = QCModelDFT.WithOptions       ( converger = _converger, orbitalBasis = "3-21g" )
_qcModelMNDO = QCModelMNDO.WithOptions      ( converger = _converger )
_qcModelORCA = None
_nbModel     = NBModelCutOff.WithDefaults   ( )

# . ORCA tests.
if not _SkipORCATests:
    try:
        _qcModelORCA = QCModelORCA.WithOptions ( keywords = [ "HF", "6-31G*", "SCFCONV10", "EXTREMESCF" ] )
        _qcModelORCA.command # . To see if ORCA is installed.
    except:
        TestScriptExit_NotInstalled ( )

# . Job data.
# . Options: doQCMM for bottom layer, QC model for upper layer and QC selection for upper layer.
_jobData = { "Water Dimer 1"     : [ ( False , _qcModelHF   , Selection.FromIterable ( range (  3     ) ) ) ,
                                     ( False , _qcModelORCA , Selection.FromIterable ( range (  3,  6 ) ) ) ,
                                     ( True  , _qcModelHF   , Selection.FromIterable ( range (  3     ) ) ) ] ,
             "Water Dimer 2"     : [ ( True  , _qcModelORCA , Selection.FromIterable ( range (  3,  6 ) ) ) ] ,
             "Cyclohexane 9"     : [ ( False , _qcModelHF   , Selection.FromIterable ( range (  6     ) ) ) ,
                                     ( True  , _qcModelORCA , Selection.FromIterable ( range (  3,  6 ) ) ) ] ,
             "Alanine Dipeptide" : [ ( False , _qcModelHF   , Selection.FromIterable ( range ( 10, 14 ) ) ) ,
                                     ( False , _qcModelORCA , Selection.FromIterable ( range ( 10, 14 ) ) ) ,
                                     ( True  , _qcModelHF   , Selection.FromIterable ( range ( 10, 14 ) ) ) ,
                                     ( True  , _qcModelORCA , Selection.FromIterable ( range ( 10, 14 ) ) ) ] }

# . Set up the systems.
_Jobs            = _jobData
_QCMMTestSystems = qcmmTestSystems

#===================================================================================================================================
# . Script.
#===================================================================================================================================
# . Header.
logFile.Header ( )

# . Initialization.
numberErrors = 0
if _TestGradients: maximumGradientDeviation = 0.0

# . Loop over systems and jobs.
for label in sorted ( _Jobs.keys ( ) ):
    testSystem = _QCMMTestSystems[label]
    for ( doQCMM, qcModelU, qcSelectionU ) in _Jobs.get ( label, [] ):

        # . Skip ORCA tests.
        if qcModelU is None: continue

        # . Get the molecule.
        molecule = testSystem.GetSystem ( doQCMM = doQCMM, nbModel = _nbModel, qcModel = _qcModelMNDO )

        # . Energy.
        try:
            energy = molecule.Energy ( doGradients = True )

            # . Define the object function.
            of = MultiLayerSystemGeometryObjectiveFunction.FromSystem ( molecule )

            # . First layer.
            of.DefineQCLayer ( qcSelectionU, qcModelU, electronicState = ElectronicState.WithOptions ( charge = 0 ) )
            of.SubsystemSummary ( )

            # . Gradient testing.
            if _TestGradients and ( len ( molecule.atoms ) < _MaximumAtoms ):
                gradientDeviation = of.TestGradients ( delta = 5.0e-04, tolerance = 1.0e-02 )
                maximumGradientDeviation = max ( maximumGradientDeviation, gradientDeviation )

        # . Error.
        except Exception as e:
            numberErrors += 1
            logFile.Paragraph ( "Error occurred> " +  e.args[0] )

# . Get the observed and reference data.
observed      = {}
referenceData = TestDataSet.WithOptions ( label = "ONIOM Energies" )
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
