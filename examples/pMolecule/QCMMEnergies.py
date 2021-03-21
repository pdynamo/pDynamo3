"""QC/MM tests."""

import math, os.path

from Definitions        import _FullVerificationSummary                       , \
                               outPath
from pBabel             import ExportSystem
from pCore              import Align                                          , \
                               logFile                                        , \
                               LogFileActive                                  , \
                               TestDataSet                                    , \
                               TestReal                                       , \
                               TestScriptExit_Fail
from pMolecule          import SystemGeometryObjectiveFunction
from pMolecule.NBModel  import NBModelCutOff                                  , \
                               NBModelFull                                    , \
                               QCMMElectrostaticModelDensityCutOffMNDO        , \
                               QCMMElectrostaticModelDensityFullGaussianBasis , \
                               QCMMElectrostaticModelDensityFullMNDO          , \
                               QCMMElectrostaticModelMultipoleCutOff          , \
                               QCMMElectrostaticModelMultipoleFull            , \
                               QCMMLennardJonesModelCutOff                    , \
                               QCMMLennardJonesModelFull
from pMolecule.QCModel  import DIISSCFConverger                               , \
                               LoewdinMultipoleEvaluator                      , \
                               MullikenMultipoleEvaluator                     , \
                               OrthogonalizationType                          , \
                               QCModelDFT                                     , \
                               QCModelMNDO
from pScientific.Arrays import ArrayPrint
from pSimulation        import ConjugateGradientMinimize_SystemGeometry
from QCMMTestSystems    import qcmmTestSystems

#===================================================================================================================================
# . Parameters for test.
#===================================================================================================================================
# . Models.
_DensityModelDFTF  = QCMMElectrostaticModelDensityFullGaussianBasis.WithDefaults ( )
_DensityModelMNDOC = QCMMElectrostaticModelDensityCutOffMNDO.WithDefaults        ( )
_DensityModelMNDOF = QCMMElectrostaticModelDensityFullMNDO.WithDefaults          ( )
_MultipoleModel0F  = QCMMElectrostaticModelMultipoleFull.WithDefaults            ( )
_MultipoleModel2C  = QCMMElectrostaticModelMultipoleCutOff.WithOptions ( multipoleOrder = 2 )
_MultipoleModel2F  = QCMMElectrostaticModelMultipoleFull.WithOptions   ( multipoleOrder = 2 )
_NBModelC          = NBModelCutOff.WithDefaults ( )
_NBModelF          = NBModelFull.WithDefaults   ( )
_QCMMLJModelC      = QCMMLennardJonesModelCutOff.WithDefaults ( )
_QCMMLJModelF      = QCMMLennardJonesModelFull.WithDefaults   ( )
_QCModels          = ( ( "am1/densityFull"       , QCModelMNDO, { "hamiltonian" : "am1" }, _NBModelF, _QCMMLJModelF, _DensityModelMNDOF ) ,
                       ( "am1/multipoleFull"     , QCModelMNDO, { "hamiltonian" : "am1" }, _NBModelF, _QCMMLJModelF, _MultipoleModel2F  ) ,
                       ( "am1/densityCutOff"     , QCModelMNDO, { "hamiltonian" : "am1" }, _NBModelC, _QCMMLJModelC, _DensityModelMNDOC ) ,
                       ( "am1/multipoleCutOff"   , QCModelMNDO, { "hamiltonian" : "am1" }, _NBModelC, _QCMMLJModelC, _MultipoleModel2C  ) ,
                       ( "pm6/densityFull"       , QCModelMNDO, { "hamiltonian" : "pm6" }, _NBModelF, _QCMMLJModelF, _DensityModelMNDOF ) ,
                       ( "pm6/multipoleFull"     , QCModelMNDO, { "hamiltonian" : "pm6" }, _NBModelF, _QCMMLJModelF, _MultipoleModel2F  ) ,
                       ( "pm6/densityCutOff"     , QCModelMNDO, { "hamiltonian" : "pm6" }, _NBModelC, _QCMMLJModelC, _DensityModelMNDOC ) ,
                       ( "pm6/multipoleCutOff"   , QCModelMNDO, { "hamiltonian" : "pm6" }, _NBModelC, _QCMMLJModelC, _MultipoleModel2C  ) ,
                       ( "hf/density"            , QCModelDFT , { "functional" : "HF",                                                    "orbitalBasis" : "svp" }, _NBModelF, _QCMMLJModelF, _DensityModelDFTF ) ,
                       ( "hf/multipole/loewdin"  , QCModelDFT , { "functional" : "HF", "multipoleEvaluator" : LoewdinMultipoleEvaluator , "orbitalBasis" : "svp" }, _NBModelF, _QCMMLJModelF, _MultipoleModel0F ) ,
                       ( "hf/multipole/mulliken" , QCModelDFT , { "functional" : "HF", "multipoleEvaluator" : MullikenMultipoleEvaluator, "orbitalBasis" : "svp" }, _NBModelF, _QCMMLJModelF, _MultipoleModel0F ) )

# . Options.
_DoMinimization = True
_MaximumAtomsG  = 10
_MaximumAtomsM  = 10
_TestGradients  = True

# . Systems.
_QCMMTestSystems = qcmmTestSystems

# . Tolerances.
_EnergyTolerance   = 0.1        # . Unused.
_GradientTolerance = 1.0e-03

# . Paths.
if _DoMinimization:
    _OutPath = os.path.join ( outPath, "QCMMStructures" )
    if not os.path.exists ( _OutPath ): os.mkdir ( _OutPath )

#===================================================================================================================================
# . Script.
#===================================================================================================================================
# . Header.
logFile.Header ( )

# . Initialization.
results      = {}
numberErrors = 0
qcLabels     = sorted ( [ label for ( label, _, _, _, _, _ ) in _QCModels ] )
systemLabels = sorted ( _QCMMTestSystems.keys ( ) )
if _TestGradients: maximumGradientDeviation = 0.0

# . Loop over QC/MM models.
for ( qcLabel, qcModelClass, qcModelOptions, nbModel, qcmmLennardJonesModel, qcmmElectrostaticModel ) in _QCModels:
    if _DoMinimization:
        outPath = os.path.join ( _OutPath, qcLabel.replace ( "/", "_" ) )
        if not os.path.exists ( outPath ): os.mkdir ( outPath )
    # . Loop over systems.
    gradientDeviations = {}
    initialEnergies    = {}
    minimizedEnergies  = {}
    results[qcLabel]   = ( gradientDeviations, initialEnergies, minimizedEnergies )
    for label in systemLabels:
        testSystem = _QCMMTestSystems[label]
        # . Get the molecule.
        qcmmModels                  = { "qcmmElectrostatic" : qcmmElectrostaticModel ,
                                        "qcmmLennardJones"  : qcmmLennardJonesModel  }
        qcModelOptions["converger"] = DIISSCFConverger.WithOptions ( densityTolerance = 1.0e-8, maximumIterations = 250 )
        qcModel                     = qcModelClass.WithOptions     ( **qcModelOptions )
        molecule                    = testSystem.GetSystem ( log = None, nbModel = nbModel, qcModel = qcModel, qcmmModels = qcmmModels )
        molecule.Summary ( )
        # . Energy.
        try:
            initialEnergies[label] = molecule.Energy ( doGradients = True )
            # . Charges.
            charges = molecule.AtomicCharges ( )
            ArrayPrint ( charges, itemFormat = "{:.3f}", title = "Charges" )
            logFile.Paragraph ( "Total Charge = {:.3f}".format ( sum ( charges ) ) )
            spins = molecule.qcModel.AtomicSpins ( molecule )
            if spins is not None:
                ArrayPrint ( spins, itemFormat = "{:.3f}", title = "Spin Densities" )
                logFile.Paragraph ( "Total Spin Density = {:.3f}".format ( sum ( spins ) ) )
            # . Gradient testing.
            if _TestGradients and ( len ( molecule.atoms ) < _MaximumAtomsG ):
                of = SystemGeometryObjectiveFunction.FromSystem ( molecule )
                gradientDeviation         = of.TestGradients ( delta = 1.0e-04, tolerance = _GradientTolerance )
                gradientDeviations[label] = gradientDeviation
                maximumGradientDeviation  = max ( maximumGradientDeviation, gradientDeviation )
            # . Minimization.
            if _DoMinimization and ( len ( molecule.atoms ) < _MaximumAtomsM ):
                ConjugateGradientMinimize_SystemGeometry ( molecule                      ,
                                                           logFrequency         = 1      ,
                                                           maximumIterations    = 10000  ,
                                                           rmsGradientTolerance = 1.0e-2 )
                minimizedEnergies[label] = molecule.Energy ( doGradients = True )
                ExportSystem ( os.path.join ( outPath, testSystem.fileOutName + ".xyz" ), molecule )
        # . Error.
        except Exception as e:
            numberErrors += 1
            logFile.Paragraph ( "Error occurred> " +  e.args[0] )

# . Output.
columns = [ 24, 24, 20 ]
if _TestGradients : columns.append (   20       ) 
if _DoMinimization: columns.extend ( [ 20, 20 ] ) 
table   = logFile.GetTable ( columns = columns )
table.Start   ( )
table.Title   ( "QC Model Energies" )
table.Heading ( "QC Model"        )
table.Heading ( "System"          )
if _DoMinimization:
    table.Heading ( "Initial Energy" )
    table.Heading ( "Final Energy"   )
    table.Heading ( "Difference"     )
else:
    table.Heading ( "Energy" )
if _TestGradients:
    table.Heading ( "Gradient Error" )
for qcLabel in qcLabels:
    ( gradientDeviations, initialEnergies, minimizedEnergies ) = results[qcLabel]
    for ( i, label ) in enumerate ( systemLabels ):
        if i == 0: table.Entry ( qcLabel, align = Align.Left )
        else:      table.Entry ( "" )
        table.Entry ( label, align = Align.Left )
        try   : tag = "{:.1f}".format ( initialEnergies[label] )
        except: tag = "-"
        table.Entry ( tag )
        if _DoMinimization:
            try   : tag = "{:.1f}".format ( minimizedEnergies[label] )
            except: tag = "-"
            table.Entry ( tag )
            try   : tag = "{:.1f}".format ( minimizedEnergies[label] - initialEnergies[label] )
            except: tag = "-"
            table.Entry ( tag )
        if _TestGradients:
            try:
                value = gradientDeviations[label]
                if value > _GradientTolerance: tag = "{:.4f}".format ( value )
                else:                          tag = "ok"
            except: tag = "-"
            table.Entry ( tag )
table.Stop ( )
# . Get the observed and reference data.
observed      = {}
referenceData = TestDataSet.WithOptions ( label = "MNDO QC/MM Energies" )
if _TestGradients:
    observed["Gradient Error"] = maximumGradientDeviation
    referenceData.AddDatum ( TestReal.WithOptions ( label                  = "Gradient Error"   , 
                                                    value                  = 0.0                , 
                                                    parent                 = referenceData      , 
                                                    absoluteErrorTolerance = _GradientTolerance ,
                                                    toleranceFormat        = "{:.4f}"           ,
                                                    valueFormat            = "{:.4f}"           ) )
# . Check for success/failure.
if len ( observed ) > 0:
    results = referenceData.VerifyAgainst ( observed )
    results.Summary ( fullSummary = _FullVerificationSummary )
    isOK    = results.WasSuccessful ( )
else:
    isOK    = True
isOK = isOK and ( numberErrors == 0 )
logFile.Footer ( )
if not isOK: TestScriptExit_Fail ( )
