"""Functions for testing of QC systems."""

import glob, math, os.path

from Definitions            import _FullVerificationSummary        , \
                                   yamlPath as _ReferenceDataPath
from pBabel                 import ImportSystem
from pCore                  import Align                           , \
                                   Clone                           , \
                                   logFile                         , \
                                   LogFileActive                   , \
                                   Selection                       , \
                                   TestDataSet                     , \
                                   TestReal                        , \
                                   TestScript_InputDataPath        , \
                                   TestScriptExit_Fail             , \
                                   YAMLMappingFile_ToObject
from pMolecule              import SystemGeometryObjectiveFunction
from pMolecule.QCModel      import DIISSCFConverger                , \
                                   ElectronicState                 , \
                                   OccupancyType
from pScientific.Arrays     import ArrayPrint
from pScientific.Statistics import Statistics

#===================================================================================================================================
# . Basic parameters.
#===================================================================================================================================
# . Parameters.
_TableDataWidth  = 20
_TableModelWidth = 45

# . Paths.
_dataPath       = TestScript_InputDataPath ( "pMolecule" )
_structuresPath = os.path.join ( os.getenv ( "PDYNAMO3_HOME" ), "structures" )

# . Tolerances.
_AbsoluteErrorTolerance         = 1.0
_BondOrderTolerance             = 0.1
_EnergyTolerance                = 0.1
_GradientAbsoluteErrorTolerance = 1.0e-03

#===================================================================================================================================
# . Class for a set of QC tests.
#===================================================================================================================================
class QCTestSystem:
    """QC test system."""

    def __init__ ( self, **options ):
        """Constructor."""
        for ( attribute, value ) in options.items ( ): setattr ( self, attribute, value )

    def GetSystem ( self, log = logFile, maximumAtoms = None ):
        """Get the system with the energy model defined."""
        # . Initialization.
        molecule  = None
        tooBig    = True
        undefined = True
        # . Get the QC model options.
        electronicStateKeywords = { "charge"           : 0                      ,
                                    "isSpinRestricted" : True                   ,
                                    "multiplicity"     : 1                      ,
                                    "occupancyType"    : OccupancyType.Cardinal }
        convergerKeywords       = getattr ( self, "convergerKeywords"       ,   {} )
        qcModelClass            = getattr ( self, "qcModelClass"            , None )
        qcModelKeywords         = getattr ( self, "qcModelKeywords"         ,   {} )
        electronicStateKeywords.update ( getattr ( self, "electronicStateKeywords" , {} ) )
        # . Basic setup.
        if   self.fileFormat == "mol": path = os.path.join ( self.dataPath, self.fileName + ".mol" )
        elif self.fileFormat == "xyz": path = os.path.join ( self.dataPath, self.fileName + ".xyz" )
        molecule        = ImportSystem  ( path )
        electronicState = ElectronicState.WithOptions ( **electronicStateKeywords )
        molecule.label  = self.label
        # . Only keep the molecule if it is not too large.
        if ( maximumAtoms is None ) or ( ( maximumAtoms is not None ) and ( len ( molecule.atoms ) <= maximumAtoms ) ):
            tooBig = False
            # . Define the QC model.
            if qcModelClass is not None:
                converger            = DIISSCFConverger.WithOptions ( **convergerKeywords )
                options              = dict ( qcModelKeywords )
                options["converger"] = converger
                try:
                    qcModel                  = qcModelClass.WithOptions ( **options )
                    molecule.electronicState = electronicState
                    molecule.DefineQCModel ( qcModel )
                    undefined = False
                except Exception as e:
                    print ( e.args[0] )
                    pass
            # . Summary.
            if LogFileActive ( log ):
                molecule.Summary ( log = log )
                log.Paragraph ( "Formula = " + molecule.atoms.FormulaString ( ) + "." )
        # . Finish up.
        return ( molecule, undefined, tooBig )

#===================================================================================================================================
# . Helper function.
#===================================================================================================================================
def _TestSystemsFromData ( data, convergerKeywords      = {}   ,
                                 electronicStateOptions = {}   ,
                                 qcModelClass           = None ,
                                 qcModelOptions         = {}   ):
    """Gather test systems from data."""
    # . Process options.
    numberZero   = 0
    if len ( electronicStateOptions ) == 0:
        esOptions   = { "" : {} }
        numberZero += 1
    else:
        esOptions   = electronicStateOptions
    if len ( qcModelOptions         ) == 0:
        qcOptions   = { "" : {} }
        numberZero += 1
    else:
        qcOptions   = qcModelOptions
    if numberZero == 2: raise QCModelError ( "There are no test system models." )
    # . Create models.
    modelOptions = {}
    for ( qcLabel, qcKeywords ) in qcOptions.items ( ):
        for ( esLabel, esKeywords ) in esOptions.items ( ):
            if   len ( esLabel ) == 0: label = qcLabel
            elif len ( qcLabel ) == 0: label = esLabel
            else:                      label = qcLabel + " " + esLabel
            modelOptions[label] = ( esKeywords, qcKeywords )
    # . Create systems.
    testSystems  = []
    for values in [ data[label] for label in sorted ( data.keys ( ) ) ]:
        for modelLabel in sorted ( modelOptions.keys ( ) ):
            ( esKeywords, qcKeywords )   = modelOptions[modelLabel]
            options                      = Clone ( values )
            options["convergerKeywords"] = convergerKeywords
            options["modelLabel"       ] = modelLabel
            options["qcModelClass"     ] = qcModelClass
            if len ( esKeywords ) > 0:
                keywords = options.get ( "electronicStateKeywords", {} )
                keywords.update ( esKeywords )
                options["electronicStateKeywords"] = keywords
            if len ( qcKeywords ) > 0:
                keywords = options.get ( "qcModelKeywords", {} )
                keywords.update ( qcKeywords )
                options["qcModelKeywords"] = keywords
            testSystems.append ( QCTestSystem ( **options ) )
    return testSystems

#===================================================================================================================================
# . Basic closed shell tests.
#===================================================================================================================================
# . Charges.
_closedShellCharges = { "cadmiumComplex"  :  2 ,
                        "chloride"        : -1 ,
                        "fch3cl"          : -1 ,
                        "histamine"       :  1 ,
                        "methylPhosphate" : -2 ,
                        "rhodiumComplex"  : -2 }

# . Sources.
_closedShellSources = ( os.path.join ( _dataPath                                      , "xyz" ) ,
                        os.path.join ( _structuresPath, "difficultSCFCases"           , "xyz" ) ,
                        os.path.join ( _structuresPath, "gaussianGeometryOptimization", "xyz" ) )

# . Return all systems.
def GetClosedShellMoleculeSystems ( convergerKeywords      = {}   ,
                                    electronicStateOptions = {}   ,
                                    qcModelClass           = None ,
                                    qcModelOptions         = {}   ):
    """Get all closed shell systems."""
    # . Get paths.
    paths = []
    for source in _closedShellSources:
        if os.path.exists ( source ):
            paths.extend ( glob.glob ( os.path.join ( source, "*.xyz" ) ) )
    # . Gather data.
    data = {}
    for path in paths:
        ( head, tail ) = os.path.split ( path )
        label    = tail[0:-4]
        fileName = label
        data[label] = { "label"      : label ,
                        "fileName"   : label ,
                        "fileFormat" : "xyz" ,
                        "dataPath"   : head  ,
                        "electronicStateKeywords" : { "charge" : _closedShellCharges.get ( label, 0 ) } }
    # . Get the systems.
    return _TestSystemsFromData ( data, convergerKeywords      = convergerKeywords      ,
                                        electronicStateOptions = electronicStateOptions ,
                                        qcModelClass           = qcModelClass           ,
                                        qcModelOptions         = qcModelOptions         )

#===================================================================================================================================
# . Basic radical tests (includes closed and open shell cases).
#===================================================================================================================================
# . Systems.
_radicalKeywordLabels = ( "fileName", "fileFormat", "dataPath", "electronicStateKeywords" )
_radicalSystems       = { "Tyrosine Dipeptide Singlet" : ( "tyrosineDipeptide", "xyz", os.path.join ( _dataPath, "xyz"      ), { "multiplicity" : 1 } ) ,
                          "Tyrosine Dipeptide Triplet" : ( "tyrosineDipeptide", "xyz", os.path.join ( _dataPath, "xyz"      ), { "multiplicity" : 3 } ) ,
                          "Allyl Radical"              : ( "allyl"            , "xyz", os.path.join ( _dataPath, "radicals" ), { "multiplicity" : 2 } ) ,
                          "Methylene Radical"          : ( "methylene"        , "xyz", os.path.join ( _dataPath, "radicals" ), { "multiplicity" : 3 } ) }

# . Extra model options.
_restrictedOptions   = { "Cardinal Restricted"              : { "isSpinRestricted" : True , "occupancyType" : OccupancyType.Cardinal           } ,
                         "Fractional Fixed Restricted"      : { "isSpinRestricted" : True , "occupancyType" : OccupancyType.FractionalFixed    } ,
                         "Fractional Variable Restricted"   : { "isSpinRestricted" : True , "occupancyType" : OccupancyType.FractionalVariable } }
_unRestrictedOptions = { "Cardinal Unrestricted"            : { "isSpinRestricted" : False, "occupancyType" : OccupancyType.Cardinal           } ,
                         "Fractional Fixed Unrestricted"    : { "isSpinRestricted" : False, "occupancyType" : OccupancyType.FractionalFixed    } ,
                         "Fractional Variable Unrestricted" : { "isSpinRestricted" : False, "occupancyType" : OccupancyType.FractionalVariable } }

# . Return all systems.
def GetRadicalMoleculeSystems ( convergerKeywords = {}   ,
                                qcModelClass      = None ,
                                qcModelOptions    = {}   ):
    """Get all radical systems."""
    # . Gather data - R and U separately.
    dataR = {}
    dataU = {}
    for label in sorted ( _radicalSystems.keys ( ) ):
        basicOptions          = { key : value for ( key, value ) in zip ( _radicalKeywordLabels, _radicalSystems[label] ) }
        basicOptions["label"] = label
        dataU[label]          = basicOptions
        if basicOptions["electronicStateKeywords"]["multiplicity"] == 1: dataR[label] = dict ( basicOptions )
    # . Get the systems.
    testSystemsR = _TestSystemsFromData ( dataR, convergerKeywords      = convergerKeywords    ,
                                                 electronicStateOptions = _restrictedOptions   ,
                                                 qcModelClass           = qcModelClass         ,
                                                 qcModelOptions         = qcModelOptions       )
    testSystemsU = _TestSystemsFromData ( dataU, convergerKeywords      = convergerKeywords    ,
                                                 electronicStateOptions = _unRestrictedOptions ,
                                                 qcModelClass           = qcModelClass         ,
                                                 qcModelOptions         = qcModelOptions       )
    return ( testSystemsR + testSystemsU )

#===================================================================================================================================
# . Functions to run tests.
#===================================================================================================================================
def RunQCTestSet ( testSystems                    ,
                   dataSetTag                     ,
                   maximumEnergyAtoms      = 30   ,
                   maximumEnergyTests      = 1000 ,
                   maximumGradientAtoms    = 10   ,
                   maximumGradientTests    = 20   ,
                   referenceDataFileName   = None ,
                   testGradients           = True ):
    """Run a set of QC tests."""
    # . Header.
    logFile.Header ( )
    # . Initialization.
    modelResults        = {}
    numberEnergyTests   = 0
    numberErrors        = 0
    numberGradientTests = 0
    if testGradients:
        maximumGradientDeviation = 0.0
    # . Energy data.
    energyResults = {}
    # . Loop over systems.
    for testSystem in testSystems:
        # . Get results (if necessary).
        modelLabel = testSystem.modelLabel
        energies   = energyResults.get ( modelLabel, None )
        if energies is None:
            energies = {}
            energyResults[modelLabel] = energies
        results = modelResults.get ( modelLabel, None )
        if results is None:
            results = {}
            modelResults[modelLabel] = results
        cycles  = results.get ( "Cycles", None )
        if cycles is None:
            cycles = []
            results["Cycles"] = cycles
        # . Get the molecule.
        ( molecule, undefined, tooBig ) = testSystem.GetSystem ( maximumAtoms = maximumEnergyAtoms )
        if tooBig: continue
        if undefined:
            results["Undefined"] = results.get ( "Undefined", 0 ) + 1
            continue
        # . Determine an energy.
        isConverged = False
        try:
            energy               = molecule.Energy ( doGradients = True )
            isConverged          = True
            results["Converged"] = results.get ( "Converged", 0 ) + 1
            cycles.append ( float ( molecule.scratch.qcEnergyReport["SCF Iterations"] ) )
        # . Error.
        except Exception as error:
            isNotConverged = False
            message        = None
            if ( len ( error.args ) > 0 ) and isinstance ( error.args[0], str ):
                message        = error.args[0]
                isNotConverged = message.endswith ( "Converger error: Too many iterations." )
            if isNotConverged:
                results["Not Converged"] = results.get ( "Not Converged", 0 ) + 1
            else:
                numberErrors            += 1
                results["Failed"]        = results.get ( "Failed"       , 0 ) + 1
            if message is None: logFile.Paragraph ( "** Unspecified error in energy calculation. **" )
            else:               logFile.Paragraph ( "** Error in energy calculation: " + message + " **" )
        # . Other properties.
        if isConverged:
            # . Accumulate observed data.
            energies[molecule.label] = energy
            numberOfAtoms            = len ( molecule.atoms )
            # . Bond orders.
            labels = []
            for i in range ( numberOfAtoms ): labels.append ( molecule.atoms[i].path )
            ( bondOrders, _, _ ) = molecule.qcModel.BondOrders ( molecule )
            ArrayPrint ( bondOrders                             ,
                         conditionFunction = ( lambda x, i, j: ( i != j ) and ( math.fabs ( x ) >= _BondOrderTolerance ) ) ,
                         itemFormat  = "{:.3f}"           ,
                         itemsPerRow = 5                  ,
                         itemWidth   = 12                 ,
                         labels      = [ labels, labels ] ,
                         title       = "Bond Orders"      )
            # . Charges.
            charges = molecule.AtomicCharges ( )
            ArrayPrint ( charges, itemFormat = "{:.3f}", itemWidth = 12, title = "Charges" )
            logFile.Paragraph ( "Total Charge = {:.3f}".format ( sum ( charges ) ) )
            spins = molecule.qcModel.AtomicSpins ( molecule )
            if spins is not None:
                ArrayPrint ( spins, itemFormat = "{:.3f}", itemWidth = 12, title = "Spin Densities" )
                logFile.Paragraph ( "Total Spin Density = {:.3f}".format ( sum ( spins ) ) )
            # . Dipole moment.
            dipole = molecule.DipoleMoment ( )
            ArrayPrint ( dipole, itemFormat = "{:.3f}", itemWidth = 12, title = "Dipole", useLabelRange = True )
            # . Gradient testing.
            if testGradients and ( numberOfAtoms < maximumGradientAtoms ) and ( numberGradientTests < maximumGradientTests ):
                of = SystemGeometryObjectiveFunction.FromSystem ( molecule )
                gradientDeviation = of.TestGradients ( )
                maximumGradientDeviation = max ( maximumGradientDeviation, gradientDeviation )
                numberGradientTests += 1
        # . Check the number of tests.
        numberEnergyTests += 1
        if numberEnergyTests >= maximumEnergyTests: break
    # . Print model results.
    models = list ( modelResults.keys ( ) )
    models.sort ( )
    table = logFile.GetTable ( columns = [ _TableModelWidth ] + 7 * [ _TableDataWidth ] )
    table.Start   ( )
    table.Title   ( "Model Results" )
    table.Heading ( "Model"         )
    table.Heading ( "Undefined"     )
    table.Heading ( "Failed"        )
    table.Heading ( "Not Converged" )
    table.Heading ( "Converged"     )
    table.Heading ( "<Cycles>"      )
    table.Heading ( "Max. Cycles"   )
    table.Heading ( "Min. Cycles"   )
    for model in models:
        results = modelResults[model]
        rCycles = results.get ( "Cycles", None )
        if ( rCycles is None ) or ( len ( rCycles ) <= 0 ): cycles = None
        else:                                               cycles = Statistics ( results["Cycles"] )
        table.Entry ( model, align = Align.Left )
        table.Entry ( "{:d}".format ( results.get ( "Undefined"    , 0 ) ) )
        table.Entry ( "{:d}".format ( results.get ( "Failed"       , 0 ) ) )
        table.Entry ( "{:d}".format ( results.get ( "Not Converged", 0 ) ) )
        table.Entry ( "{:d}".format ( results.get ( "Converged"    , 0 ) ) )
        if ( cycles is not None ) and ( cycles.size > 0 ):
            table.Entry ( "{:.1f}".format ( cycles.mean    ) )
            table.Entry ( "{:.0f}".format ( cycles.maximum ) )
            table.Entry ( "{:.0f}".format ( cycles.minimum ) )
        else:
            table.Entry ( "-" )
            table.Entry ( "-" )
            table.Entry ( "-" )
    table.Stop ( )
    # . Energy data.
    # . Verify the observed data against the reference data.
    isOK = True
    if referenceDataFileName is not None:
        referenceDataPath = os.path.join ( _ReferenceDataPath, referenceDataFileName )
        referenceData     = YAMLMappingFile_ToObject ( referenceDataPath, TestDataSet )
        results           = referenceData.VerifyAgainst ( energyResults )
        isOK              = results.WasSuccessful ( )
        results.Summary ( fullSummary = _FullVerificationSummary )
        # . Gradient data set.
        if testGradients:
            # . Get the observed and reference data.
            observed      = {}
            referenceData = TestDataSet.WithOptions ( label = dataSetTag + " Gradients" )
            observed["Gradient Error"] = maximumGradientDeviation
            referenceData.AddDatum ( TestReal.WithOptions ( label = "Gradient Error", value = 0.0, parent = referenceData, absoluteErrorTolerance = _GradientAbsoluteErrorTolerance, toleranceFormat = "{:.3f}", valueFormat = "{:.3f}" ) )
            # . Check for success/failure.
            if len ( observed ) > 0:
                results = referenceData.VerifyAgainst ( observed )
                results.Summary ( fullSummary = _FullVerificationSummary )
                isOK    = isOK and results.WasSuccessful ( )
    # . Footer.
    isOK = isOK and ( numberErrors == 0 )
    logFile.Footer ( )
    if not isOK: TestScriptExit_Fail ( )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass
