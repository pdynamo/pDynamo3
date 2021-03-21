"""pDynamo system benchmarks."""

import glob, os, os.path

from pBabel               import CHARMMParameterFileReader               , \
                                 ImportCoordinates3                      , \
                                 ImportSystem
from pCore                import Align                                   , \
                                 logFile                                 , \
                                 CPUTime                                 , \
                                 Selection                               , \
                                 TestScript_InputDataPath                , \
                                 YAMLUnpickle
from pMolecule            import EnergyModelPriority                     , \
                                 SystemWithTimings
from pMolecule.NBModel    import NBModelCutOff                           , \
                                 QCMMElectrostaticModelDensityCutOffMNDO , \
                                 QCMMLennardJonesModelCutOff             , \
                                 QCQCElectrostaticModelMultipoleCutOff   , \
                                 QCQCLennardJonesModelCutOff
from pMolecule.QCModel    import ElectronicState                         , \
                                 QCModelMNDO
from pScientific.Symmetry import CrystalSystemCubic                      , \
                                 CrystalSystemOrthorhombic               , \
                                 PeriodicBoundaryConditions
from pSimulation          import LangevinDynamics_SystemGeometry         , \
                                 VelocityVerletDynamics_SystemGeometry

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . File and path names.
_name             = "benchmarks"
_TestDataFileName = "systemData.yaml"
_TestRootPath     = TestScript_InputDataPath ( _name )

# . Other options.
_CrystalSystems       = { "Cubic"        : CrystalSystemCubic        ( ) ,
                          "Orthorhombic" : CrystalSystemOrthorhombic ( ) }
_DoDynamics           = True
_ForceNoQC            = False
_QCRegionKey          = "Large QC Region" # "Small QC Region"
_Steps                = 1000
_UseLangevin          = True
_UseSystemWithTimings = True

#===================================================================================================================================
# . Functions.
#===================================================================================================================================
def FindTests ( ):
    """Find the tests."""
    paths = glob.glob ( os.path.join ( _TestRootPath, "*" ) )
    tests = []
    for path in paths:
        if os.path.exists ( os.path.join ( path, _TestDataFileName ) ): tests.append ( path )
    tests.sort ( )
    return tests

def SetUpSystem ( path, forceNoQC = False, useSystemWithTimings = True ):
    """Set up the system."""
    # . Get system data.
    systemData = YAMLUnpickle ( os.path.join ( path, _TestDataFileName ) )
    # . Get the parameters.
    parameters = CHARMMParameterFileReader.PathsToParameters ( glob.glob ( os.path.join ( path, "*.prm" ) ) )
    # . Get the test name.
    name       = os.path.split ( path )[-1]
    # . Generate the system.
    system              = ImportSystem ( os.path.join ( path, name + ".psfx" ), isXPLOR = True, parameters = parameters )
    if useSystemWithTimings: system = SystemWithTimings.FromSystem ( system )
    system.coordinates3 = ImportCoordinates3 ( os.path.join ( path, name + ".xyz" ) )
    system.label        = systemData["Label"]
    # . Symmetry.
    if systemData.get ( "Crystal Class", None ) is not None:
        system.symmetry           = PeriodicBoundaryConditions.WithCrystalSystem ( _CrystalSystems[systemData["Crystal Class"]] )
        system.symmetryParameters = system.symmetry.MakeSymmetryParameters ( **systemData["Symmetry Parameters"] )
    # . QC model.
    hasQCModel = False
    if not forceNoQC:
        qcData = systemData.get ( _QCRegionKey, None )
        if qcData is not None:
            # . QC atoms.
            qcAtoms = set ( )
            for path in qcData["Atoms"]:
                index = system.sequence.AtomIndex ( path )
                qcAtoms.add ( index )
            # . QC model.
            qcModel                 = QCModelMNDO.WithDefaults ( )
            system.electronicState  = ElectronicState.WithOptions ( charge       = qcData.get ( "Charge"      , 0 ) ,
                                                                    multiplicity = qcData.get ( "Multiplicity", 1 ) )
            system.DefineQCModel ( qcModel, qcSelection = Selection ( qcAtoms ) )
            hasQCModel = True
    # . NB model.
    system.DefineNBModel ( NBModelCutOff.WithDefaults ( ), assignQCMMModels = False )
    # . QC/MM models.
    if hasQCModel:
        _QCMMModels = { "qcmmElectrostatic" : QCMMElectrostaticModelDensityCutOffMNDO.WithDefaults ( ) ,
                        "qcmmLennardJones"  : QCMMLennardJonesModelCutOff.WithDefaults             ( ) ,
                        "qcqcElectrostatic" : QCQCElectrostaticModelMultipoleCutOff.WithDefaults   ( ) ,
                        "qcqcLennardJones"  : QCQCLennardJonesModelCutOff.WithDefaults             ( ) }
        for key in sorted ( _QCMMModels.keys ( ) ):
            system.AddEnergyModel ( key, _QCMMModels[key], priority = EnergyModelPriority.QCMMModel )
    # . Finish up.
    return system

#===================================================================================================================================
# . Execution.
#===================================================================================================================================
# . Header.
logFile.Header ( )

# . Find the tests to run.
tests = FindTests ( )

# . Loop over tests.
results = []
for test in tests:

    # . Get the system.
    system = SetUpSystem ( test, forceNoQC = _ForceNoQC, useSystemWithTimings = _UseSystemWithTimings )
    logFile.Heading ( "Test for \"{:s}\"".format ( system.label ), includeBlankLine = True )
    system.Summary ( )

    # . Do the test.
    cpu = CPUTime ( )
    if _UseSystemWithTimings: system.TimingStart ( )
    system.Energy ( doGradients = True )
    if _DoDynamics:
        if _UseLangevin:
            LangevinDynamics_SystemGeometry       ( system                                 ,
                                                    collisionFrequency        =       25.0 ,
                                                    logFrequency              =        100 ,
                                                    steps                     =     _Steps ,
                                                    temperature               =      300.0 ,
                                                    timeStep                  =      0.001 )
        else:
            VelocityVerletDynamics_SystemGeometry ( system                                 ,
                                                    logFrequency              =        100 ,
                                                    steps                     =     _Steps ,
                                                    temperatureScaleFrequency =        100 ,
                                                    temperatureScaleOption    = "constant" ,
                                                    temperatureStart          =      300.0 ,
                                                    timeStep                  =      0.001 )
    results.append ( ( system.label, cpu.CurrentAsString ( ) ) )
    if _UseSystemWithTimings: system.TimingSummary ( orderByMagnitude = True )
    system.nbModel.StatisticsSummary ( system )

# . Finish up.
table = logFile.GetTable ( columns = [ 20, 30 ] )
table.Start   ( )
table.Title   ( "Test Results" )
table.Heading ( "Test" )
table.Heading ( "Time" )
for ( label, time ) in results:
    table.Entry ( label, align = Align.Left )
    table.Entry ( time )
table.Stop ( )

# . Footer.
logFile.Footer ( )

