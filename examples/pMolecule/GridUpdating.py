"""Grid updating tests."""

import math, os

from Definitions           import dataPath
from pBabel                import CHARMMParameterFileReader  , \
                                  ImportCoordinates3         , \
                                  ImportSystem
from pCore                 import CPUTime                    , \
                                  logFile                    , \
                                  TestScriptExit_Fail
from pMolecule.MMModel     import MMModelOPLS
from pMolecule.NBModel     import NBModelCutOff
from pScientific.Geometry3 import PairListGenerator
from pScientific.Symmetry  import CrystalSystemCubic         , \
                                  PeriodicBoundaryConditions

#===================================================================================================================================
# . Test system.
#===================================================================================================================================
class GridUpdatingTestSystem:
    """Grid updating test system."""

    def __init__ ( self, **options ):
        """Constructor."""
        self.__dict__.update ( options )

    def GetSystem ( self, log = logFile ):
        """Get the system."""
        if self.setUpType == "CHARMM":
            parameters = CHARMMParameterFileReader.PathsToParameters ( self.parameterFiles, log = log )
            system     = ImportSystem ( self.setUpFile, isXPLOR = True, log = log, parameters = parameters )
        else:
            system     = ImportSystem ( self.setUpFile )
        if self.xyzFile is not None: system.coordinates3 = ImportCoordinates3 ( self.xyzFile )
        if self.mmModel is not None: system.DefineMMModel ( self.mmModel )
        if self.hasSymmetry:
            system.symmetry           = PeriodicBoundaryConditions.WithCrystalSystem ( CrystalSystemCubic ( ) )
            system.symmetryParameters = system.symmetry.MakeSymmetryParameters ( a = self.a )
        system.label = self.label
        system.Summary ( log = log )
        return system

#===================================================================================================================================
# . Test systems.
#===================================================================================================================================
# . Test system definitions.
dataPath = os.path.join ( dataPath, "gridUpdating" )
testSystemDefinitions = ( { "hasSymmetry"    : False     ,
                            "label"          : "Crambin" ,
                            "mmModel"        : MMModelOPLS.WithParameterSet ( "protein" ) ,
                            "setUpFile"      : os.path.join ( dataPath, "crambin.mol" ) ,
                            "setUpType"      : "MOL" ,
                            "xyzFile"        : None  },
                          { "a"              : 40.0 ,
                            "hasSymmetry"    : True ,
                            "label"          : "Water Box 40x40x40" ,
                            "mmModel"        : MMModelOPLS.WithParameterSet ( "protein" ) ,
                            "setUpFile"      : os.path.join ( dataPath, "waterBox40x40x40.mol2" ) ,
                            "setUpType"      : "MOL2" ,
                            "xyzFile"        : None   },
                          { "a"              : 62.23 ,
                            "hasSymmetry"    : True  ,
                            "label"          : "DHFR Benchmark" ,
                            "mmModel"        : None ,
                            "parameterFiles" : [ os.path.join ( dataPath, "par_all22_prot.inp" ) ] ,
                            "setUpFile"      : os.path.join ( dataPath, "dhfr.psfx" ) ,
                            "setUpType"      : "CHARMM" ,
                            "xyzFile"        : os.path.join ( dataPath, "dhfr.xyz" ) } )

# . Define the test systems.
testSystems = []
for options in testSystemDefinitions:
    testSystems.append ( GridUpdatingTestSystem ( **options ) ) 

#===================================================================================================================================
# . Script.
#===================================================================================================================================
# . Header.
logFile.Header ( )

# . Set up the generators.
_Values = 11
generators = [ ( "Direct" , PairListGenerator ( cutOff = 13.5, minimumPoints = 100000 ) ) ]
for i in range ( _Values ):
    s = float ( i + 1 ) * 0.25
    generators.append ( ( "CC {:6.3f}".format ( s ), PairListGenerator ( cutOff = 13.5, cutOffCellSizeFactor = s, minimumCellExtent = 0, minimumPoints = 100, sortIndices = False, useGridByCell = True ) ) )
for i in range ( _Values ):
    s = float ( i + 1 ) * 0.25
    generators.append ( ( "PC {:6.3f}".format ( s ), PairListGenerator ( cutOff = 13.5, cutOffCellSizeFactor = s, minimumCellExtent = 0, minimumPoints = 100, sortIndices = False, useGridByCell = False ) ) )

# . Loop over test systems.
numberEnergyFailures = 0
for testSystem in testSystems:

    # . Get the system.
    system = testSystem.GetSystem ( )

    # . Energies.
    energies = []
    for ( label, generator ) in generators:
        nbModel = NBModelCutOff.WithOptions ( generator = generator )
        nbModel.Summary ( )
        system.DefineNBModel ( nbModel )
        cpuTime = CPUTime ( )
        e = system.Energy  ( doGradients = True )
        c = cpuTime.CurrentAsString ( )
        energies.append ( ( e, c, label ) )
        system.nbModel.StatisticsSummary ( system )

    # . Print summary of results.
    e0 = energies[0][0]
    n  = 0
    table = logFile.GetTable ( columns = [ 6, 10, 16, 16, 16, 20 ] )
    table.Start  ( )
    table.Title  ( "Energies and CPU Times for {:s}:".format ( system.label ) )
    for ( i, ( e, c, l ) ) in enumerate ( energies ):
        deltaE = math.fabs ( e - e0 )
        if deltaE > 1.0e-2: n += 1
        table.Entry ( "{:d}".format ( i ) )
        table.Entry (  l  )
        table.Entry ( "{:15.3f}".format ( e      ) )
        table.Entry ( "{:15.3f}".format ( e0     ) )
        table.Entry ( "{:15.3f}".format ( deltaE ) )
        table.Entry ( "    {:s}".format ( c      ) )
    table.Stop ( )
    if n == 0: logFile.Paragraph ( "All tests were successful."      )
    else:      logFile.Paragraph ( "{:d} tests failed.".format ( n ) )
    numberEnergyFailures += n

# . Footer.
logFile.Footer ( )
if numberEnergyFailures != 0: TestScriptExit_Fail ( )
