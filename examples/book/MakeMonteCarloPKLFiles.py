"""Make the pkl files for Examples 25 and 26."""

from Definitions import *

# . Header.
logFile.Header ( )

# . Molecule definitions.
methane = ImportSystem ( os.path.join ( molPath, "methane.mol" ) )
water   = ImportSystem ( os.path.join ( molPath, "water.mol"   ) )

# . Systems to create.
_ToCreate = ( ( [ methane ] + 215 * [ water ], "Methane in Water", "ch4_water215_cubicBox_mc.xyz", "ch4_water215_cubicBox_mc.pkl" ) ,
              (               216 * [ water ], "Water Box"       , "water216_cubicBox_mc.xyz"    , "water216_cubicBox_mc.pkl"     ) )

for ( molecules, label, xyzIn, pklOut ) in _ToCreate:

    system = MergeByAtom ( molecules )
    system.label = label

    xyzSystem  = ImportSystem ( os.path.join ( xyzPath, xyzIn ) )
    sideLength = float ( xyzSystem.label.split ( )[-1] )

    system.coordinates3 = xyzSystem.coordinates3

    system.DefineMMModel ( MMModelOPLS.WithParameterSet ( "bookSmallExamples" ) )
    system.DefineNBModel ( NBModelMonteCarlo.WithDefaults ( ) )

    system.symmetry           = PeriodicBoundaryConditions.WithCrystalSystem ( CrystalSystemCubic ( ) )
    system.symmetryParameters = system.symmetry.MakeSymmetryParameters ( a = sideLength )

    system.Summary ( )

    ExportSystem ( os.path.join ( pklPath, pklOut ), system )
 
    system.Energy ( )
    system.nbModel.StatisticsSummary ( system )

# . Footer.
logFile.Footer ( )
