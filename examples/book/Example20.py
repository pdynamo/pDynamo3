"""Example 20."""

from Definitions import *

# . Header.
logFile.Header ( )

# . Define the box side (in Angstroms).
_BoxSide  = 18.641

# . Define the MM and NB models.
mmModel = MMModelOPLS.WithParameterSet ( "bookSmallExamples" )
nbModel = NBModelCutOff.WithDefaults ( )

# . Generate the solvent.
solvent                    = ImportSystem ( os.path.join ( molPath, "water216_cubicBox.mol" ) )
solvent.symmetry           = PeriodicBoundaryConditions.WithCrystalSystem ( CrystalSystemCubic ( ) )
solvent.symmetryParameters = solvent.symmetry.MakeSymmetryParameters ( a = _BoxSide )
solvent.DefineMMModel ( mmModel )
solvent.DefineNBModel ( nbModel )
solvent.Summary ( )

# . Save the system for later use.
ExportSystem ( os.path.join ( scratchPath, "water216_cubicBox.pkl" ), solvent )

# . Define a normal deviate generator in a given state.
normalDeviateGenerator = NormalDeviateGenerator.WithRandomNumberGenerator ( RandomNumberGenerator.WithSeed ( 491831 ) )

# . Equilibration.
VelocityVerletDynamics_SystemGeometry ( solvent                                ,
                                        logFrequency              =        500 ,
                                        normalDeviateGenerator    = normalDeviateGenerator ,
                                        steps                     =       5000 ,
                                        timeStep                  =      0.001 ,
                                        temperatureScaleFrequency =        100 ,
                                        temperatureScaleOption    = "constant" ,
                                        temperatureStart          =      300.0 )

# . Data-collection.
trajectory = ExportTrajectory ( os.path.join ( scratchPath, "water216_cubicBox.ptGeo" ), solvent )
VelocityVerletDynamics_SystemGeometry ( solvent              ,
                                        logFrequency =   500 ,
                                        steps        = 10000 ,
                                        timeStep     = 0.001 ,
                                        trajectories = [ ( trajectory, 50 ) ] )

# . Footer.
logFile.Footer ( )
