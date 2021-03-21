"""Example 16."""

from Definitions import *

# . Header.
logFile.Header ( )

# . Define the energy models.
mmModel = MMModelOPLS.WithParameterSet ( "bookSmallExamples" )
nbModel = NBModelFull.WithDefaults ( )

# . Generate the molecule.
molecule = ImportSystem ( os.path.join ( molPath, "bala_c7eq.mol" ) )
molecule.DefineMMModel ( mmModel )
molecule.DefineNBModel ( nbModel )
molecule.Summary ( )
molecule.Energy  ( )

# . Optimization.
ConjugateGradientMinimize_SystemGeometry ( molecule                    ,
                                           maximumIterations    = 2000 ,
                                           logFrequency         =  100 ,
                                           rmsGradientTolerance =  0.1 )

# . Define a random number generator in a given state.
normalDeviateGenerator = NormalDeviateGenerator.WithRandomNumberGenerator ( RandomNumberGenerator.WithSeed ( 175189 ) )

# . Heating.
VelocityVerletDynamics_SystemGeometry ( molecule                             ,
                                        logFrequency              =      100 ,
                                        normalDeviateGenerator    = normalDeviateGenerator ,
                                        steps                     =     1000 ,
                                        timeStep                  =    0.001 ,
                                        temperatureScaleFrequency =      100 ,
                                        temperatureScaleOption    = "linear" ,
                                        temperatureStart          =     10.0 ,
                                        temperatureStop           =    300.0 )

# . Equilibration.
VelocityVerletDynamics_SystemGeometry ( molecule                               ,
                                        logFrequency              =        500 ,
                                        steps                     =       5000 ,
                                        timeStep                  =      0.001 ,
                                        temperatureScaleFrequency =        100 ,
                                        temperatureScaleOption    = "constant" ,
                                        temperatureStart          =      300.0 )

# . Data-collection.
trajectory = ExportTrajectory ( os.path.join ( scratchPath, "bala_c7eq.ptGeo" ), molecule )
VelocityVerletDynamics_SystemGeometry ( molecule             ,
                                        logFrequency =   500 ,
                                        steps        = 10000 ,
                                        timeStep     = 0.001 ,
                                        trajectories = [ ( trajectory, 100 ) ] )

# . Footer.
logFile.Footer ( )
