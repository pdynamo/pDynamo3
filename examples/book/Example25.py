"""Example 25."""

from Definitions import *

# . Header.
logFile.Header ( )

# . Read in the system.
solution = ImportSystem ( os.path.join ( pklPath, "ch4_water215_cubicBox_mc.pkl" ) )
solution.Summary ( )
solution.Energy  ( )

# . Define a random number generator.
randomNumberGenerator = RandomNumberGenerator.WithSeed ( 899311 )

# . Do a Monte Carlo simulation.
trajectory = ExportTrajectory ( os.path.join ( scratchPath, "ch4_water215_cubicBox_mc.ptGeo" ), solution )
MonteCarlo_SystemGeometry ( solution                                        ,
                            blocks                =                     20  ,
                            moves                 =                 100000  ,
                            randomNumberGenerator = randomNumberGenerator   ,
                            trajectories          = [ ( trajectory, 100 ) ] )

# . Footer.
logFile.Footer ( )
