"""Initial minimization and dynamics."""

from Definitions import *

# . Header.
logFile.Header ( )

# . Parameters.
_LogFrequency     = 1000
_NSave            = 100
_NSteps           = 100000
_TimeStep         = 0.001
_TimeBetweenSteps = _NSave * _TimeStep

# . Define the energy models.
mmModel = MMModelOPLS.WithParameterSet ( "bookSmallExamples" )
nbModel = NBModelFull.WithDefaults ( )

# . Generate the system.
system = ImportSystem ( os.path.join ( dataPath, "bala_c7eq.mol" ) )
system.DefineMMModel ( mmModel )
system.DefineNBModel ( nbModel )
system.Summary ( )
system.Energy  ( doGradients = True )

# . Optimization.
ConjugateGradientMinimize_SystemGeometry ( system                      ,
                                           maximumIterations    = 2000 ,
                                           logFrequency         =  100 ,
                                           rmsGradientTolerance =  0.1 )
system.Energy ( doGradients = True )

# . Save minimized system.
Pickle ( os.path.join ( outPath, "bAla.pkl" ), system )

# . Generator options.
randomNumberGenerator  = RandomNumberGenerator.WithSeed ( 511718 )
normalDeviateGenerator = NormalDeviateGenerator.WithRandomNumberGenerator ( randomNumberGenerator )

# . Dynamics.
trajectoryPath = os.path.join ( outPath, "bAla_md.dcd" )
trajectory     = ExportTrajectory ( trajectoryPath, system )
LangevinDynamics_SystemGeometry ( system                                 ,
                                  collisionFrequency     =          25.0 ,
                                  logFrequency           = _LogFrequency ,
                                  normalDeviateGenerator = normalDeviateGenerator ,
                                  steps                  =       _NSteps ,
                                  temperature            =         300.0 ,
                                  timeStep               =     _TimeStep ,
                                  trajectories = [ ( trajectory, _NSave ) ] )

# . Save the coordinates.
ExportSystem ( os.path.join ( outPath, "bAla_md.xyz" ), system )

# . Analysis.
trajectory = ImportTrajectory ( trajectoryPath, system )
trajectory.ReadHeader ( )

# . Loop over the frames in the trajectory.
phi = []
psi = []
while trajectory.RestoreOwnerData ( ):
    phi.append ( system.coordinates3.Dihedral ( *phiAtomIndices ) )
    psi.append ( system.coordinates3.Dihedral ( *psiAtomIndices ) )
times = [ _TimeBetweenSteps * float ( i ) for i in range ( len ( phi ) ) ]

# . Set up the statistics calculation.
phiStatistics = Statistics ( phi )
psiStatistics = Statistics ( psi )

# . Output the results.
table = logFile.GetTable ( columns = [ 20, 20, 20 ] )
table.Start   ( )
table.Title   ( "Phi/Psi Angles" )
table.Heading ( "Time" )
table.Heading ( "Phi"  )
table.Heading ( "Psi"  )
for ( t, h, s ) in zip ( times, phi, psi ):
    table.Entry ( "{:.3f}".format ( t ) )
    table.Entry ( "{:.2f}".format ( h ) )
    table.Entry ( "{:.2f}".format ( s ) )
table.Entry ( "Mean:",               align = Align.Left )
table.Entry ( "{:.2f}".format ( phiStatistics.mean ) )
table.Entry ( "{:.2f}".format ( psiStatistics.mean ) )
table.Entry ( "Standard Deviation:", align = Align.Left )
table.Entry ( "{:.2f}".format ( phiStatistics.standardDeviation ) )
table.Entry ( "{:.2f}".format ( psiStatistics.standardDeviation ) )
table.Stop ( )

# . Footer.
logFile.Footer ( )
