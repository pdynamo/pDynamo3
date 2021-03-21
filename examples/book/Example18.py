"""Example 18."""

from Definitions import *

# . Header.
logFile.Header ( )

# . Define various parameters.
_ClusterSize     = 3.0
_EnergyTolerance = 0.01
_ForceConstant   = 5.0
_Trails          = 100

# . Set up the system.
cluster = ImportSystem ( os.path.join ( molPath, "argon13.mol" ) )
cluster.DefineMMModel ( MMModelOPLS.WithParameterSet ( "lennardJones" ) )
cluster.DefineNBModel ( NBModelFull.WithDefaults ( ) )
cluster.Summary ( )

# . Define tether restraints for all the atoms.
allAtoms           = Selection.FromIterable  ( range ( len ( cluster.atoms ) ) )
origin             = Coordinates3.WithExtent (         len ( cluster.atoms )   )
origin.Set ( 0.0 )
tetherEnergyModel  = RestraintEnergyModel.FlatBottomedHarmonic ( 0.0, 0.5 * _ClusterSize, _ForceConstant )
tethers            = RestraintModel ( )
tethers["Tethers"] = RestraintMultipleTether.WithOptions ( energyModel = tetherEnergyModel ,
                                                           reference   = origin            , 
                                                           selection   = allAtoms          )

# . Define a random number generator.
randomNumberGenerator  = RandomNumberGenerator ( )
normalDeviateGenerator = NormalDeviateGenerator.WithRandomNumberGenerator ( randomNumberGenerator )

# . Initialize lists to keep energies.
pe0 = []
pe1 = []
pe2 = []
pe3 = []

# . Loop over the trials.
for i in range ( _Trails ):

    # . Reset the cluster coordinates and cluster restraints.
    cluster.coordinates3 = ImportCoordinates3 ( os.path.join ( molPath, "argon13.mol" ), log = None )
    cluster.DefineRestraintModel ( tethers )

    # . Initialize some variables for the trial.
    randomNumberGenerator.SetSeed ( 957612 + i )
    tStart = 300.0 * ( randomNumberGenerator.NextReal ( ) + 1.0 )

    # . Do a short dynamics to generate a new structure.
    VelocityVerletDynamics_SystemGeometry ( cluster                                ,
                                            log                       =       None ,
                                            normalDeviateGenerator    = normalDeviateGenerator ,
                                            steps                     =      10000 ,
                                            timeStep                  =      0.001 ,
                                            temperatureScaleFrequency =        100 ,
                                            temperatureScaleOption    = "constant" ,
                                            temperatureStart          =     tStart )

    # . Save the starting coordinates and energy.
    temporary3 = Clone ( cluster.coordinates3 )
    cluster.DefineRestraintModel ( None )
    pe0.append ( cluster.Energy ( log = None ) )

    # . Minimization.
    cluster.DefineRestraintModel ( tethers )
    ConjugateGradientMinimize_SystemGeometry ( cluster                       ,
                                               log                  =  None  ,
                                               maximumIterations    = 10000  ,
                                               rmsGradientTolerance = 1.0e-4 )
    cluster.DefineRestraintModel ( None )
    ConjugateGradientMinimize_SystemGeometry ( cluster                       ,
                                               log                  =  None  ,
                                               maximumIterations    = 10000  ,
                                               rmsGradientTolerance = 1.0e-4 )
    pe1.append ( cluster.Energy ( log = None ) )

    # . Simulated annealing from the starting coordinates.
    cluster.coordinates3 = temporary3
    cluster.DefineRestraintModel ( tethers )
    VelocityVerletDynamics_SystemGeometry ( cluster                                   ,
                                            log                       =          None ,
                                            steps                     =         40000 ,
                                            timeStep                  =         0.001 ,
                                            temperatureScaleFrequency =            10 ,
                                            temperatureScaleOption    = "exponential" ,
                                            temperatureStart          = tStart        ,
                                            temperatureStop           = tStart * math.exp ( - 10.0 ) )
    cluster.DefineRestraintModel ( None )
    pe2.append ( cluster.Energy ( log = None ) )

    # . Minimization of the annealed structure.
    ConjugateGradientMinimize_SystemGeometry ( cluster                       ,
                                               log                  =  None  ,
                                               maximumIterations    = 10000  ,
                                               rmsGradientTolerance = 1.0e-4 )
    pe3.append ( cluster.Energy ( log = None ) )

# . Prepare the energies for statistics.
stpe1 = Statistics ( pe1 )
stpe2 = Statistics ( pe2 )
stpe3 = Statistics ( pe3 )

# . Output the results.
table = logFile.GetTable ( columns = [ 10, 20, 20, 20, 20 ] )
table.Start   ( )
table.Title   ( "Optimization Results" )
table.Heading ( "Attempt"          )
table.Heading ( "Initial Energy"   )
table.Heading ( "Minimized Energy" )
table.Heading ( "Annealed Energy"  )
table.Heading ( "Final Energy"     )
for i in range ( _Trails ):
    table.Entry ( "{:d}"    .format (     i  ) )
    table.Entry ( "{:20.3f}".format ( pe0[i] ) )
    table.Entry ( "{:20.3f}".format ( pe1[i] ) )
    table.Entry ( "{:20.3f}".format ( pe2[i] ) )
    table.Entry ( "{:20.3f}".format ( pe3[i] ) )
table.Entry ( "Minimum Energies:", align = Align.Left, columnSpan = 2 )
table.Entry ( "{:20.3f}".format ( stpe1.minimum ) )
table.Entry ( "{:20.3f}".format ( stpe2.minimum ) )
table.Entry ( "{:20.3f}".format ( stpe3.minimum ) )
table.Entry ( "Frequencies:",      align = Align.Left, columnSpan = 2 )
table.Entry ( "{:d}".format ( stpe1.Count ( stpe1.minimum, tolerance = _EnergyTolerance ) ) )
table.Entry ( "{:d}".format ( stpe2.Count ( stpe2.minimum, tolerance = _EnergyTolerance ) ) )
table.Entry ( "{:d}".format ( stpe3.Count ( stpe3.minimum, tolerance = _EnergyTolerance ) ) )
table.Stop ( )

# . Footer.
logFile.Footer ( )
