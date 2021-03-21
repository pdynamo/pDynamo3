"""Example 23."""

from Definitions import *

# . Header.
logFile.Header ( )

# . Define some parameters.
_DIncrement    =   1.0
_DMinimum      =   1.5
_DName         = "dOH"
_ForceConstant =  20.0
_Windows       =     5

# . Define the atom indices.
_Oxygen   =  5
_Hydrogen = 17

# . Define the MM and NB models.
mmModel = MMModelOPLS.WithParameterSet ( "bookSmallExamples" )
nbModel = NBModelFull.WithDefaults ( )

# . Generate the molecule.
molecule = ImportSystem ( os.path.join ( molPath, "bala_c7eq.mol" ) )
molecule.DefineMMModel ( mmModel )
molecule.DefineNBModel ( nbModel )
molecule.Summary ( )

# . Read in the starting coordinates.
molecule.coordinates3 = ImportCoordinates3 ( os.path.join ( xyzPath, "bala_1pt5.xyz" ) )

# . Define an empty restraint model and assign it to the system.
restraints = RestraintModel ( )
molecule.DefineRestraintModel ( restraints )

# . Save the molecule definition.
ExportSystem ( os.path.join ( scratchPath, "bala_example23.pkl" ), molecule )

# . Define a random number generator.
randomNumberGenerator  = RandomNumberGenerator ( )
normalDeviateGenerator = NormalDeviateGenerator.WithRandomNumberGenerator ( randomNumberGenerator )

# . Loop over the values of the distance.
for i in range ( _Windows ):

    # . Reset the random number generator.
    randomNumberGenerator.SetSeed ( 291731 + i )

    # . Calculate the new restraint distance.
    distance = _DIncrement * float ( i ) + _DMinimum

    # . Define a new restraint.
    rModel    = RestraintEnergyModel.Harmonic ( distance, _ForceConstant )
    restraint = RestraintDistance.WithOptions ( energyModel = rModel, point1 = _Oxygen, point2 = _Hydrogen )
    restraints[_DName] = restraint

    # . Equilibration.
    LeapFrogDynamics_SystemGeometry ( molecule                        ,
                                      logFrequency           =   1000 ,
                                      normalDeviateGenerator = normalDeviateGenerator ,
                                      steps                  =  50000 ,
                                      temperature            =  300.0 ,
                                      temperatureControl     =   True ,
                                      temperatureCoupling    =    0.1 ,
                                      timeStep               =  0.001 )

    # . Data-collection.
    trajectory = ExportTrajectory ( os.path.join ( scratchPath, "bala_window{:d}.ptRes".format ( i ) ), molecule )
    LeapFrogDynamics_SystemGeometry ( molecule                     ,
                                      logFrequency        =   1000 ,
                                      steps               = 100000 ,
                                      temperature         =  300.0 ,
                                      temperatureControl  =   True ,
                                      temperatureCoupling =    0.1 ,
                                      timeStep            =  0.001 ,
                                      trajectories        = [ ( trajectory, 1 ) ] )

# . Footer.
logFile.Footer ( )
