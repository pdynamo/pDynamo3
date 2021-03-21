"""1-D PMFs - data generation."""

from Definitions import *

# . Header.
logFile.Header ( )

# . Dynamics parameters.
_LogFrequency     = 1000
_NSteps0          = 10000
_NSteps1          = 10000
_TimeStep         = 0.001

# . Window parameters.
_ForceConstant    = 0.1
_NumberOfWindows  = 36
_WindowIncrement  = 360.0 / float ( _NumberOfWindows )
seed0             = 511719

# . Loop over angles.
for ( tag, indices ) in ( ( "Phi", phiAtomIndices ), ( "Psi", psiAtomIndices ) ):

    # . Get the system.
    system              = Unpickle ( os.path.join ( outPath, "bAla.pkl" ) )
    system.coordinates3 = ImportCoordinates3 ( os.path.join ( outPath, "bAla_md.xyz" ) )
    system.Summary ( )
    system.Energy  ( )

    # . Define restraints.
    restraints = RestraintModel ( )
    system.DefineRestraintModel ( restraints )

    # . Initial dihedral angle.
    angle0 = system.coordinates3.Dihedral ( *indices )

    # . Generator options.
    randomNumberGenerator  = RandomNumberGenerator.WithSeed ( seed0 )
    normalDeviateGenerator = NormalDeviateGenerator.WithRandomNumberGenerator ( randomNumberGenerator )
    seed0 += 1

    # . Loop over the windows.
    for i in range ( _NumberOfWindows ):

        # . Redefine the restraint.
        angle     = angle0 + float ( i ) * _WindowIncrement
        rModel    = RestraintEnergyModel.Harmonic ( angle, _ForceConstant, period = 360.0 )
        restraint = RestraintDihedral.WithOptions ( energyModel = rModel     ,
                                                    point1      = indices[0] ,
                                                    point2      = indices[1] ,
                                                    point3      = indices[2] ,
                                                    point4      = indices[3] )
        restraints[tag] = restraint

        # . Equilibration.
        LangevinDynamics_SystemGeometry ( system                                 ,
                                          collisionFrequency     =          25.0 ,
                                          logFrequency           = _LogFrequency ,
                                          normalDeviateGenerator = normalDeviateGenerator ,
                                          steps                  =      _NSteps0 ,
                                          temperature            =         300.0 ,
                                          timeStep               =     _TimeStep )

        # . Data-collection.
        trajectoryPath = os.path.join ( outPath, "bAla_{:s}_{:d}.ptRes".format ( tag.lower ( ), i ) )
        trajectory     = ExportTrajectory ( trajectoryPath, system )
        LangevinDynamics_SystemGeometry ( system                                 ,
                                          collisionFrequency     =          25.0 ,
                                          logFrequency           = _LogFrequency ,
                                          normalDeviateGenerator = normalDeviateGenerator ,
                                          steps                  =      _NSteps1 ,
                                          temperature            =         300.0 ,
                                          timeStep               =     _TimeStep ,
                                          trajectories = [ ( trajectory, 1 ) ]   )

# . Footer.
logFile.Footer ( )
