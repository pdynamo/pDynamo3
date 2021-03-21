"""2-D PMFs - data generation."""

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
seed0             = 215689

# . Get the system.
system              = Unpickle ( os.path.join ( outPath, "bAla.pkl" ) )
system.coordinates3 = ImportCoordinates3 ( os.path.join ( outPath, "bAla_md.xyz" ) )
system.Summary ( )
system.Energy  ( )

# . Define restraints.
restraints = RestraintModel ( )
system.DefineRestraintModel ( restraints )

# . Initial dihedral angles.
phi0 = system.coordinates3.Dihedral ( *phiAtomIndices )
psi0 = system.coordinates3.Dihedral ( *psiAtomIndices )

# . Generator options.
randomNumberGenerator  = RandomNumberGenerator.WithSeed ( seed0 )
normalDeviateGenerator = NormalDeviateGenerator.WithRandomNumberGenerator ( randomNumberGenerator )
seed0 += 1

# . Loop over phi windows.
for i in range ( _NumberOfWindows ):

    # . Phi restraint.
    phi       = phi0 + float ( i ) * _WindowIncrement
    phiRModel = RestraintEnergyModel.Harmonic ( phi, _ForceConstant, period = 360.0 )
    restraint = RestraintDihedral.WithOptions ( energyModel = phiRModel         ,
                                                point1      = phiAtomIndices[0] ,
                                                point2      = phiAtomIndices[1] ,
                                                point3      = phiAtomIndices[2] ,
                                                point4      = phiAtomIndices[3] )
    restraints["Phi"] = restraint

    # . Loop over psi windows.
    for j in range ( _NumberOfWindows ):

        # . Psi restraint.
        psi       = psi0 + float ( j ) * _WindowIncrement
        psiRModel = RestraintEnergyModel.Harmonic ( psi, _ForceConstant, period = 360.0 )
        restraint = RestraintDihedral.WithOptions ( energyModel = psiRModel         ,
                                                    point1      = psiAtomIndices[0] ,
                                                    point2      = psiAtomIndices[1] ,
                                                    point3      = psiAtomIndices[2] ,
                                                    point4      = psiAtomIndices[3] )
        restraints["Psi"] = restraint

        # . Equilibration.
        LangevinDynamics_SystemGeometry ( system                                 ,
                                          collisionFrequency     =          25.0 ,
                                          logFrequency           = _LogFrequency ,
                                          normalDeviateGenerator = normalDeviateGenerator ,
                                          steps                  =      _NSteps0 ,
                                          temperature            =         300.0 ,
                                          timeStep               =     _TimeStep )

        # . Data-collection.
        trajectoryPath = os.path.join ( outPath, "bAla_phi_{:d}_psi_{:d}.ptRes".format ( i, j ) )
        trajectory     = ExportTrajectory ( trajectoryPath, system )
        LangevinDynamics_SystemGeometry ( system                                 ,
                                          collisionFrequency     =          25.0 ,
                                          logFrequency           = _LogFrequency ,
                                          normalDeviateGenerator = normalDeviateGenerator ,
                                          steps                  =      _NSteps1 ,
                                          temperature            =         300.0 ,
                                          timeStep               =     _TimeStep ,
                                          trajectories = [ ( trajectory, 1 ) ] )

# . Footer.
logFile.Footer ( )
