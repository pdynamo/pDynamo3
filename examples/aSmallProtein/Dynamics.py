"""Molecular dynamics simulation."""

from Definitions import *

# . Header.
logFile.Header ( )

# . PDB files and the number of cations to add.
_PDBPaths = (  "1UAO", "2E4E" )

# . Structures.
_Structures = ( "folded", "unfolded" )

# . Dynamics options.
_NSave  =  250
_NSteps = 1000

# . Generator options.
randomNumberGenerator  = RandomNumberGenerator.WithSeed ( 511717 )
normalDeviateGenerator = NormalDeviateGenerator.WithRandomNumberGenerator ( randomNumberGenerator )
seed                   = 171718

# . Loop over structures.
for pdbPath in _PDBPaths:

    # . Loop over folded and unfolded systems.
    for structure in _Structures:

        # . Retrieve the system.
        system = Unpickle ( os.path.join ( outPath, pdbPath + "_" + structure + "_solvated.pkl" ) )
        system.Summary ( )
        system.Energy  ( )

        # . Put the random number generator in a given state.
        randomNumberGenerator.SetSeed ( seed )
        seed += 1

        # . Do a Langevin dynamics calculation for all atoms.
        trajectory = ExportTrajectory ( os.path.join ( outPath, pdbPath + "_" + structure + "_md1.mdcrd" ), system )
        LangevinDynamics_SystemGeometry ( system                           ,
                                          collisionFrequency     =    25.0 ,
                                          logFrequency           =     100 ,
                                          normalDeviateGenerator = normalDeviateGenerator ,
                                          steps                  = _NSteps ,
                                          temperature            =   300.0 ,
                                          timeStep               =   0.001 ,
                                          trajectories = [ ( trajectory, _NSave ) ] )

        # . Save the coordinates.
        ExportSystem ( os.path.join ( outPath, pdbPath + "_" + structure + "_md1.xyz" ), system )

# . Footer.
logFile.Footer ( )
