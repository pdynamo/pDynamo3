"""Restrained minimization to unfold vacuum structures."""

from Definitions import *

# . Header.
logFile.Header ( )

# . PDB files.
_PDBPaths = ( "1UAO", "2E4E" )

# . Define parameters for the minimization and restraints.
_ForceConstant  = 100.0
_RName          = "rMD"
_TargetDistance =  20.0

# . Loop over structures.
for pdbPath in _PDBPaths:

    # . Retrieve system.
    system = Unpickle ( os.path.join ( outPath, pdbPath + ".pkl" ) )
    system.coordinates3 = ImportCoordinates3 ( os.path.join ( outPath, pdbPath + "_folded.xyz" ) )
    system.Summary ( )
    system.Energy ( )

    # . Get atom indices.
    nTer = system.sequence.AtomIndex ( "A:GLY.1:N"  )
    cTer = system.sequence.AtomIndex ( "A:GLY.10:C" )

    # . Get the starting restraint distance (rounded to the nearest tenth of an Angstrom).
    distance0 = system.coordinates3.Distance ( nTer, cTer )

    # . Set up the restraint.
    restraints = RestraintModel ( )
    system.DefineRestraintModel ( restraints )
    rModel    = RestraintEnergyModel.Harmonic ( _TargetDistance, _ForceConstant )
    restraint = RestraintDistance.WithOptions ( energyModel = rModel, point1 = nTer, point2 = cTer )
    restraints[_RName] = restraint

    # . Optimize.
    system.Energy ( doGradients = True )
    ConjugateGradientMinimize_SystemGeometry ( system,                      \
                                               logFrequency         =  100, \
                                               maximumIterations    = 5000, \
                                               rmsGradientTolerance =  0.5  )

    # . Get the energy and restraint distance.
    system.DefineRestraintModel ( None )
    energy    = system.Energy ( doGradients = True )
    distance1 = system.coordinates3.Distance ( nTer, cTer )

    # . Print the starting and stopping distances.
    print ( "\n\nStarting and stopping distances = {:.2f} {:.2f}.\n".format ( distance0, distance1 ) )

    # . Save structures.
    ExportSystem ( os.path.join ( outPath, pdbPath + "_unfolded.xyz" ), system )

# . Footer.
logFile.Footer ( )
