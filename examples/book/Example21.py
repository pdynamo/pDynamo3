"""Example 21."""

from Definitions import *

# . Header.
logFile.Header ( )

# . Read the system definition.
solvent = ImportSystem ( os.path.join ( scratchPath, "water216_cubicBox.pkl" ) )
solvent.Summary ( )

# . Select all oxygens.
indices = []
for ( i, atom ) in enumerate ( solvent.atoms ):
    if atom.atomicNumber == 8: indices.append ( i )
oxygens = Selection.FromIterable ( indices )

# . Analyse the trajectory data.
# . Radial distribution and self-diffusion functions.
trajectoryPath = os.path.join ( scratchPath, "water216_cubicBox.ptGeo" )
RadialDistributionFunction ( trajectoryPath, solvent, selection1 = oxygens )
SelfDiffusionFunction      ( trajectoryPath, solvent, selection  = oxygens )

# . Footer.
logFile.Footer ( )
