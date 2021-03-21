"""Get information about possible solvation scenarios."""

from Definitions import *

# . Header.
logFile.Header ( )

# . PDB files.
_PDBPaths = ( "1UAO", "2E4E" )

# . Structures.
_Structures = ( "folded", "unfolded" )

# . Loop over structures.
for pdbPath in _PDBPaths:

    # . Retrieve the system.
    system = Unpickle ( os.path.join ( outPath, pdbPath + ".pkl" ) )
    system.Summary ( )

    # . Loop over xyz files.
    for structure in _Structures:

        # . Get the coordinates.
        system.coordinates3 = ImportCoordinates3 ( os.path.join ( outPath, pdbPath + "_" + structure + ".xyz" ) )
        system.Energy ( )

        # . Determine solvation parameters.
        DetermineSolvationParameters ( system )

# . Footer.
logFile.Footer ( )
