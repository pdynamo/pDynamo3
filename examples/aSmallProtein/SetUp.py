"""Setup the proteins."""

from Definitions import *

# . Header.
logFile.Header ( )

# . PDB files.
_PDBPaths = ( "1UAO", "2E4E" )

# . Loop over structures.
for pdbPath in _PDBPaths:

    # . Define the energy models.
    mmModel = MMModelOPLS.WithParameterSet ( "protein" )
    nbModel = NBModelCutOff.WithDefaults ( )

    # . Set up the systems.
    system       = ImportSystem ( os.path.join ( dataPath, pdbPath + ".pdb" ), modelNumber = 1, useComponentLibrary = True )
    system.label = "Chignolin ({:s})".format ( pdbPath )
    BuildHydrogenCoordinates3FromConnectivity ( system )
    system.DefineMMModel ( mmModel )
    system.DefineNBModel ( nbModel )
    system.Summary ( )

    # . Get an initial energy.
    system.Energy ( doGradients = True )

    # . Save.
    Pickle ( os.path.join ( outPath, pdbPath + ".pkl" ), system )

# . Footer.
logFile.Footer ( )
