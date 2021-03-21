"""Generate a MM model for a system."""

import os, os.path, sys
sys.path.append ( os.path.join ( os.path.dirname ( os.path.realpath ( __file__ ) ), "..", "data" ) )

from Definitions       import outPath
from pCore             import logFile  , \
                              Pickle   , \
                              Unpickle
from pMolecule.MMModel import MMModelOPLS

# . Header.
logFile.Header ( )

# . Set up a MM model.
mmModel = MMModelOPLS.WithParameterSet ( "protein" )

# . Get the system.
system = Unpickle ( os.path.join ( outPath, "step3.pkl" ) )
system.Summary ( )

# . Add the energy model.
system.DefineMMModel ( mmModel )
system.Summary ( )

# . Save the system.
Pickle ( os.path.join ( outPath, "step4.pkl" ), system )

# . Footer.
logFile.Footer ( )
