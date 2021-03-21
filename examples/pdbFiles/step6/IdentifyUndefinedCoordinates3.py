"""Identify a system's undefined coordinates."""

import os, os.path, sys
sys.path.append ( os.path.join ( os.path.dirname ( os.path.realpath ( __file__ ) ), "..", "data" ) )

from Definitions      import outPath
from pCore            import logFile  , \
                             Unpickle
from pSimulation      import IdentifyUndefinedCoordinates3

# . Header.
logFile.Header ( )

# . Get the system.
system = Unpickle ( os.path.join ( outPath, "step4.pkl" ) )
system.Summary ( )

# . Identify undefined coordinates.
IdentifyUndefinedCoordinates3 ( system )

# . Footer.
logFile.Footer ( )
