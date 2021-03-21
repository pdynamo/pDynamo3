"""Merge systems into one."""

import os, os.path, sys
sys.path.append ( os.path.join ( os.path.dirname ( os.path.realpath ( __file__ ) ), "..", "data" ) )

from Definitions      import outPath
from pCore            import logFile  , \
                             Pickle   , \
                             Unpickle
from pSimulation      import MergeByAtom

# . Header.
logFile.Header ( )

# . Specify the file names containing the systems to merge.
pklFiles = [ "part1.pkl", "part2.pkl", "part3.pkl" ]

# . Recover the systems from the files.
parts = []
for pklFile in pklFiles:
    parts.append ( Unpickle ( os.path.join ( outPath, pklFile ) ) )

# . Merge the systems.
system = MergeByAtom ( parts )
system.Summary ( )

# . Save the system.
Pickle ( os.path.join ( outPath, "step5.pkl" ), system )

# . Footer.
logFile.Footer ( )
