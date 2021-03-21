"""Some definitions."""

import os, os.path

from pCore             import TestScript_InputDataPath    , \
                              TestScript_InputPath        , \
                              TestScript_OutputDataPath   , \
                              TestScriptExit_NotInstalled
from pMolecule.QCModel import _ORCACommand

# . Abort if an executable ORCA does not exist.
command = os.getenv ( _ORCACommand )
if  ( command is None ) or not ( os.path.isfile ( command ) and os.access ( command, os.X_OK ) ):
    TestScriptExit_NotInstalled ( )

# . For some reason, ORCA does not like dots in file names.
# . Local name.
_name  = "orca"
_nameM = "pMolecule"

# . Paths.
dataPath  = TestScript_InputDataPath  ( _name  )
dataPathM = TestScript_InputDataPath  ( _nameM )
