"""Some definitions."""

import os, os.path

from pCore             import TestScript_InputDataPath    , \
                              TestScript_InputPath        , \
                              TestScript_OutputDataPath   , \
                              TestScriptExit_NotInstalled
from pMolecule.QCModel import _DFTBCommand

# . Abort if an executable DFTB does not exist.
command = os.getenv ( _DFTBCommand )
if  ( command is None ) or not ( os.path.isfile ( command ) and os.access ( command, os.X_OK ) ):
    TestScriptExit_NotInstalled ( )

# . Local name.
_name  = "dftbPlus"
_nameM = "pMolecule"

# . Paths.
dataPath          = TestScript_InputDataPath  ( _name  )
dataPathM         = TestScript_InputDataPath  ( _nameM )
outPath           = TestScript_OutputDataPath ( _name  )
referenceDataPath = os.path.join ( dataPath, "reference" )
structuresPath    = os.path.join ( os.getenv ( "PDYNAMO3_HOME" ), "structures" )

# . Other options.
_FullVerificationSummary = False

