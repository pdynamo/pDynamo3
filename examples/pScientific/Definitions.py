"""Some definitions."""

import os, os.path

from pCore import TestScript_InputDataPath  , \
                  TestScript_InputPath      , \
                  TestScript_OutputDataPath

# . Local name.
_name  = "pScientific"

# . Paths.
dataPath       = TestScript_InputDataPath  ( _name  )
outPath        = TestScript_OutputDataPath ( _name  )
structuresPath = os.path.join ( os.getenv ( "PDYNAMO3_HOME" ), "structures" )

# . Other options.
_FullVerificationSummary = False
