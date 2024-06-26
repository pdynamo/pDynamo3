"""Some definitions."""

import os, os.path

from pCore import TestScript_InputDataPath  , \
                  TestScript_InputPath      , \
                  TestScript_OutputDataPath

# . Local name.
_name = os.path.join ( "addOns", "pyCPR" )

# . Paths.
dataPath = TestScript_InputDataPath  ( _name  )
outPath  = TestScript_OutputDataPath ( _name  )

# . Other options.
_FullVerificationSummary = False
