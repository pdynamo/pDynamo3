"""Some definitions."""

import os, os.path

from pCore import TestScript_InputPath      , \
                  TestScript_OutputDataPath

# . Local name.
_name      = "pdbFiles"

# . Input paths.
_inputPath = TestScript_InputPath ( _name )
dataPath   = os.path.join ( _inputPath, "data"  )
step2Path  = os.path.join ( _inputPath, "step2" )

# . Output paths.
outPath     = TestScript_OutputDataPath ( _name )
pdbDataPath = os.path.join ( outPath, "pdbData" )
if not os.path.exists ( pdbDataPath ): os.mkdir ( pdbDataPath )
