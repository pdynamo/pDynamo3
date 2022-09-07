"""Some definitions."""

import os, os.path

from pCore import TestScript_InputDataPath  , \
                  TestScript_InputPath      , \
                  TestScript_OutputDataPath

# . Local name.
_name  = "pSimulation"
_nameM = "pMolecule"

# . Paths.
dataPath          = TestScript_InputDataPath  ( _name  )
dataPathM         = TestScript_InputDataPath  ( _nameM )
outPath           = TestScript_OutputDataPath ( _name  )
referenceDataPath = os.path.join ( dataPath, "yaml" )
structuresPath    = os.path.join ( os.getenv ( "PDYNAMO3_HOME" ), "structures" )

# . Other options.
_FullVerificationSummary = False
