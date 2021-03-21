"""Some definitions."""

import os, os.path

from pCore import TestScript_InputDataPath  , \
                  TestScript_OutputDataPath

# . Local name.
_name = "pMolecule"

# . Paths.
dataPath          = TestScript_InputDataPath  ( _name )
outPath           = TestScript_OutputDataPath ( _name )
referenceDataPath = os.path.join ( dataPath, "reference" )
structuresPath    = os.path.join ( os.getenv ( "PDYNAMO3_HOME" ), "structures" )

# . Other options.
_FullVerificationSummary = False
