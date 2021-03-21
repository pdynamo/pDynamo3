"""Some definitions."""

import os, os.path

from addOns.pcetk import _MEADPath
from pCore        import TestScript_InputDataPath    , \
                         TestScript_InputPath        , \
                         TestScript_OutputDataPath   , \
                         TestScriptExit_NotInstalled

# . Abort if the MEAD path not found.
meadPath = os.getenv ( _MEADPath )
if ( meadPath is None ) or not ( os.path.isdir ( meadPath ) ): TestScriptExit_NotInstalled ( )

# . Local name.
_name = os.path.join ( "addOns", "pcetk" )

# . Paths.
dataPath = TestScript_InputDataPath  ( _name  )
outPath  = TestScript_OutputDataPath ( _name  )

# . Other options.
_FullVerificationSummary = False
