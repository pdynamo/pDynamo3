"""Add counterions to the system."""

import os, os.path, sys
sys.path.append ( os.path.join ( os.path.dirname ( os.path.realpath ( __file__ ) ), "..", "data" ) )

from Definitions        import dataPath       , \
                               outPath
from pBabel             import ExportSystem   , \
                               ImportSystem
from pCore              import logFile        , \
                               Pickle         , \
                               Unpickle
from pMolecule.MMModel  import MMModelOPLS
from pScientific.Arrays import Array
from pSimulation        import AddCounterIons

# . Header.
logFile.Header ( )

# . Parameters.
# . Box sizes.
_XBox = 75.0
_YBox = 60.0
_ZBox = 60.0

# . Number and type of ions to add.
_NNegative   = 27
_NPositive   = 24
_NegativeIon = "chloride"
_PositiveIon = "potassium"

# . Reorient option.
_Reorient = False

# . Get the system.
system = Unpickle ( os.path.join ( outPath, "step7.pkl" ) )
system.Summary ( )

# . Reorient the system if necessary (see the results of GetSolvationInformation.py).
if _Reorient:
    masses = Array.FromIterable ( [ atom.mass for atom in system.atoms ] )
    system.coordinates3.ToPrincipalAxes ( weights = masses )

# . Get the positive and negative ions.
if _NNegative > 0:
    anion = ImportSystem ( os.path.join ( dataPath, _NegativeIon + ".mol" ) )
    anion.DefineMMModel ( MMModelOPLS.WithParameterSet ( "bookSmallExamples" ) )
    anion.Summary ( )
if _NPositive > 0:
    cation = ImportSystem ( os.path.join ( dataPath, _PositiveIon + ".mol" ) )
    cation.DefineMMModel ( MMModelOPLS.WithParameterSet ( "bookSmallExamples" ) )
    cation.Summary ( )

# . Add the counterions.
newSystem = AddCounterIons ( system, _NNegative, anion, _NPositive, cation, ( _XBox, _YBox, _ZBox ) )

# . Save the combined system.
Pickle ( os.path.join ( outPath, "step8_a.pkl" ), newSystem )

# . Print PDB file.
ExportSystem ( os.path.join ( outPath, "step8_a.pdb" ), newSystem )

# . Footer.
logFile.Footer ( )
