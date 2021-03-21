"""A package containing miscellaneous scientific algorithms and data."""

# . To import a constant or unit use:
#   from pScientific import Constants
#   a = Constants.Atomic_Mass * b, etc.

from .Magnitudes    import Magnitude        , \
                           Magnitude_Adjust
from .PeriodicTable import IsMainGroup      , \
                           IsMetal          , \
                           IsNonMetal       , \
                           IsSemiMetal      , \
                           PeriodicTable
