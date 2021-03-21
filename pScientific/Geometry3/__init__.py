"""A package for geometrical manipulations in three dimensions."""

from .Coordinates3             import Coordinates3
from .Geometry3Error           import Geometry3Error
from .Matrix33                 import Matrix33
from .PairListGenerator        import CrossPairList_FromDoubleCoordinates3   , \
                                      CrossPairList_FromSingleCoordinates3   , \
                                      SelfPairList_FromCoordinates3          , \
                                      PairListGenerator
from .RegularGrid              import RegularGrid
from .RegularGridOccupancy     import RegularGridOccupancy
from .Transformation3          import Transformation3
from .Transformation3Container import Transformation3Container
from .Vector3                  import Vector3
