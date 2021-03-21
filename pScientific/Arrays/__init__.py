"""A package for numerical arrays."""

# . Example:
#    from pCore              import DataType
#    from pScientific.Arrays import Array, ArrayPrint, StorageType
#    x = Array.WithExtent ( 10, dataType = DataType.Real, storageType = StorageType.Antisymmetric )
#    x.Set ( 10.0 )
#    ArrayPrint ( x, title = "Test Antisymmetric Matrix" )

from .AntisymmetricMatrix   import AntisymmetricMatrix
from .Array                 import StorageType
from .ArrayError            import ArrayError
from .ArrayIterator         import ArrayIterator
from .ArrayPrint            import ArrayPrint          , \
                                   ArrayPrint2D
from .ArrayReshape          import Flatten             , \
                                   Reshape
from .BaseIterator          import BaseIterator
from .BooleanArray1D        import BooleanArray1D
from .BooleanArray2D        import BooleanArray2D
from .BooleanArrayND        import BooleanArrayND
from .BooleanIterator       import BooleanIterator
from .DoubleSymmetricMatrix import DoubleSymmetricMatrix
from .IntegerArray1D        import IntegerArray1D
from .IntegerArray2D        import IntegerArray2D
from .IntegerArrayND        import IntegerArrayND
from .IntegerIterator       import IntegerIterator
from .RealArray1D           import RealArray1D
from .RealArray2D           import RealArray2D
from .RealArrayND           import RealArrayND
from .RealIterator          import RealIterator
from .SparseSymmetricMatrix import SparseSymmetricMatrix
from .SymmetricMatrix       import SymmetricMatrix
