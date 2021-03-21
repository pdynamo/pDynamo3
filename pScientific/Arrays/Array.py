"""Numerical arrays."""

from  enum                  import Enum
from  pCore                 import DataType
from .AntisymmetricMatrix   import AntisymmetricMatrix
from .ArrayError            import ArrayError
from .BooleanArray1D        import BooleanArray1D
from .BooleanArray2D        import BooleanArray2D
from .BooleanArrayND        import BooleanArrayND
from .DoubleSymmetricMatrix import DoubleSymmetricMatrix
from .IntegerArray1D        import IntegerArray1D
from .IntegerArray2D        import IntegerArray2D
from .IntegerArrayND        import IntegerArrayND
from .RealArray1D           import RealArray1D
from .RealArray2D           import RealArray2D
from .RealArrayND           import RealArrayND
from .SymmetricMatrix       import SymmetricMatrix

#===================================================================================================================================
# . Definitions.
#===================================================================================================================================
# . Storage type.
class StorageType ( Enum ):
    """Storage types."""
    Antisymmetric   = 10
    DoubleSymmetric = 20
    Regular         = 30
    Symmetric       = 40

# . Other types?
#    SparseSymmetric - WithExtentAndSize
#    Tridiagonal

# . Properties.
_DataTypeConverters = { DataType.Boolean : bool  ,
                        DataType.Integer : int   ,
                        DataType.Real    : float }

#===================================================================================================================================
# . Methods.
#===================================================================================================================================
def _IsShapeOK ( shape ):
    """Check the shape specification."""
    isOK = isinstance ( shape, ( list, tuple ) ) and ( len ( shape ) > 0 )
    for i in shape:
        isOK = isOK and isinstance ( i, int ) and ( i >= 0 )
    if not isOK: raise ArrayError ( "Invalid array shape: {:s}.".format ( repr ( shape ) ) )
    return len ( shape )

def FromIterable ( iterable, dataType = DataType.Real, storageType = StorageType.Regular ):
    """Constructor from an iterable."""
    array = WithExtent ( len ( iterable ), dataType = dataType, storageType = storageType )
    for ( i, value ) in enumerate ( iterable ): array[i] = value
    return array

def IdentityMatrix ( extent, dataType = DataType.Real, storageType = StorageType.Regular ):
    """Identity matrix constructor."""
    array = WithShape ( ( extent, extent ), dataType = dataType, storageType = storageType )
    array.diagonal.Set    ( _DataTypeConverters[dataType] ( 1 ) )
    array.offDiagonal.Set ( _DataTypeConverters[dataType] ( 0 ) )
    return array

# . Could be made more efficient.
def WithExtent ( extent, dataType = DataType.Real, storageType = StorageType.Regular ):
    """Constructor with extent."""
    if   storageType == StorageType.Antisymmetric  : shape = 2 * [ extent ]
    elif storageType == StorageType.DoubleSymmetric: shape = 4 * [ extent ]
    elif storageType == StorageType.Symmetric      : shape = 2 * [ extent ]
    else: shape = [ extent ]
    return WithShape ( shape, dataType = dataType, storageType = storageType )

def WithExtents ( *extents, dataType = DataType.Real, storageType = StorageType.Regular ):
    """Constructor with extents."""
    return WithShape ( extents, dataType = dataType, storageType = storageType )

def WithShape ( shape, dataType = DataType.Real, storageType = StorageType.Regular ):
    """Constructor with shape."""
    rank = _IsShapeOK ( shape )
    if storageType == StorageType.Antisymmetric:
        if ( rank == 2 ) and ( shape[0] == shape[1] ) and ( dataType == DataType.Real ):
            array = AntisymmetricMatrix.WithExtent ( shape[0] )
        else:
            raise ArrayError ( "Invalid antisymmetric matrix specification: {:s}/{:s}.".format ( repr ( shape ), dataType.name ) )
    elif storageType == StorageType.DoubleSymmetric:
        if ( rank     == 4        ) and \
           ( shape[0] == shape[1] ) and \
           ( shape[0] == shape[2] ) and \
           ( shape[0] == shape[3] ) and \
           ( dataType == DataType.Real ):
            array = DoubleSymmetricMatrix.WithExtent ( shape[0] )
        else:
            raise ArrayError ( "Invalid double symmetric matrix specification: {:s}/{:s}.".format ( repr ( shape ), dataType.name ) )
    elif storageType == StorageType.Symmetric:
        if ( rank == 2 ) and ( shape[0] == shape[1] ) and ( dataType == DataType.Real ):
            array = SymmetricMatrix.WithExtent ( shape[0] )
        else:
            raise ArrayError ( "Invalid symmetric matrix specification: {:s}/{:s}.".format ( repr ( shape ), dataType.name ) )
    elif dataType == DataType.Boolean:
        if   rank == 1: array = BooleanArray1D.WithExtent  (  shape[0]           )
        elif rank == 2: array = BooleanArray2D.WithExtents (  shape[0], shape[1] )
        else:           array = BooleanArrayND.WithExtents ( *shape )
    elif dataType == DataType.Integer:
        if   rank == 1: array = IntegerArray1D.WithExtent  (  shape[0]           )
        elif rank == 2: array = IntegerArray2D.WithExtents (  shape[0], shape[1] )
        else:           array = IntegerArrayND.WithExtents ( *shape )
    elif dataType == DataType.Real:
        if   rank == 1: array = RealArray1D.WithExtent     (  shape[0]           )
        elif rank == 2: array = RealArray2D.WithExtents    (  shape[0], shape[1] )
        else:           array = RealArrayND.WithExtents    ( *shape )
    return array

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
