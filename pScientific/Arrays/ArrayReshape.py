"""Functions for reshaping arrays."""

from  pCore          import DataType
from .ArrayError     import ArrayError
from .BaseArray1D    import BaseArray1D
from .BaseArray2D    import BaseArray2D
from .BaseArrayND    import BaseArrayND
from .BooleanArray1D import BooleanArray1D
from .BooleanArray2D import BooleanArray2D
from .BooleanArrayND import BooleanArrayND
from .IntegerArray1D import IntegerArray1D
from .IntegerArray2D import IntegerArray2D
from .IntegerArrayND import IntegerArrayND
from .RealArray1D    import RealArray1D
from .RealArray2D    import RealArray2D
from .RealArrayND    import RealArrayND

#===================================================================================================================================
# . Utility methods.
#===================================================================================================================================
# . Currently done in Python - could be moved to C if necessary.
def _GetReshapeState ( array, shape ):
    """Get the state of the reshaped array."""
    # . No checking of array and shape - assumed to be OK.
    # . Get the old state and flatten as much as possible.
    oldState = array.__getstate__ ( )
    extents0 = oldState["extents"]
    rank0    = 0
    size0    = oldState["size"   ]
    strides0 = oldState["strides"]
    for d0 in range ( 1, oldState["rank"] ):
        e0 = extents0[d0]
        s0 = strides0[d0]
        if strides0[d0-1] == ( e0 * s0 ):
            extents0[rank0] *= e0
        else:
            rank0 += 1
            extents0[rank0] = e0
        strides0[rank0] = s0
    rank0 += 1
#    print ( "Old state>", oldState )
#    print ( "Flattened state>", rank0, extents0[0:rank0], strides0[0:rank0] ) 
    # . Remove unnecessary trailing extents of 1 (except the last).
    while ( extents0[rank0-1] == 1 ) and ( rank0 > 1 ):
        extents0.pop ( )
        strides0.pop ( )
        rank0 -= 1
#    print ( "State with trailing 1s removed>", rank0, extents0[0:rank0], strides0[0:rank0] ) 
    # . Determine if the old shape can be reshaped to the new one. For this to be possible, each contiguous 
    # . multiple of extents in the new shape must be equal to an extent in the (flattened) old shape.
    # . The strides need to be verified for cases where extra extents of 1 are added.
    extents1 = list ( shape )
    rank1    = len  ( shape )
    strides1 = [ size0 for d1 in range ( rank1 ) ] # . Use size0 for cases such as [12] -> [1,3,4].
    d1       = rank1-1
    for d0 in range ( rank0-1, -1, -1 ):
        e0 = extents0[d0]
        s0 = strides0[d0]
        e1 = 1
        while d1 >= 0:
            strides1[d1] = s0 * e1
            e1          *= extents1[d1]
            d1          -= 1
            if e1 >= e0: break
        if e1 != e0: raise ArrayError ( "New shape incompatible with array." )
#    print ( "New state>", rank1, extents1, strides1 ) 
    # . Finish up.
    return { "block"   : oldState["block" ] ,
             "extents" : extents1           ,
             "offset"  : oldState["offset"] ,
             "rank"    : rank1              ,
             "size"    : size0              ,
             "strides" : strides1           }

#===================================================================================================================================
# . Public methods.
#===================================================================================================================================
def Flatten ( array ):
    """Flatten an array."""
    if array.rank > 1: return Reshape ( array, [ array.size ] )
    else:              return array

def Reshape ( array, shape, resultClass = None ):
    """Reshape an array."""
    # . Check input class.
    if not ( isinstance ( array, BaseArray1D ) or \
             isinstance ( array, BaseArray2D ) or \
             isinstance ( array, BaseArrayND ) ):
        raise ArrayError ( "Unable to reshape {:s} arrays.".format ( array.__class__.__name__ ) )
    # . Check input shape.
    rank = len ( shape )
    if ( rank == 0 ) or any ( [ ( e <= 0 ) for e in shape ] ):
        raise ArrayError ( "Invalid input shape." )
    # . Get default output class.
    outClass = None
    if   array.dataType is DataType.Boolean:
        genericClass = BooleanArrayND
        if   rank == 1: outClass = BooleanArray1D
        elif rank == 2: outClass = BooleanArray2D
    elif array.dataType is DataType.Integer:
        genericClass = IntegerArrayND
        if   rank == 1: outClass = IntegerArray1D
        elif rank == 2: outClass = IntegerArray2D
    elif array.dataType is DataType.Real:
        genericClass = RealArrayND
        if   rank == 1: outClass = RealArray1D
        elif rank == 2: outClass = RealArray2D
    if outClass is None: outClass = genericClass
    # . Determine the result class.
    # . Use outClass.
    if resultClass is None:
        resultClass = outClass
    # . Use resultClass if it is a subclass of outClass or of genericClass.
    elif not ( issubclass ( resultClass, genericClass ) or \
               issubclass ( resultClass, outClass     ) ):
        raise ArrayError ( "Unable to reshape {:s} into {:s} arrays.".format ( array.__class__.__name__, resultClass.__name__ ) )
    # . Create the new array.
    newState = _GetReshapeState ( array, shape )
    new      = resultClass.Raw ( )
    new.__setstate__ ( newState )
    return new

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
