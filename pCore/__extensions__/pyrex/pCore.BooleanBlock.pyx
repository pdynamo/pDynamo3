"""Boolean arrays."""

from .CoreError     import CoreError
from .DataTypes     import DataType
from .Serialization import RawObjectConstructor

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class BooleanBlock:

    # . Methods.
    def __copy__ ( self ):
        """Copying."""
        return self.__deepcopy__ ( None )

    def __dealloc__ ( self ):
        """Finalization."""
        BooleanBlock_Deallocate ( &self.cObject )

    def __deepcopy__ ( self, memo ):
        """Copying."""
        cdef BooleanBlock new
        new = self.__class__.WithCapacity ( self.size )
        self.CopyTo ( new )
        return new

    def __getitem__ ( self, i ):
        """Get an item."""
        if isinstance ( i, int ):
            if i < 0: i += len ( self )
            if ( i < 0 ) or ( i >= len ( self ) ): raise IndexError ( "Index {:d} out of range.".format ( i ) )
            return BooleanBlock_GetItem ( self.cObject, i, NULL )
        else: raise TypeError ( "Expecting integer index not {!r}.".format ( type ( i ) ) )

    def __getstate__ ( self ):
        """Return the state."""
        cdef CInteger i
        items = []
        for i from 0 <= i < self.size:
            items.append ( BooleanBlock_GetItem ( self.cObject, i, NULL ) )
        return { "items" : items }

    def __init__ ( self, capacity ):
        """Constructor with a capacity."""
        self._Initialize ( )
        self._Allocate ( capacity )

    def __len__ ( self ):
        """Return the size of the array."""
        return self.size

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setitem__ ( self, i, CBoolean value ):
        """Set an item."""
        if isinstance ( i, int ):
            if i < 0: i += len ( self )
            if ( i < 0 ) or ( i >= len ( self ) ): raise IndexError ( "Index {:d} out of range.".format ( i ) )
            BooleanBlock_SetItem ( self.cObject, i, value, NULL )
        else: raise TypeError ( "Expecting integer index not {!r}.".format ( type ( i ) ) )

    def __setstate__ ( self, state ):
        """Set the state."""
        cdef CBoolean value
        cdef CInteger i, size
        items = state.get ( "items", None )
        if items is not None:
            size = len ( items )
            self._Allocate ( size )
            for i from 0 <= i < size:
                value = items[i]
                BooleanBlock_SetItem ( self.cObject, i, value, NULL )

    def _Allocate ( self, size ):
        """Allocation."""
        cdef CStatus status
        if self.cObject == NULL:
            status       = CStatus_OK
            self.cObject = BooleanBlock_Allocate ( size, &status )
            BooleanBlock_Reference ( self.cObject )
            if status   != CStatus_OK: raise CoreError ( "C object allocation failure." )
        else:
            raise CoreError ( "C object already exists." )

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL

    def CopyTo ( self, BooleanBlock other ):
        """Copying."""
        BooleanBlock_CopyTo ( self.cObject, other.cObject, NULL )

    def Negate ( self ):
        """Negation."""
        BooleanBlock_Negate ( self.cObject )

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def Set ( self, CBoolean value ):
        """Set all the elements."""
        BooleanBlock_Set ( self.cObject, value )

    @classmethod
    def WithCapacity ( selfClass, capacity ):
        """Constructor with a capacity."""
        return selfClass ( capacity )

    # . Properties.
    @property
    def dataType ( self ): return DataType.Boolean
    @property
    def rank ( self ): return 1
    @property
    def shape ( self ): return [ self.size ]
    @property
    def size ( self ):
        if self.cObject == NULL: return 0
        else:                    return self.cObject.capacity
