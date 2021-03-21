"""Integer arrays."""

from .CoreError     import CoreError
from .DataTypes     import DataType
from .Serialization import RawObjectConstructor

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class IntegerBlock:

    # . Methods.
    def __copy__ ( self ):
        """Copying."""
        return self.__deepcopy__ ( None )

    def __dealloc__ ( self ):
        """Finalization."""
        IntegerBlock_Deallocate ( &self.cObject )

    def __deepcopy__ ( self, memo ):
        """Copying."""
        cdef IntegerBlock new
        new = self.__class__.WithCapacity ( self.size )
        self.CopyTo ( new )
        return new

    def __getitem__ ( self, i ):
        """Get an item."""
        if isinstance ( i, int ):
            if i < 0: i += len ( self )
            if ( i < 0 ) or ( i >= len ( self ) ): raise IndexError ( "Index {:d} out of range.".format ( i ) )
            return IntegerBlock_GetItem ( self.cObject, i, NULL )
        else: raise TypeError ( "Expecting integer index not {!r}.".format ( type ( i ) ) )

    def __getstate__ ( self ):
        """Return the state."""
        cdef CInteger i
        items = []
        for i from 0 <= i < self.size:
            items.append ( IntegerBlock_GetItem ( self.cObject, i, NULL ) )
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

    def __setitem__ ( self, i, CInteger value ):
        """Set an item."""
        if isinstance ( i, int ):
            if i < 0: i += len ( self )
            if ( i < 0 ) or ( i >= len ( self ) ): raise IndexError ( "Index {:d} out of range.".format ( i ) )
            IntegerBlock_SetItem ( self.cObject, i, value, NULL )
        else: raise TypeError ( "Expecting integer index not {!r}.".format ( type ( i ) ) )

    def __setstate__ ( self, state ):
        """Set the state."""
        cdef CInteger i, size, value
        items = state.get ( "items", None )
        if items is not None:
            size = len ( items )
            self._Allocate ( size )
            for i from 0 <= i < size:
                value = items[i]
                IntegerBlock_SetItem ( self.cObject, i, value, NULL )

    def _Allocate ( self, size ):
        """Allocation."""
        cdef CStatus status
        if self.cObject == NULL:
            status       = CStatus_OK
            self.cObject = IntegerBlock_Allocate ( size, &status )
            IntegerBlock_Reference ( self.cObject )
            if status   != CStatus_OK: raise CoreError ( "C object allocation failure." )
        else:
            raise CoreError ( "C object already exists." )

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL

    def CopyTo ( self, IntegerBlock other ):
        """Copying."""
        IntegerBlock_CopyTo ( self.cObject, other.cObject, NULL )

    @classmethod
    def Merge ( selfClass, entries, information = {} ):
        """Merging."""
        n    = sum ( [ len ( entry ) for entry in entries ] )
        self = selfClass.WithCapacity ( n )
        n    = 0
        for entry in entries:
            for ( i, v ) in enumerate ( entry ): self[i+n] = entry[i]
            n += len ( entry )
        return self

    def Prune ( self, Selection selection, information = {} ):
        """Pruning."""
        new = self.__class__.WithCapacity ( len ( selection ) )
        for ( i, s ) in enumerate ( selection ): new[i] = self[s]
        return new

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def Set ( self, CInteger value ):
        """Set all the elements."""
        IntegerBlock_Set ( self.cObject, value )

    @classmethod
    def WithCapacity ( selfClass, capacity ):
        """Constructor with a capacity."""
        return selfClass ( capacity )

    # . Properties.
    @property
    def dataType ( self ): return DataType.Integer
    @property
    def rank ( self ): return 1
    @property
    def shape ( self ): return [ self.size ]
    @property
    def size ( self ):
        if self.cObject == NULL: return 0
        else:                    return self.cObject.capacity
