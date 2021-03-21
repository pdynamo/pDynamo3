"""Double symmetric matrices."""

from  pCore      import RawObjectConstructor
from .ArrayError import ArrayError
from .Slicing    import ProcessIntegerSlice

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class DoubleSymmetricMatrix:

    # . Methods.
    def __copy__ ( self ):
        """Copying."""
        cdef DoubleSymmetricMatrix clone
        cdef CStatus               cStatus = CStatus_OK
        clone         = self.__class__.Raw ( )
        clone.block   = self.block
        clone.cObject = DoubleSymmetricMatrix_CloneShallow ( self.cObject, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error cloning array." )
        return clone

    def __dealloc__ ( self ):
        """Finalization."""
        self.block     = None
        self._iterator = None
        DoubleSymmetricMatrix_Deallocate ( &self.cObject )

    def __deepcopy__ ( self, memo ):
        """Copying."""
        clone = self.__class__.WithExtent ( self.extent )
        self.CopyTo ( clone )
        return clone

    def __getattr__ ( self, name ):
        """Handle unknown attributes."""
        return getattr ( self.iterator, name )

    def __getitem__ ( self, indices ):
        """Get an item."""
        cdef CInteger i, j, k, l
        cdef CReal    value
        ( i, j, k, l ) = ProcessIntegerSlice ( self.shape, indices )
        value          = DoubleSymmetricMatrix_GetItem ( self.cObject, i, j, k, l, NULL )
        return value

    def __getstate__ ( self ):
        """Return the state."""
        return { "block"  : self.block  ,
                 "extent" : self.extent }

    def __init__ ( self, extent ):
        """Constructor with extent."""
        self._Initialize ( )
        self._Allocate ( extent )

    def __iter__ ( self ):
        return self._MakeIterator ( )

    def __len__ ( self ):
        """Return the size of the symmetricmatrix."""
        return self.size

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setitem__ ( self, indices, CReal value ):
        """Set an item."""
        cdef CInteger i, j, k, l
        ( i, j, k, l ) = ProcessIntegerSlice ( self.shape, indices )
        DoubleSymmetricMatrix_SetItem ( self.cObject, i, j, k, l, value, NULL )

    def __setstate__ ( self, state ):
        """Set the state."""
        self._Allocate ( state["extent"], block = state["block"] )

    def _Allocate ( self, CInteger extent, RealBlock block = None ):
        """Allocation."""
        cdef CRealBlock *cBlock  = NULL
        cdef CStatus     cStatus = CStatus_OK
        if block is     None:  block = RealBlock.WithCapacity ( DoubleSymmetricMatrix_ViewSize ( extent ) )
        if block is not None: cBlock = block.cObject
        self.block   = block
        self.cObject = DoubleSymmetricMatrix_FromExtentBlock ( extent, cBlock, CTrue, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error allocating double symmetric matrix." )

    def _Initialize ( self ):
        """Initialization."""
        self.cObject   = NULL
        self.block     = None
        self._iterator = None

    def _MakeIterator ( self ):
        """Make the default array iterator."""
        cdef RealIterator iterator
        cdef CIterator   *cIterator = NULL
        cdef CStatus      cStatus   = CStatus_OK
        iterator           = RealIterator.Raw ( )
        iterator.cIterator = DoubleSymmetricMatrix_MakeIterator ( self.cObject, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Unable to make iterator." )
        iterator.block     = self.block
        iterator._SetDataPointer ( )
        iterator.Reset ( )
        return iterator

    def CopyTo ( self, DoubleSymmetricMatrix other ):
        """Copying."""
        cdef CStatus cStatus = CStatus_OK
        DoubleSymmetricMatrix_CopyTo ( self.cObject, other.cObject, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error copying double symmetric matrices." )

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def Set ( self, CReal value ):
        """Set all the elements of a symmetricmatrix."""
        DoubleSymmetricMatrix_Set ( self.cObject, value )

    @classmethod
    def WithExtent ( selfClass, extent ):
        """Constructor with extent."""
        return selfClass ( extent )

    # . Properties.
    @property
    def extent ( self ):
        if self.cObject == NULL: return 0
        else:                    return self.cObject.extent
    @property
    def iterator ( self ):
        if self._iterator is None:
            self._iterator = self._MakeIterator ( )
        return self._iterator
    @property
    def rank ( self ): return 4
    @property
    def shape ( self ): return [ self.extent, self.extent, self.extent, self.extent ]
    @property
    def size ( self ):
        if self.cObject == NULL: return 0
        else:                    return self.cObject.size
