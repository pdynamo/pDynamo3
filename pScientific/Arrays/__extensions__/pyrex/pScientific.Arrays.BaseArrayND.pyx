"""Base class for N-D arrays."""

from  pCore      import RawObjectConstructor
from .ArrayError import ArrayError

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class BaseArrayND:

    # . Methods.
    def __copy__ ( self ):
        """Copying."""
        cdef BaseArrayND clone
        cdef CStatus     cStatus = CStatus_OK
        clone             = self.__class__.Raw ( )
        clone.cMultiSlice = MultiSlice_Allocate ( self.rank , &cStatus )
        clone.cView       = ViewND_Clone        ( self.cView, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error cloning array." )
        clone.block       = self.block
        clone._MakeCObject ( )
        return clone

    def __dealloc__ ( self ):
        """Finalization."""
        self.block = None
        MultiSlice_Deallocate ( &self.cMultiSlice )
        ViewND_Deallocate     ( &self.cView       )
        self._DeallocateCObject ( )

    def __deepcopy__ ( self, memo ):
        """Copying."""
        clone = self.__class__.WithShape ( self.shape )
        self.CopyTo ( clone )
        return clone

    def __getattr__ ( self, name ):
        """Handle unknown attributes."""
        return getattr ( self.iterator, name )

    def __getitem__ ( self, indices ):
        """Get an item."""
        cdef CInteger rank
        rank = CProcessMultiSlice ( self.shape, indices, self.cMultiSlice )
        if   rank == 0: return self._ScalarGet ( )
        elif rank == 1: return self._View1DGet ( )
        elif rank == 2: return self._View2DGet ( )
        else:           return self._ViewGet   ( )

    def __getstate__ ( self ):
        """Return the state."""
        return { "block"   : self.block   ,
                 "extents" : self.shape   ,
                 "offset"  : self.offset  ,
                 "rank"    : self.rank    ,
                 "size"    : self.size    ,
                 "strides" : self.strides }

    def __init__ ( self, shape ):
        """Constructor with a shape."""
        self._Initialize ( )
        self._Allocate ( shape )

    def __iter__ ( self ):
        return self._MakeIterator ( )

    def __len__ ( self ):
        """Return the size of the array."""
        return self.size

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setitem__ ( self, indices, value ):
        """Set an item."""
        cdef CInteger rank
        rank = CProcessMultiSlice ( self.shape, indices, self.cMultiSlice )
        if rank == 0: self._ScalarSet ( value )
        else:         self._ViewSet   ( value )

    def __setstate__ ( self, state ):
        """Set the state."""
        cdef CInteger  d, e, n, o, r, s
        cdef CInteger *cState  = NULL
        cdef CStatus   cStatus = CStatus_OK
        items = state.get ( "items", None )
        shape = state.get ( "shape", None )
        if ( items is not None ) and ( shape is not None ):
            self._Allocate ( shape )
            for ( i, v ) in enumerate ( items ): self.block[i] = v
        else:
            n = state["size"  ]
            r = state["rank"  ]
            o = state["offset"]
            self.cView = ViewND_FromState ( r, o, n, &cStatus )
            for ( d, ( e, s ) ) in enumerate ( zip ( state["extents"], state["strides"] ) ):
                ViewND_SetExtentStride ( self.cView, d, e, s, NULL )
            self.cMultiSlice = MultiSlice_Allocate ( self.rank, &cStatus )
            if cStatus != CStatus_OK: raise ArrayError ( "Error setting array state." )
            self.block       = state["block"]
            self._MakeCObject ( )

    def _Allocate ( self, shape ):
        """Allocation."""
        cdef CInteger  cRank
        cdef CInteger *cShape  = NULL
        cdef CStatus   cStatus = CStatus_OK
        cRank  = len ( shape )
        cShape = Integer_Allocate ( cRank, &cStatus )
        if cShape != NULL:
            for i from 0 <= i < cRank: cShape[i] = shape[i]
        self.cView       = ViewND_AllocateWithShape ( cRank, cShape, &cStatus )
        self.cMultiSlice = MultiSlice_Allocate      ( cRank,         &cStatus )
        Integer_Deallocate ( &cShape )
        if cStatus != CStatus_OK: raise ArrayError ( "Error allocating array." )
        self.block = self._GetMemoryBlock ( self.size )
        self._MakeCObject ( )

    cdef void _DeallocateCObject ( self ):
        """Deallocate the C object."""
        pass

    def _GetMemoryBlock ( self, size ):
        """Get the memory block."""
        return None

    def _GetRawArray1D ( self ):
        """Get a raw 1-D array."""
        return None

    def _GetRawArray2D ( self ):
        """Get a raw 2-D array."""
        return None

    def _GetRawIterator ( self ):
        """Get the raw iterator."""
        return None

    def _Initialize ( self ):
        """Initialization."""
        self.cMultiSlice = NULL
        self.cView       = NULL
        self.block       = None
        self._iterator   = None

    cdef void _MakeCObject ( self ):
        """Make the C object."""
        pass

    def _MakeIterator ( self ):
        """Make the default array iterator."""
        cdef BaseIterator iterator
        cdef CIterator   *cIterator = NULL
        cdef CStatus      cStatus   = CStatus_OK
        iterator           = self._GetRawIterator ( )
        iterator.cIterator = ViewND_MakeIterator ( self.cView, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Unable to make iterator." )
        iterator.block     = self.block
        iterator._SetDataPointer ( )
        iterator.Reset ( )
        return iterator

    def _ScalarGet ( self ):
        """Get a scalar."""
        return None

    def _ScalarSet ( self, value ):
        """Set a scalar."""
        pass

    def _View1DGet ( self ):
        """Get a 1-D view."""
        cdef BaseArray1D view
        cdef CStatus     cStatus = CStatus_OK
        view = self._GetRawArray1D ( )
        ViewND_View1DMultiSlice ( self.cView, self.cMultiSlice, view.cView, &cStatus )
        view._AssignBlock       ( self.block )
        if cStatus != CStatus_OK: raise ArrayError ( "Error getting 1-D view of array." )
        return view

    def _View2DGet ( self ):
        """Get a 2-D view."""
        cdef BaseArray2D view
        cdef CStatus     cStatus = CStatus_OK
        view = self._GetRawArray2D ( )
        ViewND_View2DMultiSlice ( self.cView, self.cMultiSlice, view.cView, &cStatus )
        view._AssignBlock       ( self.block )
        if cStatus != CStatus_OK: raise ArrayError ( "Error getting 2-D view of array." )
        return view

    def _ViewGet ( self ):
        """Get a view."""
        cdef BaseArrayND view
        cdef CStatus     cStatus = CStatus_OK
        view             = self.__class__.Raw ( )
        view.cView       = ViewND_ViewMultiSlice ( self.cView, self.cMultiSlice, &cStatus )
        view.cMultiSlice = MultiSlice_Allocate ( view.rank, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error getting view of array." )
        view.block       = self.block
        view._MakeCObject ( )
        return view

    # . More efficient is to get C iterator from multislice directly and then use C iterator set function.
    def _ViewSet ( self, value ):
        """Set a view."""
        view = self._ViewGet ( )
        view.Set ( value )

    def CopyTo ( self, BaseArrayND other not None ):
        """Copying."""
        self.iterator.CopyTo ( other.iterator )

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def Set ( self, value ):
        """Setting."""
        self.iterator.Set ( value )

    @classmethod
    def WithExtent ( selfClass, extent ):
        """Constructor with an extent."""
        return selfClass ( [ extent ] )

    @classmethod
    def WithExtents ( selfClass, *extents ):
        """Constructor with extents."""
        return selfClass ( extents )

    @classmethod
    def WithShape ( selfClass, shape ):
        """Constructor with shape."""
        return selfClass ( shape )

    # . Properties.
    @property
    def iterator ( self ):
        if self._iterator is None:
            self._iterator = self._MakeIterator ( )
        return self._iterator
    @property
    def offset ( self ): return ViewND_GetOffset ( self.cView )
    @property
    def rank ( self ): return ViewND_GetRank ( self.cView )
    @property
    def shape ( self ):
        cdef CInteger d
        return [ ViewND_GetExtent ( self.cView, d, NULL ) for d in range ( self.rank ) ]
    @property
    def size ( self ): return ViewND_GetSize ( self.cView )
    @property
    def strides ( self ):
        cdef CInteger d
        return [ ViewND_GetStride ( self.cView, d, NULL ) for d in range ( self.rank ) ]
