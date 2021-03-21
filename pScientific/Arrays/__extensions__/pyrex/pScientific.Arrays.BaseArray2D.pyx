"""Base class for 2-D arrays."""

from  pCore      import RawObjectConstructor
from .ArrayError import ArrayError
from .Slicing    import ProcessSlice1D

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class BaseArray2D:

    # . Methods.
    def __copy__ ( self ):
        """Copying."""
        cdef BaseArray2D clone
        clone = self.__class__.RawWithCObject ( )
        View2D_CopyTo      ( self.cView , clone.cView )
        clone._AssignBlock ( self.block )
        return clone

    def __dealloc__ ( self ):
        """Finalization."""
        self.block = None
        MultiSlice_Deallocate ( &self.cMultiSlice )
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
        else:           return self._View2DGet ( )

    def __getstate__ ( self ):
        """Return the state."""
        return { "block"   : self.block   ,
                 "extents" : self.shape   ,
                 "offset"  : self.offset  ,
                 "rank"    : self.rank    ,
                 "size"    : self.size    ,
                 "strides" : self.strides }

    def __init__ ( self, rows, columns ):
        """Constructor with extents."""
        self._Initialize ( )
        self._Allocate ( rows, columns )

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
        if   rank == 0: self._ScalarSet ( value )
        elif rank == 1: self._View1DSet ( value )
        else:           self._View2DSet ( value )

    def __setstate__ ( self, state ):
        """Set the state."""
        cdef CInteger e0, e1, n, o, s0, s1
        cdef CStatus  cStatus = CStatus_OK
        items = state.get ( "items", None )
        shape = state.get ( "shape", None )
        if ( items is not None ) and ( shape is not None ):
            self._Allocate ( shape[0], shape[1] )
            for ( i, v ) in enumerate ( items ): self.block[i] = v
        else:
            n          = state["size"   ]
            o          = state["offset" ]
            ( e0, e1 ) = state["extents"]
            ( s0, s1 ) = state["strides"]
            self._AllocateCObject ( )
            View2D_SetState ( self.cView, e0, e1, o, n, s0, s1 )
            self.cMultiSlice = MultiSlice_Allocate ( 2, &cStatus )
            if cStatus != CStatus_OK: raise ArrayError ( "Error setting array state." )
            self._AssignBlock ( state["block"] )

    def _Allocate ( self, rows, columns ):
        """Allocation."""
        cdef CStatus cStatus = CStatus_OK
        self._AllocateCObject ( )
        View2D_Initialize ( self.cView, rows, columns, &cStatus )
        self.cMultiSlice = MultiSlice_Allocate ( 2, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error allocating array." )
        block = self._GetMemoryBlock ( self.size )
        self._AssignBlock ( block )

    cdef void _AllocateCObject ( self ):
        """Allocate the C object."""
        pass

    def _AssignBlock ( self, block ):
        """Assign a block to the array."""
        pass

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
        return self.__class__.RawWithCObject ( )

    def _GetRawIterator ( self ):
        """Get the raw iterator."""
        return None

    def _Initialize ( self ):
        """Initialization."""
        self.cMultiSlice = NULL
        self.cView       = NULL
        self.block       = None
        self._iterator   = None

    def _MakeIterator ( self ):
        """Make the default array iterator."""
        cdef BaseIterator iterator
        cdef CIterator   *cIterator = NULL
        cdef CStatus      cStatus   = CStatus_OK
        iterator           = self._GetRawIterator ( )
        iterator.cIterator = View2D_MakeIterator ( self.cView, &cStatus )
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
        View2D_View1DMultiSlice ( self.cView , self.cMultiSlice, view.cView, &cStatus )
        view._AssignBlock       ( self.block )
        if cStatus != CStatus_OK: raise ArrayError ( "Error getting 1-D view of array." )
        return view

    def _View2DGet ( self ):
        """Get a 2-D view."""
        cdef BaseArray2D view
        cdef CStatus     cStatus = CStatus_OK
        view = self._GetRawArray2D ( )
        View2D_ViewMultiSlice ( self.cView , self.cMultiSlice, view.cView, &cStatus )
        view._AssignBlock     ( self.block )
        if cStatus != CStatus_OK: raise ArrayError ( "Error getting 2-D view of array." )
        return view

    def _View1DSet ( self, value ):
        """Set a view."""
        view = self._View1DGet ( )
        view.Set ( value )

    def _View2DSet ( self, value ):
        """Set a view."""
        view = self._View2DGet ( )
        view.Set ( value )

    def CopyTo ( self, BaseArray2D other not None ):
        """Copying."""
        self.iterator.CopyTo ( other.iterator )

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    @classmethod
    def RawWithCObject ( selfClass ):
        """Raw constructor with C object."""
        cdef BaseArray2D self
        cdef CStatus     cStatus = CStatus_OK
        self = selfClass.Raw ( )
        self.cMultiSlice = MultiSlice_Allocate ( 2, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error allocating array." )
        self._AllocateCObject ( )
        return self

    def RowIterator ( self, Selection selection = None ):
        """Make a row iterator."""
        cdef BaseIterator iterator
        cdef CIterator   *cIterator = NULL
        cdef CStatus      cStatus   = CStatus_OK
        if ( selection is None ) or ( len ( selection ) == 0 ):
            iterator = self.iterator
        else:
            iterator           = self._GetRawIterator ( )
            iterator.cIterator = View2D_MakeRowIterator ( self.cView, selection.cObject, &cStatus )
            if cStatus != CStatus_OK: raise ArrayError ( "Unable to make row iterator." )
            iterator.block     = self.block
            iterator._SetDataPointer ( )
        return iterator

    def Set ( self, value ):
        """Setting."""
        self.iterator.Set ( value )

    @classmethod
    def WithExtents ( selfClass, extent0, extent1 ):
        """Constructor with extents."""
        return selfClass ( extent0, extent1 )

    @classmethod
    def WithShape ( selfClass, shape ):
        """Constructor with shape."""
        return selfClass ( *shape )

    # . Properties.
    @property
    def columns ( self ): return View2D_GetColumns ( self.cView )
    @property
    def iterator ( self ):
        if self._iterator is None:
            self._iterator = self._MakeIterator ( )
        return self._iterator
    @property
    def offset ( self ): return View2D_GetOffset ( self.cView )
    @property
    def rank  ( self ): return 2
    @property
    def rows  ( self ): return View2D_GetRows ( self.cView )
    @property
    def shape ( self ): return [ self.rows, self.columns ]
    @property
    def size  ( self ): return View2D_GetSize ( self.cView )
    @property
    def strides ( self ):
        cdef CInteger d
        return [ View2D_GetStride ( self.cView, d, NULL ) for d in range ( self.rank ) ]
