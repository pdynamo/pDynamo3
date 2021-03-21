"""Base class for 1-D arrays."""

from  pCore      import RawObjectConstructor
from .ArrayError import ArrayError
from .Slicing    import ProcessSlice1D

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class BaseArray1D:

    # . Methods.
    def __copy__ ( self ):
        """Copying."""
        cdef BaseArray1D clone
        clone = self.__class__.RawWithCObject ( )
        View1D_CopyTo      ( self.cView , clone.cView )
        clone._AssignBlock ( self.block )
        return clone

    def __dealloc__ ( self ):
        """Finalization."""
        self.block = None
        self._DeallocateCObject ( )

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
        cdef CInteger start, stop, stride, extent
        ( isScalar, start, stop, stride, extent ) = ProcessSlice1D ( indices, self.extent )
        if isScalar: return self._ScalarGet ( start )
        else:        return self._ViewGet   ( start, extent, stride )

    def __getstate__ ( self ):
        """Return the state."""
        return { "block"   : self.block   ,
                 "extents" : self.shape   ,
                 "offset"  : self.offset  ,
                 "rank"    : self.rank    ,
                 "size"    : self.size    ,
                 "strides" : self.strides }

    def __init__ ( self, extent ):
        """Constructor with extent."""
        self._Initialize ( )
        self._Allocate ( extent )

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
        cdef CInteger start, stop, stride, extent
        ( isScalar, start, stop, stride, extent ) = ProcessSlice1D ( indices, self.extent )
        if isScalar: self._ScalarSet ( start, value )
        else:        self._ViewSet   ( start, extent, stride, value )

    def __setstate__ ( self, state ):
        """Set the state."""
        cdef CInteger  e, n, o, s
        items = state.get ( "items", None )
        if items is not None:
            self._Allocate ( len ( items ) )
            for ( i, v ) in enumerate ( items ): self.block[i] = v
        else:
            e = state["extents"][0]
            n = state["size"   ]
            o = state["offset" ]
            s = state["strides"][0]
            self._AllocateCObject ( )
            View1D_SetState ( self.cView, e, o, n, s )
            self._AssignBlock ( state["block"] )

    def _Allocate ( self, extent ):
        """Allocation."""
        cdef CStatus cStatus = CStatus_OK
        self._AllocateCObject ( )
        View1D_Initialize ( self.cView, extent, &cStatus )
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

    def _GetRawArray1D ( self, CInteger extent ):
        """Get a raw 1-D array."""
        return self.__class__.RawWithCObject ( )

    def _GetMemoryBlock ( self, size ):
        """Get the memory block."""
        return None

    def _GetRawIterator ( self ):
        """Get the raw iterator."""
        return None

    def _Initialize ( self ):
        """Initialization."""
        self.cView     = NULL
        self.block     = None
        self._iterator = None

    def _MakeIterator ( self ):
        """Make the default array iterator."""
        cdef BaseIterator iterator
        cdef CIterator   *cIterator = NULL
        cdef CStatus      cStatus   = CStatus_OK
        iterator           = self._GetRawIterator ( )
        iterator.cIterator = View1D_MakeIterator ( self.cView, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Unable to make iterator." )
        iterator.block     = self.block
        iterator._SetDataPointer ( )
        iterator.Reset ( )
        return iterator

    def _ScalarGet ( self, CInteger index ):
        """Get a scalar."""
        return None

    def _ScalarSet ( self, CInteger index, value ):
        """Set a scalar."""
        pass

    def _ViewGet ( self, CInteger start, CInteger extent, CInteger stride ):
        """Get a view."""
        cdef BaseArray1D view
        cdef CStatus     cStatus = CStatus_OK
        view = self._GetRawArray1D ( extent )
        View1D_View       ( self.cView , start, extent, stride, view.cView, &cStatus )
        view._AssignBlock ( self.block )
        if cStatus != CStatus_OK: raise ArrayError ( "Error getting view of array." )
        return view

    def _ViewSet ( self, CInteger start, CInteger extent, CInteger stride, value ):
        """Set a view."""
        view = self._ViewGet ( start, extent, stride )
        view.Set ( value )

    def CopyTo ( self, BaseArray1D other not None ):
        """Copying."""
        self.iterator.CopyTo ( other.iterator )

    @classmethod
    def Merge ( selfClass, items, information = {} ):
        """Merging."""
        merged  = None
        lengths = [ len ( item ) for item in items if item is not None ]
        if len ( lengths ) == len ( items ):
            merged = selfClass.WithExtent ( sum ( lengths ) )
            m      = 0
            for ( n, item ) in zip ( lengths, items ):
                item.CopyTo ( merged[m:m+n] )
                m += n
        return merged

    def Prune ( self, Selection selection, information = {} ):
        """Pruning."""
        pruned = self.__class__.WithExtent ( len ( selection ) )
        for ( i, s ) in enumerate ( selection ): pruned[i] = self[s]
        return pruned

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    @classmethod
    def RawWithCObject ( selfClass ):
        """Raw constructor with C object."""
        cdef BaseArray1D self
        cdef CStatus     cStatus = CStatus_OK
        self = selfClass.Raw  ( )
        self._AllocateCObject ( )
        return self

    def Set ( self, value ):
        """Setting."""
        self.iterator.Set ( value )

    @classmethod
    def WithExtent ( selfClass, extent ):
        """Constructor with extent."""
        return selfClass ( extent )

    @classmethod
    def WithShape ( selfClass, shape ):
        """Constructor with shape."""
        return selfClass ( *shape )

    # . Properties.
    @property
    def extent ( self ): return View1D_GetExtent ( self.cView )
    @property
    def iterator ( self ):
        if self._iterator is None:
            self._iterator = self._MakeIterator ( )
        return self._iterator
    @property
    def offset ( self ): return View1D_GetOffset ( self.cView )
    @property
    def rank ( self ): return 1
    @property
    def shape ( self ): return [ self.extent ]
    @property
    def size ( self ): return View1D_GetExtent ( self.cView )
    @property
    def strides ( self ): return [ View1D_GetStride ( self.cView ) ]
