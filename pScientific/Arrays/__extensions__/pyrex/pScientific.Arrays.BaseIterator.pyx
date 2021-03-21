"""Base class for array iterators."""

from  pCore      import Clone                , \
                        RawObjectConstructor
from .ArrayError import ArrayError

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class BaseIterator:

    # . Methods.
    def __copy__ ( self ):
        """Copying."""
        cdef BaseIterator new
        cdef CStatus      cStatus = CStatus_OK
        new           = self.__class__.Raw ( )
        new.cIterator = Iterator_Clone ( self.cIterator, &cStatus )
        if cStatus   != CStatus_OK: raise ArrayError ( "Error cloning iterator." )
        new.block     = self.block
        new._SetDataPointer ( )
        return new

    def __dealloc__ ( self ):
        """Finalization."""
        self.block = None
        if self.isOwner:
            Iterator_Deallocate ( &self.cIterator )

    def __deepcopy__ ( self, memo ):
        """Copying."""
        cdef BaseIterator new
        cdef CStatus      cStatus = CStatus_OK
        new           = self.__class__.Raw ( )
        new.block     = Clone ( self.block )
        new.cIterator = Iterator_Clone ( self.cIterator, &cStatus )
        if cStatus   != CStatus_OK: raise ArrayError ( "Error cloning iterator." )
        new._SetDataPointer ( )
        return new

    def __getstate__ ( self ):
        """Return the state."""
        cdef CInteger  n = 0
        cdef CInteger *cState  = NULL
        cdef CStatus   cStatus = CStatus_OK
        cState = Iterator_Dump ( self.cIterator, &n, &cStatus )
        pState = [ cState[i] for i in range ( n ) ]
        Integer_Deallocate ( &cState )
        if cStatus != CStatus_OK: raise ArrayError ( "Error getting iterator state." )
        return { "block" : self.block ,
                 "state" : pState     }

    def __init__ ( self ):
        """Constructor."""
        self._Initialize ( )
        self._Allocate   ( )

    def __iter__ ( self ):
        self.Reset ( )
        return self

    def __len__ ( self ):
        """Return the size of the array."""
        return self.size

    def __next__ ( self ):
        raise StopIteration

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ):
        """Set the state."""
        cdef CInteger  i, n = 0
        cdef CInteger *cState  = NULL
        cdef CStatus   cStatus = CStatus_OK
        pState = state["state"]
        n      = len ( pState )
        cState = Integer_Allocate ( n, &cStatus )
        for i from 0 <= i < n: cState[i] = pState[i]
        self.cIterator = Iterator_Load ( n, cState, &cStatus )
        Integer_Deallocate ( &cState )
        if cStatus != CStatus_OK: raise ArrayError ( "Error setting iterator state." )
        self.block = state["block"]
        self._SetDataPointer ( )

    # . Not very useful by itself so do something else here?
    def _Allocate ( self ):
        """Allocation."""
        cdef CStatus cStatus = CStatus_OK
        self.cIterator = Iterator_Allocate ( &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error allocating iterator." )
        self.isOwner = True

    def _Initialize ( self ):
        """Initialization."""
        self.cIterator = NULL
        self.block     = None
        self.isOwner   = True # . The default is always owner.

    cdef _SetDataPointer ( self ):
        """Set the data pointer."""
        pass

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def Reset ( self ):
        """Reset the iterator."""
        Iterator_Reset ( self.cIterator )

    @property
    def dataType ( self ): return self.block.dataType

    @property
    def isFinished ( self ):
        """Is the iterator finished?"""
        return ( Iterator_CurrentIndex ( self.cIterator ) < 0 )

    @property
    def size ( self ):
        """How many items are there in the iterator?"""
        return Iterator_GetSize ( self.cIterator )

