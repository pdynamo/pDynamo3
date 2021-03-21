"""Real N-D arrays."""

from .ArrayError   import ArrayError
from .RealIterator import RealIterator

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class RealArrayND ( BaseArrayND ):

    cdef void _DeallocateCObject ( self ):
        """Deallocate the C object."""
        if self.cObject != NULL:
            self.cObject.data = NULL
            self.cObject.view = NULL
            RealArrayND_Deallocate ( &self.cObject )

    def _GetMemoryBlock ( self, size ):
        """Get the memory block."""
        return RealBlock.WithCapacity ( size )

    def _GetRawArray1D ( self ):
        """Get a raw 1-D array."""
        cdef RealArray1D array
        array = RealArray1D.Raw ( )
        array._AllocateCObject  ( )
        return array

    def _GetRawArray2D ( self ):
        """Get a raw 2-D array."""
        cdef RealArray2D array
        cdef CStatus     cStatus = CStatus_OK
        array             = RealArray2D.Raw  ( )
        array.cMultiSlice = MultiSlice_Allocate ( 2, &cStatus )
        array._AllocateCObject ( )
        if cStatus != CStatus_OK: raise ArrayError ( "Error getting raw 2-D array." )
        return array

    def _GetRawIterator ( self ):
        """Get the raw iterator."""
        return RealIterator.Raw ( )

    def _Initialize ( self ):
        """Initialization."""
        super ( RealArrayND, self )._Initialize ( )
        self.cObject = NULL

    cdef void _MakeCObject ( self ):
        """Make the C object."""
        cdef RealBlock block
        cdef CStatus   cStatus = CStatus_OK
        if self.cObject == NULL:
            block        = self.block
            self.cObject = RealArrayND_FromViewBlock ( self.cView, block.cObject, CTrue, &cStatus )
            if cStatus != CStatus_OK: raise ArrayError ( "Error making array C object." )

    def _ScalarGet ( self ):
        """Get a scalar."""
        cdef CReal   value
        cdef CStatus cStatus = CStatus_OK
        value = RealArrayND_GetItemMultiSlice ( self.cObject, self.cMultiSlice, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error getting item." )
        return value

    def _ScalarSet ( self, CReal value ):
        """Set a scalar."""
        cdef CStatus cStatus = CStatus_OK
        RealArrayND_SetItemMultiSlice ( self.cObject, self.cMultiSlice, value, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error setting item." )

    def CopyTo ( RealArrayND self, RealArrayND other not None ):
        """Copying."""
        cdef CStatus cStatus = CStatus_OK
        RealArrayND_CopyTo ( self.cObject, other.cObject, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error copying array." )

    def Set ( RealArrayND self, CReal value ):
        """Setting."""
        cdef CStatus cStatus = CStatus_OK
        RealArrayND_Set ( self.cObject, value, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error setting array." )
