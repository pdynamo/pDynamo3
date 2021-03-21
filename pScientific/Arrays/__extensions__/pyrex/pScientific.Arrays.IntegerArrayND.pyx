"""Integer N-D arrays."""

from .ArrayError      import ArrayError
from .IntegerIterator import IntegerIterator

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class IntegerArrayND ( BaseArrayND ):

    cdef void _DeallocateCObject ( self ):
        """Deallocate the C object."""
        if self.cObject != NULL:
            self.cObject.data = NULL
            self.cObject.view = NULL
            IntegerArrayND_Deallocate ( &self.cObject )

    def _GetMemoryBlock ( self, size ):
        """Get the memory block."""
        return IntegerBlock.WithCapacity ( size )

    def _GetRawArray1D ( self ):
        """Get a raw 1-D array."""
        cdef IntegerArray1D array
        array = IntegerArray1D.Raw ( )
        array._AllocateCObject  ( )
        return array

    def _GetRawArray2D ( self ):
        """Get a raw 2-D array."""
        cdef IntegerArray2D array
        cdef CStatus        cStatus = CStatus_OK
        array             = IntegerArray2D.Raw  ( )
        array.cMultiSlice = MultiSlice_Allocate ( 2, &cStatus )
        array._AllocateCObject ( )
        if cStatus != CStatus_OK: raise ArrayError ( "Error getting raw 2-D array." )
        return array

    def _GetRawIterator ( self ):
        """Get the raw iterator."""
        return IntegerIterator.Raw ( )

    def _Initialize ( self ):
        """Initialization."""
        super ( IntegerArrayND, self )._Initialize ( )
        self.cObject = NULL

    cdef void _MakeCObject ( self ):
        """Make the C object."""
        cdef IntegerBlock block
        cdef CStatus      cStatus = CStatus_OK
        if self.cObject == NULL:
            block        = self.block
            self.cObject = IntegerArrayND_FromViewBlock ( self.cView, block.cObject, CTrue, &cStatus )
            if cStatus != CStatus_OK: raise ArrayError ( "Error making array C object." )

    def _ScalarGet ( self ):
        """Get a scalar."""
        cdef CInteger value
        cdef CStatus  cStatus = CStatus_OK
        value = IntegerArrayND_GetItemMultiSlice ( self.cObject, self.cMultiSlice, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error getting item." )
        return value

    def _ScalarSet ( self, CInteger value ):
        """Set a scalar."""
        cdef CStatus cStatus = CStatus_OK
        IntegerArrayND_SetItemMultiSlice ( self.cObject, self.cMultiSlice, value, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error setting item." )

    def CopyTo ( IntegerArrayND self, IntegerArrayND other not None ):
        """Copying."""
        cdef CStatus cStatus = CStatus_OK
        IntegerArrayND_CopyTo ( self.cObject, other.cObject, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error copying array." )

    def Set ( IntegerArrayND self, CInteger value ):
        """Setting."""
        cdef CStatus cStatus = CStatus_OK
        IntegerArrayND_Set ( self.cObject, value, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error setting array." )
