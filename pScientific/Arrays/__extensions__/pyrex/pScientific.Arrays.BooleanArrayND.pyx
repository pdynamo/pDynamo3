"""Boolean N-D arrays."""

from .ArrayError      import ArrayError
from .BooleanIterator import BooleanIterator

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class BooleanArrayND ( BaseArrayND ):

    cdef void _DeallocateCObject ( self ):
        """Deallocate the C object."""
        if self.cObject != NULL:
            self.cObject.data = NULL
            self.cObject.view = NULL
            BooleanArrayND_Deallocate ( &self.cObject )

    def _GetMemoryBlock ( self, size ):
        """Get the memory block."""
        return BooleanBlock.WithCapacity ( size )

    def _GetRawArray1D ( self ):
        """Get a raw 1-D array."""
        cdef BooleanArray1D array
        array = BooleanArray1D.Raw ( )
        array._AllocateCObject  ( )
        return array

    def _GetRawArray2D ( self ):
        """Get a raw 2-D array."""
        cdef BooleanArray2D array
        cdef CStatus        cStatus = CStatus_OK
        array             = BooleanArray2D.Raw  ( )
        array.cMultiSlice = MultiSlice_Allocate ( 2, &cStatus )
        array._AllocateCObject ( )
        if cStatus != CStatus_OK: raise ArrayError ( "Error getting raw 2-D array." )
        return array

    def _GetRawIterator ( self ):
        """Get the raw iterator."""
        return BooleanIterator.Raw ( )

    def _Initialize ( self ):
        """Initialization."""
        super ( BooleanArrayND, self )._Initialize ( )
        self.cObject = NULL

    cdef void _MakeCObject ( self ):
        """Make the C object."""
        cdef BooleanBlock block
        cdef CStatus      cStatus = CStatus_OK
        if self.cObject == NULL:
            block        = self.block
            self.cObject = BooleanArrayND_FromViewBlock ( self.cView, block.cObject, CTrue, &cStatus )
            if cStatus != CStatus_OK: raise ArrayError ( "Error making array C object." )

    def _ScalarGet ( self ):
        """Get a scalar."""
        cdef CBoolean value
        cdef CStatus  cStatus = CStatus_OK
        value = BooleanArrayND_GetItemMultiSlice ( self.cObject, self.cMultiSlice, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error getting item." )
        return value

    def _ScalarSet ( self, CBoolean value ):
        """Set a scalar."""
        cdef CStatus cStatus = CStatus_OK
        BooleanArrayND_SetItemMultiSlice ( self.cObject, self.cMultiSlice, value, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error setting item." )

    def CopyTo ( BooleanArrayND self, BooleanArrayND other not None ):
        """Copying."""
        cdef CStatus cStatus = CStatus_OK
        BooleanArrayND_CopyTo ( self.cObject, other.cObject, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error copying array." )

    def Set ( BooleanArrayND self, CBoolean value ):
        """Setting."""
        cdef CStatus cStatus = CStatus_OK
        BooleanArrayND_Set ( self.cObject, value, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error setting array." )
