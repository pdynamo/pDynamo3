"""Integer 2-D arrays."""

from .ArrayError      import ArrayError
from .IntegerIterator import IntegerIterator

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class IntegerArray2D ( BaseArray2D ):

    cdef void _AllocateCObject ( self ):
        """Allocate the C object."""
        cdef CStatus cStatus = CStatus_OK
        if self.cObject == NULL:
            self.cObject = IntegerArray2D_Allocate ( &cStatus )
            if cStatus != CStatus_OK: raise ArrayError ( "Error allocating C object." )
            self.cView = < CView2D * > self.cObject

    def _AssignBlock ( self, IntegerBlock block not None ):
        """Assign a block to the array."""
        cdef CStatus cStatus = CStatus_OK
        self.block = block
        IntegerArray2D_AssignBlock ( self.cObject, block.cObject, CTrue, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error assigning block to array." )

    cdef void _DeallocateCObject ( self ):
        """Deallocate the C object."""
        IntegerArray2D_Deallocate ( &self.cObject )

    def _GetMemoryBlock ( self, size ):
        """Get the memory block."""
        return IntegerBlock.WithCapacity ( size )

    def _GetRawArray1D ( self ):
        """Get a raw 1-D array."""
        return IntegerArray1D.RawWithCObject ( )

    def _GetRawIterator ( self ):
        """Get the raw iterator."""
        return IntegerIterator.Raw ( )

    def _Initialize ( self ):
        """Initialization."""
        super ( IntegerArray2D, self )._Initialize ( )
        self.cObject = NULL

    def _ScalarGet ( self ):
        """Get a scalar."""
        cdef CInteger cValue
        cdef CStatus  cStatus = CStatus_OK
        cValue = IntegerArray2D_GetItemMultiSlice ( self.cObject, self.cMultiSlice, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error getting item." )
        return cValue

    def _ScalarSet ( self, value ):
        """Set a scalar."""
        cdef CInteger cValue  = value
        cdef CStatus  cStatus = CStatus_OK
        IntegerArray2D_SetItemMultiSlice ( self.cObject, self.cMultiSlice, cValue, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error setting item." )

    def CopyTo ( IntegerArray2D self, IntegerArray2D other not None ):
        """Copying."""
        cdef CStatus cStatus = CStatus_OK
        IntegerArray2D_CopyTo ( self.cObject, other.cObject, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error copying array." )

    def Set ( IntegerArray2D self, value ):
        """Setting."""
        cdef CInteger cValue = value
        IntegerArray2D_Set ( self.cObject, cValue )
