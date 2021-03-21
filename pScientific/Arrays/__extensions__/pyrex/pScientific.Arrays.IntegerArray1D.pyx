"""Integer 1-D arrays."""

from .ArrayError      import ArrayError
from .IntegerIterator import IntegerIterator

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class IntegerArray1D ( BaseArray1D ):

    cdef void _AllocateCObject ( self ):
        """Allocate the C object."""
        cdef CStatus cStatus = CStatus_OK
        if self.cObject == NULL:
            self.cObject = IntegerArray1D_Allocate ( &cStatus )
            if cStatus != CStatus_OK: raise ArrayError ( "Error allocating C object." )
            self.cView = < CView1D * > self.cObject

    def _AssignBlock ( self, IntegerBlock block not None ):
        """Assign a block to the array."""
        cdef CStatus cStatus = CStatus_OK
        self.block = block
        IntegerArray1D_AssignBlock ( self.cObject, block.cObject, CTrue, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error assigning block to array." )

    cdef void _DeallocateCObject ( self ):
        """Deallocate the C object."""
        IntegerArray1D_Deallocate ( &self.cObject )

    def _GetMemoryBlock ( self, size ):
        """Get the memory block."""
        return IntegerBlock.WithCapacity ( size )

    def _GetRawIterator ( self ):
        """Get the raw iterator."""
        return IntegerIterator.Raw ( )

    def _Initialize ( self ):
        """Initialization."""
        super ( IntegerArray1D, self )._Initialize ( )
        self.cObject = NULL

    def _ScalarGet ( self, CInteger index ):
        """Get a scalar."""
        cdef CInteger cValue
        cdef CStatus  cStatus = CStatus_OK
        cValue = IntegerArray1D_GetItem ( self.cObject, index, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error getting item." )
        return cValue

    def _ScalarSet ( self, CInteger index, value ):
        """Set a scalar."""
        cdef CInteger cValue  = value
        cdef CStatus  cStatus = CStatus_OK
        IntegerArray1D_SetItem ( self.cObject, index, cValue, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error setting item." )

    def CopyTo ( self, IntegerArray1D other not None ):
        """Copying."""
        cdef CStatus cStatus = CStatus_OK
        IntegerArray1D_CopyTo ( self.cObject, other.cObject, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error copying array." )

    def Set ( self, value ):
        """Setting."""
        cdef CInteger cValue = value
        IntegerArray1D_Set ( self.cObject, cValue )

    def Sort ( self ):
        """Sorting."""
        IntegerArray1D_Sort ( self.cObject )
