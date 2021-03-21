"""Real 1-D arrays."""

from .ArrayError   import ArrayError
from .RealIterator import RealIterator

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class RealArray1D ( BaseArray1D ):

    cdef void _AllocateCObject ( self ):
        """Allocate the C object."""
        cdef CStatus cStatus = CStatus_OK
        if self.cObject == NULL:
            self.cObject = RealArray1D_Allocate ( &cStatus )
            if cStatus != CStatus_OK: raise ArrayError ( "Error allocating C object." )
            self.cView = < CView1D * > self.cObject

    def _AssignBlock ( self, RealBlock block not None ):
        """Assign a block to the array."""
        cdef CStatus cStatus = CStatus_OK
        self.block = block
        RealArray1D_AssignBlock ( self.cObject, block.cObject, CTrue, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error assigning block to array." )

    cdef void _DeallocateCObject ( self ):
        """Deallocate the C object."""
        RealArray1D_Deallocate ( &self.cObject )

    def _GetMemoryBlock ( self, size ):
        """Get the memory block."""
        return RealBlock.WithCapacity ( size )

    def _GetRawIterator ( self ):
        """Get the raw iterator."""
        return RealIterator.Raw ( )

    def _Initialize ( self ):
        """Initialization."""
        super ( RealArray1D, self )._Initialize ( )
        self.cObject = NULL

    def _ScalarGet ( self, CInteger index ):
        """Get a scalar."""
        cdef CReal    cValue
        cdef CStatus  cStatus = CStatus_OK
        cValue = RealArray1D_GetItem ( self.cObject, index, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error getting item." )
        return cValue

    def _ScalarSet ( self, CInteger index, value ):
        """Set a scalar."""
        cdef CStatus cStatus = CStatus_OK
        RealArray1D_SetItem ( self.cObject, index, value, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error setting item." )

    def CopyTo ( self, RealArray1D other not None ):
        """Copying."""
        cdef CStatus cStatus = CStatus_OK
        RealArray1D_CopyTo ( self.cObject, other.cObject, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error copying array." )

    def Set ( self, value ):
        """Setting."""
        RealArray1D_Set ( self.cObject, value )

    def Sort ( self ):
        """Sorting."""
        RealArray1D_Sort ( self.cObject )
