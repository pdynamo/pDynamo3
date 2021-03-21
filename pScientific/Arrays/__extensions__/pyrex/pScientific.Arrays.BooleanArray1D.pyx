"""Boolean 1-D arrays."""

from .ArrayError      import ArrayError
from .BooleanIterator import BooleanIterator

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class BooleanArray1D ( BaseArray1D ):

    cdef void _AllocateCObject ( self ):
        """Allocate the C object."""
        cdef CStatus cStatus = CStatus_OK
        if self.cObject == NULL:
            self.cObject = BooleanArray1D_Allocate ( &cStatus )
            if cStatus != CStatus_OK: raise ArrayError ( "Error allocating C object." )
            self.cView = < CView1D * > self.cObject

    def _AssignBlock ( self, BooleanBlock block not None ):
        """Assign a block to the array."""
        cdef CStatus cStatus = CStatus_OK
        self.block = block
        BooleanArray1D_AssignBlock ( self.cObject, block.cObject, CTrue, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error assigning block to array." )

    cdef void _DeallocateCObject ( self ):
        """Deallocate the C object."""
        BooleanArray1D_Deallocate ( &self.cObject )

    def _GetMemoryBlock ( self, size ):
        """Get the memory block."""
        return BooleanBlock.WithCapacity ( size )

    def _GetRawIterator ( self ):
        """Get the raw iterator."""
        return BooleanIterator.Raw ( )

    def _Initialize ( self ):
        """Initialization."""
        super ( BooleanArray1D, self )._Initialize ( )
        self.cObject = NULL

    def _ScalarGet ( self, CInteger index ):
        """Get a scalar."""
        cdef CBoolean cValue
        cdef CStatus  cStatus = CStatus_OK
        cValue = BooleanArray1D_GetItem ( self.cObject, index, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error getting item." )
        if cValue == CTrue: return True
        else:               return False

    def _ScalarSet ( self, CInteger index, value ):
        """Set a scalar."""
        cdef CBoolean cValue
        cdef CStatus  cStatus = CStatus_OK
        if value: cValue = CTrue
        else:     cValue = CFalse
        BooleanArray1D_SetItem ( self.cObject, index, cValue, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error setting item." )

    def CopyTo ( BooleanArray1D self, BooleanArray1D other not None ):
        """Copying."""
        cdef CStatus cStatus = CStatus_OK
        BooleanArray1D_CopyTo ( self.cObject, other.cObject, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error copying array." )

    def Set ( BooleanArray1D self, value ):
        """Setting."""
        cdef CBoolean cValue
        if value: cValue = CTrue
        else:     cValue = CFalse
        BooleanArray1D_Set ( self.cObject, cValue )
