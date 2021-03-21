"""Boolean 2-D arrays."""

from .ArrayError      import ArrayError
from .BooleanIterator import BooleanIterator

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class BooleanArray2D ( BaseArray2D ):

    cdef void _AllocateCObject ( self ):
        """Allocate the C object."""
        cdef CStatus cStatus = CStatus_OK
        if self.cObject == NULL:
            self.cObject = BooleanArray2D_Allocate ( &cStatus )
            if cStatus != CStatus_OK: raise ArrayError ( "Error allocating C object." )
            self.cView = < CView2D * > self.cObject

    def _AssignBlock ( self, BooleanBlock block not None ):
        """Assign a block to the array."""
        cdef CStatus cStatus = CStatus_OK
        self.block = block
        BooleanArray2D_AssignBlock ( self.cObject, block.cObject, CTrue, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error assigning block to array." )

    cdef void _DeallocateCObject ( self ):
        """Deallocate the C object."""
        BooleanArray2D_Deallocate ( &self.cObject )

    def _GetMemoryBlock ( self, size ):
        """Get the memory block."""
        return BooleanBlock.WithCapacity ( size )

    def _GetRawArray1D ( self ):
        """Get a raw 1-D array."""
        return BooleanArray1D.RawWithCObject ( )

    def _GetRawIterator ( self ):
        """Get the raw iterator."""
        return BooleanIterator.Raw ( )

    def _Initialize ( self ):
        """Initialization."""
        super ( BooleanArray2D, self )._Initialize ( )
        self.cObject = NULL

    def _ScalarGet ( self ):
        """Get a scalar."""
        cdef CBoolean cValue
        cdef CStatus  cStatus = CStatus_OK
        cValue = BooleanArray2D_GetItemMultiSlice ( self.cObject, self.cMultiSlice, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error getting item." )
        if cValue == CTrue: return True
        else:               return False

    def _ScalarSet ( self, value ):
        """Set a scalar."""
        cdef CBoolean cValue
        cdef CStatus  cStatus = CStatus_OK
        if value: cValue = CTrue
        else:     cValue = CFalse
        BooleanArray2D_SetItemMultiSlice ( self.cObject, self.cMultiSlice, cValue, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error setting item." )

    def CopyTo ( BooleanArray2D self, BooleanArray2D other not None ):
        """Copying."""
        cdef CStatus cStatus = CStatus_OK
        BooleanArray2D_CopyTo ( self.cObject, other.cObject, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error copying array." )

    def Set ( BooleanArray2D self, value ):
        """Setting."""
        cdef CBoolean cValue
        if value: cValue = CTrue
        else:     cValue = CFalse
        BooleanArray2D_Set ( self.cObject, cValue )
