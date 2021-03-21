"""Integer iterator."""

from .ArrayError import ArrayError

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class IntegerIterator ( BaseIterator ):
    """Integer iterator."""

    def __next__ ( self ):
        cdef CInteger i
        i = Iterator_NextIndex ( self.cIterator )
        if i >= 0: return self.cData[i]
        else:      raise StopIteration

    def _Initialize ( self ):
        """Initialization."""
        super ( IntegerIterator, self )._Initialize ( )
        self.cData = NULL

    cdef _SetDataPointer ( self ):
        """Set the data pointer."""
        cdef IntegerBlock block = self.block
        if ( block is not None ) and ( block.cObject != NULL ):
            self.cData = &(block.cObject.items[Iterator_DataOffSet ( self.cIterator )])

    # . Unary functions.
    def AbsoluteMaximum ( self ):
        """Absolute maximum."""
        return IntegerIterator_AbsoluteMaximum ( self.cIterator, self.cData )

    def CountSmall ( self, CInteger tolerance ):
        """Count the number of values in the array whose absolute value is less than or equal to tolerance."""
        return IntegerIterator_CountSmall ( self.cIterator, self.cData, tolerance )

    def DotSelf ( self ):
        """Dot product of self."""
        return IntegerIterator_DotSelf ( self.cIterator, self.cData )

    def FilterGreaterThan ( self, CInteger tolerance, CInteger value = 0 ):
        """Filter values in the array whose value is greater than or equal to tolerance."""
        cdef CStatus cStatus = CStatus_OK
        IntegerIterator_FilterGreaterThan ( self.cIterator, self.cData, tolerance, value, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error in filter greater than operation." )

    def FilterLessThan ( self, CInteger tolerance, CInteger value = 0 ):
        """Filter values in the array whose value is less than or equal to tolerance."""
        cdef CStatus cStatus = CStatus_OK
        IntegerIterator_FilterLessThan ( self.cIterator, self.cData, tolerance, value, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error in filter less than operation." )

    def FilterSmall ( self, CInteger tolerance, CInteger value = 0 ):
        """Filter values in the array whose absolute value is less than or equal to tolerance."""
        cdef CStatus cStatus = CStatus_OK
        IntegerIterator_FilterSmall ( self.cIterator, self.cData, tolerance, value, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error in filter small operation." )

    def Increment ( self, CInteger value ):
        """Increment the values of an array."""
        cdef CStatus cStatus = CStatus_OK
        IntegerIterator_Increment ( self.cIterator, self.cData, value, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error in increment operation." )

    def Maximum ( self ):
        """Maximum of values."""
        return IntegerIterator_Maximum ( self.cIterator, self.cData )

    def Minimum ( self ):
        """Minimum of values."""
        return IntegerIterator_Minimum ( self.cIterator, self.cData )

    def Product ( self ):
        """Product of values."""
        return IntegerIterator_Product ( self.cIterator, self.cData )

    def Scale ( self, CInteger value ):
        """Scale the values of an array."""
        cdef CStatus cStatus = CStatus_OK
        IntegerIterator_Scale ( self.cIterator, self.cData, value, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error in scale operation." )

    def Set ( self, CInteger value ):
        """Set the values of an array."""
        cdef CStatus cStatus = CStatus_OK
        IntegerIterator_Set ( self.cIterator, self.cData, value, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error in set operation." )

    def Sparsity ( self, CInteger tolerance = 1 ): # . Smallest non-zero value!
        """Sparsity of array."""
        return IntegerIterator_Sparsity ( self.cIterator, self.cData, tolerance )

    def Sum ( self ):
        """Sum of values."""
        return IntegerIterator_Sum ( self.cIterator, self.cData )

    # . Binary functions.
    def Add ( self, other not None, CInteger scale = 1 ):
        """Add two arrays with optional scaling."""
        cdef IntegerIterator iOther
        cdef CStatus         cStatus = CStatus_OK
        if isinstance ( other, IntegerIterator ): iOther = other
        else:                                     iOther = other.iterator
        IntegerIterator_Add ( self.cIterator, self.cData, iOther.cIterator, iOther.cData, scale, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error in add operation." )

    def CopyTo ( self, other not None ):
        """Copy |self| to |other|."""
        cdef IntegerIterator iOther
        cdef CStatus         cStatus = CStatus_OK
        if isinstance ( other, IntegerIterator ): iOther = other
        else:                                     iOther = other.iterator
        IntegerIterator_CopyTo ( self.cIterator, self.cData, iOther.cIterator, iOther.cData, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error in copy operation." )

    def Dot ( self, other not None ):
        """Dot two arrays."""
        cdef IntegerIterator iOther
        cdef CInteger        result  = 0
        cdef CStatus         cStatus = CStatus_OK
        if isinstance ( other, IntegerIterator ): iOther = other
        else:                                     iOther = other.iterator
        result = IntegerIterator_Dot ( self.cIterator, self.cData, iOther.cIterator, iOther.cData, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error in dot operation." )
        return result

    def Multiply ( self, other not None ):
        """Multiply two arrays."""
        cdef IntegerIterator iOther
        cdef CStatus         cStatus = CStatus_OK
        if isinstance ( other, IntegerIterator ): iOther = other
        else:                                     iOther = other.iterator
        IntegerIterator_Multiply ( self.cIterator, self.cData, iOther.cIterator, iOther.cData, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error in multiply operation." )

    def Swap ( self, other not None ):
        """Swap two arrays."""
        cdef IntegerIterator iOther
        cdef CStatus         cStatus = CStatus_OK
        if isinstance ( other, IntegerIterator ): iOther = other
        else:                                     iOther = other.iterator
        IntegerIterator_Swap ( self.cIterator, self.cData, iOther.cIterator, iOther.cData, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error in swap operation." )
