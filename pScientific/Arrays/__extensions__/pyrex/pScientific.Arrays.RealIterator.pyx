"""Real iterator."""

from .ArrayError import ArrayError

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Defaults for some tolerances, etc.
_DefaultSmall = 1.0e-300

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class RealIterator ( BaseIterator ):
    """Integer iterator."""

    def __next__ ( self ):
        cdef CInteger i
        i = Iterator_NextIndex ( self.cIterator )
        if i >= 0: return self.cData[i]
        else:      raise StopIteration

    def _Initialize ( self ):
        """Initialization."""
        super ( RealIterator, self )._Initialize ( )
        self.cData = NULL

    cdef _SetDataPointer ( self ):
        """Set the data pointer."""
        cdef RealBlock block = self.block
        if ( block is not None ) and ( block.cObject != NULL ):
            self.cData = &(block.cObject.items[Iterator_DataOffSet ( self.cIterator )])

    # . Unary functions.
    def Absolute ( self ):
        """Absolute values."""
        RealIterator_Absolute ( self.cIterator, self.cData )

    def AbsoluteMaximum ( self ):
        """Absolute maximum."""
        return RealIterator_AbsoluteMaximum ( self.cIterator, self.cData )

    def CountSmall ( self, CReal tolerance ):
        """Count the number of values in the array whose absolute value is less than or equal to tolerance."""
        return RealIterator_CountSmall ( self.cIterator, self.cData, tolerance )

    def DotSelf ( self ):
        """Dot product of self."""
        return RealIterator_DotSelf ( self.cIterator, self.cData )

    def Exponential ( self ):
        """Exponential."""
        cdef CStatus cStatus = CStatus_OK
        RealIterator_Exponential ( self.cIterator, self.cData, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error in exponential operation." )

    def FilterGreaterThan ( self, CReal tolerance, CReal value = 0.0 ):
        """Filter values in the array whose value is greater than or equal to tolerance."""
        cdef CStatus cStatus = CStatus_OK
        RealIterator_FilterGreaterThan ( self.cIterator, self.cData, tolerance, value, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error in filter greater than operation." )

    def FilterLessThan ( self, CReal tolerance, CReal value = 0.0 ):
        """Filter values in the array whose value is less than or equal to tolerance."""
        cdef CStatus cStatus = CStatus_OK
        RealIterator_FilterLessThan ( self.cIterator, self.cData, tolerance, value, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error in filter less than operation." )

    def FilterSmall ( self, CReal tolerance, CReal value = 0.0 ):
        """Filter values in the array whose absolute value is less than or equal to tolerance."""
        cdef CStatus cStatus = CStatus_OK
        RealIterator_FilterSmall ( self.cIterator, self.cData, tolerance, value, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error in filter small operation." )

    def Increment ( self, CReal value ):
        """Increment the values of an array."""
        cdef CStatus cStatus = CStatus_OK
        RealIterator_Increment ( self.cIterator, self.cData, value, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error in increment operation." )

    def Maximum ( self ):
        """Maximum of values."""
        return RealIterator_Maximum ( self.cIterator, self.cData )

    def Minimum ( self ):
        """Minimum of values."""
        return RealIterator_Minimum ( self.cIterator, self.cData )

    def NaturalLogarithm ( self ):
        """Natural logarithm."""
        cdef CStatus cStatus = CStatus_OK
        RealIterator_NaturalLogarithm ( self.cIterator, self.cData, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error in natural logarithm operation." )

    def Norm2 ( self ):
        """Norm-2."""
        return RealIterator_Norm2 ( self.cIterator, self.cData )

    def Normalize ( self, CReal tolerance = _DefaultSmall ):
        """Normalization."""
        cdef CStatus cStatus = CStatus_OK
        RealIterator_Normalize ( self.cIterator, self.cData, tolerance, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error in normalize operation." )

    def Power ( self, CReal power ):
        """Power."""
        cdef CStatus cStatus = CStatus_OK
        RealIterator_Power ( self.cIterator, self.cData, power, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error in power operation." )

    def Product ( self ):
        """Product of values."""
        return RealIterator_Product ( self.cIterator, self.cData )

    def Reciprocate ( self, CReal tolerance = _DefaultSmall, CReal value = 0.0 ):
        """Reciprocate."""
        cdef CStatus cStatus = CStatus_OK
        RealIterator_Reciprocate ( self.cIterator, self.cData, tolerance, value, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error in reciprocate operation." )

    # . Read as "reciprocate and then power" so power is positive here!
    def ReciprocatePower ( self, CReal power, CReal tolerance = _DefaultSmall, CReal value = 0.0 ):
        """Reciprocate."""
        cdef CStatus cStatus = CStatus_OK
        RealIterator_ReciprocatePower ( self.cIterator, self.cData, power, tolerance, value, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error in reciprocate power operation." )

    def RootMeanSquare ( self ):
        """Root mean square."""
        return RealIterator_RootMeanSquare ( self.cIterator, self.cData )

    def Scale ( self, CReal value ):
        """Scale the values of an array."""
        cdef CStatus cStatus = CStatus_OK
        RealIterator_Scale ( self.cIterator, self.cData, value, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error in scale operation." )

    def Set ( self, CReal value ):
        """Set the values of an array."""
        cdef CStatus cStatus = CStatus_OK
        RealIterator_Set ( self.cIterator, self.cData, value, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error in set operation." )

    def Sparsity ( self, CReal tolerance = _DefaultSmall ):
        """Sparsity of array."""
        return RealIterator_Sparsity ( self.cIterator, self.cData, tolerance )

    def Square ( self ):
        """Square."""
        cdef CStatus cStatus = CStatus_OK
        RealIterator_Square ( self.cIterator, self.cData, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error in square operation." )

    def SquareRoot ( self ):
        """Square root."""
        cdef CStatus cStatus = CStatus_OK
        RealIterator_SquareRoot ( self.cIterator, self.cData, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error in square root operation." )

    def Sum ( self ):
        """Sum of values."""
        return RealIterator_Sum ( self.cIterator, self.cData )

    # . Binary functions.
    def Add ( self, other not None, CReal scale = 1.0 ):
        """Add two arrays with optional scaling."""
        cdef RealIterator iOther
        cdef CStatus      cStatus = CStatus_OK
        if isinstance ( other, RealIterator ): iOther = other
        else:                                  iOther = other.iterator
        RealIterator_Add ( self.cIterator, self.cData, iOther.cIterator, iOther.cData, scale, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error in add operation." )

    def CopyTo ( self, other not None ):
        """Copy |self| to |other|."""
        cdef RealIterator iOther
        cdef CStatus      cStatus = CStatus_OK
        if isinstance ( other, RealIterator ): iOther = other
        else:                                  iOther = other.iterator
        RealIterator_CopyTo ( self.cIterator, self.cData, iOther.cIterator, iOther.cData, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error in copy operation." )

    def Divide ( self, other not None, CReal tolerance = _DefaultSmall, CReal value = 0.0 ):
        """Divide self by other."""
        cdef RealIterator iOther
        cdef CStatus      cStatus = CStatus_OK
        if isinstance ( other, RealIterator ): iOther = other
        else:                                  iOther = other.iterator
        RealIterator_Divide ( self.cIterator, self.cData, iOther.cIterator, iOther.cData, tolerance, value, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error in divide operation." )

    def Dot ( self, other not None ):
        """Dot two arrays."""
        cdef CReal        result  = 0.0
        cdef RealIterator iOther
        cdef CStatus      cStatus = CStatus_OK
        if isinstance ( other, RealIterator ): iOther = other
        else:                                  iOther = other.iterator
        result = RealIterator_Dot ( self.cIterator, self.cData, iOther.cIterator, iOther.cData, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error in dot operation." )
        return result

    def Multiply ( self, other not None ):
        """Multiply two arrays."""
        cdef RealIterator iOther
        cdef CStatus      cStatus = CStatus_OK
        if isinstance ( other, RealIterator ): iOther = other
        else:                                  iOther = other.iterator
        RealIterator_Multiply ( self.cIterator, self.cData, iOther.cIterator, iOther.cData, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error in multiply operation." )

    def Swap ( self, other not None ):
        """Swap two arrays."""
        cdef RealIterator iOther
        cdef CStatus      cStatus = CStatus_OK
        if isinstance ( other, RealIterator ): iOther = other
        else:                                  iOther = other.iterator
        RealIterator_Swap ( self.cIterator, self.cData, iOther.cIterator, iOther.cData, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error in swap operation." )
