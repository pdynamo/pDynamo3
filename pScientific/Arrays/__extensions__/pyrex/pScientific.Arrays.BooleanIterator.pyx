"""Boolean iterator."""

from .ArrayError import ArrayError

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class BooleanIterator ( BaseIterator ):
    """Boolean iterator."""

    def __next__ ( self ):
        cdef CInteger i
        i = Iterator_NextIndex ( self.cIterator )
        if i >= 0: return self.cData[i]
        else:      raise StopIteration

    def _Initialize ( self ):
        """Initialization."""
        super ( BooleanIterator, self )._Initialize ( )
        self.cData = NULL

    cdef _SetDataPointer ( self ):
        """Set the data pointer."""
        cdef BooleanBlock block = self.block
        if ( block is not None ) and ( block.cObject != NULL ):
            self.cData = &(block.cObject.items[Iterator_DataOffSet ( self.cIterator )])

    # . Unary functions.
    def All ( self ):
        """Are all values of the array true?"""
        return ( BooleanIterator_All ( self.cIterator, self.cData ) == CTrue )

    def Any ( self ):
        """Are any values of the array true?"""
        return ( BooleanIterator_Any ( self.cIterator, self.cData ) == CTrue )

    def CountFalse ( self ):
        """Count the number of false values in the array."""
        return BooleanIterator_CountFalse ( self.cIterator, self.cData )

    def CountTrue ( self ):
        """Count the number of true values in the array."""
        return BooleanIterator_CountTrue ( self.cIterator, self.cData )

    def Not ( self ):
        """Negation."""
        cdef CStatus cStatus = CStatus_OK
        BooleanIterator_Not ( self.cIterator, self.cData, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error in NOT operation." )

    def Set ( self, value ):
        """Set the values of an array."""
        cdef CBoolean cValue
        cdef CStatus cStatus = CStatus_OK
        if value: cValue = CTrue
        else:     cValue = CFalse
        BooleanIterator_Set ( self.cIterator, self.cData, cValue, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error in set operation." )

    # . Binary functions.
    def And ( self, other not None ):
        """AND two arrays."""
        cdef BooleanIterator iOther
        cdef CStatus         cStatus = CStatus_OK
        if isinstance ( other, BooleanIterator ): iOther = other
        else:                                     iOther = other.iterator
        BooleanIterator_And ( self.cIterator, self.cData, iOther.cIterator, iOther.cData, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error in AND operation." )

    def CopyTo ( self, other not None ):
        """Copy |self| to |other|."""
        cdef BooleanIterator iOther
        cdef CStatus         cStatus = CStatus_OK
        if isinstance ( other, BooleanIterator ): iOther = other
        else:                                     iOther = other.iterator
        BooleanIterator_CopyTo ( self.cIterator, self.cData, iOther.cIterator, iOther.cData, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error in copy operation." )

    def Or ( self, other not None ):
        """OR two arrays."""
        cdef BooleanIterator iOther
        cdef CStatus         cStatus = CStatus_OK
        if isinstance ( other, BooleanIterator ): iOther = other
        else:                                     iOther = other.iterator
        BooleanIterator_Or ( self.cIterator, self.cData, iOther.cIterator, iOther.cData, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error in OR operation." )

    def Swap ( self, other not None ):
        """Swap two arrays."""
        cdef BooleanIterator iOther
        cdef CStatus         cStatus = CStatus_OK
        if isinstance ( other, BooleanIterator ): iOther = other
        else:                                     iOther = other.iterator
        BooleanIterator_Swap ( self.cIterator, self.cData, iOther.cIterator, iOther.cData, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error in swap operation." )

    def Xor ( self, other not None ):
        """XOR two arrays."""
        cdef BooleanIterator iOther
        cdef CStatus         cStatus = CStatus_OK
        if isinstance ( other, BooleanIterator ): iOther = other
        else:                                     iOther = other.iterator
        BooleanIterator_Xor ( self.cIterator, self.cData, iOther.cIterator, iOther.cData, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error in XOR operation." )
