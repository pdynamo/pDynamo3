"""Regular grid occupancy."""

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class RegularGridOccupancy:
    """Regular grid occupancy."""

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            RegularGridOccupancy_Deallocate ( &self.cObject )
            self.isOwner = False

    def __init__ ( self ):
        """Constructor."""
        self._Initialize ( )

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = False

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self
