"""Slicing functions."""

from .ArrayError import ArrayError

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class MultiSlice:

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            MultiSlice_Deallocate ( &self.cObject )
            self.isOwner = False

    def __init__ ( self, capacity ):
        """Constructor with capacity."""
        self._Initialize ( )
        self._Allocate ( capacity )

    def __len__ ( MultiSlice self ):
        """Return the capacity of the symmetricmatrix."""
        return self.capacity

    def _Allocate ( self, CInteger capacity ):
        """Allocation."""
        cdef CStatus cStatus = CStatus_OK
        self.cObject = MultiSlice_Allocate ( capacity, &cStatus )
        self.isOwner = True
        if cStatus != CStatus_OK: raise ArrayError ( "Error allocating multislice." )

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = False

    @classmethod
    def WithCapacity ( selfClass, capacity ):
        """Constructor with capacity."""
        return selfClass ( capacity )

    # . Properties.
    @property
    def capacity ( self ):
        return MultiSlice_GetCapacity ( self.cObject )
    @property
    def rank ( self ):
        return MultiSlice_GetRank ( self.cObject )
    @property
    def size ( self ):
        return MultiSlice_GetSize ( self.cObject )

#===================================================================================================================================
# . Functions for slice processing.
#===================================================================================================================================
# . No checking done when setting multiSlice as already done.
# . Omitting trailing slice indices is allowed - in this case they are set to the full rank.
cdef CInteger CProcessMultiSlice ( shape, indices, CMultiSlice *multiSlice ):
    """Process a general N-D slice."""
    if not isinstance ( indices, tuple ): indices = ( indices, )
    nI = len ( indices )
    nS = len ( shape   )
    if nS < nI: raise ArrayError ( "Invalid N-D slice rank." )
    for ( d, ( index, extent ) ) in enumerate ( zip ( indices, shape ) ):
        ( isScalar, start, stop, stride, n ) = ProcessSlice1D ( index, extent )
        MultiSlice_SetSliceIndices ( multiSlice, d, ( CTrue if isScalar else CFalse ), start, stop, stride, n )
    if nS > nI:
        for ( d, extent ) in enumerate ( shape[nI:nS] ):
            MultiSlice_SetSliceIndices ( multiSlice, d+nI, CFalse, 0, extent, 1, extent )
    MultiSlice_SetRank ( multiSlice )
    return MultiSlice_GetRank ( multiSlice )

def ProcessIntegerSlice ( shape, indices ):
    """Process an integer N-D slice."""
    if not isinstance ( indices, tuple ): indices = ( indices, )
    if len ( shape ) != len ( indices ): raise ArrayError ( "Invalid N-D slice rank." ) # . Trailing slices not allowed.
    indexes = []
    for ( d, ( index, extent ) ) in enumerate ( zip ( indices, shape ) ):
        if isinstance ( index, int ):
            if index < 0: index += extent
            if ( index < 0 ) or ( index >= extent ): raise IndexError ( "Index out of range." ) # . Needed for Python iteration.
            indexes.append ( index )
        else: raise ArrayError ( "Invalid integer slice." )
    return indexes

def ProcessSlice1D ( index, extent ):
    """Process an arbitrary single slice argument."""
    isOK = True
    if isinstance ( index, int ):
        isScalar = True
        if index < 0: index += extent
        if ( index < 0 ) or ( index >= extent ): raise IndexError ( "Index out of range." ) # . Needed for Python iteration.
        start    = index
        stop     = index + 1
        stride   = 1
        n        = 1
    elif isinstance ( index, slice ):
        isScalar = False
        try:
            ( start, stop, stride ) = index.indices ( extent )
            if   ( stride > 0 ) and ( stop > start ): n = ( stop - start - 1 ) // stride + 1
            elif ( stride < 0 ) and ( stop < start ): n = ( stop - start + 1 ) // stride + 1
            else: n = 0
        except:
            isOK = False
    if not isOK: raise ArrayError ( "Invalid slice argument." )
    return ( isScalar, start, stop, stride, n )
