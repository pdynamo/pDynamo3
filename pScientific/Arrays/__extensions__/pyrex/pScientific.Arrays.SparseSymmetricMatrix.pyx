"""Sparse symmetric matrices."""

from  pCore      import DataType      , \
                        logFile       , \
                        LogFileActive
from .ArrayError import ArrayError

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class SparseSymmetricMatrix:

    # . Public methods.
    def __copy__ ( self ):
        """Copying."""
        cdef SparseSymmetricMatrix new
        new = self.__class__.WithExtentAndSize ( self.extent, self.size )
        self.CopyTo ( new )
        return new

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            SparseSymmetricMatrix_Deallocate ( &self.cObject )
            self.isOwner = False

    def __deepcopy__ ( self, memo ):
        """Copying."""
        return self.__copy__ ( )

#   def __getstate__ ( self ):
#       """Return the state."""
#       cdef CInteger  i, n
#       n     = self.extent
#       items = []
#       for i from 0 <= i < self.size:
#           items.append ( self.cObject.data[i] )
#       return { "items" : items, "shape" : [ n, n, n, n ], "storage" : "LowerTriangleRowMajor" }

    def __init__ ( self, extent, size ):
        """Constructor with extent and size."""
        self._Initialize ( )
        self._Allocate ( extent, size )

    def __len__ ( SparseSymmetricMatrix self ):
        """Return the size of the symmetricmatrix."""
        return self.size

#   def __reduce_ex__ ( self, protocol ):
#       """Pickling protocol."""
#       return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

#   def __setstate__ ( self, state ):
#       """Set the state."""
#       items = state["items"]
#       self._Allocate ( state["shape"][0] )
#       for i from 0 <= i < self.size:
#           self.cObject.data[i] = items[i]

    def _Allocate ( self, CInteger extent, CInteger size ):
        """Allocation."""
        cdef CStatus cStatus = CStatus_OK
        self.cObject = SparseSymmetricMatrix_Allocate ( extent, size, &cStatus )
        self.isOwner = True
        if cStatus != CStatus_OK: raise ArrayError ( "Error allocating sparse symmetric matrix." )

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = False

    def CopyTo ( self, SparseSymmetricMatrix other ):
        """Copying."""
        cdef CStatus cStatus = CStatus_OK
        SparseSymmetricMatrix_CopyTo ( self.cObject, other.cObject, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error copying sparse symmetric matrices." )

    def MakeDiagonalPreconditioner ( self, RealArray1D preconditioner not None, tolerance = None ):
        """Make a diagonal preconditioner given the matrix."""
        cdef CReal   cTolerance
        cdef CReal  *pTolerance = NULL
        cdef CStatus cStatus    = CStatus_OK
        if tolerance is not None:
            cTolerance =   tolerance
            pTolerance = &cTolerance
        SparseSymmetricMatrix_MakeDiagonalPreconditioner ( self.cObject, preconditioner.cObject, pTolerance, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error making diagonal preconditioner." )

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

#   def Set ( SparseSymmetricMatrix self, CReal value ):
#       """Set all the elements of a symmetricmatrix."""
#       SparseSymmetricMatrix_Set ( self.cObject, value )

    @classmethod
    def WithExtentAndSize ( selfClass, extent, size ):
        """Constructor with extent and size."""
        return selfClass ( extent, size )

    # . Properties.
    @property
    def columns ( self ): return self.extent
    @property
    def dataType ( self ): return DataType.Real
    @property
    def extent ( self ):
        if self.cObject == NULL: return 0
        else:                    return self.cObject.extent
    @property
    def isSquare ( self ): return True
    @property
    def rank ( self ): return 2
    @property
    def rows ( self ): return self.extent
    @property
    def shape ( self ): return [ self.extent, self.extent ]
    @property
    def size ( self ):
        if self.cObject == NULL: return 0
        else:                    return self.cObject.size
