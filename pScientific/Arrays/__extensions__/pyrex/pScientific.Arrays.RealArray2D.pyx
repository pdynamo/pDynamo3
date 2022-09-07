"""Real 2-D arrays."""

from .ArrayError   import ArrayError
from .RealIterator import RealIterator

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class RealArray2D ( BaseArray2D ):

    cdef void _AllocateCObject ( self ):
        """Allocate the C object."""
        cdef CStatus cStatus = CStatus_OK
        if self.cObject == NULL:
            self.cObject = RealArray2D_Allocate ( &cStatus )
            if cStatus != CStatus_OK: raise ArrayError ( "Error allocating C object." )
            self.cView = < CView2D * > self.cObject

    def _AssignBlock ( self, RealBlock block not None ):
        """Assign a block to the array."""
        cdef CStatus cStatus = CStatus_OK
        self.block = block
        RealArray2D_AssignBlock ( self.cObject, block.cObject, CTrue, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error assigning block to array." )

    cdef void _DeallocateCObject ( self ):
        """Deallocate the C object."""
        RealArray2D_Deallocate ( &self.cObject )

    def _GetMemoryBlock ( self, size ):
        """Get the memory block."""
        return RealBlock.WithCapacity ( size )

    def _GetRawArray1D ( self ):
        """Get a raw 1-D array."""
        return RealArray1D.RawWithCObject ( )

    def _GetRawIterator ( self ):
        """Get the raw iterator."""
        return RealIterator.Raw ( )

    def _Initialize ( self ):
        """Initialization."""
        super ( RealArray2D, self )._Initialize ( )
        self.cObject = NULL

    def _ScalarGet ( self ):
        """Get a scalar."""
        cdef CReal   cValue
        cdef CStatus cStatus = CStatus_OK
        cValue = RealArray2D_GetItemMultiSlice ( self.cObject, self.cMultiSlice, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error getting item." )
        return cValue

    def _ScalarSet ( self, value ):
        """Set a scalar."""
        cdef CStatus cStatus = CStatus_OK
        RealArray2D_SetItemMultiSlice ( self.cObject, self.cMultiSlice, value, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error setting item." )

    def CopyTo ( RealArray2D self, RealArray2D other not None ):
        """Copying."""
        cdef CStatus cStatus = CStatus_OK
        RealArray2D_CopyTo ( self.cObject, other.cObject, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error copying array." )

    def DiagonalOfProduct ( self, RealArray2D other    not None  ,
                                  RealArray1D diagonal not None  ,
                                              aTranspose = False ,
                                              bTranspose = False ,
                                              initialize = True  ):
        """The diagonal of the product of A and B with optional transposition."""
        cdef CBoolean cATranspose
        cdef CBoolean cBTranspose
        cdef CStatus  cStatus = CStatus_OK
        if aTranspose: cATranspose = CTrue
        else:          cATranspose = CFalse
        if bTranspose: cBTranspose = CTrue
        else:          cBTranspose = CFalse
        if initialize: diagonal.Set ( 0.0 )
        RealArray2D_DiagonalOfProduct ( self.cObject     ,
                                        cATranspose      ,
                                        other.cObject    ,
                                        cBTranspose      ,
                                        diagonal.cObject ,
                                        &cStatus         )
        if cStatus != CStatus_OK: raise ArrayError ( "Non-conformable arrays in diagonal of product." )

    def MatrixMultiply ( self, RealArray2D x, RealArray2D y, alpha = 1.0, beta = 0.0, xTranspose = False, yTranspose = False ):
        """Calculate self -> beta * self + alpha * x * y. x and y can be transposed."""
        cdef CBoolean cXTranspose
        cdef CBoolean cYTranspose
        cdef CStatus  cStatus = CStatus_OK
        if xTranspose: cXTranspose = CTrue
        else:          cXTranspose = CFalse
        if yTranspose: cYTranspose = CTrue
        else:          cYTranspose = CFalse
        RealArray2D_MatrixMultiply ( cXTranspose, cYTranspose, alpha, x.cObject, y.cObject, beta, self.cObject, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error multiplying matrix by matrix." )

    def OrthogonalizeColumns ( self, RealArray1D constants = None ):
        """Orthogonalize the columns of an array."""
        cdef CInteger      n
        cdef CRealArray1D *cConstants
        cdef CStatus       cStatus = CStatus_OK
        if constants is None: cConstants = NULL
        else:                 cConstants = constants.cObject
        n = RealArray2D_GramSchmidtOrthogonalize ( self.cObject, cConstants, NULL, NULL, NULL, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error orthogonalizing array columns." )
        return n

    def ProjectOutOfArray ( self, RealArray1D vector ):
        """Project the array out of another one."""
        cdef CStatus cStatus = CStatus_OK
        RealArray2D_ProjectOutOfArray1D ( self.cObject, vector.cObject, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error projecting out array." )

    def Set ( RealArray2D self, value ):
        """Setting."""
        RealArray2D_Set ( self.cObject, value )

    def TraceOfProduct ( self, RealArray2D other not None ):
        """Trace of product."""
        cdef CReal   value
        cdef CStatus cStatus = CStatus_OK
        value = RealArray2D_TraceOfProduct ( self.cObject, other.cObject, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error determining trace of product." )
        return value

    def Transpose ( self ):
        """Transpose a square array."""
        cdef CStatus cStatus = CStatus_OK
        RealArray2D_TransposeSquare ( self.cObject, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Unable to transpose a non-square array in-place." )

    def VectorMultiply ( self, RealArray1D x not None        ,
                               RealArray1D y not None        ,
                               CReal       alpha     = 1.0   ,
                               CReal       beta      = 0.0   ,
                                           transpose = False ):
        """Calculate y -> beta * y + alpha * self * x. self can be transposed."""
        cdef CBoolean aTranpose
        cdef CStatus  cStatus = CStatus_OK
        if transpose: aTranspose = CTrue
        else:         aTranspose = CFalse
        RealArray2D_VectorMultiply ( aTranspose, alpha, self.cObject, x.cObject, beta, y.cObject, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error multiplying matrix by vector." )
