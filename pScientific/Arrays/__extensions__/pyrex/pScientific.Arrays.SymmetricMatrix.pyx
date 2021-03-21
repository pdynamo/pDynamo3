"""Symmetric matrices."""

from  pCore      import RawObjectConstructor
from .ArrayError import ArrayError
from .Slicing    import ProcessIntegerSlice

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . These data must correspond to the enums in the C source.
# . Updating options.
_UpdatingOptions = { "BFGS"   : SymmetricMatrixUpdating_BFGS   ,
                     "BOFILL" : SymmetricMatrixUpdating_Bofill ,
                     "MS"     : SymmetricMatrixUpdating_MS     ,
                     "POWELL" : SymmetricMatrixUpdating_Powell }

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class SymmetricMatrix:

    # . Methods.
    def __copy__ ( self ):
        """Copying."""
        cdef SymmetricMatrix clone
        cdef CStatus         cStatus = CStatus_OK
        clone         = self.__class__.Raw ( )
        clone.block   = self.block
        clone.cObject = SymmetricMatrix_CloneShallow ( self.cObject, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error cloning array." )
        return clone

    def __dealloc__ ( self ):
        """Finalization."""
        self.block     = None
        self._iterator = None
        SymmetricMatrix_Deallocate ( &self.cObject )

    def __deepcopy__ ( self, memo ):
        """Copying."""
        clone = self.__class__.WithExtent ( self.extent )
        self.CopyTo ( clone )
        return clone

    def __getattr__ ( self, name ):
        """Handle unknown attributes."""
        return getattr ( self.iterator, name )

    def __getitem__ ( self, indices ):
        """Get an item."""
        cdef CInteger i, j
        cdef CReal    value
        ( i, j ) = ProcessIntegerSlice ( self.shape, indices )
        value    = SymmetricMatrix_GetItem ( self.cObject, i, j, NULL )
        return value

    def __getstate__ ( self ):
        """Return the state."""
        return { "block"  : self.block  ,
                 "extent" : self.extent }

    def __init__ ( self, extent ):
        """Constructor with extent."""
        self._Initialize ( )
        self._Allocate ( extent )

    def __iter__ ( self ):
        return self._MakeIterator ( )

    def __len__ ( self ):
        """Return the size of the symmetricmatrix."""
        return self.size

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setitem__ ( self, indices, CReal value ):
        """Set an item."""
        cdef CInteger i, j
        ( i, j ) = ProcessIntegerSlice ( self.shape, indices )
        SymmetricMatrix_SetItem ( self.cObject, i, j, value, NULL )

    def __setstate__ ( self, state ):
        """Set the state."""
        self._Allocate ( state["extent"], block = state["block"] )

    def _Allocate ( self, CInteger extent, RealBlock block = None ):
        """Allocation."""
        cdef CRealBlock *cBlock  = NULL
        cdef CStatus     cStatus = CStatus_OK
        if block is     None:  block = RealBlock.WithCapacity ( SymmetricMatrix_ViewSize ( extent ) )
        if block is not None: cBlock = block.cObject
        self.block   = block
        self.cObject = SymmetricMatrix_FromExtentBlock ( extent, cBlock, CTrue, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error allocating symmetric matrix." )

    def _Initialize ( self ):
        """Initialization."""
        self.cObject   = NULL
        self.block     = None
        self._iterator = None

    def _MakeIterator ( self ):
        """Make the default array iterator."""
        cdef RealIterator iterator
        cdef CIterator   *cIterator = NULL
        cdef CStatus      cStatus   = CStatus_OK
        iterator           = RealIterator.Raw ( )
        iterator.cIterator = SymmetricMatrix_MakeIterator ( self.cObject, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Unable to make iterator." )
        iterator.block     = self.block
        iterator._SetDataPointer ( )
        iterator.Reset ( )
        return iterator

    def AbsoluteMaximum ( self ):
        """Return the maximum absolute value in the matrix."""
        return SymmetricMatrix_AbsoluteMaximum ( self.cObject )

    def Add ( self, SymmetricMatrix other, CReal scale = 1.0 ):
        """Add a scaled matrix."""
        cdef CStatus cStatus = CStatus_OK
        SymmetricMatrix_Add ( self.cObject, scale, other.cObject, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error adding symmetric matrices." )

    def CopyTo ( self, SymmetricMatrix other ):
        """Copying."""
        cdef CStatus cStatus = CStatus_OK
        SymmetricMatrix_CopyTo ( self.cObject, other.cObject, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error copying symmetric matrices." )

    def DiagonalOfProduct ( self, SymmetricMatrix other  not None ,
                                  RealArray1D     result not None ):
        """Diagonal of the product of two symmetric matrices."""
        cdef CStatus cStatus = CStatus_OK
        SymmetricMatrix_DiagonalOfProduct ( self.cObject, other.cObject, result.cObject, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error taking diagonal of two symmetric matrices." )

    def DiagonalOfTransform ( self, RealArray2D matrix not None ,
                                    RealArray1D result not None , useTranspose = False ):
        """Take the diagonal of the transform of the matrix, either X^T * S * X or X * S * X^T."""
        cdef CBoolean cUseTranspose
        cdef CStatus  cStatus = CStatus_OK
        if useTranspose: cUseTranspose = CTrue
        else:            cUseTranspose = CFalse
        SymmetricMatrix_DiagonalOfTransform ( self.cObject, matrix.cObject, cUseTranspose, result.cObject, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error transforming symmetric matrix." )

    def FromSquare ( self, RealArray2D square not None ):
        """Copy from a square matrix."""
        cdef CStatus cStatus = CStatus_OK
        SymmetricMatrix_CopyFromRealArray2D ( self.cObject, square.cObject, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error copying from square matrix." )

    def MakeFromEigenSystem ( self, n, RealArray1D d not None, RealArray2D v not None, initialize = True ):
        """Make given a set of diagonal values and vectors."""
        cdef CBoolean cInitialize
        cdef CStatus  cStatus = CStatus_OK
        if initialize: cInitialize = CTrue
        else:          cInitialize = CFalse
        SymmetricMatrix_MakeFromEigensystem ( self.cObject, cInitialize, n, d.cObject, v.cObject, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error making symmetric matrix from eigensystem." )

    def MatrixMultiply ( self, SymmetricMatrix other not None, RealArray2D result not None ):
        """Multiply two symmetric matrices."""
        cdef CStatus cStatus = CStatus_OK
        SymmetricMatrix_SymmetricMatrixMultiply ( self.cObject, other.cObject, result.cObject, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error multiplying two symmetric matrices." )

    def PostMatrixMultiply ( self, RealArray2D x not None, RealArray2D y not None, xTranspose = False ):
        """Calculate y = self * x. x can be transposed."""
        cdef CBoolean cXTranspose
        cdef CStatus  cStatus      = CStatus_OK
        if xTranspose: cXTranspose = CTrue
        else:          cXTranspose = CFalse
        SymmetricMatrix_PostMatrixMultiply ( self.cObject, x.cObject, cXTranspose, y.cObject, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error multiplying symmetric matrix by matrix." )

    def PreMatrixMultiply ( self, RealArray2D x not None, RealArray2D y not None, xTranspose = False ):
        """Calculate y = x * self. x can be transposed."""
        cdef CBoolean cXTranspose
        cdef CStatus  cStatus      = CStatus_OK
        if xTranspose: cXTranspose = CTrue
        else:          cXTranspose = CFalse
        SymmetricMatrix_PreMatrixMultiply ( self.cObject, x.cObject, cXTranspose, y.cObject, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error multiplying matrix by symmetric matrix." )

    def ProjectOutVectors ( self, RealArray2D vectors ):
        """Project a set of vectors from the matrix."""
        cdef SymmetricMatrix projected
        cdef CStatus         cStatus = CStatus_OK
        projected = self.__class__.WithExtent ( self.extent )
        SymmetricMatrix_ProjectOut ( self.cObject, vectors.cObject, projected.cObject, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error projecting out vectors from symmetric matrix." )
        return projected

    def Raise ( self, RealArray2D vectors, CReal value ):
        """Make 'vectors' eigenvectors of the matrix with eigenvalues 'value'."""
        cdef CStatus cStatus = CStatus_OK
        SymmetricMatrix_Raise ( self.cObject, vectors.cObject, value, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error raising vectors from symmetric matrix." )

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def RootMeanSquare ( self ):
        """Determine the RMS value of the elements."""
        return SymmetricMatrix_RootMeanSquare ( self.cObject )

    def Scale ( self, CReal value ):
        """Scale all the elements of a symmetricmatrix."""
        SymmetricMatrix_Scale ( self.cObject, value )

    def Set ( self, CReal value ):
        """Set all the elements of a symmetricmatrix."""
        SymmetricMatrix_Set ( self.cObject, value )

    def SetColumn ( self, CInteger  column, RealArray1D vector ):
        """Set a column from a vector."""
        cdef CStatus cStatus = CStatus_OK
        SymmetricMatrix_SetColumn ( self.cObject, column, vector.cObject, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error setting symmetric matrix column." )

    def SumDifference ( self, SymmetricMatrix other not None ):
        """In-place sum and difference of two matrices."""
        cdef CStatus cStatus = CStatus_OK
        SymmetricMatrix_SumDifference ( self.cObject, other.cObject, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error in symmetric matrix sum/difference." )

    def SymmetricTransform ( self, SymmetricMatrix matrix not None ,
                                   SymmetricMatrix result not None ):
        """Transform the matrix by X * S * X."""
        cdef CStatus  cStatus = CStatus_OK
        SymmetricMatrix_SymmetricTransform ( self.cObject, matrix.cObject, result.cObject, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error transforming symmetric matrix." )

    def ToSquare ( self, RealArray2D square not None ):
        """Copy to a square matrix."""
        cdef CStatus cStatus = CStatus_OK
        SymmetricMatrix_CopyToRealArray2D ( self.cObject, square.cObject, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error copying to square matrix." )

    def Trace ( self ):
        """Trace."""
        return SymmetricMatrix_Trace ( self.cObject )

    def TraceOfProduct ( self, SymmetricMatrix other not None ):
        """Trace of product."""
        cdef CReal   value
        cdef CStatus cStatus = CStatus_OK
        value = SymmetricMatrix_TraceOfProduct ( self.cObject, other.cObject, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error determining trace of product." )
        return value

    def Transform ( self, RealArray2D     matrix not None ,
                          SymmetricMatrix result not None , useTranspose = False ):
        """Transform the matrix by X^T * S * X or by X * S * X^T."""
        cdef CBoolean cUseTranspose
        cdef CStatus  cStatus = CStatus_OK
        if useTranspose: cUseTranspose = CTrue
        else:            cUseTranspose = CFalse
        SymmetricMatrix_Transform ( self.cObject, matrix.cObject, cUseTranspose, result.cObject, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error transforming symmetric matrix." )

    def Update ( self, RealArray1D dx not None, RealArray1D dg not None, option = "BFGS", tolerance = None ):
        """Update the matrix."""
        cdef CReal    cTolerance
        cdef CReal   *pTolerance = NULL
        cdef CStatus  cStatus    = CStatus_OK
        cdef CSymmetricMatrixUpdating_Option formula
        formula = _UpdatingOptions.get ( option.upper ( ), SymmetricMatrixUpdating_BFGS )
        if tolerance is not None:
            cTolerance =   tolerance
            pTolerance = &cTolerance
        SymmetricMatrix_Update ( self.cObject, dx.cObject, dg.cObject, formula, pTolerance, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error updating symmetric matrix." )

    def VectorMultiply ( self, RealArray1D other, RealArray1D result ):
        """Multiply by a vector."""
        cdef CStatus cStatus = CStatus_OK
        SymmetricMatrix_VectorMultiply ( self.cObject, other.cObject, result.cObject, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error multiplying symmetric matrix by a vector." )

    @classmethod
    def WithExtent ( selfClass, extent ):
        """Constructor with extent."""
        return selfClass ( extent )

    # . Properties.
    @property
    def columns ( self ): return self.extent
    @property
    def extent ( self ):
        if self.cObject == NULL: return 0
        else:                    return self.cObject.extent
    @property
    def isSquare ( self ): return True
    @property
    def iterator ( self ):
        if self._iterator is None:
            self._iterator = self._MakeIterator ( )
        return self._iterator
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
