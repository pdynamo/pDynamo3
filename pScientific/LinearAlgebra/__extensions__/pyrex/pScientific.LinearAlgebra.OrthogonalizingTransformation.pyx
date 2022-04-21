"""Orthogonalizing transformation functions."""

from  enum               import Enum
from .LinearAlgebraError import LinearAlgebraError

# . Inverse is just Y = S * X.

#===================================================================================================================================
# . Definitions.
#===================================================================================================================================
class OrthogonalizationMethod ( Enum ):
    """The orthogonalization method."""
    Canonical = 1
    Diagonal  = 2
    Symmetric = 3

#===================================================================================================================================
# . Functions.
#===================================================================================================================================
# . The eigenvalues and eigenvectors are also returned.
def OrthogonalizingTransformation_Make ( SymmetricMatrix integrals not None, orthogonalizationMethod     ,
                                                                             diagonalTolerance   = None  ,
                                                                             eigenValueTolerance = None  ,
                                                                             preserveInput       = True  ):
    """Calculate an orthogonalizing transformation and its inverse given a symmetric metric matrix."""
    cdef RealArray1D              eigenValues
    cdef RealArray2D              eigenVectors
    cdef RealArray2D              X
    cdef CBoolean                 cPreserveInput
    cdef CInteger                 d
    cdef CInteger                 N
    cdef COrthogonalizationMethod cMethod
    cdef CReal                    cDTolerance    = 0.0
    cdef CReal                    cETolerance    = 0.0
    cdef CReal                   *pDTolerance    = NULL
    cdef CReal                   *pETolerance    = NULL
    cdef CStatus                  cStatus        = CStatus_OK
    cMethod = orthogonalizationMethod.value
    if diagonalTolerance is not None:
        cDTolerance = diagonalTolerance
        pDTolerance = &cDTolerance
    if eigenValueTolerance is not None:
        cETolerance = eigenValueTolerance
        pETolerance = &cETolerance
    if preserveInput: cPreserveInput = CTrue
    else:             cPreserveInput = CFalse
    d            = integrals.rows
    eigenValues  = RealArray1D.WithExtent  ( d    )
    eigenVectors = RealArray2D.WithExtents ( d, d )
    X            = RealArray2D.WithExtents ( d, d )
    N = COrthogonalizingTransformation ( integrals.cObject    ,
                                         cMethod              ,
                                         cPreserveInput       ,
                                         pDTolerance          ,
                                         pETolerance          ,
                                         eigenValues.cObject  ,
                                         eigenVectors.cObject ,
                                         X.cObject            ,
                                         &cStatus             )
    if cStatus != CStatus_OK: raise LinearAlgebraError ( "Error making the orthogonalizing transformation ({:d}/{:d},{:d}).".format ( N, d, cStatus ) )
    if N < d: return ( X[0:d,0:N], eigenValues, eigenVectors )
    else:     return ( X         , eigenValues, eigenVectors )
