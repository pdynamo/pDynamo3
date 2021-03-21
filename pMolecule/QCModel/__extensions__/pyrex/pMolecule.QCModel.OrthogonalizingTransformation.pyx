"""Orthogonalizing transformation functions."""

from .QCModelError import QCModelError

# . Inverse is just Y = S * X.

#===================================================================================================================================
# . Functions.
#===================================================================================================================================
# . The eigenvalues and eigenvectors are also returned.
def OrthogonalizingTransformation_Make ( SymmetricMatrix integrals not None, doCanonical         = False ,
                                                                             preserveInput       = True  ,
                                                                             eigenValueTolerance = None  ):
    """Calculate an orthogonalizing transformation and its inverse given a symmetric matrix of two-body integrals."""
    cdef RealArray1D eigenValues
    cdef RealArray2D eigenVectors
    cdef RealArray2D X
    cdef CBoolean    cDoCanonical
    cdef CBoolean    cPreserveInput
    cdef CInteger    d
    cdef CInteger    N
    cdef CReal       cTolerance    = 0.0
    cdef CReal      *pTolerance    = NULL
    cdef CStatus     cStatus       = CStatus_OK
    if doCanonical  : cDoCanonical   = CTrue
    else:             cDoCanonical   = CFalse
    if preserveInput: cPreserveInput = CTrue
    else:             cPreserveInput = CFalse
    if eigenValueTolerance is not None:
        cTolerance = eigenValueTolerance
        pTolerance = &cTolerance
    d            = integrals.rows
    eigenValues  = RealArray1D.WithExtent  ( d    )
    eigenVectors = RealArray2D.WithExtents ( d, d )
    X            = RealArray2D.WithExtents ( d, d )
    N = COrthogonalizingTransformation ( integrals.cObject    ,
                                         cDoCanonical         ,
                                         cPreserveInput       ,
                                         pTolerance           ,
                                         eigenValues.cObject  ,
                                         eigenVectors.cObject ,
                                         X.cObject            ,
                                         &cStatus             )
    if cStatus != CStatus_OK: raise QCModelError ( "Error making the orthogonalizing transformation ({:d}/{:d},{:d}).".format ( N, d, cStatus ) )
    if N < d: return ( X[0:d,0:N], eigenValues, eigenVectors )
    else:     return ( X         , eigenValues, eigenVectors )
