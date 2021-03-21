from pCore.CPrimitiveTypes          cimport CBoolean     , \
                                            CFalse       , \
                                            CInteger     , \
                                            CReal        , \
                                            CTrue
from pCore.Status                   cimport CStatus      , \
                                            CStatus_OK
from pScientific.Arrays.RealArray1D cimport CRealArray1D , \
                                            RealArray1D

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "SparseSymmetricMatrix.h":

    ctypedef struct CSparseSymmetricMatrix "SparseSymmetricMatrix":
        CInteger extent
        CInteger size

    ctypedef struct CSparseSymmetricMatrixRowItemIterator "SparseSymmetricMatrixRowItemIterator":
        pass

    cdef CSparseSymmetricMatrix *SparseSymmetricMatrix_Allocate                               ( CInteger                 extent                 ,
                                                                                                CInteger                 size                   ,
                                                                                                CStatus                 *status                 )
    cdef void                    SparseSymmetricMatrix_AppendItem                             ( CSparseSymmetricMatrix  *self                   ,
                                                                                                CInteger                 i                      ,
                                                                                                CInteger                 j                      ,
                                                                                                CReal                    value                  ,
                                                                                                CStatus                 *status                 )
    cdef void                    SparseSymmetricMatrix_ApplyIncompleteCholeskyDecomposition   ( CSparseSymmetricMatrix  *self                   ,
                                                                                                CRealArray1D            *b                      ,
                                                                                                CRealArray1D            *x                      )
    cdef void                    SparseSymmetricMatrix_Canonicalize                           ( CSparseSymmetricMatrix  *self                   ,
                                                                                                CStatus                 *status                 )
    cdef void                    SparseSymmetricMatrix_Clear                                  ( CSparseSymmetricMatrix  *self                   )
    cdef CSparseSymmetricMatrix *SparseSymmetricMatrix_Clone                                  ( CSparseSymmetricMatrix  *self                   ,
                                                                                                CStatus                 *status                 )
    cdef void                    SparseSymmetricMatrix_ComputeIncompleteCholeskyDecomposition ( CSparseSymmetricMatrix  *self                   ,
                                                                                                CReal                    alpha                  ,
                                                                                                CInteger                *numberOfModifiedPivots ,
                                                                                                CStatus                 *status                 )
    cdef void                    SparseSymmetricMatrix_CopyTo                                 ( CSparseSymmetricMatrix  *self                   ,
                                                                                                CSparseSymmetricMatrix  *other                  ,
                                                                                                CStatus                 *status                 )
    cdef void                    SparseSymmetricMatrix_Deallocate                             ( CSparseSymmetricMatrix **self                   )
    cdef void                    SparseSymmetricMatrix_GetDiagonal                            ( CSparseSymmetricMatrix  *self                   ,
                                                                                                CRealArray1D            *diagonal               ,
                                                                                                CStatus                 *status                 )
    cdef void                    SparseSymmetricMatrix_MakeDiagonalPreconditioner             ( CSparseSymmetricMatrix  *self                   ,
                                                                                                CRealArray1D            *preconditioner         ,
                                                                                                CReal                   *tolerance              ,
                                                                                                CStatus                 *status                 )
    cdef void                    SparseSymmetricMatrix_Print                                  ( CSparseSymmetricMatrix  *self                   )
    cdef void                    SparseSymmetricMatrix_VectorMultiply                         ( CSparseSymmetricMatrix  *self                   ,
                                                                                                CRealArray1D            *x                      ,
                                                                                                CRealArray1D            *y                      ,
                                                                                                CStatus                 *status                 )

    cdef void                    SparseSymmetricMatrixRowItemIterator_Initialize              ( CSparseSymmetricMatrixRowItemIterator *self     ,
                                                                                                CSparseSymmetricMatrix                *target   ,
                                                                                                CInteger                               row      ,
                                                                                                CStatus                               *status   )
    cdef CInteger                SparseSymmetricMatrixRowItemIterator_Next                    ( CSparseSymmetricMatrixRowItemIterator *self     )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class SparseSymmetricMatrix:

    cdef CSparseSymmetricMatrix *cObject
    cdef public object           isOwner
