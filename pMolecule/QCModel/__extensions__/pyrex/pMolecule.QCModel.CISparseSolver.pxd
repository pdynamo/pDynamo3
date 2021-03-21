from pCore.CPrimitiveTypes                    cimport CBoolean                     , \
                                                      CFalse                       , \
                                                      CInteger                     , \
                                                      CReal                        , \
                                                      CTrue
from pScientific.Arrays.IntegerArray1D        cimport IntegerArray1D               , \
                                                      IntegerArray1D_PointerToData
from pScientific.Arrays.RealArray1D           cimport CRealArray1D                 , \
                                                      RealArray1D                  , \
                                                      RealArray1D_PointerToData
from pScientific.Arrays.RealArray2D           cimport RealArray2D                  , \
                                                      RealArray2D_PointerToData
from pScientific.Arrays.SparseSymmetricMatrix cimport SparseSymmetricMatrix        , \
                                                      SparseSymmetricMatrix_Print
from pScientific.LinearAlgebra.PrimmeSolver   cimport CPrimmeParams                , \
                                                      dprimme                      , \
                                                      DYNAMIC                      , \
                                                      primme_display_params        , \
                                                      primme_initialize            , \
                                                      primme_set_method

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "CISparseSolver.h":

    cdef void CISparseSolver_ApplyMatrix         ( void *xVoid, void *yVoid, CInteger *blockSize, CPrimmeParams *primme )
    cdef void CISparseSolver_ApplyPreconditioner ( void *xVoid, void *yVoid, CInteger *blockSize, CPrimmeParams *primme )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class CISparseSolver:

    cdef CPrimmeParams                cObject
    cdef public RealArray1D           eigenValues
    cdef public RealArray2D           eigenVectors
    cdef public IntegerArray1D        iWork
    cdef public SparseSymmetricMatrix matrix
    cdef public RealArray1D           preconditioner
    cdef public RealArray1D           residualNorms
    cdef public RealArray1D           rWork
