from pCore.CPrimitiveTypes              cimport CBoolean         , \
                                                CFalse           , \
                                                CInteger         , \
                                                CReal            , \
                                                CTrue
from pCore.Status                       cimport CStatus          , \
                                                CStatus_OK
from pScientific.Arrays.IntegerArray1D  cimport CIntegerArray1D  , \
                                                IntegerArray1D
from pScientific.Arrays.RealArray1D     cimport CRealArray1D     , \
                                                RealArray1D
from pScientific.Arrays.RealArray2D     cimport CRealArray2D     , \
                                                RealArray2D
from pScientific.Arrays.SymmetricMatrix cimport CSymmetricMatrix , \
                                                SymmetricMatrix

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "DenseDeterminants.h":

    # . Functions.
    cdef void  Matrix_LUPFactorizationInPlace       ( CRealArray2D     *self              ,
                                                      CIntegerArray1D  *pivots            ,
                                                      CStatus          *status            )
    cdef CReal SquareMatrix_Determinant             ( CRealArray2D     *self              ,
                                                      CStatus          *status            )

cdef extern from "DenseEigenvalueSolvers.h":

    # . Functions.
    cdef void  SymmetricMatrix_EigenvaluesSolve     ( CSymmetricMatrix *self              ,
                                                      CBoolean          preserveInput     ,
                                                      CInteger          lower             ,
                                                      CInteger          upper             ,
                                                      CRealArray1D     *eigenValues       ,
                                                      CRealArray2D     *eigenVectors      ,
                                                      CBoolean          isColumnMajor     ,
                                                      CStatus          *status            )

cdef extern from "DenseLinearEquationSolvers.h":

    # . Functions.
    cdef void  SquareMatrix_LinearEquationsSolve    ( CRealArray2D     *self              ,
                                                      CRealArray1D     *rhs               ,
                                                      CBoolean          preserveInput     ,
                                                      CRealArray1D     *solution          ,
                                                      CStatus          *status            )
    cdef void  SymmetricMatrix_LinearEquationsSolve ( CSymmetricMatrix *self              ,
                                                      CRealArray1D     *rhs               ,
                                                      CRealArray1D     *solution          ,
                                                      CStatus          *status            )

cdef extern from "DenseMatrixPower.h":

    # . Functions.
    cdef void Matrix_PseudoInverse                  ( CRealArray2D     *self              ,
                                                      CBoolean          preserveInput     ,
                                                      CReal            *relativeTolerance ,
                                                      CRealArray2D     *inverse           ,
                                                      CInteger         *rank              ,
                                                      CReal            *condition         ,
                                                      CStatus          *status            )
    cdef void SymmetricMatrix_InversePower          ( CSymmetricMatrix *self              ,
                                                      CBoolean          preserveInput     ,
                                                      CReal             power             ,
                                                      CReal             tolerance         ,
                                                      CSymmetricMatrix *result            ,
                                                      CStatus          *status            )
    cdef void SymmetricMatrix_Power                 ( CSymmetricMatrix *self              ,
                                                      CBoolean          preserveInput     ,
                                                      CReal             power             ,
                                                      CSymmetricMatrix *result            ,
                                                      CStatus          *status            )

#===================================================================================================================================
