"""Miscellaneous dense linear algebra operations."""

from  pCore              import logFile            , \
                                LogFileActive
from .LinearAlgebraError import LinearAlgebraError

#===================================================================================================================================
# . Functions.
#===================================================================================================================================
def Determinant ( RealArray2D matrix not None ):
    """Find the determinant of a matrix."""
    cdef CReal    determinant
    cdef CStatus  cStatus = CStatus_OK
    determinant = SquareMatrix_Determinant ( matrix.cObject, &cStatus )
    if cStatus != CStatus_OK: raise LinearAlgebraError ( "Determinant calculation error." )
    return determinant

#-----------------------------------------------------------------------------------------------------------------------------------
def EigenPairs ( SymmetricMatrix matrix, RealArray1D eigenValues, RealArray2D eigenVectors, columnMajor = False, lower = -1, preserveInput = True, upper = -1 ):
    """Find the eigenvalues and eigenvectors of a matrix."""
    cdef CBoolean cColumnMajor
    cdef CBoolean cPreserveInput
    cdef CInteger cLower, cUpper
    cdef CStatus  cStatus = CStatus_OK
    if columnMajor  : cColumnMajor   = CTrue
    else:             cColumnMajor   = CFalse
    if preserveInput: cPreserveInput = CTrue
    else:             cPreserveInput = CFalse
    if lower < 0: cLower = 0
    else:         cLower = lower
    if upper < 0: cUpper = matrix.rows
    else:         cUpper = upper
    SymmetricMatrix_EigenvaluesSolve ( matrix.cObject       ,
                                       cPreserveInput       ,
                                       cLower               ,
                                       cUpper               ,
                                       eigenValues.cObject  ,
                                       eigenVectors.cObject ,
                                       cColumnMajor         ,
                                       &cStatus             )
    if cStatus != CStatus_OK: raise LinearAlgebraError ( "Eigenpairs solution error." )
    return { "Eigenvalues" : eigenValues, "Eigenvectors" : eigenVectors }

#-----------------------------------------------------------------------------------------------------------------------------------
def EigenValues ( SymmetricMatrix matrix, RealArray1D eigenValues, columnMajor = False, lower = -1, preserveInput = True, upper = -1 ):
    """Find the eigenvalues of a matrix."""
    cdef CBoolean cColumnMajor
    cdef CBoolean cPreserveInput
    cdef CInteger cLower, cUpper
    cdef CStatus  cStatus = CStatus_OK
    if columnMajor  : cColumnMajor   = CTrue
    else:             cColumnMajor   = CFalse
    if preserveInput: cPreserveInput = CTrue
    else:             cPreserveInput = CFalse
    if lower < 0: cLower = 0
    else:         cLower = lower
    if upper < 0: cUpper = matrix.rows
    else:         cUpper = upper
    SymmetricMatrix_EigenvaluesSolve ( matrix.cObject       ,
                                       cPreserveInput       ,
                                       cLower               ,
                                       cUpper               ,
                                       eigenValues.cObject  ,
                                       NULL                 ,
                                       cColumnMajor         ,
                                       &cStatus             )
    if cStatus != CStatus_OK: raise LinearAlgebraError ( "Eigenvalues solution error." )
    return eigenValues

#-----------------------------------------------------------------------------------------------------------------------------------
def LinearEquations ( matrix, RealArray1D rhs, preserveInput = True, RealArray1D solution = None ):
    """Solve a set of linear equations."""
    cdef RealArray2D     rArray2D
    cdef SymmetricMatrix rSymmetric
    cdef CBoolean        cPreserveInput
    cdef CStatus         cStatus        = CStatus_OK
    if preserveInput   : cPreserveInput = CTrue
    else:                cPreserveInput = CFalse
    if solution is None: solution       = rhs
    if isinstance ( matrix, RealArray2D ):
        rArray2D = matrix
        SquareMatrix_LinearEquationsSolve    ( rArray2D.cObject   ,
                                               rhs.cObject        ,
                                               cPreserveInput     ,
                                               solution.cObject   ,
                                               &cStatus           )
    elif isinstance ( matrix, SymmetricMatrix ):
        rSymmetric = matrix
        SymmetricMatrix_LinearEquationsSolve ( rSymmetric.cObject ,
                                               rhs.cObject        ,
                                               solution.cObject   ,
                                               &cStatus           )
    else: raise LinearAlgebraError ( "Invalid linear equations matrix type." )
    if cStatus != CStatus_OK: raise LinearAlgebraError ( "Linear equations solution error." )
    return solution

#-----------------------------------------------------------------------------------------------------------------------------------
def LinearLeastSquaresBySVD ( RealArray2D matrix, RealArray1D rhs, preserveInput = True, relativeTolerance = 1.0e-12, RealArray1D solution = None ):
    """Solve a linear least squares problem using SVD."""
    cdef CBoolean cPreserveInput
    cdef CInteger cRank
    cdef CReal    cCondition
    cdef CReal    cRelativeTolerance
    cdef CStatus  cStatus = CStatus_OK
    if preserveInput   : cPreserveInput = CTrue
    else:                cPreserveInput = CFalse
    if solution is None: solution       = rhs
    cRelativeTolerance = relativeTolerance
    LinearLeastSquaresSVDSolve ( matrix.cObject      ,
                                 rhs.cObject         ,
                                 cPreserveInput      ,
                                 &cRelativeTolerance ,
                                 solution.cObject    ,
                                 &cRank              ,
                                 &cCondition         ,
                                 &cStatus            )
    if cStatus != CStatus_OK: raise LinearAlgebraError ( "Linear least squares by SVD solution error." )
    return { "Solution" : solution, "Rank" : cRank, "Condition Number" : cCondition }

#-----------------------------------------------------------------------------------------------------------------------------------
def MatrixPower ( SymmetricMatrix matrix not None  ,
                  CReal           power            ,
                  SymmetricMatrix result not None  ,
                              preserveInput = True ):
    """Find the power of a matrix."""
    cdef CBoolean cPreserveInput
    cdef CStatus  cStatus = CStatus_OK
    if preserveInput: cPreserveInput = CTrue
    else:             cPreserveInput = CFalse
    SymmetricMatrix_Power ( matrix.cObject ,
                            cPreserveInput ,
                            power          ,
                            result.cObject ,
                            &cStatus       )
    if cStatus != CStatus_OK: raise LinearAlgebraError ( "Error making matrix power." )

#-----------------------------------------------------------------------------------------------------------------------------------
def MatrixPowerInverse ( SymmetricMatrix matrix not None  ,
                         CReal           power            , # . Must be positive!
                         CReal           tolerance        ,
                         SymmetricMatrix result not None  ,
                                     preserveInput = True ):
    """Find the power of a matrix."""
    cdef CBoolean cPreserveInput
    cdef CStatus  cStatus = CStatus_OK
    if preserveInput: cPreserveInput = CTrue
    else:             cPreserveInput = CFalse
    SymmetricMatrix_InversePower ( matrix.cObject ,
                                   cPreserveInput ,
                                   power          ,
                                   tolerance      ,
                                   result.cObject ,
                                   &cStatus       )
    if cStatus != CStatus_OK: raise LinearAlgebraError ( "Error making matrix inverse power." )

#===================================================================================================================================
