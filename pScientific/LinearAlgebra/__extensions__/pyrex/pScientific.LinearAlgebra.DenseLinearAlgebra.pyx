"""Miscellaneous dense linear algebra operations."""

from   pCore              import logFile            , \
                                 LogFileActive
from ..Arrays             import Array
from  .LinearAlgebraError import LinearAlgebraError

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . SVD relative tolerance.
_SVDTolerance = 1.0e-12

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
def LinearLeastSquaresBySVD ( RealArray2D matrix not None                   ,
                                          rhs    not None                   ,
                                          preserveInput     = True          ,
                                          relativeTolerance = _SVDTolerance ,
                                          solution          = None          ):
    """Solve a linear least squares problem using SVD with a vector or matrix RHS."""
    # . No check is made for overdetermined or underdetermined systems.
    # . Get the inverse.
    inverse = Array.WithExtents ( matrix.columns, matrix.rows )
    report  = MatrixPseudoInverse ( matrix, inverse, preserveInput = True, relativeTolerance = relativeTolerance )
    # . Set the solution array.
    if isinstance ( rhs, RealArray1D ):
        if solution is None: solution = Array.WithExtent ( matrix.columns )
        inverse.VectorMultiply ( rhs, solution )
    else:
        if solution is None: solution = Array.WithExtents ( matrix.columns, rhs.columns )
        solution.MatrixMultiply ( inverse, rhs )
    report["Solution"] = solution
    return report

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

#-----------------------------------------------------------------------------------------------------------------------------------
def MatrixPseudoInverse ( RealArray2D matrix  not None            ,
                          RealArray2D inverse not None            ,
                                      preserveInput     = True    ,
                                      relativeTolerance = 1.0e-12 ):
    """Find the pseudo-inverse of a matrix."""
    cdef CBoolean cPreserveInput
    cdef CInteger cRank
    cdef CReal    cCondition
    cdef CReal    cRelativeTolerance
    cdef CStatus  cStatus = CStatus_OK
    if preserveInput   : cPreserveInput = CTrue
    else:                cPreserveInput = CFalse
    cRelativeTolerance = relativeTolerance
    Matrix_PseudoInverse ( matrix.cObject      ,
                          cPreserveInput      ,
                          &cRelativeTolerance ,
                          inverse.cObject     ,
                          &cRank              ,
                          &cCondition         ,
                          &cStatus            )
    if cStatus != CStatus_OK: raise LinearAlgebraError ( "Error finding matrix pseudo-inverse." )
    return { "Rank" : cRank, "Condition Number" : cCondition }

#===================================================================================================================================
