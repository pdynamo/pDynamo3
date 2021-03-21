/*==================================================================================================================================
! . Functions for canonical and symmetric orthogonalization.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>

# include "DenseEigenvalueSolvers.h"
# include "Integer.h"
# include "OrthogonalizingTransformation.h"
# include "RealArray1D.h"

/*
!
! . Definitions: integrals = S, forward transformation = X, inverse = Y = S * X, such that Y^T * X = X^T * Y = I.
!

Code to do this is:

Y = RealArray2D_AllocateWithExtents ( SymmetricMatrix_Extent ( S ), View2D_Columns ( X ), status ) ;
SymmetricMatrix_PostMatrixMultiply ( S, X, False, Y, status ) ;

*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check eigenvectors.
!---------------------------------------------------------------------------------------------------------------------------------*/
/*# define _CheckEigenvectors*/
# ifdef _CheckEigenvectors
static void CheckEigenvectors ( const Integer d, RealArray1D *eigenValues, RealArray2D *eigenVectors )
{
    auto Integer      i ;
    auto Real         m ;
    auto RealArray2D *A = NULL ;
    printf ( "\nEigenvalues:\n" ) ;
    RealArray1D_Print ( eigenValues ) ;
    A = RealArray2D_AllocateWithExtents ( d, d, NULL ) ;
    RealArray2D_MatrixMultiply ( True, False, 1.0e+00, eigenVectors, eigenVectors, 0.0e+00, A, NULL ) ;
    for ( i = 0 ; i < d ; i++ ) Array2D_Item ( A, i, i ) = 0.0e+00 ;
    m = RealArray2D_AbsoluteMaximum ( A ) ;
    if ( m > 1.0e-06 )
    {
        printf ( "\nU^T * U:" ) ;
        RealArray2D_Print ( A ) ;
    }
    RealArray2D_MatrixMultiply ( False, True, 1.0e+00, eigenVectors, eigenVectors, 0.0e+00, A, NULL ) ;
    for ( i = 0 ; i < d ; i++ ) Array2D_Item ( A, i, i ) = 0.0e+00 ;
    m = RealArray2D_AbsoluteMaximum ( A ) ;
    if ( m > 1.0e-06 )
    {
        printf ( "\nU * U^T:" ) ;
        RealArray2D_Print ( A ) ;
    }
    RealArray2D_Deallocate ( &A ) ;
}
# endif

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check orthogonalization.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real CheckOrthogonalization ( RealArray2D *transformation ,
                              RealArray2D *inverse        ,
                              Status      *status         )
{
    Real deviation = 0.0e+00 ;
    if ( ( transformation != NULL ) && ( inverse != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer      i, i0, i1, t0, t1 ;
        auto RealArray2D *m ;
        i0 = View2D_Rows ( inverse        ) ; i1 = View2D_Columns ( inverse        ) ;
        t0 = View2D_Rows ( transformation ) ; t1 = View2D_Columns ( transformation ) ;
        if ( ( i0 != t0 ) || ( i1 != t1 ) ) { Status_Set ( status, Status_NonConformableArrays ) ; }
        else
        {
            m = RealArray2D_AllocateWithExtents ( i1, i1, status ) ;
            RealArray2D_MatrixMultiply ( True , False, 1.0e+00, inverse, transformation, 0.0e+00, m, NULL ) ; /* Y^T * X */
            for ( i = 0 ; i < i1 ; i++ ) Array2D_Item ( m, i, i ) -= 1.0e+00 ;
            deviation = RealArray2D_AbsoluteMaximum ( m ) ;
            RealArray2D_Deallocate ( &m ) ;
        }
    }
    return deviation ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . A canonical or symmetric orthogonalizing transformation, such that X^T * S * X = I.
! . X must be allocated on entry, and the number of linearly-independent vectors is returned.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _EigenValueTolerance 1.0e-10
Integer OrthogonalizingTransformation (       SymmetricMatrix *S                   ,
                                        const Boolean          doCanonical         ,
                                        const Boolean          preserveInput       ,
                                        const Real            *eigenValueTolerance ,
                                              RealArray1D     *eigenValues         ,
                                              RealArray2D     *eigenVectors        ,
                                              RealArray2D     *X                   ,
                                              Status          *status              )
{
    Integer numberOfVectors = -1 ;
    if ( ( S != NULL ) && ( X != NULL ) && ( Status_IsOK ( status ) ) )
    {
        auto Boolean isOK ;
        auto Integer d = SymmetricMatrix_Extent ( S ) ;
        isOK = ( View2D_Columns ( X ) == d ) && ( View2D_Rows ( X ) == d ) ;
        if ( eigenValues  != NULL ) isOK = isOK && ( View1D_Extent  ( eigenValues  ) == d ) ;
        if ( eigenVectors != NULL ) isOK = isOK && ( View2D_Columns ( eigenVectors ) == d )
                                                && ( View2D_Rows    ( eigenVectors ) == d ) ;
        if ( ! isOK ) { Status_Set ( status, Status_NonConformableArrays ) ; }
        else
        {
            auto Real         eTolerance ;
            auto RealArray1D *eigenValuesL  = NULL, *inverseEigenValues = NULL ;
            auto RealArray2D *eigenVectorsL = NULL, *scaledEigenVectors = NULL ;
            if ( eigenValueTolerance == NULL ) eTolerance = _EigenValueTolerance   ;
            else                               eTolerance = (*eigenValueTolerance) ;
            inverseEigenValues = RealArray1D_AllocateWithExtent ( d, status ) ;
            if ( eigenValues  == NULL ) eigenValuesL       = RealArray1D_AllocateWithExtent  ( d   , status ) ;
            else                        eigenValuesL       = eigenValues ;
            if ( eigenVectors == NULL ) eigenVectorsL      = RealArray2D_AllocateWithExtents ( d, d, status ) ;
            else                        eigenVectorsL      = eigenVectors ;
            if ( ! doCanonical        ) scaledEigenVectors = RealArray2D_AllocateWithExtents ( d, d, status ) ;
            if ( Status_IsOK ( status ) )
            {
                auto Integer      i, n ;
                auto Real         eValue , result  ;
                auto RealArray1D  iColumn, nColumn ;
                /* . Diagonalization and find number of linearly-independent vectors. */
                SymmetricMatrix_EigenvaluesSolve ( S, preserveInput, 0, d, eigenValuesL, eigenVectorsL, False, status ) ;
# ifdef _CheckEigenvectors
printf ( "\nIntegrals:\n" ) ; fflush ( stdout ) ;
SymmetricMatrix_Print ( S ) ; fflush ( stdout ) ;
CheckEigenvectors ( d, eigenValuesL, eigenVectorsL ) ; fflush ( stdout ) ;
printf ( "\nEigenvalues:\n" ) ; fflush ( stdout ) ;
RealArray1D_Print ( eigenValuesL ) ; fflush ( stdout ) ;
printf ( "\nEigenvectors:\n" ) ; fflush ( stdout ) ;
RealArray2D_Print ( eigenVectorsL ) ; fflush ( stdout ) ;
# endif
                for ( i = 0, n = 0 ; i < d ; i++ )
                {
                    eValue = Array1D_Item ( eigenValuesL, i ) ;
                    if ( eValue > eTolerance ) { result =  1.0e+00 / sqrt ( eValue ) ; n++ ; }
                    else                       { result = -1.0e+00 ; }
                    Array1D_Item ( inverseEigenValues, i ) = result ;
                }
                numberOfVectors = n ;
                /* . Check for error. */
                if ( ( n == 0 ) || ( ( ! doCanonical ) && ( n != d ) ) ) { Status_Set ( status, Status_AlgorithmError ) ; }
                /* . Canonical orthogonalization. */
                else if ( doCanonical )
                {
                    for ( i = 0, n = 0 ; i < d ; i++ )
                    {
                        eValue = Array1D_Item ( inverseEigenValues, i ) ;
                        if ( eValue > 0.0e+00 )
                        {
                            RealArray2D_ColumnView ( eigenVectorsL, i, False, &iColumn, NULL ) ;
                            RealArray2D_ColumnView ( X            , n, False, &nColumn, NULL ) ;
                            RealArray1D_CopyTo     ( &iColumn     ,           &nColumn, NULL ) ;
                            RealArray1D_Scale      ( &nColumn     , eValue                   ) ;
                            n++ ;
                        }
                    }
                }
                /* . Symmetric orthogonalization. */
                else
                {
                    RealArray2D_CopyTo ( eigenVectorsL, scaledEigenVectors, NULL ) ;
                    for ( i = 0 ; i < d ; i++ )
                    {
                        RealArray2D_ColumnView ( scaledEigenVectors, i, False, &iColumn, NULL ) ;
                        RealArray1D_Scale ( &iColumn, Array1D_Item ( inverseEigenValues, i ) ) ;
                    }
                    RealArray2D_MatrixMultiply ( False, True, 1.0e+00, eigenVectorsL, scaledEigenVectors, 0.0e+00, X, NULL ) ;
                }
            }
            /* . Finish up. */
            RealArray1D_Deallocate ( &inverseEigenValues ) ;
            RealArray2D_Deallocate ( &scaledEigenVectors ) ;
            if ( eigenValues  == NULL ) RealArray1D_Deallocate ( &eigenValuesL  ) ;
            if ( eigenVectors == NULL ) RealArray2D_Deallocate ( &eigenVectorsL ) ;
        }
    }
    return numberOfVectors ;
}
# undef _EigenValueTolerance
