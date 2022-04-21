/*==================================================================================================================================
! . Gaussian basis set orthonormalization.
!=================================================================================================================================*/

# include <math.h>
# include <stdarg.h>
# include <stdlib.h>

# include "GaussianBasisIntegrals_f1Xg1.h"
# include "GaussianBasisOrthonormalize.h"
# include "IntegerUtilities.h"
# include "NumericalMacros.h"
# include "RealUtilities.h"
# include "SymmetricMatrix.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Orthonormalize the basis.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _DiagonalTolerance   1.0e-10
# define _EigenValueTolerance 1.0e-30
void GaussianBasis_Orthonormalize (       GaussianBasis          *self         ,
                                    const GaussianBasisOperator   operator     ,
                                    const OrthogonalizationMethod method       ,
                                          Integer                *nIndependent ,
                                          Real                   *deviation    ,
                                          RealArray2D            *MOut         ,
                                          RealArray2D            *XOut         ,
                                          RealArray2D            *YOut         ,
                                          Status                 *status       )
{
    if ( deviation    != NULL ) (*deviation)    = -1.0e+00 ;
    if ( nIndependent != NULL ) (*nIndependent) = -1 ;
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer          d, n, nI = 0, nR = 0, s2 ;
        auto Integer         *iWork = NULL ;
        auto Real            *rWork = NULL ;
        auto RealArray2D     *M = NULL, *X = NULL, *Y = NULL ;
        auto SymmetricMatrix *D = NULL ;
        /* . Fill the primitive CCBF for the basis. */
        GaussianBasis_Finalize ( self, status ) ;
        /* . Allocate space. */
        d  = self->nBasis ;
        D  = SymmetricMatrix_AllocateWithExtent ( d, status ) ;
        if ( MOut == NULL ) M = RealArray2D_AllocateWithExtents ( d, d, status ) ;
        else                M = MOut ;
        if ( XOut == NULL ) X = RealArray2D_AllocateWithExtents ( d, d, status ) ;
        else                X = XOut ;
        if ( YOut == NULL ) Y = RealArray2D_AllocateWithExtents ( d, d, status ) ;
        else                Y = YOut ;
        n  = GaussianBasis_LargestShell ( self, True ) ;
        s2 = n*n ;
             if ( operator == GaussianBasisOperator_AntiCoulomb ) { nI = 6*s2 ; nR = 3*s2 ; }
        else if ( operator == GaussianBasisOperator_Coulomb     ) { nI = 3*s2 ; nR = 3*s2 ; }
        else if ( operator == GaussianBasisOperator_Overlap     ) { nI = 0    ; nR = 2*s2 ; }
        if ( nI > 0 ) iWork = Integer_Allocate ( nI, status ) ;
        if ( nR > 0 ) rWork = Real_Allocate    ( nR, status ) ;
        if ( ( View2D_Columns ( X ) != d ) || ( View2D_Rows ( X ) != d ) ||
             ( View2D_Columns ( Y ) != d ) || ( View2D_Rows ( Y ) != d ) ) Status_Set ( status, Status_NonConformableArrays ) ;
        if ( Status_IsOK ( status ) )
        {
            auto Integer N ;
            auto Real    dTolerance = _DiagonalTolerance   ,
                         eTolerance = _EigenValueTolerance ;
            auto Real    r[3] = { 0.0e+00, 0.0e+00, 0.0e+00 } , /* . Fix to ensure full matrix is calculated. */
                         s[3] = { 0.0e+00, 0.0e+00, 0.0e+00 } ;
            /* . Calculate metric matrix. */
                 if ( operator == GaussianBasisOperator_AntiCoulomb ) GaussianBasisIntegrals_f1Ag1i ( self, r, self, s, s2, iWork, rWork, M ) ;
            else if ( operator == GaussianBasisOperator_Coulomb     ) GaussianBasisIntegrals_f1Cg1i ( self, r, self, s, s2, iWork, rWork, M ) ;
            else if ( operator == GaussianBasisOperator_Overlap     ) GaussianBasisIntegrals_f1Og1i ( self, r, self, s, s2,        rWork, M ) ;
            else Status_Set ( status, Status_AlgorithmError ) ; /* . Not currently available! */
            SymmetricMatrix_CopyFromRealArray2D ( D, M, status ) ;
            /* . Determine X and Y - latter wasteful if N < d. */
            N = OrthogonalizingTransformation ( D, method, True, &dTolerance, &eTolerance, NULL, NULL, X, status ) ;
            SymmetricMatrix_PostMatrixMultiply ( D, X, False, Y, status ) ;
            /* . Reporting. */
            if ( nIndependent != NULL ) (*nIndependent) = N ;
            if ( deviation != NULL )
            {
                auto RealArray2D XView, YView ; /* . Views unnecessary if N = d. */
                RealArray2D_View ( X, 0, 0, d, N, 1, 1, False, &XView, status ) ;
                RealArray2D_View ( Y, 0, 0, d, N, 1, 1, False, &YView, status ) ;
                (*deviation) = CheckOrthogonalization ( &XView, &YView, status ) ;
            }
        }
        /* . Finish up. */
        Integer_Deallocate         ( &iWork ) ;
        Real_Deallocate            ( &rWork ) ;
        SymmetricMatrix_Deallocate ( &D     ) ;
        if ( MOut == NULL ) RealArray2D_Deallocate ( &M ) ;
        if ( XOut == NULL ) RealArray2D_Deallocate ( &X ) ;
        if ( YOut == NULL ) RealArray2D_Deallocate ( &Y ) ;
    }
}
# undef _DiagonalTolerance
# undef _EigenValueTolerance
