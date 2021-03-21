/*==================================================================================================================================
! . Functions for the transformation of Gaussian integrals.
!=================================================================================================================================*/

# include "GaussianBasisTransform.h"
# include "Integer.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Transform a matrix of integrals corresponding to the bases i and j.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasis_TransformIntegrals2 ( RealArray2D **integrals, const RealArray2D *ic2o, const RealArray2D *jc2o )
{

    if ( ( (*integrals) != NULL ) && ( ic2o != NULL ) && ( jc2o != NULL ) )
    {
        auto Integer      nCi, nOi, nOj ;
        auto RealArray2D *temp1 = NULL, *temp2 = NULL ;
        nCi   = View2D_Rows    ( ic2o ) ;
        nOi   = View2D_Columns ( ic2o ) ;
        nOj   = View2D_Columns ( jc2o ) ;
        temp1 = RealArray2D_AllocateWithExtents ( nCi, nOj, NULL ) ;
        temp2 = RealArray2D_AllocateWithExtents ( nOi, nOj, NULL ) ;
        RealArray2D_MatrixMultiply  ( False, False, 1.0e+00, (*integrals), jc2o , 0.0e+00, temp1, NULL ) ;
        RealArray2D_MatrixMultiply  ( True,  False, 1.0e+00, ic2o        , temp1, 0.0e+00, temp2, NULL ) ;
        RealArray2D_Deallocate ( integrals ) ;
        RealArray2D_Deallocate ( &temp1    ) ;
        (*integrals) = temp2 ;
    }
}
