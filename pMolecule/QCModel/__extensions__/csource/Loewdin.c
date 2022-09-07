/*==================================================================================================================================
! . Loewdin charge analysis.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>

# include "Integer.h"
# include "Loewdin.h"
# include "NumericalMacros.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Atomic charges.
! . Charges incremented here.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Loewdin_AtomicCharges ( const IntegerArray1D  *basisIndices ,
                             const SymmetricMatrix *loewdinT     ,
                             const SymmetricMatrix *density      ,
                                   RealArray1D     *charges      ,
                                   Status          *status       )
{
    if ( ( basisIndices != NULL ) &&
         ( charges      != NULL ) &&
         ( density      != NULL ) &&
         ( loewdinT     != NULL ) )
    {
        auto SymmetricMatrix *ps = NULL ;
        ps = SymmetricMatrix_AllocateWithExtent ( SymmetricMatrix_Extent ( loewdinT ), status ) ;
        SymmetricMatrix_SymmetricTransform ( density, loewdinT, ps, status ) ;
        if ( Status_IsOK ( status ) )
        {
            auto Integer i, u, u0, u1 ;
            auto Real    f ;
            for ( i = 0 ; i < View1D_Extent ( charges ) ; i++ )
            {
                u0 = Array1D_Item ( basisIndices, i   ) ;
                u1 = Array1D_Item ( basisIndices, i+1 ) ;
                for ( u = u0, f = 0.0e+00 ; u < u1 ; u++ ) f += SymmetricMatrix_Item ( ps, u, u ) ;
                Array1D_Item ( charges, i ) -= f ;
            }
        }
        SymmetricMatrix_Deallocate ( &ps ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Bond orders.
! . Bond orders incremented here.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Loewdin_BondOrders ( const IntegerArray1D  *basisIndices ,
                          const SymmetricMatrix *loewdinT     ,
                          const SymmetricMatrix *density      ,
                                SymmetricMatrix *bondOrders   ,
                                Status          *status       )
{
    if ( ( basisIndices != NULL ) &&
         ( bondOrders   != NULL ) &&
         ( density      != NULL ) &&
         ( loewdinT     != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto SymmetricMatrix *ps = NULL ;
        ps = SymmetricMatrix_AllocateWithExtent ( SymmetricMatrix_Extent ( loewdinT ), status ) ;
        SymmetricMatrix_SymmetricTransform ( density, loewdinT, ps, status ) ;
        if ( Status_IsOK ( status ) )
        {
            auto Integer i, j, u, u0, u1, v, v0, v1 ;
            auto Real    sum ;
            for ( i = 0 ; i < SymmetricMatrix_Extent ( bondOrders ) ; i++ )
            {
                u0 = Array1D_Item ( basisIndices, i   ) ;
                u1 = Array1D_Item ( basisIndices, i+1 ) ;
                /* . Off-diagonal. */
                for ( j = 0 ; j < i ; j++ )
                {
                    v0 = Array1D_Item ( basisIndices, j   ) ;
                    v1 = Array1D_Item ( basisIndices, j+1 ) ;
                    for ( u = u0, sum = 0.0e+00 ; u < u1 ; u++ )
                    {
                        for ( v = v0 ; v < v1 ; v++ ) sum += pow ( SymmetricMatrix_Item ( ps, u, v ), 2 ) ;
                    }
                    SymmetricMatrix_Item ( bondOrders, i, j ) += sum ;
                }
                /* . Diagonal. */
                for ( u = u0, sum = 0.0e+00 ; u < u1 ; u++ )
                {
                    for ( v = u0 ; v < u ; v++ ) sum += ( 2.0e+00 * pow ( SymmetricMatrix_Item ( ps, u, v ), 2 ) ) ;
                    sum += pow ( SymmetricMatrix_Item ( ps, u, u ), 2 ) ;
                }
                SymmetricMatrix_Item ( bondOrders, i, i ) += sum ;
            }
        }
        SymmetricMatrix_Deallocate ( &ps ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Charge density derivatives.
! . Fock incremented here only.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Loewdin_ChargeDensityDerivatives ( const IntegerArray1D  *basisIndices   ,
                                        const RealArray1D     *potentials     , /* . = dXdQ. */
                                        const SymmetricMatrix *loewdinT       ,
                                              SymmetricMatrix *fock           )
{
    if ( ( basisIndices != NULL ) &&
         ( fock         != NULL ) &&
         ( loewdinT     != NULL ) &&
         ( potentials   != NULL ) )
    {
        auto Integer i, u, u0, u1, v, w ;
        auto Real    f, p ;
        for ( v = 0 ; v < SymmetricMatrix_Extent ( fock ) ; v++ )
        {
            for ( w = 0 ; w <= v ; w++ )
            {
                for ( i = 0, f = 0.0e+00 ; i < View1D_Extent ( potentials ) ; i++ )
                {
                    p  = Array1D_Item ( potentials  , i   ) ;
                    u0 = Array1D_Item ( basisIndices, i   ) ;
                    u1 = Array1D_Item ( basisIndices, i+1 ) ;
                    for ( u = u0 ; u < u1 ; u++ ) f += ( p * SymmetricMatrix_GetItem ( loewdinT, v, u, NULL ) * SymmetricMatrix_GetItem ( loewdinT, w, u, NULL ) ) ;
                }
                SymmetricMatrix_Item ( fock, v, w ) -= f ;
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Charge restraint W-matrix and core term.
! . Only basic checking is done.
! . The input W matrix is initialized on entry.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real Loewdin_ChargeRestraintMatrix ( const IntegerArray1D  *basisIndices   ,
                                     const RealArray1D     *nuclearCharges ,
                                     const IntegerArray1D  *crIndices      ,
                                     const RealArray1D     *crWeights      ,
                                     const Boolean          isSpin         ,
                                     const SymmetricMatrix *loewdinT       ,
                                           SymmetricMatrix *W              )
{
    Real core = 0.0e+00 ;
    if ( ( basisIndices   != NULL ) &&
         ( crIndices      != NULL ) &&
         ( crWeights      != NULL ) &&
         ( loewdinT       != NULL ) &&
         ( nuclearCharges != NULL ) &&
         ( W              != NULL ) )
    {
        auto Integer a, i, r, s, u, u0, u1 ;
        auto Real    f, w ;
        SymmetricMatrix_Set ( W, 0.0e+00 ) ;
        for ( i = 0 ; i < View1D_Extent ( crIndices ) ; i++ )
        {
            a = Array1D_Item ( crIndices, i ) ;
            w = Array1D_Item ( crWeights, i ) ;
            if ( ! isSpin )
            {    
                 core += ( w * Array1D_Item ( nuclearCharges, a ) ) ;
                 w    *= -1.0e+00 ; /* . w is -1.0 for electrons in this case. */
            }
            u0 = Array1D_Item ( basisIndices, a   ) ;
            u1 = Array1D_Item ( basisIndices, a+1 ) ;
            for ( r = 0 ; r < SymmetricMatrix_Extent ( W ); r++ )
            {
                for ( s = 0 ; s <= r ; s++ )
                {
                    for ( u = u0, f = 0.0e+00 ; u < u1 ; u++ )
                    {
                        /* . Use GetItem here as r, s, and u unordered! */
                        f += ( SymmetricMatrix_GetItem ( loewdinT, r, u, NULL ) * SymmetricMatrix_GetItem ( loewdinT, s, u, NULL ) ) ;
                    }
                    SymmetricMatrix_Item ( W, r, s ) += ( w * f ) ;
                }
            }
        }
    }
    return core ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Charge restraint weighted density.
! . This method needs to be called for each restraint separately together with the derivative of the restraint energy model
! . with respect to the restraint, dRdL, the overlap eigenvectors and the Z matrix, which is the appropriate density post-multiplied
! . by the Loewdin transformation.
! . Only a partial weighted density is formed. The full matrix needs to transformed afterwards by the overlap factors.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _Factor 2.0e+00 /* . The appropriate weight factor. */
void Loewdin_ChargeRestraintWeightedDensity ( const IntegerArray1D  *basisIndices ,
                                              const IntegerArray1D  *crIndices    ,
                                              const RealArray1D     *crWeights    ,
                                              const Boolean          isSpin       ,
                                              const Real             dRdL         ,
                                              const RealArray2D     *eigenVectors ,
                                              const RealArray2D     *Z            ,
                                                    SymmetricMatrix *A            )
{
    if ( ( basisIndices != NULL    ) &&
         ( crIndices    != NULL    ) &&
         ( crWeights    != NULL    ) &&
         ( eigenVectors != NULL    ) &&
         ( Z            != NULL    ) &&
         ( A            != NULL    ) &&
         ( dRdL         != 0.0e+00 ) ) /* . Do nothing if dRdL is zero! */
    {
        auto Integer a, i, j, n, r, u, u0, u1 ;
        auto Real    f, sum, w ;
        f = _Factor * dRdL ;
        if ( ! isSpin ) f *= -1.0e+00 ;
        n = View2D_Columns ( eigenVectors ) ;
        for ( r = 0 ; r < View1D_Extent ( crIndices ) ; r++ )
        {
            a =     Array1D_Item ( crIndices, r ) ;
            w = f * Array1D_Item ( crWeights, r ) ;
            u0 = Array1D_Item ( basisIndices, a   ) ;
            u1 = Array1D_Item ( basisIndices, a+1 ) ;
            for ( i = 0 ; i < n ; i++ )
            {
                for ( j = 0 ; j <= i ; j++ )
                {
                    sum = 0.0e+00 ;
                    for ( u = u0 ; u < u1 ; u++ )
                    {
                        sum += ( Array2D_Item ( eigenVectors, u, i ) * Array2D_Item ( Z, u, j ) +
                                 Array2D_Item ( eigenVectors, u, j ) * Array2D_Item ( Z, u, i ) ) ;
                    }
                    SymmetricMatrix_Item ( A, i, j ) += ( w * sum ) ;
                }
            }
        }
    }
}
# undef _Factor

/*----------------------------------------------------------------------------------------------------------------------------------
! . Weighted density matrix.
! . wDensity incremented here only.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _EigenValueTolerance 1.0e-10
void Loewdin_WeightedDensity ( const IntegerArray1D  *basisIndices        ,
                               const RealArray1D     *potentials          , /* . = dXdQ. */
                               const RealArray1D     *eigenValues         ,
                               const RealArray2D     *eigenVectors        ,
                               const SymmetricMatrix *loewdinT            ,
                               const SymmetricMatrix *density             ,
                                     Real            *eigenValueTolerance ,
                                     SymmetricMatrix *wDensity            ,
                                     Status          *status              )
{
    if ( ( basisIndices != NULL ) &&
         ( density      != NULL ) &&
         ( eigenValues  != NULL ) &&
         ( eigenVectors != NULL ) &&
         ( loewdinT     != NULL ) &&
         ( potentials   != NULL ) &&
         ( wDensity     != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto Boolean isOK ;
        auto Integer nA, nE, nO ;
        /* . Get and check dimensions. */
        nA   = View1D_Extent          ( potentials  ) ;
        nE   = View1D_Extent          ( eigenValues ) ;
        nO   = SymmetricMatrix_Extent ( density     ) ;
        isOK = ( View1D_Extent          ( basisIndices     ) == ( nA + 1 ) ) &&
               ( View2D_Columns         ( eigenVectors     ) ==   nE       ) &&
               ( View2D_Rows            ( eigenVectors     ) ==   nO       ) &&
               ( SymmetricMatrix_Extent ( loewdinT         ) ==   nO       ) &&
               ( Array1D_Item           ( basisIndices, nA ) ==   nO       ) &&
               ( SymmetricMatrix_Extent ( wDensity         ) ==   nO       ) ;
        if ( ! isOK ) Status_Set ( status, Status_NonConformableArrays ) ;
        if ( Status_IsOK ( status ) )
        {
            auto RealArray1D     *tempEV ;
            auto SymmetricMatrix *tempNE, *tempNO ;
            /* . Allocate space. */
            tempEV = RealArray1D_AllocateWithExtent     ( nE, status ) ;
            tempNE = SymmetricMatrix_AllocateWithExtent ( nE, status ) ;
            tempNO = SymmetricMatrix_AllocateWithExtent ( nO, status ) ;
            if ( Status_IsOK ( status ) )
            {
                auto Integer i, j, u, u0, u1, v, v0, v1, w ;
                auto Real    a, ab, b, f, tolerance ;
                if ( eigenValueTolerance == NULL ) tolerance =  _EigenValueTolerance  ;
                else                               tolerance = (*eigenValueTolerance) ;
                /* . Symmetrized core matrix. */
                for ( i = 0 ; i < View1D_Extent ( potentials ) ; i++ )
                {
                    a  = Array1D_Item ( potentials  , i   ) ;
                    u0 = Array1D_Item ( basisIndices, i   ) ;
                    u1 = Array1D_Item ( basisIndices, i+1 ) ;
                    for ( u = u0 ; u < u1 ; u++ )
                    {
                         for ( j = 0 ; j <= i ; j++ )
                        {
                            b  = Array1D_Item ( potentials  , j ) ;
                            v0 = Array1D_Item ( basisIndices, j ) ;
                            v1 = Minimum ( u, Array1D_Item ( basisIndices, j+1 ) - 1 ) ;
                            for ( v = v0 ; v <= v1 ; v++ )
                            {
                                for ( w = 0, f = 0.0e+00 ; w < nO ; w++ )
                                {
                                    f += b * SymmetricMatrix_GetItem ( density, w, u, NULL ) * SymmetricMatrix_GetItem ( loewdinT, w, v, NULL ) +
                                         a * SymmetricMatrix_GetItem ( density, w, v, NULL ) * SymmetricMatrix_GetItem ( loewdinT, w, u, NULL ) ;
                                }
                                SymmetricMatrix_Item ( tempNO, u, v ) = f ;
                            }
                        }
                    }
                }
                /* . First transformation. */
                SymmetricMatrix_Transform ( tempNO, eigenVectors, False, tempNE, status ) ;
                /* . Get the square roots of the eigenvalues. */
                for ( u = 0 ; u < nE ; u++ )
                {
                    a = Array1D_Item ( eigenValues, u ) ;
                    Array1D_Item ( tempEV, u ) = ( a < 0.0e+00 ? 0.0e+00 : sqrt ( a ) ) ;
                }
                /* . Scale by the eigenvalue factors. */
                for ( u = 0 ; u < nE ; u++ )
                {
                    a = Array1D_Item ( tempEV, u ) ;
                    for ( v = 0 ; v <= u ; v++ )
                    {
                        b  = Array1D_Item ( tempEV, v ) ;
                        ab = a + b ;
                        if ( ab > tolerance ) SymmetricMatrix_Item ( tempNE, u, v ) /= ab      ;
                        else                  SymmetricMatrix_Item ( tempNE, u, v )  = 0.0e+00 ;
                    }
                }
                /* . Second transformation. */
                /* . When rationalize Arrays can use instead:
                !    SymmetricMatrix_Transform ( wDM, tempNE, eigenVectors, True, -2.0e+00, 1.0e+00, status ) ;
                !    In this case doesn't save much as tempNO already exists.
                */
                SymmetricMatrix_Transform ( tempNE, eigenVectors, True, tempNO, status ) ;
                /* . Add in the contributions to the density. */
                SymmetricMatrix_Add ( wDensity, -2.0e+00, tempNO, status ) ;
            }
            /* . Deallocate space. */
            RealArray1D_Deallocate     ( &tempEV ) ;
            SymmetricMatrix_Deallocate ( &tempNE ) ;
            SymmetricMatrix_Deallocate ( &tempNO ) ;
        }
    }
}
# undef _EigenValueTolerance
