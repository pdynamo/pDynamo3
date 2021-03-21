/*==================================================================================================================================
! . Mulliken charge analysis.
!=================================================================================================================================*/

# include "Integer.h"
# include "Mulliken.h"
# include "Real.h"

/* . Basis indices refer to whatever representation the density, overlap and Fock matrices are in. */

/*----------------------------------------------------------------------------------------------------------------------------------
! . Atomic charges.
! . Charges incremented here.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Mulliken_AtomicCharges ( const IntegerArray1D  *basisIndices ,
                              const SymmetricMatrix *density      ,
                              const SymmetricMatrix *overlap      ,
                                    RealArray1D     *charges      )
{
    if ( ( basisIndices != NULL ) &&
         ( charges      != NULL ) &&
         ( density      != NULL ) &&
         ( overlap      != NULL ) )
    {
        auto Integer i, n, u, u0, u1, v ;
        auto Real    ps ;
        n = SymmetricMatrix_Extent ( density ) ;
        for ( i = 0 ; i < View1D_Extent ( charges ) ; i++ )
        {
            u0 = Array1D_Item ( basisIndices, i   ) ;
            u1 = Array1D_Item ( basisIndices, i+1 ) ;
            for ( u = u0, ps = 0.0e+00 ; u < u1 ; u++ )
            {
                for ( v = 0 ; v < u ; v++ ) ps += ( SymmetricMatrix_Item ( density, u, v ) * SymmetricMatrix_Item ( overlap, u, v ) ) ;
                for ( v = u ; v < n ; v++ ) ps += ( SymmetricMatrix_Item ( density, v, u ) * SymmetricMatrix_Item ( overlap, v, u ) ) ;
            }
            Array1D_Item ( charges, i ) -= ps ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Bond orders.
! . Bond orders incremented here.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Mulliken_BondOrders ( const IntegerArray1D  *basisIndices ,
                           const SymmetricMatrix *density      ,
                           const SymmetricMatrix *overlap      ,
                                 SymmetricMatrix *bondOrders   ,
                                 Status          *status       )
{
    if ( ( basisIndices != NULL ) &&
         ( bondOrders   != NULL ) &&
         ( density      != NULL ) &&
         ( overlap      != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto Integer      n  = SymmetricMatrix_Extent ( overlap ) ;
        auto RealArray2D *ps = NULL ;
        ps = RealArray2D_AllocateWithExtents    ( n, n, status ) ;
        SymmetricMatrix_SymmetricMatrixMultiply ( density, overlap, ps, status ) ;
        if ( Status_IsOK ( status ) )
        {
            auto Integer i, j, u, u0, u1, v, v0, v1 ;
            auto Real    sum ;
            for ( i = 0 ; i < SymmetricMatrix_Extent ( bondOrders ) ; i++ )
            {
                u0 = Array1D_Item ( basisIndices, i   ) ;
                u1 = Array1D_Item ( basisIndices, i+1 ) ;
                /* . Off-diagonal and diagonal. */
                for ( j = 0 ; j <= i ; j++ )
                {
                    v0 = Array1D_Item ( basisIndices, j   ) ;
                    v1 = Array1D_Item ( basisIndices, j+1 ) ;
                    for ( u = u0, sum = 0.0e+00 ; u < u1 ; u++ )
                    {
                        for ( v = v0 ; v < v1 ; v++ ) sum += ( Array2D_Item ( ps, u, v ) * Array2D_Item ( ps, v, u ) ) ;
                    }
                    SymmetricMatrix_Item ( bondOrders, i, j ) += sum ;
                }
            }
        }
        RealArray2D_Deallocate ( &ps ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Charge density derivatives.
! . Fock incremented here only.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Mulliken_ChargeDensityDerivatives ( const IntegerArray1D  *basisIndices   ,
                                         const RealArray1D     *potentials     , /* . = dXdQ. */
                                         const SymmetricMatrix *overlap        ,
                                               SymmetricMatrix *fock           )
{
    if ( ( basisIndices != NULL ) &&
         ( fock         != NULL ) &&
         ( overlap      != NULL ) &&
         ( potentials   != NULL ) )
    {
        auto Integer i, j, u, u0, u1, v, v0, v1 ;
        auto Real    p, pq, q ;
        for ( i = 0 ; i < View1D_Extent ( potentials ) ; i++ )
        {
            p  = Array1D_Item ( potentials  , i   ) ;
            u0 = Array1D_Item ( basisIndices, i   ) ;
            u1 = Array1D_Item ( basisIndices, i+1 ) ;
            /* . Off-diagonal. */
            for ( j = 0 ; j < i ; j++ )
            {
                q  = Array1D_Item ( potentials  , j   ) ;
                v0 = Array1D_Item ( basisIndices, j   ) ;
                v1 = Array1D_Item ( basisIndices, j+1 ) ;
                pq = 0.5e+00 * ( p + q ) ;
                for ( u = u0 ; u < u1 ; u++ )
                {
                    for ( v = v0 ; v < v1 ; v++ ) SymmetricMatrix_Item ( fock, u, v ) -= ( pq * SymmetricMatrix_Item ( overlap, u, v ) ) ;
                }
            }
            /* . Diagonal. */
            for ( u = u0 ; u < u1 ; u++ )
            {
                for ( v = u0 ; v <= u ; v++ ) SymmetricMatrix_Item ( fock, u, v ) -= ( p * SymmetricMatrix_Item ( overlap, u, v ) ) ;
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Weighted density matrix.
! . wDensity incremented here only.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Mulliken_WeightedDensity ( const IntegerArray1D  *basisIndices   ,
                                const RealArray1D     *potentials     , /* . = dXdQ. */
                                const SymmetricMatrix *density        ,
                                      SymmetricMatrix *wDensity       )
{
    if ( ( basisIndices != NULL ) &&
         ( density      != NULL ) &&
         ( potentials   != NULL ) &&
         ( wDensity     != NULL ) )
    {
        auto Integer i, j, u, u0, u1, v, v0, v1 ;
        auto Real    p, pq, q ;
        for ( i = 0 ; i < View1D_Extent ( potentials ) ; i++ )
        {
            p  = Array1D_Item ( potentials  , i   ) ;
            u0 = Array1D_Item ( basisIndices, i   ) ;
            u1 = Array1D_Item ( basisIndices, i+1 ) ;
            /* . Off-diagonal only. */
            for ( j = 0 ; j < i ; j++ )
            {
                q  = Array1D_Item ( potentials  , j   ) ;
                v0 = Array1D_Item ( basisIndices, j   ) ;
                v1 = Array1D_Item ( basisIndices, j+1 ) ;
                pq = - ( p + q ) ;
                for ( u = u0 ; u < u1 ; u++ )
                {
                    for ( v = v0 ; v < v1 ; v++ ) SymmetricMatrix_Item ( wDensity, u, v ) += ( pq * SymmetricMatrix_Item ( density, u, v ) ) ;
                }
            }
        }
    }
}
